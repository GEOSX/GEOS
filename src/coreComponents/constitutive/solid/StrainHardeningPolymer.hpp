/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StrainHardeningPolymer.hpp
 * @brief Simple flow stress model for strain hardening polymer
 *
 * This damage model is intended for use with damage-field partitioning (DFG) within the
 * MPM solver, but can also be used without DFG by any solver. It is only appropriate for
 * schemes implementing explicit time integration. CC: TODO add description of material model
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_STRAINHARDENINGPOLYMER_HPP_
#define GEOS_CONSTITUTIVE_SOLID_STRAINHARDENINGPOLYMER_HPP_

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "PropertyConversions.hpp"
#include "SolidModelDiscretizationOpsFullyAnisotroipic.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class StrainHardeningPolymerUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class StrainHardeningPolymerUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] deformationGradient The ArrayView holding the deformation gradient for each element/particle.
   * @param[in] plasticStrain The ArrayView holding the plastic strain for each quadrature point
   * @param[in] damage The ArrayView holding the damage for each quadrature point.
   * @param[in] jacobian The ArrayView holding the jacobian for each quadrature point.
   * @param[in] yieldStrength The ArrayView holding the current yield strength
   * @param[in] strainHardeningSlope The strain hardening slope
   * @param[in] shearSofteningMagnitude The shear softening magnitude
   * @param[in] shearSofteningShapeParameter1 The shear softening shape parameter 1
   * @param[in] shearSofteningShapeParameter2 The shear softening shape parameter 2
   * @param[in] maximumStretch The maximum stretch
   * @param[in] thermalSoftening not currently implemented (CC: TODO)
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  StrainHardeningPolymerUpdates( arrayView3d< real64 > const & deformationGradient,
                                 arrayView3d< real64 > const & plasticStrain,
                                 arrayView2d< real64 > const & damage,
                                 arrayView2d< real64 > const & jacobian,
                                 arrayView1d< real64 > const & yieldStrength,
                                 real64 const & strainHardeningSlope,
                                 real64 const & shearSofteningMagnitude,
                                 real64 const & shearSofteningShapeParameter1,
                                 real64 const & shearSofteningShapeParameter2,
                                 real64 const & maximumStretch,
                                 // arrayView2d< real64 > const & thermalSoftening,
                                 arrayView1d< real64 const > const & bulkModulus,
                                 arrayView1d< real64 const > const & shearModulus,
                                 arrayView1d< real64 const > const & thermalExpansionCoefficient,
                                 arrayView3d< real64, solid::STRESS_USD > const & newStress,
                                 arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                                 arrayView2d< real64 > const & density,
                                 arrayView2d< real64 > const & wavespeed,
                                 bool const & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus,
                             shearModulus,
                             thermalExpansionCoefficient,
                             newStress,
                             oldStress,
                             density,
                             wavespeed,
                             disableInelasticity ),
    m_deformationGradient( deformationGradient ),
    m_plasticStrain( plasticStrain ),
    m_damage( damage ),
    m_jacobian( jacobian ),
    m_yieldStrength( yieldStrength ),
    m_strainHardeningSlope( strainHardeningSlope ),
    m_shearSofteningMagnitude( shearSofteningMagnitude ),
    m_shearSofteningShapeParameter1( shearSofteningShapeParameter1 ),
    m_shearSofteningShapeParameter2( shearSofteningShapeParameter2 ),
    m_maximumStretch( maximumStretch )
  {}

  /// Default copy constructor
  StrainHardeningPolymerUpdates( StrainHardeningPolymerUpdates const & ) = default;

  /// Default move constructor
  StrainHardeningPolymerUpdates( StrainHardeningPolymerUpdates && ) = default;

  /// Deleted default constructor
  StrainHardeningPolymerUpdates() = delete;

  /// Deleted copy assignment operator
  StrainHardeningPolymerUpdates & operator=( StrainHardeningPolymerUpdates const & ) = delete;

  /// Deleted move assignment operator
  StrainHardeningPolymerUpdates & operator=( StrainHardeningPolymerUpdates && ) =  delete;

  /// Use the uncompressed version of the stiffness bilinear form
  using DiscretizationOps = SolidModelDiscretizationOpsFullyAnisotroipic; // TODO: typo in anistropic (fix in DiscOps PR)

  // Bring in base implementations to prevent hiding warnings
  using ElasticIsotropicUpdates::smallStrainUpdate;

  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & timeIncrement,
                          real64 const ( &strainIncrement )[6],
                          real64 ( &stress )[6],
                          real64 ( &stiffness )[6][6] ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( &stress )[6],
                                  DiscretizationOps & stiffness ) const;

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;


  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate_StressOnly( localIndex const k,
                                             localIndex const q,
                                             real64 const & timeIncrement,
                                             real64 const ( & beginningRotation )[3][3],
                                             real64 const ( & endRotation )[3][3],
                                             real64 const ( &strainIncrement )[6],
                                             real64 ( &stress )[6] ) const override;

  GEOS_HOST_DEVICE
  void smallStrainUpdateHelper( localIndex const k,
                                localIndex const q,
                                real64 const timeIncrement,
                                real64 const ( & beginningRotation )[3][3],
                                real64 const ( & endRotation )[3][3],
                                real64 const ( &strainIncrement )[6],
                                real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  void computePlasticStrainIncrement ( localIndex const k,
                                       localIndex const q,
                                       const real64 timeIncrement,
                                       real64 const ( &strainIncrement )[6],
                                       real64 const ( &stressIncrement )[6],
                                       real64 ( & plasticStrainIncrement )[6] ) const;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the deformation gradient for each element/particle.
  arrayView3d< real64 > const m_deformationGradient;

  /// A reference to the ArrayView holding the plastic strain for each quadrature point.
  arrayView3d< real64 > const m_plasticStrain;

  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the yield strength for each element/particle
  arrayView1d< real64 > const m_yieldStrength;

  /// The strain hardening slope
  real64 const m_strainHardeningSlope;

  /// The shear softening magnitude
  real64 m_shearSofteningMagnitude;

  /// The shear softening shape parameter 1
  real64 m_shearSofteningShapeParameter1;

  /// The shear softening shape parameter 2
  real64 m_shearSofteningShapeParameter2;

  /// The compressive strength
  real64 const m_maximumStretch;

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & timeIncrement,
                                                       real64 const ( &strainIncrement )[6],
                                                       real64 ( & stress )[6],
                                                       real64 ( & stiffness )[6][6] ) const
{
  // CC: we don't call this for the MPM solver but other solvers may want to use this
  // Need to resolve if rotation matrix is always passed and update accordingly
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for StrainHardeningPolymer" );

  // // elastic predictor (assume strainIncrement is all elastic)
  // ElasticIsotropicUpdates::smallStrainUpdate( k, 
  //                                             q, 
  //                                             timeIncrement, 
  //                                             strainIncrement, 
  //                                             stress, 
  //                                             stiffness );
  // m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  // if( m_disableInelasticity )
  // {
  //   return;
  // }

  // // call the constitutive model
  // StrainHardeningPolymerUpdates::smallStrainUpdateHelper( k, 
  //                                                         q, 
  //                                                         timeIncrement,
  //                                                         strainIncrement, 
  //                                                         stress );

  // // It doesn't make sense to modify stiffness with this model

  // // save new stress and return
  // saveStress( k, q, stress );
  // return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate( localIndex const k,
                                                       localIndex const q,
                                                       real64 const & timeIncrement,
                                                       real64 const ( &strainIncrement )[6],
                                                       real64 ( & stress )[6],
                                                       DiscretizationOps & stiffness ) const
{
  // CC: we don't call this for the MPM solver but other solvers may want to use this
  // Need to resolve if rotation matrix is always passed and update accordingly
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "smallStrainUpdate not implemented for StrainHardeningPolymer" );

  // smallStrainUpdate( k, 
  //                    q, 
  //                    timeIncrement, 
  //                    strainIncrement, 
  //                    stress, 
  //                    stiffness.m_c );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                         localIndex const q,
                                                         real64 const & timeIncrement,
                                                         real64 const ( & strainIncrement )[6],
                                                         real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "smallStrainUpdateStressOnly overload not implemented for CeramicDamage" );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                                  localIndex const q,
                                                                  real64 const & timeIncrement,
                                                                  real64 const ( & beginningRotation )[3][3],
                                                                  real64 const ( & endRotation )[3][3],
                                                                  real64 const ( & strainIncrement )[6],
                                                                  real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k, 
                                                         q, 
                                                         timeIncrement,
                                                         strainIncrement, 
                                                         stress );
  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // call the constitutive model
  StrainHardeningPolymerUpdates::smallStrainUpdateHelper( k, 
                                                          q, 
                                                          timeIncrement,
                                                          beginningRotation,
                                                          endRotation,
                                                          strainIncrement, 
                                                          stress );

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::smallStrainUpdateHelper( localIndex const k,
                                                             localIndex const q,
                                                             real64 const timeIncrement,
                                                             real64 const ( & beginningRotation )[3][3],
                                                             real64 const ( & endRotation )[3][3],
                                                             real64 const ( & strainIncrement )[6],
                                                             real64 ( & stress )[6] ) const
{
  // real64 yieldStrength = m_yieldStrength[k];

  // Store old stress for plastic strain increment
  real64 oldStress[6] = { 0 };
  LvArray::tensorOps::copy< 6 >( oldStress, stress);

  // decompose into mean (P) and von Mises (Q) stress invariants
  real64 trialP;
  real64 trialQ;
  real64 deviator[6] = { 0 };
  twoInvariant::stressDecomposition( oldStress,
                                     trialP,
                                     trialQ,
                                     deviator );

    // CC: model needs the unrotated deformation gradient
    // Right stretch tensor
    real64 rotationTranspose[3][3];
    LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation );

    real64 oldPlasticStrain[6] = { 0 };
    LvArray::tensorOps::copy< 6 >(oldPlasticStrain, m_plasticStrain[k][q]);
    oldPlasticStrain[3] *= 0.5;
    oldPlasticStrain[4] *= 0.5;
    oldPlasticStrain[5] *= 0.5;

    real64 unrotatedOldPlasticStrain[6] = { 0 };
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >(unrotatedOldPlasticStrain, rotationTranspose, oldPlasticStrain);

    unrotatedOldPlasticStrain[3] *= 2.0;
    unrotatedOldPlasticStrain[4] *= 2.0;
    unrotatedOldPlasticStrain[5] *= 2.0;

    real64 unrotatedDeformationGradient[3][3] = { { 0 } };
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3>( unrotatedDeformationGradient, rotationTranspose, m_deformationGradient[k] );

    real64 U[6] = { 0.0 };
    LvArray::tensorOps::denseToSymmetric< 3 >( U, unrotatedDeformationGradient );

    real64 stretch[3] = { 0 };
    real64 eigenVectors[3][3] = { { 0 } };
    LvArray::tensorOps::symEigenvectors< 3 >( stretch, eigenVectors, U );

    // Find the largest eigenvalues
    real64 maximumStretch = 0.0;
    for( localIndex i = 0; i < 3; ++i )
    {
        maximumStretch = LvArray::math::max( stretch[i], maximumStretch );
    }

    if(maximumStretch > m_maximumStretch)
    {
        m_damage[k][q] = 1.0;
    }

    // Return to yield surface requires iterative solution
    // Implemented fixed points, however a newton solver may be more efficient and applicable
    real64 tol = 1e-10; // CC: need to experiment with these for the best options
    int maxEvals = 100; // Same cas above

    real64 yieldStrength = m_yieldStrength[k];
    real64 oldYieldStrength = yieldStrength;
    real64 unrotatedTempPlasticStrain[6] = { 0 };
    real64 plasticStrainIncrement[6] = { 0 };
    for(int iter=0; iter < maxEvals; ++iter)
    {
      LvArray::tensorOps::copy< 6 >(unrotatedTempPlasticStrain, unrotatedOldPlasticStrain);
      LvArray::tensorOps::add< 6 >(unrotatedTempPlasticStrain, plasticStrainIncrement);
      
      // Compute magnitude of plastic strain tensor
      real64 gamma_p = 0.0;
      for( int i = 0; i < 6; i++ )
      {
        gamma_p += 0.5*( 1 + (i < 3) ) * unrotatedTempPlasticStrain[i] * unrotatedTempPlasticStrain[i];
      }
      gamma_p = LvArray::math::sqrt( gamma_p );
      
      // This term starts at value r0 and decays with plastic shear strain to give plastic softening.
      // Put in a check to prevent roundoff error.
      real64 gamma_by_r1_to_r2 = powf( gamma_p / m_shearSofteningShapeParameter1, m_shearSofteningShapeParameter2 );

      // Compute change in yield strength
      real64 plasticSoftening = m_shearSofteningMagnitude * LvArray::math::exp( LvArray::math::max( -1.0 * gamma_by_r1_to_r2, -16.0 ) );
      real64 stretchHardening = m_strainHardeningSlope * ( maximumStretch * maximumStretch - 1.0 / maximumStretch );
      yieldStrength = m_yieldStrength[k] + plasticSoftening + stretchHardening; // CC: debugging disabling change in yield strength

      // // CC: need to add this later
      // real64 thermalStrengthReduction = 1.0;
      // if(m_thermalSoftening)
      // {
      //   thermalStrengthReduction = computeThermalStrengthReduction();
      // }
      // m_yieldStrength[k] *= thermalStrengthReduction;

      // check yield function
      if( trialQ > yieldStrength || iter > 0 ){

        // re-construct stress = P*eye + sqrt(2/3)*Q*nhat
        real64 stressTemp[6] = {0};
        twoInvariant::stressRecomposition( trialP,
                                           yieldStrength,
                                           deviator,
                                           stressTemp );

        // Increment plastic strain
        real64 stressIncrement[6] = {0};
        LvArray::tensorOps::copy< 6 >(stressIncrement, stressTemp);
        LvArray::tensorOps::subtract< 6 >(stressIncrement, oldStress);

        // increment plastic strain
        computePlasticStrainIncrement( k,
                                       q,
                                       timeIncrement,           
                                       strainIncrement,
                                       stressIncrement,
                                       plasticStrainIncrement );

        real64 unrotatedNewPlasticStrain[6] = { 0 };
        LvArray::tensorOps::copy< 6 >(unrotatedNewPlasticStrain, unrotatedOldPlasticStrain);
        LvArray::tensorOps::add< 6 >(unrotatedNewPlasticStrain, plasticStrainIncrement);

        if(abs(yieldStrength - oldYieldStrength) < tol)
        {
          unrotatedNewPlasticStrain[3] *= 0.5;
          unrotatedNewPlasticStrain[4] *= 0.5;
          unrotatedNewPlasticStrain[5] *= 0.5;
          real64 newPlasticStrain[6] = { 0 };
          LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >(newPlasticStrain, endRotation, unrotatedNewPlasticStrain);
          newPlasticStrain[3] *= 2.0;
          newPlasticStrain[4] *= 2.0;
          newPlasticStrain[5] *= 2.0;

          LvArray::tensorOps::copy< 6 >(m_plasticStrain[k][q], newPlasticStrain);
          LvArray::tensorOps::copy< 6 >(stress, stressTemp);
          return;
        }
        else
        {
          oldYieldStrength = yieldStrength;
        }
      } 
      else 
      {
        return;
      }
    }

    GEOS_ERROR("Plastic strain of StrainHardeningPolymer model did not converge within max evals.");
}


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void StrainHardeningPolymerUpdates::computePlasticStrainIncrement ( localIndex const k,
                                                                    localIndex const q,
                                                                    const real64 timeIncrement,
                                                                    real64 const ( & strainIncrement )[6],
                                                                    real64 const ( & stressIncrement )[6],
                                                                    real64 ( & plasticStrainIncrement )[6] ) const
{ 
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( timeIncrement );
  
  // For hypo-elastic models we compute the increment in plastic strain assuming
  // for some increment in total strain and stress and elastic properties.

  // Isotroptic-deviatoric decomposition;
  real64 trialP;
  real64 trialQ;
  real64 stressIncrementDeviator[6];
  twoInvariant::stressDecomposition( stressIncrement,
                                     trialP,
                                     trialQ,
                                     stressIncrementDeviator );

  real64 stressIncrementIsostatic[6] = {0};
  stressIncrementIsostatic[0] = trialP;
  stressIncrementIsostatic[1] = trialP;
  stressIncrementIsostatic[2] = trialP;

  // For damage or softening it there may be cases where bulk or shear are approx 0, 
  // so we need to be careful that we compute this
  real64 elasticStrainIncrement[6] = {0};
  for( int i = 0; i < 6; ++i )
  {
    if (m_bulkModulus[k] > 1.0e-12)
    {
      // CC: off diagonal elements need x2 for strain
      elasticStrainIncrement[i] += ( 1 + (i >= 3) ) * stressIncrementIsostatic[i] * 1.0/3.0/m_bulkModulus[k];
    }
    if (m_shearModulus[k] > 1.0e-12)
    {
      elasticStrainIncrement[i] += ( 1 + (i >= 3) ) * LvArray::math::sqrt(2/3) * trialQ * stressIncrementDeviator[i] * 1.0/2.0/m_shearModulus[k];
    }
  }

  LvArray::tensorOps::copy< 6 >( plasticStrainIncrement, strainIncrement);
  LvArray::tensorOps::subtract< 6 >( plasticStrainIncrement, elasticStrainIncrement);
}


/**
 * @class StrainHardeningPolymer
 *
 * Strain hardening polymer material model.
 */
class StrainHardeningPolymer : public ElasticIsotropic
{
public:

  /// @typedef Alias for StrainHardeningPolymerUpdates
  using KernelWrapper = StrainHardeningPolymerUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  StrainHardeningPolymer( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~StrainHardeningPolymer() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "StrainHardeningPolymer";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static string catalogName() { return m_catalogNameString; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /**
   * Keys for data specified in this class.
   */
  struct viewKeyStruct : public SolidBase::viewKeyStruct
  {
    // string/key for element/particle deformation gradient value
    static constexpr char const * deformationGradientString() { return "deformationGradient"; }

    // string/key for quadrature point plastic strain value
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }

    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for yield strength
    static constexpr char const * yieldStrengthString() { return "yieldStrength"; }

    /// string/key for strain hardening slope
    static constexpr char const * strainHardeningSlopeString() { return "strainHardeningSlope"; }

    /// string/key for shear softening magnitude
    static constexpr char const * shearSofteningMagnitudeString() { return "shearSofteningMagnitude"; }

    /// string/key for shear softening shape parameter 1
    static constexpr char const * shearSofteningShapeParameter1String() { return "shearSofteningShapeParameter1"; }

    /// string/key for shear softening shape parameter 2
    static constexpr char const * shearSofteningShapeParameter2String() { return "shearSofteningShapeParameter2"; }

    /// string/key for default yield strength
    static constexpr char const * defaultYieldStrengthString() { return "defaultYieldStrength"; }

    /// string/key for maximum stretch
    static constexpr char const * maximumStretchString() { return "maximumStretch"; }
  };

  /**
   * @brief Create a instantiation of the StrainHardeningPolymerUpdates class that refers to the data in this.
   * @return An instantiation of StrainHardeningPolymerUpdates.
   */
  StrainHardeningPolymerUpdates createKernelUpdates() const
  {
    return StrainHardeningPolymerUpdates( m_deformationGradient,
                                          m_plasticStrain,
                                          m_damage,
                                          m_jacobian,
                                          m_yieldStrength,
                                          m_strainHardeningSlope,
                                          m_shearSofteningMagnitude,
                                          m_shearSofteningShapeParameter1,
                                          m_shearSofteningShapeParameter2,
                                          m_maximumStretch,
                                          m_bulkModulus,
                                          m_shearModulus,
                                          m_thermalExpansionCoefficient,
                                          m_newStress,
                                          m_oldStress,
                                          m_density,
                                          m_wavespeed,
                                          m_disableInelasticity );
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams )
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_deformationGradient,
                          m_plasticStrain,
                          m_damage,
                          m_jacobian,
                          m_yieldStrength,
                          m_strainHardeningSlope,
                          m_shearSofteningMagnitude,
                          m_shearSofteningShapeParameter1,
                          m_shearSofteningShapeParameter2,
                          m_maximumStretch,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_density,
                          m_wavespeed,
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;
  /// State variable: The deformation gradient values for each element/particle.
  array3d< real64 > m_deformationGradient;

  /// State variable: The plastic strain values for each quadrature point
  array3d< real64 > m_plasticStrain;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The jacobian of the deformation gradient for each quadrature point
  array2d< real64 > m_jacobian;

  /// State variable: The yield strength 
  array1d< real64 > m_yieldStrength;

  /// Material parameter: The value of strain hardening slope
  real64 m_strainHardeningSlope;

  /// Material parameter: The value of shear softening magnitude
  real64 m_shearSofteningMagnitude;

  /// Material parameter: The value of shear softening shape parameter 1
  real64 m_shearSofteningShapeParameter1;

  /// Material parameter: The value of shear softening shape parameter 2
  real64 m_shearSofteningShapeParameter2;

  /// Material parameter: The value of default yield strength
  real64 m_defaultYieldStrength;

  /// Material parameter: The value of maximum theoretical strength
  real64 m_maximumStretch;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_STRAINHARDENINGPOLYMER_HPP_ */
