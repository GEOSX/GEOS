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
 *  @file VonMisesJ.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_VONMISESJ_HPP_
#define GEOS_CONSTITUTIVE_SOLID_VONMISESJ_HPP_

#include "ElasticIsotropic.hpp"
#include "InvariantDecompositions.hpp"
#include "SolidModelDiscretizationOpsIsotropic.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class VonMisesJUpdates
 *
 * Class to provide elastic isotropic material updates that may be
 * called from a kernel function.
 */
class VonMisesJUpdates : public ElasticIsotropicUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] bulkModulus  The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress    The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress    The ArrayView holding the old stress data for each quadrature point.
   * @param[in] disableInelasticity Flag to disable plasticity for inelastic models
   */
  VonMisesJUpdates( arrayView1d< real64 > const & yieldStrength,
                    arrayView3d< real64 const > const & deformationGradient,
                    arrayView3d< real64 const > const & velocityGradient,
                    arrayView3d< real64 > const & plasticStrain,
                    arrayView1d< real64 const > const & bulkModulus,
                    arrayView1d< real64 const > const & shearModulus,
                    arrayView1d< real64 const > const & thermalExpansionCoefficient,
                    arrayView3d< real64, solid::STRESS_USD > const & newStress,
                    arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                    arrayView2d< real64 > const & density,
                    arrayView2d< real64 > const & wavespeed,
                    const bool & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus, 
                             shearModulus, 
                             thermalExpansionCoefficient,
                             newStress, 
                             oldStress,
                             density,
                             wavespeed, 
                             disableInelasticity ),
    m_yieldStrength( yieldStrength ),
    m_deformationGradient( deformationGradient ),
    m_velocityGradient( velocityGradient ),
    m_plasticStrain( plasticStrain )
  {}

  /// Deleted default constructor
  VonMisesJUpdates() = delete;

  /// Default copy constructor
  VonMisesJUpdates( VonMisesJUpdates const & ) = default;

  /// Default move constructor
  VonMisesJUpdates( VonMisesJUpdates && ) = default;

  /// Deleted copy assignment operator
  VonMisesJUpdates & operator=( VonMisesJUpdates const & ) = delete;

  /// Deleted move assignment operator
  VonMisesJUpdates & operator=( VonMisesJUpdates && ) =  delete;

  /// Use the "isotropic" form of inner product compression
  using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /// Use base version of saveConvergedState
  using ElasticIsotropicUpdates::saveConvergedState;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         real64 ( &stiffness )[6][6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate( localIndex const k,
                                         localIndex const q,
                                         real64 const ( &totalStrain )[6],
                                         real64 ( &stress )[6],
                                         DiscretizationOps & stiffness ) const final;

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
                                  DiscretizationOps & stiffness ) const override;

  GEOS_HOST_DEVICE
  virtual void getElasticStiffness( localIndex const k,
                                    localIndex const q,
                                    real64 ( &stiffness )[6][6] ) const override;

  GEOS_HOST_DEVICE
  virtual void getElasticStrain( localIndex const k,
                                 localIndex const q,
                                 real64 ( &elasticStrain )[6] ) const override final;

  GEOS_HOST_DEVICE
  void computePlasticStrainIncrement ( localIndex const k,
                                       localIndex const q,
                                       const real64 timeIncrement,
                                       real64 const ( &strainIncrement )[6],
                                       real64 const ( &stressIncrement )[6],
                                       real64 ( & plasticStrainIncrement )[6] ) const;

protected:

  /// A reference to the ArrayView holding the yield strength for each element.
  arrayView1d< real64 > const m_yieldStrength;

  /// A reference to the ArrayView holding the deformation gradient for each element.
  arrayView3d< real64 const > const m_deformationGradient;
  
  /// A reference to the ArrayView holding the velocity gradient for each element.
  arrayView3d< real64 const > const m_velocityGradient;

  /// A reference to the ArrayView holding the plasticStrain for each element.
  arrayView3d< real64 > const m_plasticStrain;

};


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::getElasticStiffness( localIndex const k,
                                            localIndex const q,
                                            real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( q );
  real64 const G = m_shearModulus[k];
  real64 const lambda = conversions::bulkModAndShearMod::toFirstLame( m_bulkModulus[k], G );

  LvArray::tensorOps::fill< 6, 6 >( stiffness, 0 );

  stiffness[0][0] = lambda + 2*G;
  stiffness[0][1] = lambda;
  stiffness[0][2] = lambda;

  stiffness[1][0] = lambda;
  stiffness[1][1] = lambda + 2*G;
  stiffness[1][2] = lambda;

  stiffness[2][0] = lambda;
  stiffness[2][1] = lambda;
  stiffness[2][2] = lambda + 2*G;

  stiffness[3][3] = G;
  stiffness[4][4] = G;
  stiffness[5][5] = G;
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::getElasticStrain( localIndex const k,
                                         localIndex const q,
                                         real64 ( & elasticStrain)[6] ) const
{
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( elasticStrain );
    GEOS_ERROR( "getElasticStrain not implemented for VonMisesJ model");
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                            localIndex const q,
                                                            real64 const ( &totalStrain )[6],
                                                            real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k ); 
  GEOS_UNUSED_VAR( q );                                          
  GEOS_UNUSED_VAR( totalStrain );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "smallStrainNoStateUpdate_StressOnly is not implemented for VonMisesJ model" );
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &totalStrain )[6],
                                                 real64 ( & stress )[6],
                                                 real64 ( & stiffness )[6][6] ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  getElasticStiffness( k, q, stiffness );
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::smallStrainNoStateUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const ( &totalStrain )[6],
                                                 real64 ( & stress )[6],
                                                 DiscretizationOps & stiffness ) const
{
  smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
  stiffness.m_bulkModulus = m_bulkModulus[k];
  stiffness.m_shearModulus = m_shearModulus[k];
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & timeIncrement,
                                                     real64 const ( &strainIncrement )[6],
                                                     real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k ); 
  GEOS_UNUSED_VAR( q );                                          
  GEOS_UNUSED_VAR( timeIncrement );
  GEOS_UNUSED_VAR( strainIncrement );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "SmallStrainUpdate_StressOnly is not implemented for VonMisesJ model" );
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                     localIndex const q,
                                                     real64 const & timeIncrement,
                                                     real64 const ( & beginningRotation )[3][3],
                                                     real64 const ( & endRotation )[3][3],
                                                     real64 const ( & strainIncrement )[6],
                                                     real64 ( & stress )[6] ) const
{
  real64 previousStress[6] = { 0 };
  LvArray::tensorOps::copy< 6 >( previousStress, m_oldStress[k][q] );

  real64 trialP;
  real64 trialQ;
  real64 oldDeviatoricStress[6] = { 0 };
  twoInvariant::stressDecomposition( previousStress,
                                     trialP,
                                     trialQ,
                                     oldDeviatoricStress );
  LvArray::tensorOps::scale< 6 >( oldDeviatoricStress, LvArray::math::sqrt(2.0 / 3.0)*trialQ );

  // Exactly compute pressure 
  real64 J = LvArray::tensorOps::determinant< 3 >( m_deformationGradient[k] );
  real64 pressure = -m_bulkModulus[k] * std::log( J );

  // Hypoelastically compute deviatoric stress
  real64 rotationTranspose[3][3] = { { 0 } };
  LvArray::tensorOps::transpose< 3, 3 >( rotationTranspose, beginningRotation ); 

  real64 unrotatedVelocityGradient[3][3]  = { { 0 } };
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( unrotatedVelocityGradient, rotationTranspose, m_velocityGradient[k] );

  real64 unrotatedVelocityGradientTranspose[3][3]  = { { 0 } };
  LvArray::tensorOps::transpose< 3, 3 >( unrotatedVelocityGradientTranspose, unrotatedVelocityGradient );

  // CC: Is there an LvArray operation to get the symmetric part of a matrix?
  real64 denseD[3][3] = { { 0 } };
  LvArray::tensorOps::copy< 3, 3 >( denseD, unrotatedVelocityGradient );
  LvArray::tensorOps::add< 3, 3 >( denseD, unrotatedVelocityGradientTranspose );
  LvArray::tensorOps::scale< 3, 3 >( denseD, 0.5 );

  real64 D[6] = { 0 };
  LvArray::tensorOps::denseToSymmetric<3>( D, denseD );

  // Get deviatoric part of D
  real64 trD = LvArray::tensorOps::symTrace< 3 >( D ) / 3.0;
  D[0] -= trD;
  D[1] -= trD;
  D[2] -= trD;

  real64 deviatoricStressIncrement[6] = { 0 };
  LvArray::tensorOps::scaledCopy< 6 >( deviatoricStressIncrement, D, 2.0 * m_shearModulus[k] * timeIncrement );

  LvArray::tensorOps::copy< 6 >( stress, oldDeviatoricStress );
  LvArray::tensorOps::add< 6 >( stress, deviatoricStressIncrement );
  LvArray::tensorOps::symAddIdentity< 3 >( stress, -pressure);

  if( m_disableInelasticity )
  {
    saveStress( k, q, stress );
    return;
  }

  real64 deviator[6] = { 0 };
  twoInvariant::stressDecomposition( stress,
                                     trialP,
                                     trialQ,
                                     deviator );

  real64 yieldStrength = m_yieldStrength[k];
  if( trialQ > yieldStrength )
  {
    real64 oldStress[6] = { 0 };
    LvArray::tensorOps::copy< 6 >( oldStress, stress );

    // re-construct stress = P*eye + sqrt(2/3)*Q*nhat
    twoInvariant::stressRecomposition( trialP,
                                       yieldStrength,
                                       deviator,
                                       stress );

    real64 stressIncrement[6] = {0};
    LvArray::tensorOps::copy< 6 >( stressIncrement, stress );
    LvArray::tensorOps::subtract< 6 >( stressIncrement, oldStress );

    // Compute plastic strain increment
    real64 plasticStrainIncrement[6] = {0};
    computePlasticStrainIncrement( k,
                                   q,
                                   timeIncrement,           
                                   strainIncrement,
                                   stressIncrement,
                                   plasticStrainIncrement );

    // Increment plastic strain
    real64 oldPlasticStrain[6] = { 0 };
    LvArray::tensorOps::copy< 6 >( oldPlasticStrain, m_plasticStrain[k][q] );
    oldPlasticStrain[3] *= 0.5;
    oldPlasticStrain[4] *= 0.5;
    oldPlasticStrain[5] *= 0.5;

    real64 unrotatedOldPlasticStrain[6] = { 0 };
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( unrotatedOldPlasticStrain, rotationTranspose, oldPlasticStrain );

    unrotatedOldPlasticStrain[3] *= 2.0;
    unrotatedOldPlasticStrain[4] *= 2.0;
    unrotatedOldPlasticStrain[5] *= 2.0;

    real64 unrotatedNewPlasticStrain[6] = { 0 };
    LvArray::tensorOps::copy< 6 >( unrotatedNewPlasticStrain, unrotatedOldPlasticStrain );
    LvArray::tensorOps::add< 6 >( unrotatedNewPlasticStrain, plasticStrainIncrement );

    unrotatedNewPlasticStrain[3] *= 0.5;
    unrotatedNewPlasticStrain[4] *= 0.5;
    unrotatedNewPlasticStrain[5] *= 0.5;
    real64 newPlasticStrain[6] = { 0 };
    LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( newPlasticStrain, endRotation, unrotatedNewPlasticStrain );
    newPlasticStrain[3] *= 2.0;
    newPlasticStrain[4] *= 2.0;
    newPlasticStrain[5] *= 2.0;

    LvArray::tensorOps::copy< 6 >( m_plasticStrain[k][q], newPlasticStrain );
  }

  saveStress( k, q, stress );
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::smallStrainUpdate( localIndex const k,
                                          localIndex const q,
                                          real64 const & timeIncrement,
                                          real64 const ( &strainIncrement )[6],
                                          real64 ( & stress )[6],
                                          real64 ( & stiffness )[6][6] ) const
{
  smallStrainUpdate_StressOnly( k, 
                                q, 
                                timeIncrement,
                                strainIncrement, 
                                stress );
  getElasticStiffness( k, q, stiffness );
}


GEOS_HOST_DEVICE
inline
void VonMisesJUpdates::smallStrainUpdate( localIndex const k,
                                                 localIndex const q,
                                                 real64 const & timeIncrement,
                                                 real64 const ( &strainIncrement )[6],
                                                 real64 ( & stress )[6],
                                                 DiscretizationOps & stiffness ) const
{
  smallStrainUpdate_StressOnly( k,
                                q, 
                                timeIncrement,
                                strainIncrement, 
                                stress );
  stiffness.m_bulkModulus = m_bulkModulus[k];
  stiffness.m_shearModulus = m_shearModulus[k];
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void VonMisesJUpdates::computePlasticStrainIncrement ( localIndex const k,
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
 * @class VonMisesJ
 *
 * Class to provide an elastic isotropic material response.
 */
class VonMisesJ : public ElasticIsotropic
{
public:

  /// Alias for VonMisesJUpdates
  using KernelWrapper = VonMisesJUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  VonMisesJ( string const & name, 
                    Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~VonMisesJ() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "VonMisesJ";

  /**
   * @brief Static catalog string
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string catalogName() { return m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  ///@}

  /// Keys for data specified in this class.
  struct viewKeyStruct : public ElasticIsotropic::viewKeyStruct
  {
    /// string/key for default yield strength
    static constexpr char const * defaultYieldStrengthString() { return "defaultYieldStrength"; }

    /// string/key for yield strength
    static constexpr char const * yieldStrengthString() { return "yieldStrength"; }

    /// string/key for deformation gradient
    static constexpr char const * deformationGradientString() { return "deformationGradient"; }

    /// string/key for velocity gradient
    static constexpr char const * velocityGradientString() { return "velocityGradient"; }

    /// string/key for plastic strain
    static constexpr char const * plasticStrainString() { return "plasticStrain"; }
  };

 GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getYieldStrength() const final
  {
    return m_yieldStrength;
  }

  GEOS_HOST_DEVICE
  virtual arrayView3d< real64 const > getDeformationGradient() const final
  {
    return m_deformationGradient;
  }

  GEOS_HOST_DEVICE
  virtual arrayView3d< real64 const > getVelocityGradient() const final
  {
    return m_velocityGradient;
  }

  GEOS_HOST_DEVICE
  virtual arrayView3d< real64 const > getPlasticStrain() const final
  {
    return m_plasticStrain;
  }
  /**
   * @brief Create a instantiation of the VonMisesJUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of VonMisesJUpdate.
   */
  VonMisesJUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return VonMisesJUpdates( m_yieldStrength,
                               m_deformationGradient,
                               m_velocityGradient,
                               m_plasticStrain,
                               m_bulkModulus,
                               m_shearModulus,
                               m_thermalExpansionCoefficient,
                               m_newStress,
                               m_oldStress,
                               m_density,
                               m_wavespeed,
                               m_disableInelasticity );
    }
    else // for "no state" updates, pass empty views to avoid transfer of stress data to device
    {
      return VonMisesJUpdates( m_yieldStrength,
                               m_deformationGradient,
                               m_velocityGradient,
                               m_plasticStrain,
                               m_bulkModulus,
                               m_shearModulus,
                               m_thermalExpansionCoefficient,
                               arrayView3d< real64, solid::STRESS_USD >(),
                               arrayView3d< real64, solid::STRESS_USD >(),
                               m_density,
                               m_wavespeed,
                               m_disableInelasticity );
    }
  }

  /**
   * @brief Construct an update kernel for a derived type.
   * @tparam UPDATE_KERNEL The type of update kernel from the derived type.
   * @tparam PARAMS The parameter pack to hold the constructor parameters for
   *   the derived update kernel.
   * @param constructorParams The constructor parameter for the derived type.
   * @return An @p UPDATE_KERNEL object.
   */
  template< typename UPDATE_KERNEL, typename ... PARAMS >
  UPDATE_KERNEL createDerivedKernelUpdates( PARAMS && ... constructorParams ) const
  {
    return UPDATE_KERNEL( std::forward< PARAMS >( constructorParams )...,
                          m_yieldStrength,
                          m_deformationGradient,
                          m_velocityGradient,
                          m_plasticStrain,
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

  /// Post-process XML data
  virtual void postInputInitialization() override;

  /// The default value of the yield strength for any new allocations
  real64 m_defaultYieldStrength;

  /// The yield strength for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_yieldStrength;

  /// The deformation gradient for each upper level dimension (i.e. cell) of *this
  array3d< real64 > m_deformationGradient;

  /// The velocity gradient for each upper level dimension (i.e. cell) of *this
  array3d< real64 > m_velocityGradient;

  /// The plastic strain for each upper level dimension (i.e. cell) of *this
  array3d< real64 > m_plasticStrain;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_VONMISESJ_HPP_ */
