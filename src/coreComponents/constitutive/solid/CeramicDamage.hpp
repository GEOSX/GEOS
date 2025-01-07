/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CeramicDamage.hpp
 * @brief Simple damage model for modeling material failure in brittle materials.
 *
 * This damage model is intended for use with damage-field partitioning (DFG) within the
 * MPM solver, but can also be used without DFG by any solver. It is only appropriate for
 * schemes implementing explicit time integration. The model is really a hybrid plasticity/
 * damage model in the sense that we assume damaged material behaves like granular material
 * and hence follows a modified Mohr-Coulomb law. The modifications are that at low pressures,
 * the shape of the yield surface is modified to resemble a maximum principal stress criterion,
 * and at high pressures the shape converges on the von Mises yield surface. The damage
 * parameter results in softening of the deviatoric response i.e. causes the yield surface to
 * contract. Furthermore, damage is used to scale back tensile pressure: p = (1 - d) * pTrial.
 * pTrial is calculatd as pTrial = -k * log(J), where the Jacobian J of the material motion is
 * integrated and tracked by this model.
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP
#define GEOS_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP

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
 * @class CeramicDamageUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class CeramicDamageUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] damage The ArrayView holding the damage for each quardrature point.
   * @param[in] jacobian The ArrayView holding the jacobian for each quardrature point.
   * @param[in] lengthScale The ArrayView holding the length scale for each element.
   * @param[in] tensileStrength The unconfined tensile strength.
   * @param[in] compressiveStrength The unconfined compressive strength.
   * @param[in] maximumStrength The theoretical maximum strength.
   * @param[in] crackSpeed The crack speed.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  CeramicDamageUpdates( arrayView2d< real64 > const & damage,
                        arrayView2d< real64 > const & jacobian,
                        arrayView1d< real64 > const & lengthScale,
                        real64 const & tensileStrength,
                        real64 const & compressiveStrength,
                        real64 const & maximumStrength,
                        real64 const & crackSpeed,
                        arrayView1d< real64 const > const & bulkModulus,
                        arrayView1d< real64 const > const & shearModulus,
                        arrayView1d< real64 const > const & thermalExpansionCoefficient,
                        arrayView3d< real64, solid::STRESS_USD > const & newStress,
                        arrayView3d< real64, solid::STRESS_USD > const & oldStress,
                        bool const & disableInelasticity ):
    ElasticIsotropicUpdates( bulkModulus, shearModulus, thermalExpansionCoefficient, newStress, oldStress, disableInelasticity ),
    m_damage( damage ),
    m_jacobian( jacobian ),
    m_lengthScale( lengthScale ),
    m_tensileStrength( tensileStrength ),
    m_compressiveStrength( compressiveStrength ),
    m_maximumStrength( maximumStrength ),
    m_crackSpeed( crackSpeed )
  {}

  /// Default copy constructor
  CeramicDamageUpdates( CeramicDamageUpdates const & ) = default;

  /// Default move constructor
  CeramicDamageUpdates( CeramicDamageUpdates && ) = default;

  /// Deleted default constructor
  CeramicDamageUpdates() = delete;

  /// Deleted copy assignment operator
  CeramicDamageUpdates & operator=( CeramicDamageUpdates const & ) = delete;

  /// Deleted move assignment operator
  CeramicDamageUpdates & operator=( CeramicDamageUpdates && ) =  delete;

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
  void smallStrainUpdateHelper( localIndex const k,
                                localIndex const q,
                                real64 const dt,
                                real64 ( &stress )[6] ) const;

  GEOS_HOST_DEVICE
  real64 getStrength( const real64 damage,      // damage
                      const real64 pressure,    // pressure
                      const real64 J2,          // J2 invariant of stress
                      const real64 J3,          // J3 invariant of stress
                      const real64 mu,          // friction slope
                      const real64 Yt0 ) const; // strength parameter

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the damage for each quadrature point.
  arrayView2d< real64 > const m_damage;

  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

  /// A reference to the ArrayView holding the length scale for each element/particle.
  arrayView1d< real64 > const m_lengthScale;

  /// The tensile strength
  real64 const m_tensileStrength;

  /// The compressive strength
  real64 const m_compressiveStrength;

  /// The maximum theoretical strength
  real64 const m_maximumStrength;

  // The crack speed
  real64 const m_crackSpeed;
};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const & timeIncrement,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              real64 ( & stiffness )[6][6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );
  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // call the constitutive model
  CeramicDamageUpdates::smallStrainUpdateHelper( k, q, timeIncrement, stress );

  // It doesn't make sense to modify stiffness with this model

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate( localIndex const k,
                                              localIndex const q,
                                              real64 const & timeIncrement,
                                              real64 const ( &strainIncrement )[6],
                                              real64 ( & stress )[6],
                                              DiscretizationOps & stiffness ) const
{
  smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness.m_c );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                                         localIndex const q,
                                                         real64 const & timeIncrement,
                                                         real64 const ( &strainIncrement )[6],
                                                         real64 ( & stress )[6] ) const
{
  // elastic predictor (assume strainIncrement is all elastic)
  ElasticIsotropicUpdates::smallStrainUpdate_StressOnly( k, q, timeIncrement, strainIncrement, stress );
  m_jacobian[k][q] *= exp( strainIncrement[0] + strainIncrement[1] + strainIncrement[2] );

  if( m_disableInelasticity )
  {
    return;
  }

  // call the constitutive model
  CeramicDamageUpdates::smallStrainUpdateHelper( k, q, timeIncrement, stress );

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void CeramicDamageUpdates::smallStrainUpdateHelper( localIndex const k,
                                                    localIndex const q,
                                                    real64 const dt,
                                                    real64 ( & stress )[6] ) const
{
  // get failure time
  real64 tFail = m_lengthScale[k] / m_crackSpeed;

  // get trial pressure
  real64 trialPressure = -m_bulkModulus[k] * log( m_jacobian[k][q] );
  real64 pressure = trialPressure;

  // cohesion slope
  real64 mu = 0.5773502691896258;

  // Intermediate strength parameter
  real64 Yt0 = fmax( 0.5 * m_tensileStrength, fmin( 2.0 * m_tensileStrength, (3.0 * m_compressiveStrength * m_tensileStrength) / (2.0 * m_compressiveStrength + m_tensileStrength) ) );
  Yt0 = fmin( Yt0, ( 3.0 * m_compressiveStrength - m_compressiveStrength * mu ) / ( 3 + mu ) );

  // Compute the vertex pressure (should be pmin0 < 0) for the undamaged yield surface.
  real64 pmin0 = -( 2.0 * m_compressiveStrength * Yt0 ) / ( 3.0 * ( m_compressiveStrength - Yt0 ) );
  pmin0 = fmin( pmin0, -1.0e-12 );
  real64 pmin = ( 1.0 - m_damage[k][q] ) * pmin0;

  // Enforce vertex solution
  if( trialPressure < 0 )
  {
    pressure = trialPressure * ( 1.0 - m_damage[k][q] ); // Tensile cutoff pressure (negative value in tension), scaled by damage. Goes to 0
                                                         // as damage -> 1.
  }
  if( pressure < pmin ) // TODO: pressure or trial pressure?
  {
    // Increment damage
    m_damage[k][q] = fmin( m_damage[k][q] + dt / tFail, 1.0 );

    // Pressure is on the vertex
    pressure = ( 1.0 - m_damage[k][q] ) * pmin0;

    // updated stress is isotropic, at the vertex:
    stress[0] = -pressure;
    stress[1] = -pressure;
    stress[2] = -pressure;
    stress[3] = 0.0;
    stress[4] = 0.0;
    stress[5] = 0.0;
  }
  else // Enforce strength solution
  {
    real64 meanStress;    // negative of pressure
    real64 vonMises;      // von Mises stress
    real64 deviator[6];   // direction of stress deviator
    twoInvariant::stressDecomposition( stress,
                                       meanStress,
                                       vonMises,
                                       deviator );

    real64 brittleDuctileTransitionPressure = m_maximumStrength / mu;
    real64 J2 = vonMises * vonMises / 3.0;
    real64 J3 = vonMises * vonMises * vonMises *
                ( deviator[0] * deviator[1] * deviator[2] +
                  2.0 * deviator[3] * deviator[4] * deviator[5] -
                  deviator[0] * deviator[3] * deviator[3] -
                  deviator[1] * deviator[4] * deviator[4] -
                  deviator[2] * deviator[5] * deviator[5] );

    // Find the strength
    real64 strength = CeramicDamageUpdates::getStrength( m_damage[k][q], pressure, J2, J3, mu, Yt0 );

    // Increment damage and get new associated yield surface
    real64 newDeviatorMagnitude = vonMises;
    if( vonMises > strength )
    {
      if( pressure <= brittleDuctileTransitionPressure )
      {
        m_damage[k][q] = fmin( m_damage[k][q] + dt / tFail, 1.0 );
        strength = CeramicDamageUpdates::getStrength( m_damage[k][q], pressure, J2, J3, mu, Yt0 );
      }
      newDeviatorMagnitude = strength;
    }

    // Radial return
    twoInvariant::stressRecomposition( -pressure,
                                       newDeviatorMagnitude,
                                       deviator,
                                       stress );
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 CeramicDamageUpdates::getStrength( const real64 damage,     // damage
                                          const real64 pressure,   // pressure
                                          const real64 J2,         // J2 invariant of stress
                                          const real64 J3,         // J3 invariant of stress
                                          const real64 mu,         // friction slope
                                          const real64 Yt0 ) const // strength parameter
{
  auto ceramicY10 = [&]( const real64 pLocal,   // pressure
                         const real64 dLocal,   // damage,
                         const real64 muLocal,  // friction slope
                         const real64 Yt0Local, // strength parameter
                         real64 & f )      // OUTPUT: Yield surface
  {
    f =
      (((3.0 + dLocal * (-3.0 + muLocal)) * m_compressiveStrength + (-3.0 + dLocal * (3.0 + muLocal)) * Yt0Local) *
       (pLocal - (2.0 * (dLocal - 1.0) * m_compressiveStrength * Yt0Local) / (3.0 * (m_compressiveStrength - Yt0Local)))) / (m_compressiveStrength + Yt0Local);
  };

  // Bounding pressure values
  real64 p1 = m_compressiveStrength / 3.0;
  real64 p2 = m_maximumStrength / mu;

  // Get the scaling associated with p1
  real64 dfdp1 = ((3.0 + damage * (-3.0 + mu)) * m_compressiveStrength + (-3.0 + damage * (3.0 + mu)) * Yt0) / (m_compressiveStrength + Yt0);
  real64 gamma1 = 1.0;
  if( J2 > 1e-32 )
  {
    real64 psi = fmin( 2.0, fmax( 0.5, 1.0 / (1.0 + dfdp1 / 3.0)));
    real64 theta = (1.0 / 3.0) * asin( fmin( 1.0, fmax( -1.0, -0.5 * J3 * pow( 3.0 / J2, 1.5 ))));
    real64 cosPi6plusTheta = cos( 0.5235987755982989 + theta );
    real64 num = (2.0 * psi - 1.0) * (2.0 * psi - 1.0) + 4.0 * (1.0 - psi * psi) * cosPi6plusTheta * cosPi6plusTheta;
    real64 denom = 2.0 * (1.0 - psi * psi) * cosPi6plusTheta + (2.0 * psi - 1.0) * sqrt( fmax( 0.0, -4.0 * psi + 5.0 * psi * psi + 4.0 * (1.0 - psi * psi) * cosPi6plusTheta * cosPi6plusTheta ));
    if( denom > 1.0e-12 )
    {
      gamma1 = (1 - damage) * num / denom + damage;
    }
  }

  // Determine scaled strength
  if( pressure < p1 )
  {
    real64 f;
    ceramicY10( pressure, damage, mu, Yt0, f );
    return( f / gamma1 );
  }
  else if( pressure < p2 )
  {
    real64 f, m1, f1, y1, y2;
    m1 = dfdp1 / gamma1;
    ceramicY10( p1, damage, mu, Yt0, f1 );
    y1 = f1 / gamma1;
    y2 = m_maximumStrength;
    f = pow((pressure - p2) / (p1 - p2), m1 * (p1 - p2) / (y1 - y2)) * (y1 - y2) + y2;
    return( f );
  }
  else
  {
    return( m_maximumStrength );
  }
}



/**
 * @class CeramicDamage
 *
 * Ceramic damage material model.
 */
class CeramicDamage : public ElasticIsotropic
{
public:

  /// @typedef Alias for CeramicDamageUpdates
  using KernelWrapper = CeramicDamageUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  CeramicDamage( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~CeramicDamage() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "CeramicDamage";

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
    /// string/key for quadrature point damage value
    static constexpr char const * damageString() { return "damage"; }

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for element/particle length scale
    static constexpr char const * lengthScaleString() { return "lengthScale"; }

    /// string/key for tensile strength
    static constexpr char const * tensileStrengthString() { return "tensileStrength"; }

    /// string/key for compressive strength
    static constexpr char const * compressiveStrengthString() { return "compressiveStrength"; }

    /// string/key for maximum strength
    static constexpr char const * maximumStrengthString() { return "maximumStrength"; }

    /// string/key for crack speed
    static constexpr char const * crackSpeedString() { return "crackSpeed"; }
  };

  /**
   * @brief Create a instantiation of the CeramicDamageUpdate class that refers to the data in this.
   * @return An instantiation of CeramicDamageUpdate.
   */
  CeramicDamageUpdates createKernelUpdates() const
  {
    return CeramicDamageUpdates( m_damage,
                                 m_jacobian,
                                 m_lengthScale,
                                 m_tensileStrength,
                                 m_compressiveStrength,
                                 m_maximumStrength,
                                 m_crackSpeed,
                                 m_bulkModulus,
                                 m_shearModulus,
                                 m_thermalExpansionCoefficient,
                                 m_newStress,
                                 m_oldStress,
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
                          m_damage,
                          m_jacobian,
                          m_lengthScale,
                          m_tensileStrength,
                          m_compressiveStrength,
                          m_maximumStrength,
                          m_crackSpeed,
                          m_bulkModulus,
                          m_shearModulus,
                          m_thermalExpansionCoefficient,
                          m_newStress,
                          m_oldStress,
                          m_disableInelasticity );
  }


protected:
  virtual void postInputInitialization() override;

  /// State variable: The damage values for each quadrature point
  array2d< real64 > m_damage;

  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

  /// Discretization-sized variable: The length scale for each element/particle
  array1d< real64 > m_lengthScale;

  /// Material parameter: The value of unconfined tensile strength
  real64 m_tensileStrength;

  /// Material parameter: The value of unconfined compressive strength
  real64 m_compressiveStrength;

  /// Material parameter: The value of maximum theoretical strength
  real64 m_maximumStrength;

  /// Material parameter: The value of crack speed
  real64 m_crackSpeed;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_KINEMATICDAMAGE_HPP_ */
