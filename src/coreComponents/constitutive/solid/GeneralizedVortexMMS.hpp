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
 * @file GeneralizedVortexMMS.hpp
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

#ifndef GEOSX_GENERALIZEDVORTEXMMS_HPP
#define GEOSX_GENERALIZEDVORTEXMMS_HPP

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
 * @class GeneralizedVortexMMSUpdates
 *
 * Class to provide material updates that may be
 * called from a kernel function.
 */
class GeneralizedVortexMMSUpdates : public ElasticIsotropicUpdates
{
public:

  /**
   * @brief Constructor
   * @param[in] jacobian The ArrayView holding the jacobian for each quardrature point.
   * @param[in] bulkModulus The ArrayView holding the bulk modulus data for each element.
   * @param[in] shearModulus The ArrayView holding the shear modulus data for each element.
   * @param[in] thermalExpansionCoefficient The ArrayView holding the thermal expansion coefficient data for each element.
   * @param[in] newStress The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress The ArrayView holding the old stress data for each quadrature point.
   */
  GeneralizedVortexMMSUpdates(arrayView2d< real64 > const & jacobian,
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
                             wavepseed,
                             disableInelasticity ),
    m_jacobian( jacobian )
  {}

  /// Default copy constructor
  GeneralizedVortexMMSUpdates( GeneralizedVortexMMSUpdates const & ) = default;

  /// Default move constructor
  GeneralizedVortexMMSUpdates( GeneralizedVortexMMSUpdates && ) = default;

  /// Deleted default constructor
  GeneralizedVortexMMSUpdates() = delete;

  /// Deleted copy assignment operator
  GeneralizedVortexMMSUpdates & operator=( GeneralizedVortexMMSUpdates const & ) = delete;

  /// Deleted move assignment operator
  GeneralizedVortexMMSUpdates & operator=( GeneralizedVortexMMSUpdates && ) =  delete;

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
  GEOS_FORCE_INLINE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
  }

private:
  /// A reference to the ArrayView holding the jacobian for each quadrature point.
  arrayView2d< real64 > const m_jacobian;

};


GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeneralizedVortexMMSUpdates::smallStrainUpdate( localIndex const k,
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
  GeneralizedVortexMMSUpdates::smallStrainUpdateHelper( k, q, timeIncrement, stress );

  // It doesn't make sense to modify stiffness with this model

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeneralizedVortexMMSUpdates::smallStrainUpdate( localIndex const k,
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
void GeneralizedVortexMMSUpdates::smallStrainUpdate_StressOnly( localIndex const k,
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
  GeneralizedVortexMMSUpdates::smallStrainUpdateHelper( k, q, timeIncrement, stress );

  // save new stress and return
  saveStress( k, q, stress );
  return;
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GeneralizedVortexMMSUpdates::smallStrainUpdateHelper( localIndex const k,
                                                    localIndex const q,
                                                    real64 const dt,
                                                    real64 ( & stress )[6] ) const
{
  // get trial pressure
  real64 trialPressure = -m_bulkModulus[k] * log( m_jacobian[k][q] );
  real64 pressure = trialPressure;


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
    real64 strength = GeneralizedVortexMMSUpdates::getStrength( m_damage[k][q], pressure, J2, J3, mu, Yt0 );

    // Increment damage and get new associated yield surface
    real64 newDeviatorMagnitude = vonMises;
    if( vonMises > strength )
    {
      if( pressure <= brittleDuctileTransitionPressure )
      {
        m_damage[k][q] = fmin( m_damage[k][q] + dt / tFail, 1.0 );
        strength = GeneralizedVortexMMSUpdates::getStrength( m_damage[k][q], pressure, J2, J3, mu, Yt0 );
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

/**
 * @class GeneralizedVortexMMS
 *
 * Ceramic damage material model.
 */
class GeneralizedVortexMMS : public ElasticIsotropic
{
public:

  /// @typedef Alias for GeneralizedVortexMMSUpdates
  using KernelWrapper = GeneralizedVortexMMSUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  GeneralizedVortexMMS( string const & name, Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~GeneralizedVortexMMS() override;


  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "GeneralizedVortexMMS";

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

    /// string/key for quadrature point jacobian value
    static constexpr char const * jacobianString() { return "jacobian"; }

  };

  /**
   * @brief Create a instantiation of the GeneralizedVortexMMSUpdate class that refers to the data in this.
   * @return An instantiation of GeneralizedVortexMMSUpdate.
   */
  GeneralizedVortexMMSUpdates createKernelUpdates() const
  {
    return GeneralizedVortexMMSUpdates(m_jacobian,
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
                          m_jacobian,
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

  /// State variable: The jacobian of the deformation
  array2d< real64 > m_jacobian;

};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOSX_CONSTITUTIVE_SOLID_GENERALIZEDVORTEXMMS_HPP_ */
