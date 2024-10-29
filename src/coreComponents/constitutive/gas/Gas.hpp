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
 *  @file Gas.hpp
 */

#ifndef GEOS_CONSTITUTIVE_GAS_GAS_HPP_
#define GEOS_CONSTITUTIVE_GAS_GAS_HPP_

#include "constitutive/ContinuumBase.hpp"
// #include "SolidModelDiscretizationOpsIsotropic.hpp"
// #include "constitutive/ExponentialRelation.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class GasUpdates
 *
 * Class to provide elastic isotropic material updates that may be
 * called from a kernel function.
 */
class GasUpdates : public ContinuumBaseUpdates
{
public:
  /**
   * @brief Constructor
   * @param[in] referencePressure    The value of the reference pressure 
   * @param[in] referenceTemperature The value of the reference temperature
   * @param[in] jacobian             The ArrayView holding the jacobian data for each element.
   * @param[in] temperature          The ArrayView holding the temperature data for each element.
   * @param[in] newStress            The ArrayView holding the new stress data for each quadrature point.
   * @param[in] oldStress            The ArrayView holding the old stress data for each quadrature point.
   */
  GasUpdates( real64 const referencePressure,
              real64 const referenceTemperature,
              arrayView1d< real64 const > const & jacobian,
              arrayView1d< real64 const > const & temperature,
              arrayView3d< real64, solid::STRESS_USD > const & newStress,
              arrayView3d< real64, solid::STRESS_USD > const & oldStress,
              arrayView2d< real64 > const & density,
              arrayView2d< real64 > const & wavespeed ):
    ContinuumBaseUpdates( newStress, 
                          oldStress,
                          density,
                          wavespeed ),
    m_referencePressure( referencePressure ),
    m_referenceTemperature( referenceTemperature ),
    m_jacobian( jacobian ),
    m_temperature( temperature ) 
  {}

  /// Deleted default constructor
  GasUpdates() = delete;

  /// Default copy constructor
  GasUpdates( GasUpdates const & ) = default;

  /// Default move constructor
  GasUpdates( GasUpdates && ) = default;

  /// Deleted copy assignment operator
  GasUpdates & operator=( GasUpdates const & ) = delete;

  /// Deleted move assignment operator
  GasUpdates & operator=( GasUpdates && ) =  delete;

  // /// Use the "isotropic" form of inner product compression
  // using DiscretizationOps = SolidModelDiscretizationOpsIsotropic;

  /// Use base version of saveConvergedState
  using ContinuumBaseUpdates::saveConvergedState;

  GEOS_HOST_DEVICE
  virtual void smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                    localIndex const q,
                                                    real64 const ( &totalStrain )[6],
                                                    real64 ( &stress )[6] ) const override;

  // GEOS_HOST_DEVICE
  // virtual void smallStrainNoStateUpdate( localIndex const k,
  //                                        localIndex const q,
  //                                        real64 const ( &totalStrain )[6],
  //                                        real64 ( &stress )[6],
  //                                        real64 ( &stiffness )[6][6] ) const override;

  // GEOS_HOST_DEVICE
  // virtual void smallStrainNoStateUpdate( localIndex const k,
  //                                        localIndex const q,
  //                                        real64 const ( &totalStrain )[6],
  //                                        real64 ( &stress )[6],
  //                                        DiscretizationOps & stiffness ) const;

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

  // GEOS_HOST_DEVICE
  // virtual void smallStrainUpdate( localIndex const k,
  //                                 localIndex const q,
  //                                 real64 const & timeIncrement,
  //                                 real64 const ( &strainIncrement )[6],
  //                                 real64 ( &stress )[6],
  //                                 DiscretizationOps & stiffness ) const;

  // TODO: confirm hyper stress/strain measures before activatiing
  // CC: need hyperelastic model for hyperelasticMMS constitutive model in MPM
  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6] ) const override final;

  GEOS_HOST_DEVICE
  virtual void hyperUpdate( localIndex const k,
                            localIndex const q,
                            real64 const ( & FminusI )[3][3],
                            real64 ( & stress )[6],
                            real64 ( & stiffness )[6][6] ) const override final;

protected:

  // The value of the reference pressure
  real64 const m_referencePressure;

  // The value of the reference temperature
  real64 const m_referenceTemperature;

  /// A reference to the ArrayView holding the jacobian for each element.
  arrayView1d< real64 const > const m_jacobian;

  /// A reference to the ArrayView holding the temperature for each element.
  arrayView1d< real64 const > const m_temperature;

};


GEOS_HOST_DEVICE
inline
void GasUpdates::smallStrainNoStateUpdate_StressOnly( localIndex const k,
                                                      localIndex const q,
                                                      real64 const ( &totalStrain )[6],
                                                      real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( totalStrain );

  real64 pressure =  m_referencePressure / m_jacobian[k] * ( m_temperature[k] / m_referenceTemperature );

  stress[0] = -pressure;
  stress[1] = -pressure;
  stress[2] = -pressure;

  stress[3] = 0.0;
  stress[4] = 0.0;
  stress[5] = 0.0;

  m_wavespeed[k][0] = sqrt( 1.4*pressure / m_density[k][0]);
}


// GEOS_HOST_DEVICE
// inline
// void GasUpdates::smallStrainNoStateUpdate( localIndex const k,
//                                                         localIndex const q,
//                                                         real64 const ( &totalStrain )[6],
//                                                         real64 ( & stress )[6],
//                                                         real64 ( & stiffness )[6][6] ) const
// {
//   smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
//   getElasticStiffness( k, q, stiffness );
// }


// GEOS_HOST_DEVICE
// inline
// void GasUpdates::smallStrainNoStateUpdate( localIndex const k,
//                                            localIndex const q,
//                                            real64 const ( &totalStrain )[6],
//                                            real64 ( & stress )[6],
//                                            DiscretizationOps & stiffness ) const
// {
//   smallStrainNoStateUpdate_StressOnly( k, q, totalStrain, stress );
//   stiffness.m_bulkModulus = m_bulkModulus[k];
//   stiffness.m_shearModulus = m_shearModulus[k];
// }


GEOS_HOST_DEVICE
inline
void GasUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                               localIndex const q,
                                               real64 const & timeIncrement,
                                               real64 const ( &strainIncrement )[6],
                                               real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( timeIncrement );
  smallStrainNoStateUpdate_StressOnly( k, q, strainIncrement, stress ); // stress  = incrementalStress
  saveStress( k, q, stress );                                           // m_newStress = stress
}


GEOS_HOST_DEVICE
inline
void GasUpdates::smallStrainUpdate_StressOnly( localIndex const k,
                                               localIndex const q,
                                               real64 const & timeIncrement,
                                               real64 const ( & beginningRotation )[3][3],
                                               real64 const ( & endRotation )[3][3],
                                               real64 const ( & strainIncrement )[6],
                                               real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( beginningRotation );
  GEOS_UNUSED_VAR( endRotation );

  smallStrainUpdate_StressOnly( k,
                                q,
                                timeIncrement,
                                strainIncrement,
                                stress );
}


// GEOS_HOST_DEVICE
// inline
// void GasUpdates::smallStrainUpdate( localIndex const k,
//                                                  localIndex const q,
//                                                  real64 const & timeIncrement,
//                                                  real64 const ( &strainIncrement )[6],
//                                                  real64 ( & stress )[6],
//                                                  real64 ( & stiffness )[6][6] ) const
// {
//   smallStrainUpdate_StressOnly( k, 
//                                 q, 
//                                 timeIncrement,
//                                 strainIncrement, 
//                                 stress );
//   getElasticStiffness( k, q, stiffness );
// }


// GEOS_HOST_DEVICE
// inline
// void GasUpdates::smallStrainUpdate( localIndex const k,
//                                                  localIndex const q,
//                                                  real64 const & timeIncrement,
//                                                  real64 const ( &strainIncrement )[6],
//                                                  real64 ( & stress )[6],
//                                                  DiscretizationOps & stiffness ) const
// {
//   smallStrainUpdate_StressOnly( k,
//                                 q, 
//                                 timeIncrement,
//                                 strainIncrement, 
//                                 stress );
//   stiffness.m_bulkModulus = m_bulkModulus[k];
//   stiffness.m_shearModulus = m_shearModulus[k];
// }

// GEOS_HOST_DEVICE
// GEOS_FORCE_INLINE
// void GasUpdates::viscousStateUpdate( localIndex const k,
//                                                   localIndex const q,
//                                                   real64 beta ) const
// {
//   GEOS_UNUSED_VAR( k );
//   GEOS_UNUSED_VAR( q );
//   GEOS_UNUSED_VAR( beta );
// }

// CC: place filler for MPM hyperelasticMMS model
// There is a commented out implementation below from another person
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GasUpdates::hyperUpdate( localIndex const k,
                          localIndex const q,
                          real64 const ( & FminusI )[3][3],
                          real64 ( & stress )[6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( FminusI );
  GEOS_UNUSED_VAR( stress );
  GEOS_ERROR( "hyperUpdate() not implemented for this model" );
}

// CC: place filler for MPM hyperelasticMMS model
// There is a commented out implementation below from another person
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void GasUpdates::hyperUpdate( localIndex const k,
                          localIndex const q,
                          real64 const ( & FminusI )[3][3],
                          real64 ( & stress )[6],
                          real64 ( & stiffness )[6][6] ) const
{
  GEOS_UNUSED_VAR( k );
  GEOS_UNUSED_VAR( q );
  GEOS_UNUSED_VAR( FminusI );
  GEOS_UNUSED_VAR( stress );
  GEOS_UNUSED_VAR( stiffness );
  GEOS_ERROR( "hyperUpdate() not implemented for this model" );
}

/**
 * @class Gas
 *
 * Class to provide an elastic isotropic material response.
 */
class Gas : public ContinuumBase
{
public:

  /// Alias for GasUpdates
  using KernelWrapper = GasUpdates;

  /**
   * constructor
   * @param[in] name name of the instance in the catalog
   * @param[in] parent the group which contains this instance
   */
  Gas( string const & name, 
                    Group * const parent );

  /**
   * Default Destructor
   */
  virtual ~Gas() override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "Gas";

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
  struct viewKeyStruct : public ContinuumBase::viewKeyStruct
  {
    /// string/key for reference pressure
    static constexpr char const * referencePressureString() { return "referencePressure"; }

    /// string/key for default reference volume
    static constexpr char const * referenceVolumeString() { return "referenceVolume"; }

    /// string/key for default reference temperature
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }

    /// string/key for volume
    static constexpr char const * jacobianString() { return "jacobian"; }

    /// string/key for temperature
    static constexpr char const * temperatureString() { return "temperature"; }
  };

  /**
   * @brief Accessor for jacobian
   * @return A const reference to arrayView1d<real64> containing the jacobian (at every element).
   */
  arrayView1d< real64 > const jacobian() { return m_jacobian; }

  /**
   * @brief Const accessor for jacobian
   * @return A const reference to arrayView1d<real64 const> containing jacobian (at every element).
   */
  arrayView1d< real64 const > const jacobian() const { return m_jacobian; }

  /**
   * @brief Accessor for temperature
   * @return A const reference to arrayView1d<real64> containing the temperature (at every element).
   */
  arrayView1d< real64 > const temperature() { return m_temperature; }

  /**
   * @brief Const accessor for temperature
   * @return A const reference to arrayView1d<real64 const> containing the temperature (at every element).
   */
  arrayView1d< real64 const > const temperature() const { return m_temperature; }

  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getJacobian() const final
  {
    return m_jacobian;
  }
  
  GEOS_HOST_DEVICE
  virtual arrayView1d< real64 const > getTemperature() const final
  {
    return m_temperature;
  }

  /**
   * @brief Create a instantiation of the GasUpdate class
   *        that refers to the data in this.
   * @param includeState Flag whether to pass state arrays that may not be needed for "no-state" updates
   * @return An instantiation of GasUpdate.
   */
  GasUpdates createKernelUpdates( bool const includeState = true ) const
  {
    if( includeState )
    {
      return GasUpdates( m_referencePressure,
                         m_referenceTemperature,
                         m_jacobian,
                         m_temperature,
                         m_newStress,
                         m_oldStress,
                         m_density,
                         m_wavespeed );
    }
    else // for "no state" updates, pass empty views to avoid transfer of stress data to device
    {
      return GasUpdates( m_referencePressure,
                         m_referenceTemperature,
                         m_jacobian,
                         m_temperature,
                         arrayView3d< real64, solid::STRESS_USD >(),
                         arrayView3d< real64, solid::STRESS_USD >(),
                         m_density,
                         m_wavespeed );
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
                          m_referencePressure,
                          m_referenceTemperature,
                          m_jacobian,
                          m_temperature,
                          m_newStress,
                          m_oldStress,
                          m_wavespeed );
  }

  /**
   * @brief Allocate constitutive arrays
   * @param parent Object's parent group (element subregion)
   * @param numConstitutivePointsPerParentIndex Number of quadrature points per element
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

protected:

  /// Post-process XML data
  virtual void postInputInitialization() override;

  /// The reference value of the pressure.
  real64 m_referencePressure;

  /// The reference value of the temperature.
  real64 m_referenceTemperature;

  /// The jacobian for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_jacobian;

  /// The temperature for each upper level dimension (i.e. cell) of *this
  array1d< real64 > m_temperature;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_GAS_GAS_HPP_ */
