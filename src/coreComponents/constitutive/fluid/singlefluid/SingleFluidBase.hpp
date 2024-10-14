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
 * @file SingleFluidBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SINGLEFLUIDBASE_HPP
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SINGLEFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{

namespace constitutive
{
//START_SPHINX_INCLUDE_01
/**
 * @brief Base class for single-phase fluid model kernel wrappers.
 */
class SingleFluidBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOS_HOST_DEVICE
  localIndex numElems() const { return m_density.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOS_HOST_DEVICE
  localIndex numGauss() const { return m_density.size( 1 ); };

protected:

  /**
   * @brief Constructor.
   * @param density     fluid density
   * @param dDens_dPres derivative of density w.r.t. pressure
   * @param viscosity   fluid viscosity
   * @param dVisc_dPres derivative of viscosity w.r.t. pressure
   */
  SingleFluidBaseUpdate( arrayView2d< real64 > const & density,
                         arrayView2d< real64 > const & dDens_dPres,
                         arrayView2d< real64 > const & viscosity,
                         arrayView2d< real64 > const & dVisc_dPres )
    : m_density( density ),
    m_dDens_dPres( dDens_dPres ),
    m_viscosity( viscosity ),
    m_dVisc_dPres( dVisc_dPres )
  {}

  /**
   * @brief Copy constructor.
   */
  SingleFluidBaseUpdate( SingleFluidBaseUpdate const & ) = default;

  /**
   * @brief Move constructor.
   */
  SingleFluidBaseUpdate( SingleFluidBaseUpdate && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return reference to this object
   */
  SingleFluidBaseUpdate & operator=( SingleFluidBaseUpdate const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return reference to this object
   */
  SingleFluidBaseUpdate & operator=( SingleFluidBaseUpdate && ) = delete;


  /// Fluid density
  arrayView2d< real64 > m_density;

  /// Derivative of density w.r.t. pressure
  arrayView2d< real64 > m_dDens_dPres;

  /// Fluid viscosity
  arrayView2d< real64 > m_viscosity;

  /// Derivative of viscosity w.r.t. pressure
  arrayView2d< real64 > m_dVisc_dPres;

//END_SPHINX_INCLUDE_01

public:
  /**
   * @brief Compute function to update properties in a cell without returning derivatives.
   * @details This delegates the call to the fluid wrapper using the value and derivative function.
   *          This is used for initialisation and boundary conditions.
   * @param[in] fluidWrapper the actual fluid kernel
   * @param[in]  pressure the target pressure value
   * @param[out] density fluid density
   * @param[out] viscosity fluid viscosity
   */
  template< typename FLUIDWRAPPER >
  GEOS_HOST_DEVICE
  static void computeValues( FLUIDWRAPPER const fluidWrapper,
                             real64 const pressure,
                             real64 & density,
                             real64 & viscosity )
  {
    real64 dDensity_dPressure = 0.0;
    real64 dViscosity_dPressure = 0.0;
    fluidWrapper.compute( pressure,
                          density,
                          dDensity_dPressure,
                          viscosity,
                          dViscosity_dPressure );
  }

//START_SPHINX_INCLUDE_02
private:
  /**
   * @brief Compute fluid properties and derivatives at a single point.
   * @param[in]  pressure the target pressure value
   * @param[out] density fluid density
   * @param[out] dDensity_dPressure fluid density derivative w.r.t. pressure
   * @param[out] viscosity fluid viscosity
   * @param[out] dViscosity_dPressure fluid viscosity derivative w.r.t. pressure
   */
  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const = 0;

  /**
   * @brief Compute fluid properties and derivatives at a single point.
   * @param[in]  pressure the target pressure value
   * @param[in]  temperature the target temperature value
   * @param[out] density fluid density
   * @param[out] dDensity_dPressure fluid density derivative w.r.t. pressure
   * @param[out] dDensity_dTemperature fluid density derivative w.r.t. temperature
   * @param[out] viscosity fluid viscosity
   * @param[out] dViscosity_dPressure fluid viscosity derivative w.r.t. pressure
   * @param[out] dViscosity_dTemperature fluid viscosity derivative w.r.t. temperature
   * @param[out] internalEnergy fluid internal energy
   * @param[out] dInternalEnergy_dPressure fluid internal energy derivative w.r.t. pressure
   * @param[out] dInternalEnergy_dTemperature fluid internal energy derivative w.r.t. temperature
   * @param[out] enthalpy fluid enthalpy
   * @param[out] dEnthalpy_dPressure fluid enthalpy derivative w.r.t. pressure
   * @param[out] dEnthalpy_dTemperature fluid enthalpy derivative w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & dDensity_dTemperature,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure,
                        real64 & dViscosity_dTemperature,
                        real64 & internalEnergy,
                        real64 & dInternalEnergy_dPressure,
                        real64 & dInternalEnergy_dTemperature,
                        real64 & enthalpy,
                        real64 & dEnthalpy_dPressure,
                        real64 & dEnthalpy_dTemperature ) const = 0;

  /**
   * @brief Update fluid state at a single point.
   * @param[in] k        element index
   * @param[in] q        gauss point index
   * @param[in] pressure the target pressure value
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure ) const = 0;

  /**
   * @brief Update fluid state at a single point.
   * @param[in] k           element index
   * @param[in] q           gauss point index
   * @param[in] pressure    the target pressure value
   * @param[in] temperature the target temperature value
   */
  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature ) const = 0;

};
//END_SPHINX_INCLUDE_02

/**
 * @brief Base class for single-phase fluid models.
 */
class SingleFluidBase : public ConstitutiveBase
{
public:

  /**
   * @brief Constructor.
   * @param name name of the group
   * @param parent pointer to parent group
   */
  SingleFluidBase( string const & name, Group * const parent );

  /**
   * @brief Initialize the model
   */
  void initializeState() const;

  virtual void saveConvergedState() const override;

  // *** ConstitutiveBase interface

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** SingleFluid-specific interface

  arrayView2d< real64 const > density() const { return m_density; }

  arrayView2d< real64 const > dDensity_dPressure() const { return m_dDensity_dPressure; }

  arrayView2d< real64 > dDensity_dTemperature() { return m_dDensity_dTemperature; }
  arrayView2d< real64 const > dDensity_dTemperature() const { return m_dDensity_dTemperature; }

  arrayView2d< real64 const > density_n() const { return m_density_n; }

  arrayView2d< real64 const > viscosity() const { return m_viscosity; }

  arrayView2d< real64 const > dViscosity_dPressure() const { return m_dViscosity_dPressure; }

  arrayView2d< real64 > dViscosity_dTemperature() { return m_dViscosity_dTemperature; }
  arrayView2d< real64 const > dViscosity_dTemperature() const { return m_dViscosity_dTemperature; }

  arrayView2d< real64 > internalEnergy() { return m_internalEnergy; }
  arrayView2d< real64 const > internalEnergy() const { return m_internalEnergy; }

  arrayView2d< real64 > internalEnergy_n() { return m_internalEnergy_n; }
  arrayView2d< real64 const > internalEnergy_n() const { return m_internalEnergy_n; }

  arrayView2d< real64 > dInternalEnergy_dPressure() { return m_dInternalEnergy_dPressure; }
  arrayView2d< real64 const > dInternalEnergy_dPressure() const { return m_dInternalEnergy_dPressure; }

  arrayView2d< real64 > dInternalEnergy_dTemperature() { return m_dInternalEnergy_dTemperature; }
  arrayView2d< real64 const > dInternalEnergy_dTemperature() const { return m_dInternalEnergy_dTemperature; }

  arrayView2d< real64 > enthalpy() { return m_enthalpy; }
  arrayView2d< real64 const > enthalpy() const { return m_enthalpy; }

  arrayView2d< real64 > dEnthalpy_dPressure() { return m_dEnthalpy_dPressure; }
  arrayView2d< real64 const > dEnthalpy_dPressure() const { return m_dEnthalpy_dPressure; }

  arrayView2d< real64 > dEnthalpy_dTemperature() { return m_dEnthalpy_dTemperature; }
  arrayView2d< real64 const > dEnthalpy_dTemperature() const { return m_dEnthalpy_dTemperature; }

  virtual real64 defaultDensity() const = 0;
  virtual real64 defaultViscosity() const = 0;

protected:

  virtual void postInputInitialization() override;

  //START_SPHINX_INCLUDE_00
  array2d< real64 > m_density;
  array2d< real64 > m_dDensity_dPressure;
  array2d< real64 > m_dDensity_dTemperature;

  array2d< real64 > m_density_n;

  array2d< real64 > m_viscosity;
  array2d< real64 > m_dViscosity_dPressure;
  array2d< real64 > m_dViscosity_dTemperature;

  array2d< real64 > m_internalEnergy;
  array2d< real64 > m_internalEnergy_n;
  array2d< real64 > m_dInternalEnergy_dPressure;
  array2d< real64 > m_dInternalEnergy_dTemperature;

  array2d< real64 > m_enthalpy;
  array2d< real64 > m_dEnthalpy_dPressure;
  array2d< real64 > m_dEnthalpy_dTemperature;
  //END_SPHINX_INCLUDE_00
};

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SINGLEFLUIDBASE_HPP
