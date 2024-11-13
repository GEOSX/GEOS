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
 * @file FugacityCalculator.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FUGACITYCALCULATOR_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FUGACITYCALCULATOR_HPP_

#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct FugacityCalculator
{
  /**
   * @brief Calculate the log fugacity for a phase
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] equationOfState The equation of state
   * @param[out] logFugacity the calculated log fugacity
   */
  template< int USD >
  GEOS_HOST_DEVICE
  static void computeLogFugacity( integer const numComps,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  EquationOfStateType const equationOfState,
                                  arraySlice1d< real64 > const & logFugacity );

  /**
   * @brief Calculate the derivatives for the log fugacity for a phase
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] equationOfState The equation of state
   * @param[in] logFugacity the calculated log fugacity
   * @param[out] logFugacityDerivs the calculated derivatives of the log fugacity
   */
  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  static void computeLogFugacityDerivatives( integer const numComps,
                                             real64 const pressure,
                                             real64 const temperature,
                                             arraySlice1d< real64 const, USD1 > const & composition,
                                             ComponentProperties::KernelWrapper const & componentProperties,
                                             EquationOfStateType const equationOfState,
                                             arraySlice1d< real64 const > const & logFugacity,
                                             arraySlice2d< real64, USD2 > const & logFugacityDerivs );
};

template< int USD >
GEOS_HOST_DEVICE
void FugacityCalculator::computeLogFugacity( integer const numComps,
                                             real64 const pressure,
                                             real64 const temperature,
                                             arraySlice1d< real64 const, USD > const & composition,
                                             ComponentProperties::KernelWrapper const & componentProperties,
                                             EquationOfStateType const equationOfState,
                                             arraySlice1d< real64 > const & logFugacity )
{
  if( equationOfState == EquationOfStateType::PengRobinson )
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::
    computeLogFugacityCoefficients( numComps,
                                    pressure,
                                    temperature,
                                    composition,
                                    componentProperties,
                                    logFugacity );
  }
  else if( equationOfState == EquationOfStateType::SoaveRedlichKwong )
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
    computeLogFugacityCoefficients( numComps,
                                    pressure,
                                    temperature,
                                    composition,
                                    componentProperties,
                                    logFugacity );
  }
}

template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void FugacityCalculator::computeLogFugacityDerivatives( integer const numComps,
                                                        real64 const pressure,
                                                        real64 const temperature,
                                                        arraySlice1d< real64 const, USD1 > const & composition,
                                                        ComponentProperties::KernelWrapper const & componentProperties,
                                                        EquationOfStateType const equationOfState,
                                                        arraySlice1d< real64 const > const & logFugacity,
                                                        arraySlice2d< real64, USD2 > const & logFugacityDerivs )
{
  if( equationOfState == EquationOfStateType::PengRobinson )
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::
    computeLogFugacityCoefficients( numComps,
                                    pressure,
                                    temperature,
                                    composition,
                                    componentProperties,
                                    logFugacity,
                                    logFugacityDerivs );
  }
  else if( equationOfState == EquationOfStateType::SoaveRedlichKwong )
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
    computeLogFugacityCoefficients( numComps,
                                    pressure,
                                    temperature,
                                    composition,
                                    componentProperties,
                                    logFugacity,
                                    logFugacityDerivs );
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos


#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FUGACITYCALCULATOR_HPP_
