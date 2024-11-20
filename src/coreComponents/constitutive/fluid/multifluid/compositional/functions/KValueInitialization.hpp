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
 * @file KValueInitialization.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_KVALUEINITIALIZATION_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_KVALUEINITIALIZATION_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct KValueInitialization
{
public:
  /**
   * @brief Calculate gas-liquid k-values based on the Wilson caorrelation
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[out] kValues the calculated k-values
   **/
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeWilsonGasLiquidKvalue( integer const numComps,
                                real64 const pressure,
                                real64 const temperature,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                arraySlice1d< real64, USD > const & kValues )
  {
    arrayView1d< real64 const > const & criticalPressure = componentProperties.m_componentCriticalPressure;
    arrayView1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
    arrayView1d< real64 const > const & acentricFactor = componentProperties.m_componentAcentricFactor;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const pr = criticalPressure[ic] / pressure;
      real64 const tr = criticalTemperature[ic] / temperature;
      kValues[ic] = pr * exp( 5.37 * ( 1.0 + acentricFactor[ic] ) * ( 1.0 - tr ) );
    }
  }

/**
 * @brief Calculate water-gas k-value
 * @param[in] pressure pressure
 * @param[in] temperature temperature
 * @return The water component k-value
 **/
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static double
  computeWaterGasKvalue( double pressure,
                         double temperature )
  {
    return exp( -4844.168051 / temperature + 12.93022442 ) * 1.0e5 / pressure;
  }

};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_KVALUEINITIALIZATION_HPP_
