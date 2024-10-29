/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalProperties.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct CompositionalProperties
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Compute the molar density and derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] volumeShift dimensional volume shift parameters
   * @param[in] compressibilityFactor compressibility factor (z-factor)
   * @param[in] compressibilityFactorDerivs derivatives of the compressibility factor (z-factor)
   * @param[out] molarDensity the calculated molar density
   * @param[out] molarDensityDerivs derivatives of the molar density
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static void computeMolarDensity( integer const numComps,
                                   real64 const pressure,
                                   real64 const temperature,
                                   arraySlice1d< real64 const, USD1 > const & composition,
                                   arraySlice1d< real64 const > const & volumeShift,
                                   real64 const compressibilityFactor,
                                   arraySlice1d< real64 const > const & compressibilityFactorDerivs,
                                   real64 & molarDensity,
                                   arraySlice1d< real64, USD2 > const & molarDensityDerivs );

  /**
   * @brief Compute the mass density and derivatives
   * @param[in] numComps number of components
   * @param[in] composition composition of the mixture
   * @param[in] molecularWeight the component molecular weights
   * @param[in] molarDensity the mixture molar density
   * @param[in] molarDensityDerivs derivatives of the molar density
   * @param[out] massDensity mass density
   * @param[out] massDensityDerivs derivatives of the mass density
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static void computeMassDensity( integer const numComps,
                                  arraySlice1d< real64 const, USD1 > const & composition,
                                  arraySlice1d< real64 const > const & molecularWeight,
                                  real64 const molarDensity,
                                  arraySlice1d< real64 const, USD2 > const & molarDensityDerivs,
                                  real64 & massDensity,
                                  arraySlice1d< real64, USD2 > const & massDensityDerivs );
};

} //namespace compositional

} // namespace constitutive

} // namespace geos

// Implementation
#include "CompositionalPropertiesImpl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALPROPERTIES_HPP_
