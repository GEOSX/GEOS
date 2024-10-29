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
 * @file CO2SolubilityDuanSun.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITYDUANSUN_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITYDUANSUN_HPP_

#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"

namespace geos
{
namespace constitutive
{
namespace PVTProps
{

struct CO2SolubilityDuanSun
{
/**
 * @brief Create CO2 and H2O solubility table based on Duan and Sun (2003)
 * @details Each generated table is a 2D table with lookup properties pressure (in Pa) and
 *          temperature (in degC). The returned CO2 solubility is in mole of CO2 per kg of
 *          H2O and the returned water vapourisation is in moles of H2O per kg of CO2.
 * @param[in] functionName The name of the model
 * @param[in] tableCoords The values of pressure and temperature
 * @param[in] salinity The salinity of the brine
 * @param[in] tolerance Tolerance to be used in solving for the solubility
 * @param[out] co2SolubilityValues The CO2 solubility values (mol/kg) at the given pressures and temperatures
 * @param[out] h2oSolubilityValues The H2O solubility values (mol/kg) at the given pressures and temperatures
 */
  static void populateSolubilityTables( string const & functionName,
                                        PTTableCoordinates const & tableCoords,
                                        real64 const & salinity,
                                        real64 const & tolerance,
                                        array1d< real64 > const & co2SolubilityValues,
                                        array1d< real64 > const & h2oSolubilityValues );
};

} // end namespace PVTProps
} // end namespace constitutive
} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2SOLUBILITYDUANSUN_HPP_
