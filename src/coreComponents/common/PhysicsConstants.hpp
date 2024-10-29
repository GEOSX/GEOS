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
 * @file PhysicsConstants.hpp
 * @brief Regroups useful constants that are globally used for math and physics computations.
 */
#ifndef GEOS_COMMON_PHYSICSCONSTANTS_HPP_
#define GEOS_COMMON_PHYSICSCONSTANTS_HPP_

namespace geos
{

namespace constants
{

/**
 * @brief Zero degree Celsius in Kelvin
 */
constexpr double zeroDegreesCelsiusInKelvin = 273.15;

/**
 * @brief Shorthand for pi
 */
constexpr double pi = 3.141592653589793238;

/**
 * @brief Universal gas constant
 */
constexpr double gasConstant = 8.31446261815324;

} // end namespace constants

} // end namespace geos

#endif //GEOS_COMMON_PHYSICSCONSTANTS_HPP_
