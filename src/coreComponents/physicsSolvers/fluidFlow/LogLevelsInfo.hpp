/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LogLevelsInfo.hpp
 * This file contains common log level informations for physics solvers
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_LOGLEVELSINFO_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_LOGLEVELSINFO_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace logInfo
{

/**
 * @name Common LogLevels info structures. They must comply with the `is_log_level_info` trait.
 */
///@{

/// @cond DO_NOT_DOCUMENT

struct AggregatedSourceFluxStats
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Print aggregated statistics of all source fluxes in a mesh"; }
};

struct CFL
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "CFL information"; }
};

struct Crossflow
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Crossflow information"; }
};

struct DetailedRegionsSourceFluxStats
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription() { return "Print statistics for each source flux in each regions"; }
};


struct DetailedSourceFluxStats
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription() { return "Print statistics for each source flux in a mesh"; }
};

struct StencilConnection
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Log stencil stored connection"; }
};

struct StencilInitialization
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on stencil Initialization"; }
};

/// @endcond
///@}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_LOGLEVELSINFO_HPP
