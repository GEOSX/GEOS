/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LogLevelsInfo.hpp
 * This file contains common log level informations for physics solvers
 */

#ifndef GEOS_PHYSICSSOLVERS_LOGLEVELSINFO_HPP
#define GEOS_PHYSICSSOLVERS_LOGLEVELSINFO_HPP

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

struct Fields
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription() { return "The summary of declared fields and coupling"; }
};

struct LineSearch
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Line search information"; }
};

struct Solution
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Solution information (scaling, maximum changes, quality check)"; }
};

struct Convergence
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Convergence information"; }
};

struct TimeStep
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Time step information"; }
};

struct LinearSolver
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Linear solver information"; }
};

struct NonlinearSolver
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Nonlinear solver information"; }
};

struct Timers
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Solver timers information"; }
};

struct Initialization
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Initialization information"; }
};

/// @endcond
///@}

}

}

#endif // GEOS_PHYSICSSOLVERS_LOGLEVELSINFO_HPP
