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

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_LOGLEVELSINFO_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_LOGLEVELSINFO_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace logInfo
{

struct DASType
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "DAS type information"; }
};

struct PMLParameters
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on Perfectly match layers parameters"; }
};

/// @endcond
///@}

}

}

#endif // GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_LOGLEVELSINFO_HPP
