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

#ifndef GEOS_FILEIO_LOGLEVELSINFO_HPP
#define GEOS_FILEIO_LOGLEVELSINFO_HPP

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

struct DataCollectorInitialization
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription() { return "Information on Time history Initialization"; }
};

struct ChomboIOInitialization
{
  static constexpr int getMinLogLevel() { return 1; }
  static constexpr std::string_view getDescription() { return "Information on chomboIO coupling Initialization"; }
};

struct OutputEvents
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription() { return "Information on output events (VTK/ChomboIO/HDF5)"; }
};


struct HDF5Writing
{
  static constexpr int getMinLogLevel() { return 3; }
  static constexpr std::string_view getDescription() { return "Information on buffered data in an HDF5 file "; }
};

/// @endcond
///@}

}

}

#endif // GEOS_FILEIO_LOGLEVELSINFO_HPP
