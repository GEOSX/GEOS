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
 * @file LogLevelsInfo.hpp
 * This file contains log level informations for contact solvers
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_LOGLEVELSINFO_HPP
#define GEOS_PHYSICSSOLVERS_CONTACT_LOGLEVELSINFO_HPP

#include "common/DataTypes.hpp"

namespace geos
{

namespace logInfo
{

/**
 * @name Wells LogLevels info structures. They must comply with the `is_log_level_info` trait.
 */
///@{

/// @cond DO_NOT_DOCUMENT

struct Configuration
{
  static constexpr int getMinLogLevel() { return 2; }
  static constexpr std::string_view getDescription() { return "Configuration information"; }
};

/// @endcond
///@}

}

}

#endif // GEOS_PHYSICSSOLVERS_CONTACT_LOGLEVELSINFO_HPP
