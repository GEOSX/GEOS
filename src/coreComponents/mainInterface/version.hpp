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
 * @file version.hpp
 *
 * Utilities for getting and printing version info of GEOS,
 * its build environment and dependencies.
 */

#include <string>

namespace geos
{
/**
 * @brief Get GEOSX version.
 * @return The full version string.
 */
std::string getVersion();

/**
 * @brief output version info for dependencies to log
 */
void outputVersionInfo();

}
