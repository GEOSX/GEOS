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

#ifndef GEOS_MAININTERFACE_INITIALIZATION_HPP_
#define GEOS_MAININTERFACE_INITIALIZATION_HPP_

// Source includes
#include "common/initializeEnvironment.hpp"

namespace geos
{

/**
 * @brief Parse the command line options and populate @p commandLineOptions with the results.
 * @param argc The number of command line arguments.
 * @param argv The command line arguments.
 * @return The command line options.
 */
std::unique_ptr< CommandLineOptions > parseCommandLineOptions( int argc, char * * argv );

/**
 * @brief Perform the basic GEOSX initialization and optionally parse the command line input.
 * @param [in] argc The number of command line arguments.
 * @param [in,out] argv The command line arguments.
 * @param [in] parseCommandLine True iff the command line options should be parsed.
 * @return The command line options, if @c parseCommandLine is @c false then the returned value
 *   is default constructed (empty).
 */
std::unique_ptr< CommandLineOptions > basicSetup( int argc, char * argv[], bool const parseCommandLine=false );

/**
 * @brief Perform the basic GEOSX cleanup.
 */
void basicCleanup();



} // namespace geos

#endif // GEOS_MAININTERFACE_INITIALIZATION_HPP_
