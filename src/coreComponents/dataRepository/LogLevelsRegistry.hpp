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
 * @file LogLevelsRegistry.hpp
 */

#ifndef GEOS_COMMON_LOGLEVELSREGISTRY_HPP
#define GEOS_COMMON_LOGLEVELSREGISTRY_HPP

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"

namespace geos
{

/**
 * @brief Keep track of log level documention for a group
 */
class LogLevelsRegistry
{
public:

  /**
   * @brief Add a log description for a wrapper
   * @param level The minimum log level
   * @param description The description for the log level
   */
  void addEntry( integer level, std::string_view description );

  /**
   * @brief Construct the log level string description for a wrapper
   * @return The log level string description
   */
  string buildLogLevelDescription() const;

private:

  /**
   * @brief Map for building the log level string for each wrapper.
   *        key : a logLevel condition, values : a set of description for a corresponding loglevel.
   */
  std::map< integer, std::vector< std::string > > m_logLevelsDescriptions;

};

}

#endif
