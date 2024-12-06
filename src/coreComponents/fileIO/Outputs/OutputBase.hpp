/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file OutputBase.hpp
 */
#ifndef GEOS_FILEIO_OUTPUTS_OUTPUTBASE_HPP_
#define GEOS_FILEIO_OUTPUTS_OUTPUTBASE_HPP_

#include "dataRepository/Group.hpp"
#include "dataRepository/ExecutableGroup.hpp"
#include "dataRepository/LogLevelsInfo.hpp"  // For logInfo namespace
#include "common/Timer.hpp"

namespace geos
{

namespace logInfo
{
/**
 * @brief Log information structure for output timing data
 * @details Provides configuration for logging output operation timing information
 */
struct OutputTimers
{
  /**
   * @brief Get the minimum log level for output timing information
   * @return The minimum log level
   */
  static constexpr int getMinLogLevel() { return 1; }

  /**
   * @brief Get the description of this timing category
   * @return Description string for output timing information
   */
  static constexpr std::string_view getDescription() { return "Output timers information"; }
};
}

/**
 * @class OutputBase
 *
 * A base class for output types.
 */
class OutputBase : public ExecutableGroup
{
public:
  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  explicit OutputBase( string const & name, Group * const parent );

  /// Destructor
  virtual ~OutputBase() override;

  /**
   * @brief Setter for the output directory
   * @param  outputDir The output directory
   **/
  static void setOutputDirectory( string const & outputDir );

  /**
   * @brief Getter for the output directory
   * @return The output directory
   **/
  static string const & getOutputDirectory();

  /**
   * @brief Setter for the file name root
   * @param root The file name root
   **/
  static void setFileNameRoot( string const & root );

  /**
   * @brief Getter for the file name root
   * @return The file name root
   **/
  static string const & getFileNameRoot();

  /// Method for setting up output directories.
  virtual void setupDirectoryStructure();

  // Catalog interface
  /// @cond DO_NOT_DOCUMENT
  using CatalogInterface = dataRepository::CatalogInterface< OutputBase, string const &, Group * const >;
  static CatalogInterface::CatalogType & getCatalog();

  // Catalog view keys
  struct viewKeysStruct
  {
    static constexpr auto childDirectoryString = "childDirectory";
    static constexpr auto parallelThreadsString = "parallelThreads";
  } outputBaseViewKeys;
  /// @endcond

  /**
   * @brief Get the path of the child directory where output will be written
   * @return The directory path
   **/
  string childDirectory() const { return m_childDirectory; }

  /**
   * @brief Get the number of parallel threads to use to write plotfiles
   * @return The number of threads
   **/
  integer parallelThreads() const { return m_parallelThreads; }



protected:
  /**
   * @brief Do initialization prior to calling initialization operations
   *        on the subgroups.
   * @param group The root group
   **/
  virtual void initializePreSubGroups() override;

  /// Timer used to track duration of file writing operations
  std::chrono::system_clock::duration m_outputTimer;

  /// @copydoc geos::ExecutableGroup::cleanup
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

private:
  string m_childDirectory;
  integer m_parallelThreads;

};


} /* namespace geos */

#endif /* GEOS_FILEIO_OUTPUTS_OUTPUTBASE_HPP_ */
