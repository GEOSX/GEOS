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
 * @file StatOutputController.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_STATOUTPUTCONTROLLER_HPP_
#define GEOS_PHYSICSSOLVERS_STATOUTPUTCONTROLLER_HPP_
#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseStatistics.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseStatistics.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include "dataRepository/HistoryDataSpec.hpp"
#include "events/tasks/TaskBase.hpp"
#include "events/tasks/TasksManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "fileIO/Outputs/TimeHistoryOutput.hpp"
#include "fileIO/Outputs/OutputManager.hpp"
#include "fileIO/timeHistory/PackCollection.hpp"
#include "events/PeriodicEvent.hpp"

#include <functional>

namespace geos
{
/**
 * @brief StatOutputController
 * Class responsible for creating components to output regions statistics in hdf file
 * Class in charge of creating component
 */
class StatOutputController : public TaskBase
{

public:
/**
 * @brief Construct a new Stat Output Controller object
 *
 * @param name
 * @param parent
 */
  StatOutputController( string const & name,
                        Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "StatOutputController"; }

  /**
   * @copydoc dataRepository::Group::createChild()
   *
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.:
   * @code{.unparsed}
   * <Events>
   *   <StatOutputController  name="testStats" ...>
   *     <CompositionalMultiphaseStatistics name="sub_event" .../>
   *   </StatOutputController>
   * </Events>
   * @endcode
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @defgroup Tasks Interface Functions
   * This function implements the interface defined by the abstract TaskBase class
   * Execute the computation of aggregate statistics / packCollection / TimeHistoryOutput
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

private:

  string m_outputDir;

/**
 * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
 */
  struct viewKeyStruct
  {
    /// String for the flag deciding the computation of the CFL numbers
    constexpr static char const * outputDirString() { return "outputDirectory"; }
  };

  void initializePreSubGroups() override;

  /**
   * @brief Generate a packCollection component
   * @param taskManager The task manager
   * @param regionName The region
   * @param path The region path
   * @param fieldName The region name statistic
   * @param packCollectionPaths Array containing all pack collection paths
   */
  void generatePackCollection( TasksManager & taskManager,
                               string_view regionName,
                               string_view path,
                               string_view fieldName,
                               string_array & packCollectionPaths );

  /**
   * @brief Generate a TimeHistory component and automatically add the sources (in this case packCollection path) and the output file
   * @param outputManager The OutputManager
   * @param regionName The region name
   * @param packCollectionPaths Array containing all pack collection paths
   */
  void generateTimeHistory( OutputManager & outputManager,
                            string_view regionName,
                            string_array const & packCollectionPaths );

  /// The statistics component class.
  TaskBase * m_statistics;
  /// Store all generated time histories
  std::vector< TimeHistoryOutput * > m_timeHistories;
  /// Store all generated pack collection
  std::vector< PackCollection * >  m_packCollections;

  /**
   * @brief Apply the lambda expression to the supported types
   * @tparam LAMBDA The lambda type
   * @param lambda The lambda to be evaluated
   */
  template< typename LAMBDA >
  void forSubStats( LAMBDA lambda );

};

} /* namespace geos */

#endif
