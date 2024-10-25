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
 * @file TimeHistoryOutput.hpp
 */

#ifndef GEOS_FILEIO_OUTPUTS_HISTORYOUTPUT_HPP_
#define GEOS_FILEIO_OUTPUTS_HISTORYOUTPUT_HPP_

#include "OutputBase.hpp"
#include "fileIO/timeHistory/HistoryCollection.hpp"
#include "fileIO/timeHistory/BufferedHistoryIO.hpp"
#include "fileIO/timeHistory/HDFHistoryIO.hpp"

#include "LvArray/src/Array.hpp" // just for collector

namespace geos
{

/**
 * @class TimeHistoryOutput
 *
 * A class for creating time history output
 */
class TimeHistoryOutput : public OutputBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group(string const & name, Group * const parent)
  TimeHistoryOutput( string const & name,
                     Group * const parent );

  /// Destructor
  virtual ~TimeHistoryOutput() override
  { }

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "TimeHistory"; }

  /**
   * @brief Perform initalization after all subgroups have been initialized.
   *   Check for existing files and data spaces/sets on restart, else create the
   *   file and data spaces/sets and output any set index metatata for the time history.
   * @note There are operations in this function that are collective on the GEOSX comm.
   */
  virtual void initializePostInitialConditionsPostSubGroups() override;

  /**
   * @brief Performs re-initialization of certain variable depending on the solver being used.
   */
  virtual void reinit() override;

  /**
   * @brief Set the output filename (This is usefull for pygeosx user)
   * @param root The string name of the output file
   */
  void setFileName( string const & root );

  /**
   * @brief Writes out a time history file.
   * @copydoc EventBase::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;
  /**
   * @brief Writes out a time history file at the end of the simulation.
   * @copydoc ExecutableGroup::cleanup()
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeys
  {
    static constexpr char const * timeHistoryOutputTargetString() { return "sources"; }
    static constexpr char const * timeHistoryOutputFilenameString() { return "filename"; }
    static constexpr char const * timeHistoryOutputFormatString() { return "format"; }
    static constexpr char const * timeHistoryRestartString() { return "restart"; }

    dataRepository::ViewKey timeHistoryOutputTarget = { "sources" };
    dataRepository::ViewKey timeHistoryOutputFilename = { "filename" };
    dataRepository::ViewKey timeHistoryOutputFormat = { "format" };
    dataRepository::ViewKey timeHistoryRestart = { "restart" };
  } timeHistoryOutputViewKeys;
  /// @endcond

  /**
   * @brief Return PyHistoryOutput type.
   * @return Return PyHistoryOutput type.
   */
#if defined(GEOS_USE_PYGEOSX)
  virtual PyTypeObject * getPythonType() const override;
#endif

private:

  /**
   * @brief Initialize a time history collector to write to an MPI comm-specific file collectively.
   * @param group The ProblemManager cast to a Group
   * @param collector The HistoryCollector to intialize
   */
  void initCollectorParallel( DomainPartition const & domain, HistoryCollection & collector );

  /// The paths of the collectors to collect history from.
  string_array m_collectorPaths;
  /// The file format of the time history file.
  string m_format;
  /// The filename of the time history file.
  string m_filename;
  /// The discrete number of time history states expected to be written to the file
  integer m_recordCount;
  /// The buffered time history output objects for each collector to collect data into and to use to configure/write to file.
  std::vector< std::unique_ptr< BufferedHistoryIO > > m_io;
};
}

#endif
