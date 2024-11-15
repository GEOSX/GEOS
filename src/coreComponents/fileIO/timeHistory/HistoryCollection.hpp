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

#ifndef GEOS_FILEIO_TIMEHISTORY_HISTORYCOLLECTION_HPP_
#define GEOS_FILEIO_TIMEHISTORY_HISTORYCOLLECTION_HPP_

#include "dataRepository/BufferOpsDevice.hpp"
#include "dataRepository/HistoryDataSpec.hpp"
#include "events/tasks/TaskBase.hpp"
#include "mesh/DomainPartition.hpp"

#include <functional>

namespace geos
{

/**
 * @class HistoryCollection
 *
 * A task class for serializing time history data into a buffer for later I/O.
 */
class HistoryCollection : public TaskBase
{
public:
  /// @copydoc geos::dataRepository::Group::Group(string const & name, Group * const parent)
  HistoryCollection( string const & name, Group * parent ):
    TaskBase( name, parent )
  {}

  /// Type of buffer provider
  using BufferProvider = std::function< buffer_unit_type * ( localIndex ) >;
  /// Type of time buffer provider
  using TimeBufferProvider = std::function< buffer_unit_type *() >;

  /// Forwarding public initializing function...
  void initializePostSubGroups() override {};

  /**
   * @brief Get the number of discrete collection operations this collector conducts.
   * @return The number of collection operations for this collector.
   */
  virtual localIndex numCollectors() const = 0;

  /**
   * @brief Get the metadata for what this collector collects.
   * @param[in] domain The DomainPartition.
   * @param[in] collectionIdx Which collected item to get metadata for.
   * @return A HistoryMetadata object describing  the history data being collected by this collector.
   */
  virtual HistoryMetadata getMetaData( DomainPartition const & domain, localIndex collectionIdx ) const = 0;

  /**
   * @brief Get the name of the object being targeted for collection.
   * @return The collection target's name
   */
  virtual const string & getTargetName() const = 0;

  /**
   * @brief Register a callback that provides the current head of the time history data buffer.
   * @param[in] collectionIdx Which collection item to register the buffer callback for.
   * @param[in] bufferProvider A functional that when invoked returns a pointer to the head of a buffer at least large enough to
   *                           serialize one time step of history data into.
   * @note This is typically meant to callback to BufferedHistoryIO::getBufferHead()
   */
  virtual void registerBufferProvider( localIndex collectionIdx, BufferProvider bufferProvider ) = 0;

  /**
   * @brief Get a metadata object relating the the Time variable itself.
   * @return A HistroyMetadata object describing the Time variable.
   */
  virtual HistoryMetadata getTimeMetaData() const = 0;

  /**
   * @brief Register a callback that gives the current head of the time data buffer.
   * @param[in] timeBufferProvider A functional that when invoked returns a pointer to the head of a buffer at least large enough to
   *                               serialize one instance of the Time variable into.
   * @note This is typically meant to callback to BufferedHistoryIO::GetBufferHead()
   */
  virtual void registerTimeBufferProvider( TimeBufferProvider timeBufferProvider ) = 0;

  /**
   * @brief Get the number of collectors of meta-information (set indices, etc) writing time-independent information during initialization.
   * @return The number of collectors of meta-information for this collector.
   */
  virtual localIndex numMetaDataCollectors() const = 0;

  /**
   * @brief Get a collector of meta-information for this collector.
   * @param[in] metaIdx Which of the meta-info collectors to return. (see HistoryCollection::numMetaDataCollectors()).
   * @return A mutable reference.
   */
  virtual HistoryCollection & getMetaDataCollector( localIndex metaIdx ) = 0;
};

}

#endif
