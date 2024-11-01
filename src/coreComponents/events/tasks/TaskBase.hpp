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
 * @file TaskBase.hpp
 */

#ifndef GEOS_EVENTS_TASKS_TASKBASE_HPP_
#define GEOS_EVENTS_TASKS_TASKBASE_HPP_


#include <limits>

#include "dataRepository/ExecutableGroup.hpp"
#include "common/DataTypes.hpp"
namespace geos
{

/**
 * @class TaskBase
 * A class to execute a task when directed by the TasksManager
 */
class TaskBase : public ExecutableGroup
{
public:
  /// @copydoc geos::dataRepository::Group::Group(string const & name, Group * const parent)
  explicit TaskBase( string const & name,
                     Group * const parent );
  virtual ~TaskBase( ) override;

  /// The catalog interface type for TaskBase
  using CatalogInterface = dataRepository::CatalogInterface< TaskBase, string const &, Group * const >;
  /**
   * @brief Get the catalog interface for the TaskBase
   * @return the Catalog for TaskBase.
   */
  static CatalogInterface::CatalogType & getCatalog();

  /// @copydoc geos::ExecutableGroup::execute
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    GEOS_UNUSED_VAR( time_n );
    GEOS_UNUSED_VAR( dt );
    GEOS_UNUSED_VAR( cycleNumber );
    GEOS_UNUSED_VAR( eventCounter );
    GEOS_UNUSED_VAR( eventProgress );
    GEOS_UNUSED_VAR( domain );
    GEOS_ERROR( "NOT IMPLEMENTED" );
    return false;
  }

  /// @copydoc geos::dataRepository::Group::postInputInitialization( )
  void postInputInitialization() override;
};

} /* namespace */

#endif /* GEOS_EVENTS_TASKS_TASKBASE_HPP_ */
