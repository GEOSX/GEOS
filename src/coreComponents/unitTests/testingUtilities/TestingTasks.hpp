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

/*
 * @file TestingTasks.hpp
 */

#ifndef GEOS_EVENTS_TASKS_TESTINGTASKS_HPP_
#define GEOS_EVENTS_TASKS_TESTINGTASKS_HPP_

#include "events/tasks/TaskBase.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

namespace geos
{
namespace testing
{


/**
 * @brief This Task allows to do checks each timestep during the simulation by calling a provided functor.
 * To be executed, it must be - as every task - declared in the provided xmlInput within the Tasks node and called by an even.
 * As this Group is designed for developpers and not user oriented, it should not appear in the documentation nor in the xsd.
 * Question: could this task be used in the integratedTests ?
 */
class TimeStepChecker : public TaskBase
{
public:
  /**
   * @brief Construct a new TimeStepChecker Task
   * @param name name in xsd
   * @param parent parent group in hierarchy
   */
  TimeStepChecker( string const & name, Group * const parent );

  /**
   * @brief Set the functor that must be called each time this Task is executed.
   * @param func the functor to execute
   */
  void setTimeStepCheckingFunction( std::function< void(real64) > func )
  { m_checkTimeStepFunction = std::function< void(real64) >( func ); }

  /**
   * @return Get the tested time-step count
   */
  integer getTestedTimeStepCount() const
  { return m_timestepId; }

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "TimeStepChecker"; }

  virtual bool execute( real64 time_n, real64 dt, integer cycleNumber,
                        integer eventCounter, real64 eventProgress,
                        DomainPartition & domain );

private:
  std::function< void(real64) > m_checkTimeStepFunction;
  int m_timestepId = 0;
};


} // namespace testing

} // namespace geos

#endif //GEOS_EVENTS_TASKS_TESTINGTASKS_HPP_
