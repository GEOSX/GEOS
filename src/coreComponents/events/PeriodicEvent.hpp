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
 * @file PeriodicEvent.hpp
 */

#ifndef GEOS_EVENTS_PERIODICEVENT_HPP_
#define GEOS_EVENTS_PERIODICEVENT_HPP_

#include "events/EventBase.hpp"

namespace geos
{

/**
 * @class PeriodicEvent
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class PeriodicEvent : public EventBase
{
public:

  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  PeriodicEvent( const string & name,
                 Group * const parent );

  /// Destructor
  virtual ~PeriodicEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "PeriodicEvent"; }

  /**
   * @copydoc EventBase::estimateEventTiming()
   * @note Estimate the expected number of cycles until an event is expected to trigger.
   *       The event frequency can be specified in terms of:
   *        - time (timeFrequency > 0, units = seconds)
   *        - or cycle (cycleFrequency >= 0, units = cycles)
   * @note In addition, there is an optional function input that will be called if the
   * the nominal forecast (based on timing) is zero.
   */
  virtual void estimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    DomainPartition & domain ) override;

  /**
   * @brief Determine if an optional function f should be called, and call it if so.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle The current simulation cycle.
   * @param domain The DomainPartition upcast to a Group.
   * If the event forecast is zero, and an optional function (f) is specified, then
   * this method will be called to see if the event should be triggered or ignored.
   * For example, this could be used to periodically check the condition of the mesh,
   * and trigger a cleanup if necessary.
   *
   * If functionInputObject is not specified:
   *   - The argument to the function will be the current time
   *   - The event will be executed if f(t) >= eventThreshold
   *
   * If functionInputObject is specified:
   *   - The function will be called on the object, with the arguments given by functionInputSetname
   *   - The function manager will return a set of statistics
   *   - functionStatOption selects the statistic to compare against the eventThreshold (0 = min, 1 = average, 2 = max)
   *   - The event will be executed if f(object, arguments)[stat] >= eventThreshold
   */
  void checkOptionalFunctionThreshold( real64 const time,
                                       real64 const dt,
                                       integer const cycle,
                                       DomainPartition & domain );

  /**
   * @copydoc EventBase::getEventTypeDtRequest()
   * Grab the next time-step.  If requested, then limit the requested
   * dt to exactly match the time frequency
   */
  virtual real64 getEventTypeDtRequest( real64 const time ) override;

  /**
   * @copydoc ExecutableGroup::cleanup()
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @copydoc EventBase::validate
   */
  virtual void validate() const override;

  /// A pointer to an optional function
  dataRepository::Group * m_functionTarget;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * timeFrequencyString() { return "timeFrequency"; }
    static constexpr char const * cycleFrequencyString() { return "cycleFrequency"; }
    static constexpr char const * targetExactTimestepString() { return "targetExactTimestep"; }
    static constexpr char const * functionNameString() { return "function"; }
    static constexpr char const * functionInputObjectString() { return "object"; }
    static constexpr char const * functionInputSetnameString() { return "set"; }
    static constexpr char const * functionSetNamesString() { return "setNames"; }
    static constexpr char const * functionStatOptionString() { return "stat"; }
    static constexpr char const * eventThresholdString() { return "threshold"; }

    dataRepository::ViewKey timeFrequency = { timeFrequencyString() };
    dataRepository::ViewKey cycleFrequency = { cycleFrequencyString() };
    dataRepository::ViewKey targetExactTimestep = { targetExactTimestepString() };
    dataRepository::ViewKey functionName = { functionNameString() };
    dataRepository::ViewKey functionInputObject = { functionInputObjectString() };
    dataRepository::ViewKey functionInputSetname = { functionInputSetnameString() };
    dataRepository::ViewKey functionSetNames = { functionSetNamesString() };
    dataRepository::ViewKey functionStatOption = { functionStatOptionString() };
    dataRepository::ViewKey eventThreshold = { eventThresholdString() };
  } periodicEventViewKeys;
  /// @endcond

  /// The event time frequency
  real64 m_timeFrequency;
  /// The event cycle frequency
  integer m_cycleFrequency;
  /// Whether to target the exact timestep
  integer m_targetExactTimestep;
  /// The optional function's name
  string m_functionName;
  /// The name of the optional function input object
  string m_functionInputObject;
  /// The name of the optional function input set
  string m_functionInputSetname;
  /// The optional funciton's statistic option
  integer m_functionStatOption;
  /// The event threshold
  real64 m_eventThreshold;

};

} /* namespace geos */

#endif /* GEOS_EVENTS_PERIODICEVENT_HPP_ */
