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
 * @file HaltEvent.hpp
 */


#ifndef GEOS_EVENTS_HALTEVENT_HPP_
#define GEOS_EVENTS_HALTEVENT_HPP_

#include "events/EventBase.hpp"

namespace geos
{

/**
 * @class HaltEvent
 * An event type that is designed to look at the external clock.
 * This is useful for managing wall time limitations.
 */
class HaltEvent : public EventBase
{
public:
  /**
   * @brief Main constructor.
   * @param name The name of the object in the data repository.
   * @param parent The parent of this object in the data repository.
   **/
  HaltEvent( const string & name,
             Group * const parent );

  /// Destructor
  virtual ~HaltEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "HaltEvent"; }

  /**
   * @copydoc EventBase::estimateEventTiming()
   * @note This event is designed to look at the external clock. Currently,
   * if the event is triggered it will set a flag, which will
   * instruct the code to exit.  This is useful for managing walltime
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

  /// External start time
  real64 m_externalStartTime;
  /// External last time
  real64 m_externalLastTime;
  /// External time increment
  real64 m_externalDt;
  /// Max runtime
  real64 m_maxRuntime;
  /// A pointer to an optional function
  dataRepository::Group * m_functionTarget;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * maxRuntimeString() { return "maxRuntime"; }
    static constexpr char const * functionNameString() { return "function"; }
    static constexpr char const * functionInputObjectString() { return "object"; }
    static constexpr char const * functionInputSetnameString() { return "set"; }
    static constexpr char const * functionSetNamesString() { return "setNames"; }
    static constexpr char const * functionStatOptionString() { return "stat"; }
    static constexpr char const * eventThresholdString() { return "threshold"; }

    dataRepository::ViewKey maxRuntime = { maxRuntimeString() };
    dataRepository::ViewKey functionName = { functionNameString() };
    dataRepository::ViewKey functionInputObject = { functionInputObjectString() };
    dataRepository::ViewKey functionInputSetname = { functionInputSetnameString() };
    dataRepository::ViewKey functionSetNames = { functionSetNamesString() };
    dataRepository::ViewKey functionStatOption = { functionStatOptionString() };
    dataRepository::ViewKey eventThreshold = { eventThresholdString() };
  } haltEventViewKeys;
  /// @endcond

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

#endif /* GEOS_EVENTS_HALTEVENT_HPP_ */
