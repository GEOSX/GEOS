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
 * @file EventBase.hpp
 */

#ifndef GEOS_EVENTS_EVENTSBASE_HPP_
#define GEOS_EVENTS_EVENTSBASE_HPP_

#include "dataRepository/Group.hpp"
#include "dataRepository/ExecutableGroup.hpp"


namespace geos
{

/**
 * @class EventBase
 * A base class for managing code event targets (solver applications, etc.)
 */
class EventBase : public ExecutableGroup
{
public:
  /**
   * @brief Main constructor.
   * @param name The name of the object in the data repository.
   * @param parent The parent of this object in the data repository.
   **/
  explicit EventBase( string const & name,
                      Group * const parent );

  /// Destructor
  virtual ~EventBase() override;

  /**
   * @brief If the event forecast is equal to 1, then signal the targets to prepare for execution
   *        during the next cycle.
   * @param time The current simulation time.
   * @param dt The current time increment.
   * @param cycle The current cycle.
   * @param domain The DomainPartition the event is occuring on up-casted to a Group.
   */
  virtual void signalToPrepareForExecution( real64 const time,
                                            real64 const dt,
                                            integer const cycle,
                                            DomainPartition & domain ) override;
  /**
   * @brief If the event forecast is equal to 0, then call the step function on its target and/or children.
   * @copydoc ExecutableGroup::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Call the execute method on the target and/or children if present.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle The current simulation cycle.
   * @param domain The DomainPartition up-casted to a Group.
   */
  void step( real64 const time,
             real64 const dt,
             integer const cycle,
             dataRepository::Group * domain );

  /**
   * @copydoc dataRepository::Group::createChild()
   *
   * An event may have an arbitrary number of sub-events defined as children in the input xml.
   * e.g.:
   * @code{.unparsed}
   * <Events>
   *   <PeriodicEvent name="base_event" ...>
   *     <PeriodicEvent name="sub_event" .../>
   *     ...
   *   </PeriodicEvent>
   * </Events>
   * @endcode
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @brief Process input data to retrieve targeted objects internally.
   * The target object for an event may be specified via the keyword "target" in the input xml.
   * This string is empty by default and uses getGroupByPath() method in Group, which returns
   * a pointer to the target using a unix-style path as an input (both absolute and relative paths work).
   * This involves a lot of string parsing, so we do it once during initialization.
   */
  void getTargetReferences();

  /**
   * @brief Events are triggered based upon their forecast values, which are defined
   *        as the expected number of code cycles before they are executed.  This method
   *        will call EstimateEventTiming (defined in each subclass) on this event and
   *        its children.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle the current simulation cycle.
   * @param domain The problem domain up-cast to a Group.
   */
  virtual void checkEvents( real64 const time,
                            real64 const dt,
                            integer const cycle,
                            DomainPartition & domain );

  /**
   * @brief Perform the calculations to estimate the timing of the event.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle the current simulation cycle.
   * @param domain The problem domain up-cast to a Group.
   */
  virtual void estimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    DomainPartition & domain ) = 0;

  /**
   * @brief Collect time-step size requests from targets and/or children.
   * @param time The current simulation time.
   * @return The requested time step.
   */
  virtual real64 getTimestepRequest( real64 const time ) override;

  /**
   * @brief Get event-specifit dt requests.
   * @param time The current simulation time.
   * @return The requested time step.
   */
  virtual real64 getEventTypeDtRequest( real64 const time )
  {
    GEOS_UNUSED_VAR( time );
    return std::numeric_limits< real64 >::max();
  }

  /**
   * @brief Helper function to validate the consistency of the event input
   * @note We cannot use postInputInitialization here because we can perform the validation only after the m_target pointer is set
   */
  virtual void validate() const {};

  /**
   * @brief Count the number of events/sub-events
   * @param[out] eventCounters The event count for each event/sub-event.
   */
  void getExecutionOrder( array1d< integer > & eventCounters );

  /**
   * @brief Update the event progress for the event/sub-events.
   * @note This method is used to determine how to handle the timestamp for an event
   * @note If the event occurs after anything targeting a SolverBase object, then
   *       set the m_isPostSolverEvent flag.  If set, then the time passed to the target
   *       will be time + dt.
   * @param eventCounters The event count for each event/sub-event.
   */
  void setProgressIndicator( array1d< integer > & eventCounters );

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * eventTargetString() { return "target"; }
    static constexpr char const * beginTimeString() { return "beginTime"; }
    static constexpr char const * endTimeString() { return "endTime"; }
    static constexpr char const * forceDtString() { return "forceDt"; }
    static constexpr char const * maxEventDtString() { return "maxEventDt"; }
    static constexpr char const * lastTimeString() { return "lastTime"; }
    static constexpr char const * lastCycleString() { return "lastCycle"; }
    static constexpr char const * eventForecastString() { return "eventForecast"; }
    static constexpr char const * targetExactStartStopString() { return "targetExactStartStop"; }
    static constexpr char const * currentSubEventString() { return "currentSubEvent"; }
    static constexpr char const * isTargetExecutingString() { return "isTargetExecuting"; }
    static constexpr char const * finalDtStretchString() { return "finalDtStretch"; }

    dataRepository::ViewKey eventTarget = { eventTargetString() };
    dataRepository::ViewKey beginTime = { beginTimeString() };
    dataRepository::ViewKey endTime = { endTimeString() };
    dataRepository::ViewKey forceDt = { forceDtString() };
    dataRepository::ViewKey maxEventDt = { maxEventDtString() };
    dataRepository::ViewKey lastTime = { lastTimeString() };
    dataRepository::ViewKey lastCycle = { lastCycleString() };
    dataRepository::ViewKey eventForecast = { eventForecastString() };
    dataRepository::ViewKey targetExactStartStop = { targetExactStartStopString() };
    dataRepository::ViewKey currentSubEvent = { currentSubEventString() };
    dataRepository::ViewKey isTargetExecuting = { isTargetExecutingString() };
  } viewKeys;
  /// @endcond

  /// Catalog interface
  using CatalogInterface = dataRepository::CatalogInterface< EventBase, string const &, Group * const >;
  /// @copydoc dataRepository::Group::getCatalog()
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Get the sum of the exit flags for the event/sub-events from the last execution.
   * @return The sum of the exit flags for the event/sub-events.
   */
  integer getExitFlag();

  /**
   * @brief Set this event objects exit flag.
   * @param flag The exit flag value.
   */
  void setExitFlag( integer flag ){ m_exitFlag = flag; }

  /**
   * @brief Get the current time increment request for this event.
   * @return The current time increment request.
   */
  real64  getCurrentEventDtRequest() const { return m_currentEventDtRequest; }

  /**
   * @brief Get the forecast of the current event.
   * @return The forecast.
   *
   * The `getForecast` getter only exists for debugging purpose.
   * Prefer the predicate versions instead (isReadyForExec(), hasToPrepareForExec(), isIdle()).
   */
  integer getForecast() const
  { return m_eventForecast; }

  /**
   * @brief Check if the event is ready for execution.
   * @return @p true if ready, @p false otherwise.
   */
  bool isReadyForExec() const
  { return m_eventForecast <= 0; }

  /**
   * @brief Check if the event must be preparing for execution.
   * @return @p true if it must prepare, @p false otherwise.
   */
  bool hasToPrepareForExec() const
  { return m_eventForecast == 1; }

  /**
   * @brief Check if the event is idle.
   * @return @p true if it is idle, @p false otherwise.
   */
  bool isIdle() const
  { return m_eventForecast > 1; }

  /**
   * @brief Get the string name of the target.
   * @return @p string name of the target.
   */
  string getEventName() const
  {
    return m_eventTarget;
  }

protected:

  /**
   * @brief Define the event as ready for execution.
   */
  void setReadyForExec()
  { m_eventForecast = 0; }

  /**
   * @brief Define that the event should prepare for execution.
   */
  void setPrepareForExec()
  { m_eventForecast = 1; }

  /**
   * @brief Define the event as idle.
   */
  void setIdle()
  { m_eventForecast = std::numeric_limits< decltype( m_eventForecast ) >::max(); }

  /**
   * @brief Sets the forecast
   * @param forecast The forecast provided as an integer.
   *
   * If the forecast is 0 or below, the event is considered being "ready for exec".
   * If it equals 1, it is in "prepare for exec" status. Above, the event is "idle".
   * If you can, you may prefer the dedicated setters (setReadyForExec(), setPrepareForExec(), setIdle()).
   */
  void setForecast( integer forecast )
  { m_eventForecast = forecast; }

  /**
   * @brief Get the target of this event.
   * @return The target of this event.
   */
  ExecutableGroup * getEventTarget() const
  { return m_target; }

  /**
   * @brief Is the event active?
   * @param time The time at which we want to check if the event is active.
   * @return @p true if active, @p false otherwise.
   */
  bool isActive( real64 const time ) const
  { return ( time >= m_beginTime ) && ( time < m_endTime ); }

  /**
   * @brief Is the event ready for cleanup?
   * @param time The time at which we want to check if the event is cleanup.
   * @return @p true if is ready for cleanup, @p false otherwise.
   */
  bool isReadyForCleanup( real64 const time ) const
  { return isActive( time ) || isZero( time - m_endTime ); }


  /// The last time the event occurred.
  real64 m_lastTime;
  /// The last cycle the event occurred.
  integer m_lastCycle;

private:
  string m_eventTarget;
  real64 m_beginTime;
  real64 m_endTime;
  real64 m_forceDt;
  real64 m_maxEventDt;
  real64 m_finalDtStretch;
  integer m_targetExactStartStop;
  integer m_currentSubEvent;
  integer m_targetExecFlag;
  integer m_eventForecast;
  integer m_exitFlag;
  integer m_eventCount;
  integer m_timeStepEventCount;
  real64 m_eventProgress;
  real64 m_currentEventDtRequest;

  /// A pointer to the optional event target
  ExecutableGroup * m_target;
};

} /* namespace geos */

#endif /* GEOS_EVENTS_EVENTSBASE_HPP_ */
