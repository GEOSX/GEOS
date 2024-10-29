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


#ifndef GEOS_EVENTS_EVENTMANAGER_HPP_
#define GEOS_EVENTS_EVENTMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "EventBase.hpp"

namespace geos
{

class DomainPartition;

namespace dataRepository
{
namespace keys
{
string const Events( "Events" );
}
}

/**
 * @class EventManager
 *
 * A class for managing code events.
 */
class EventManager : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for the EventManager
   * @param[in] name the name of the EventManager
   * @param[in] parent group this EventManager
   */
  EventManager( string const & name,
                Group * const parent );

  /**
   * @brief Default destructor for the EventManager
   */
  virtual ~EventManager() override;

  /**
   * @brief Create a child EventManager
   * @param[in] childKey the key of the Event to be added
   * @param[in] childName the name of the child to be added
   * @return the Event
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief This method is used to expand any catalogs in the data structure
   */
  virtual void expandObjectCatalogs() override;

  /**
   * @brief The main execution loop for the code.
   * @param[in] domain the current DomainPartition on which the Event will be ran.
   * @return True iff the simulation exited early, and needs to be run again to completion.
   * @details During each cycle, it will:
   *   - Calculate the event forecast (number of cycles until its expected execution)
   *   - Signal an event to prepare (forecast == 1)
   *   - Execute an event (forecast == 0)
   *   - Determine dt for the next cycle
   *   - Advance time, cycle, etc.
   */
  bool run( DomainPartition & domain );

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{
  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * minTimeString() { return "minTime"; }
    static constexpr char const * maxTimeString() { return "maxTime"; }
    static constexpr char const * maxCycleString() { return "maxCycle"; }

    static constexpr char const * timeString() { return "time"; }
    static constexpr char const * dtString() { return "dt"; }
    static constexpr char const * cycleString() { return "cycle"; }
    static constexpr char const * currentSubEventString() { return "currentSubEvent"; }

    static constexpr char const * timeOutputFormat() { return "timeOutputFormat"; }


    dataRepository::ViewKey time = { "time" };
    dataRepository::ViewKey dt = { "dt" };
    dataRepository::ViewKey cycle = { "cycle" };
    dataRepository::ViewKey minTime = { "minTime" };
    dataRepository::ViewKey maxTime = { "maxTime" };
    dataRepository::ViewKey maxCycle = { "maxCycle" };
    dataRepository::ViewKey currentSubEvent = { "currentSubEvent" };
  } viewKeys;
  /// @endcond
  ///@}

  /// Alias to access the object catalog for EventBase derived types.
  using CatalogInterface = dataRepository::CatalogInterface< EventBase, string const &, Group * const >;

  /// @copydoc dataRepository::Group::getCatalog()
  static CatalogInterface::CatalogType & getCatalog();

  /// enum class defining the format of the time output in the log
  enum class TimeOutputFormat : integer
  {
    seconds,
    minutes,
    hours,
    days,
    years,
    full
  };

private:


  /**
   * @brief ouput time information to the log
   *
   */
  void outputTime() const;

  /// Min time for a simulation
  real64 m_minTime;

  /// Max time for a simulation
  real64 m_maxTime;

  /// Maximum number of cycles for a simulation
  integer m_maxCycle;

  /// Simulation timestamp at the beginning of the cycle
  real64 m_time;

  /// Current timestep request
  real64 m_dt;

  /// Current cycle
  integer m_cycle;

  /// Current subevent index
  integer m_currentSubEvent;

  /// time output type
  TimeOutputFormat m_timeOutputFormat;
};

/// valid strings fort the time output enum.
ENUM_STRINGS( EventManager::TimeOutputFormat,
              "seconds",
              "minutes",
              "hours",
              "days",
              "years",
              "full" );

} /* namespace geos */

#endif /* GEOS_EVENTS_EVENTMANAGER_HPP_ */
