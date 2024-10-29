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
 * @file SoloEvent.hpp
 */

#ifndef GEOS_EVENTS_SOLOEVENT_HPP_
#define GEOS_EVENTS_SOLOEVENT_HPP_

#include "events/EventBase.hpp"

namespace geos
{

/**
 * @class SoloEvent
 *
 * An event type for events that occur only once.
 */
class SoloEvent : public EventBase
{
public:

  /// @copydoc geos::dataRepository::Group::Group( string const & name, Group * const parent )
  SoloEvent( const string & name,
             Group * const parent );

  /// Destructor
  virtual ~SoloEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string catalogName() { return "SoloEvent"; }

  /**
   * @copydoc EventBase::estimateEventTiming()
   */
  virtual void estimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    DomainPartition & domain ) override;

  /**
   * @copydoc EventBase::getEventTypeDtRequest()
   */
  virtual real64 getEventTypeDtRequest( real64 const time ) override;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * targetTimeString() { return "targetTime"; }
    static constexpr char const * targetCycleString() { return "targetCycle"; }
    static constexpr char const * targetExactTimestepString() { return "targetExactTimestep"; }

    dataRepository::ViewKey targetTime = { targetTimeString() };
    dataRepository::ViewKey targetCycle = { targetCycleString() };
    dataRepository::ViewKey targetExactTimestep = { targetExactTimestepString() };
  } SoloEventViewKeys;
  /// @endcond

  /// The target time
  real64 m_targetTime;
  /// The target cycle
  integer m_targetCycle;
  /// Whether to target the exact time step
  integer m_targetExactTimestep;

};

} /* namespace geos */

#endif /* GEOS_EVENTS_SOLOEVENT_HPP_ */
