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
  /// External start time
  real64 m_externalStartTime;
  /// External last time
  real64 m_externalLastTime;
  /// External time increment
  real64 m_externalDt;
  /// Max runtime
  real64 m_maxRuntime;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr char const * maxRuntimeString() { return "maxRuntime"; }

    dataRepository::ViewKey maxRuntime = { maxRuntimeString() };
  } haltEventViewKeys;
  /// @endcond

};

} /* namespace geos */

#endif /* GEOS_EVENTS_HALTEVENT_HPP_ */
