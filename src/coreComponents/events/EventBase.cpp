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
 * @file EventBase.cpp
 */

#include "EventBase.hpp"
#include <cstring>

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

namespace geos
{

using namespace dataRepository;


EventBase::EventBase( const string & name,
                      Group * const parent ):
  ExecutableGroup( name, parent ),
  m_lastTime( -1.0e100 ),
  m_lastCycle( -1.0e9 ),
  m_eventTarget( "" ),
  m_beginTime( 0.0 ),
  m_endTime( 1e100 ),
  m_forceDt( -1.0 ),
  m_maxEventDt( -1.0 ),
  m_finalDtStretch( 1e-3 ),
  m_targetExactStartStop( 0 ),
  m_currentSubEvent( 0 ),
  m_targetExecFlag( 0 ),
  m_eventForecast( 0 ),
  m_exitFlag( 0 ),
  m_eventCount( 0 ),
  m_timeStepEventCount( 0 ),
  m_eventProgress( 0 ),
  m_currentEventDtRequest( 0.0 ),
  m_target( nullptr )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  // This enables logLevel filtering
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::eventTargetString(), &m_eventTarget ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the object to be executed when the event criteria are met." );

  registerWrapper( viewKeyStruct::beginTimeString(), &m_beginTime ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Start time of this event." );

  registerWrapper( viewKeyStruct::endTimeString(), &m_endTime ).
    setApplyDefaultValue( 1e100 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "End time of this event." );

  registerWrapper( viewKeyStruct::forceDtString(), &m_forceDt ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "While active, this event will request this timestep value (ignoring any children/targets requests)." );

  registerWrapper( viewKeyStruct::maxEventDtString(), &m_maxEventDt ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "While active, this event will request a timestep <= this value (depending upon any child/target requests)." );

  registerWrapper( viewKeyStruct::finalDtStretchString(), &m_finalDtStretch ).
    setApplyDefaultValue( 1e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Allow the final dt request for this event to grow by this percentage to match the endTime exactly." );

  registerWrapper( viewKeyStruct::targetExactStartStopString(), &m_targetExactStartStop ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "If this option is set, the event will reduce its timestep requests to match any specified beginTime/endTimes exactly." );

  registerWrapper( viewKeyStruct::lastTimeString(), &m_lastTime ).
    setApplyDefaultValue( -1.0e100 ).
    setDescription( "Last event occurrence (time)" );

  registerWrapper( viewKeyStruct::lastCycleString(), &m_lastCycle ).
    setApplyDefaultValue( -1.0e9 ).
    setDescription( "Last event occurrence (cycle)" );

  registerWrapper( viewKeyStruct::eventForecastString(), &m_eventForecast ).
    setDescription( "Indicates when the event is expected to execute" );

  registerWrapper( viewKeyStruct::currentSubEventString(), &m_currentSubEvent ).
    setDescription( "Index of the current subevent" );

  registerWrapper( viewKeyStruct::isTargetExecutingString(), &m_targetExecFlag ).
    setDescription( "Index of the current subevent" );
}


EventBase::~EventBase()
{}


EventBase::CatalogInterface::CatalogType & EventBase::getCatalog()
{
  static EventBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

Group * EventBase::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding Event: " << childKey << ", " << childName );
  std::unique_ptr< EventBase > event = EventBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< EventBase >( childName, std::move( event ) );
}


void EventBase::expandObjectCatalogs()
{
  // Only add children if the parent is of type EventManager
  // otherwise, this would fall into a loop
  if( strcmp( this->getParent().getName().c_str(), "Events" ) == 0 )
  {
    for( auto & catalogIter: EventBase::getCatalog())
    {
      createChild( catalogIter.first, catalogIter.first );
    }
  }
}


void EventBase::getTargetReferences()
{
  if( !m_eventTarget.empty())
  {
    try
    {
      m_target = &this->getGroupByPath< ExecutableGroup >( m_eventTarget );
    }
    catch( std::exception const & e )
    {
      throw InputError( e, GEOS_FMT( "Error while reading {}:\n",
                                     getWrapperDataContext( viewKeyStruct::eventTargetString() ) ) );
    }
  }

  this->forSubGroups< EventBase >( []( EventBase & subEvent )
  {
    subEvent.getTargetReferences();
  } );
}


void EventBase::checkEvents( real64 const time,
                             real64 const dt,
                             integer const cycle,
                             DomainPartition & domain )
{
  // Check event status
  if( time < m_beginTime )
  {
    if( dt <= 0 )
    {
      this->setIdle();
    }
    else
    {
      this->setForecast( int( ( m_beginTime - time ) / dt ) );
    }
  }
  else if( time >= m_endTime )
  {
    this->setIdle();
  }
  else
  {
    this->estimateEventTiming( time, dt, cycle, domain );

    // Check sub-events
    this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
    {
      subEvent.checkEvents( time, dt, cycle, domain );
    } );
  }
}


void EventBase::signalToPrepareForExecution( real64 const time,
                                             real64 const dt,
                                             integer const cycle,
                                             DomainPartition & domain )
{
  if( m_target != nullptr )
  {
    m_target->signalToPrepareForExecution( time, dt, cycle, domain );
  }

  this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
  {
    if( subEvent.hasToPrepareForExec() )
    {
      subEvent.signalToPrepareForExecution( time, dt, cycle, domain );
    }
  } );
}


bool EventBase::execute( real64 const time_n,
                         real64 const dt,
                         const integer cycleNumber,
                         integer const,
                         real64 const,
                         DomainPartition & domain )
{
  bool earlyReturn = false;

  // If m_targetExecFlag is set, then the code has resumed at a point
  // after the target has executed.
  if((m_target != nullptr) && (m_targetExecFlag == 0))
  {
    m_targetExecFlag = 1;
    earlyReturn = earlyReturn ||
                  m_target->execute( time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain );
  }

  // Iterate through the sub-event list using the managed integer m_currentSubEvent
  // This allows for  restart runs to pick up where they left off.
  for(; m_currentSubEvent < this->numSubGroups(); ++m_currentSubEvent )
  {
    EventBase * subEvent = static_cast< EventBase * >( this->getSubGroups()[m_currentSubEvent] );

    // Print debug information for logLevel >= 1
    GEOS_LOG_LEVEL_RANK_0( 1,
                           "          SubEvent: " << m_currentSubEvent << " (" << subEvent->getName() << "), dt_request=" << subEvent->getCurrentEventDtRequest() << ", forecast=" <<
                           subEvent->getForecast() );

    if( subEvent->isReadyForExec() )
    {
      earlyReturn = earlyReturn ||
                    subEvent->execute( time_n, dt, cycleNumber, m_eventCount, m_eventProgress, domain );
    }
  }

  // Update the event status
  m_targetExecFlag = 0;
  m_currentSubEvent = 0;
  m_lastTime = time_n;
  m_lastCycle = cycleNumber;

  return earlyReturn;
}


real64 EventBase::getTimestepRequest( real64 const time )
{
  m_currentEventDtRequest = std::numeric_limits< real64 >::max() / 2.0;

  // Events and their targets may request a max dt when active
  if( isActive( time ) )
  {
    if( m_forceDt > 0 )
    {
      // Override the event dt request
      m_currentEventDtRequest = m_forceDt;
    }
    else
    {
      if( m_maxEventDt > 0 )
      {
        // Limit the event dt request
        m_currentEventDtRequest = std::min( m_maxEventDt, m_currentEventDtRequest );
      }

      // Get the event-specific dt request
      m_currentEventDtRequest = std::min( m_currentEventDtRequest, getEventTypeDtRequest( time ));

      // Get the target's dt request if the event has the potential to execute this cycle
      if( ( !this->isIdle() ) && ( m_target != nullptr ) )
      {
        m_currentEventDtRequest = std::min( m_currentEventDtRequest, m_target->getTimestepRequest( time ));
      }

      // Get the sub-event dt requests
      this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
      {
        m_currentEventDtRequest = std::min( m_currentEventDtRequest, subEvent.getTimestepRequest( time ));
      } );
    }
  }

  // Try to respect the start/stop times of the event window
  if( m_targetExactStartStop == 1 )
  {
    // Using this instead of the raw time will prevent falling into dt = 0 loops:
    real64 tmp_t = std::nextafter( time, time + 1.0 );

    if( tmp_t < m_beginTime )
    {
      m_currentEventDtRequest = std::min( m_beginTime - time, m_currentEventDtRequest );
    }
    else if( tmp_t < m_endTime )
    {
      // If the current dt request exceeds the end time, cut it
      // Otherwise, if it falls just short of the end time, grow it.
      if( time + m_currentEventDtRequest * (1.0 + m_finalDtStretch) > m_endTime )
      {
        m_currentEventDtRequest = m_endTime - time;
      }
    }
  }

  return m_currentEventDtRequest;
}


integer EventBase::getExitFlag()
{
  this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
  {
    m_exitFlag += subEvent.getExitFlag();
  } );

  return m_exitFlag;
}



void EventBase::getExecutionOrder( array1d< integer > & eventCounters )
{
  // The first entry counts all events, the second tracks solver events
  m_eventCount = eventCounters[0];
  m_timeStepEventCount = eventCounters[1];

  // Increment counters
  ++eventCounters[0];
  if( m_target != nullptr )
  {
    if( m_target->getTimesteppingBehavior() == ExecutableGroup::TimesteppingBehavior::DeterminesTimeStepSize )
    {
      ++eventCounters[1];
    }
  }

  this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
  {
    subEvent.getExecutionOrder( eventCounters );
  } );
}


void EventBase::setProgressIndicator( array1d< integer > & eventCounters )
{
  // Calculate the event progress indicator
  // This is defined as the percent completion through the executaion loop
  // with respect to the beginning of the event.
  if( eventCounters[1] > 0 )
  {
    m_eventProgress = static_cast< real64 >(m_timeStepEventCount) / static_cast< real64 >(eventCounters[1]);
  }

  // Do this for child events
  this->forSubGroups< EventBase >( [&]( EventBase & subEvent )
  {
    subEvent.setProgressIndicator( eventCounters );
  } );
}

} /* namespace geos */
