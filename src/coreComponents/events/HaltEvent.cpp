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

#include "HaltEvent.hpp"
#include <sys/time.h>

/**
 * @file HaltEvent.cpp
 */

namespace geos
{

using namespace dataRepository;


HaltEvent::HaltEvent( const string & name,
                      Group * const parent ):
  EventBase( name, parent ),
  m_externalStartTime( 0.0 ),
  m_externalLastTime( 0.0 ),
  m_externalDt( 0.0 ),
  m_maxRuntime( 0.0 )
{
  timeval tim;
  gettimeofday( &tim, nullptr );
  m_externalStartTime = tim.tv_sec + (tim.tv_usec / 1000000.0);
  m_externalLastTime = m_externalStartTime;

  registerWrapper( viewKeyStruct::maxRuntimeString(), &m_maxRuntime ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The maximum allowable runtime for the job." );
}


HaltEvent::~HaltEvent()
{}


void HaltEvent::estimateEventTiming( real64 const GEOS_UNUSED_PARAM( time ),
                                     real64 const GEOS_UNUSED_PARAM( dt ),
                                     integer const GEOS_UNUSED_PARAM( cycle ),
                                     DomainPartition & GEOS_UNUSED_PARAM( domain ))
{
  // Check run time
  timeval tim;
  gettimeofday( &tim, nullptr );
  real64 currentTime = tim.tv_sec + (tim.tv_usec / 1000000.0);

  // Update values
  m_externalDt = currentTime - m_externalLastTime;
  m_externalLastTime = currentTime;
  integer forecast = static_cast< integer >((m_maxRuntime - (currentTime - m_externalStartTime)) / m_externalDt);

  // The timing for the ranks may differ slightly, so synchronize
  // TODO: Only do the communication when you are close to the end?
#ifdef GEOS_USE_MPI
  integer forecast_global;
  MPI_Allreduce( &forecast, &forecast_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
  forecast = forecast_global;
#endif

  setForecast( forecast );

  if( this->isReadyForExec() )
  {
    this->setExitFlag( 1 );
  }
}


REGISTER_CATALOG_ENTRY( EventBase, HaltEvent, string const &, Group * const )
} /* namespace geos */
