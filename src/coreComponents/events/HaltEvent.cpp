/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "HaltEvent.hpp"
#include <sys/time.h>

#include "common/Format.hpp"
#include "functions/FunctionManager.hpp"

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
  m_maxRuntime( 0.0 ),
  m_functionTarget( nullptr ),
  m_functionName(),
  m_functionInputObject(),
  m_functionInputSetname(),
  m_functionStatOption( 0 ),
  m_eventThreshold( 0.0 )
{
  timeval tim;
  gettimeofday( &tim, nullptr );
  m_externalStartTime = tim.tv_sec + (tim.tv_usec / 1000000.0);
  m_externalLastTime = m_externalStartTime;

  registerWrapper( viewKeyStruct::maxRuntimeString(), &m_maxRuntime ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The maximum allowable runtime for the job." );

  registerWrapper( viewKeyStruct::functionNameString(), &m_functionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of an optional function to evaluate when the time/cycle criteria are met."
                    "If the result is greater than the specified eventThreshold, the function will continue to execute." );

  registerWrapper( viewKeyStruct::functionInputObjectString(), &m_functionInputObject ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "If the optional function requires an object as an input, specify its path here." );

  registerWrapper( viewKeyStruct::functionInputSetnameString(), &m_functionInputSetname ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "If the optional function is applied to an object, specify the setname to evaluate (default = everything)." );

  registerWrapper( viewKeyStruct::functionStatOptionString(), &m_functionStatOption ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "If the optional function is applied to an object, specify the statistic to compare to the eventThreshold."
                    "The current options include: min, avg, and max." );

  registerWrapper( viewKeyStruct::eventThresholdString(), &m_eventThreshold ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "If the optional function is used, the event will execute if the value returned by the function exceeds this threshold." );

}


HaltEvent::~HaltEvent()
{}


void HaltEvent::estimateEventTiming( real64 const time,
                                     real64 const dt,
                                     integer const cycle,
                                     DomainPartition & domain )
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

  if (this->isReadyForExec())
  {
    this->setExitFlag(1);
  }
  else if (!m_functionName.empty())
  {
    checkOptionalFunctionThreshold(time, dt, cycle, domain);
  }
}

void HaltEvent::checkOptionalFunctionThreshold(real64 const time,
                                               real64 const GEOS_UNUSED_PARAM(dt),
                                               integer const GEOS_UNUSED_PARAM(cycle),
                                               DomainPartition &GEOS_UNUSED_PARAM(domain))
{
  // Grab the function
  FunctionManager & functionManager = FunctionManager::getInstance();
  FunctionBase & function = functionManager.getGroup< FunctionBase >( m_functionName );

  real64 result = 0.0;
  if( m_functionInputObject.empty())
  {
    // This is a time-only function
    result = function.evaluate( &time );
  }
  else
  {
    // Link the target object
    if( m_functionTarget == nullptr )
    {
      m_functionTarget = &this->getGroupByPath( m_functionInputObject );
    }

    // Get the set
    SortedArray< localIndex > mySet;
    if( m_functionInputSetname.empty())
    {
      for( localIndex ii=0; ii<m_functionTarget->size(); ++ii )
      {
        mySet.insert( ii );
      }
    }
    else
    {
      dataRepository::Group const & sets = m_functionTarget->getGroup( haltEventViewKeys.functionSetNames );
      SortedArrayView< localIndex const > const &
      functionSet = sets.getReference< SortedArray< localIndex > >( m_functionInputSetname );

      for( localIndex const index : functionSet )
      {
        mySet.insert( index );
      }
    }

    // Find the function (min, average, max)
    real64_array stats = function.evaluateStats( *m_functionTarget, time, mySet );
    result = stats[m_functionStatOption];

    //GEOS_LOG(GEOS_FMT("-----\n SGAS: {}, Threshold: {}",result, m_eventThreshold));

    // Because the function applied to an object may differ by rank, synchronize
    // (Note: this shouldn't occur very often, since it is only called if the base forecast <= 0)
#ifdef GEOSX_USE_MPI
    real64 result_global;
    switch (m_functionStatOption)
    {
    case 0:
    {
      MPI_Allreduce(&result, &result_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      result = result_global;
      break;
    }
    case 1:
    {
      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      MPI_Allreduce(&result, &result_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      result = result_global / nprocs;
      break;
    }
    case 2:
    {
      MPI_Allreduce(&result, &result_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      result = result_global;
    }
    }
#endif   
  }

  // Forcast event
  if( result > m_eventThreshold )
  {
    this->setExitFlag( 1 );
  }
}


REGISTER_CATALOG_ENTRY( EventBase, HaltEvent, string const &, Group * const )
} /* namespace geos */
