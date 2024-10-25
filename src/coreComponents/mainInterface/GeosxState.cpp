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

// Source includes
#include "GeosxState.hpp"
#include "dataRepository/Utilities.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "common/Timer.hpp"

// TPL includes
#include <conduit.hpp>

#if defined( GEOS_USE_CALIPER )
  #include <caliper/cali-manager.h>
#endif

// System includes
#include <ostream>

namespace geos
{

GeosxState * currentGlobalState = nullptr;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState & getGlobalState()
{
  GEOS_ERROR_IF( currentGlobalState == nullptr,
                 "The state has not been created." );

  return *currentGlobalState;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string durationToString( std::chrono::system_clock::duration const duration )
{
  // If we want to print HH::MM::SS (maybe in addition to seconds-only):
  // return GEOS_FMT( "{:%T}", duration );
  double const seconds = std::chrono::duration_cast< std::chrono::milliseconds >( duration ).count() / 1000.0;
  return GEOS_FMT( "{:>20.3f}s", seconds );
}

std::ostream & operator<<( std::ostream & os, State const state )
{
  if( state == State::UNINITIALIZED )
  {
    return os << "State::UNINITIALIZED";
  }
  if( state == State::INITIALIZED )
  {
    return os << "State::INITIALIZED";
  }
  if( state == State::READY_TO_RUN )
  {
    return os << "State::READY_TO_RUN";
  }
  if( state == State::COMPLETED )
  {
    return os << "State::COMPLETED";
  }

  GEOS_ERROR( "Unrecognized state. The integral value is: " << static_cast< int >( state ) );
  return os;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::GeosxState( std::unique_ptr< CommandLineOptions > && commandLineOptions ):
  m_state( State::UNINITIALIZED ),
  m_commandLineOptions( std::move( commandLineOptions ) ),
  m_rootNode( std::make_unique< conduit::Node >() ),
  m_problemManager( nullptr ),
  m_commTools( std::make_unique< CommunicationTools >() ),
#if defined( GEOS_USE_CALIPER )
  m_caliperManager( std::make_unique< cali::ConfigManager >() ),
#endif
  m_initTime(),
  m_runTime()
{
  Timer timer( m_initTime );

#if defined( GEOS_USE_CALIPER )
  setupCaliper( *m_caliperManager, getCommandLineOptions() );
#endif

  string restartFileName;
  if( ProblemManager::parseRestart( restartFileName, getCommandLineOptions() ) )
  {
    GEOS_LOG_RANK_0( "Loading restart file " << restartFileName );
    dataRepository::loadTree( restartFileName, getRootConduitNode() );
  }

  m_problemManager = std::make_unique< ProblemManager >( getRootConduitNode() );

  GEOS_ERROR_IF( currentGlobalState != nullptr, "Only one state can exist at a time." );
  currentGlobalState = this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::~GeosxState()
{
#if defined( GEOS_USE_CALIPER )
  m_caliperManager->flush();
#endif

  GEOS_ERROR_IF( currentGlobalState != this, "This shouldn't be possible." );
  currentGlobalState = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GeosxState::initializeDataRepository()
{
  GEOS_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOS_THROW_IF_NE( m_state, State::UNINITIALIZED, std::logic_error );

  getProblemManager().parseCommandLineInput();

  if( !getProblemManager().getSchemaFileName().empty() )
  {
    getProblemManager().generateDocumentation();
    m_state = State::INITIALIZED;
    return false;
  }

  getProblemManager().parseInputFile();
  getProblemManager().problemSetup();

  m_state = State::INITIALIZED;

  if( m_commandLineOptions->printMemoryUsage >= 0.0 )
  {
    dataRepository::printMemoryAllocation( getProblemManager(), 0, m_commandLineOptions->printMemoryUsage );
  }

  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::applyInitialConditions()
{
  GEOS_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOS_THROW_IF_NE( m_state, State::INITIALIZED, std::logic_error );

  getProblemManager().applyInitialConditions();

  if( getCommandLineOptions().beginFromRestart )
  {
    getProblemManager().readRestartOverwrite();
  }

  m_state = State::READY_TO_RUN;
  MpiWrapper::barrier( MPI_COMM_GEOS );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::run()
{
  GEOS_MARK_FUNCTION;
  Timer timer( m_runTime );

  GEOS_THROW_IF_NE( m_state, State::READY_TO_RUN, std::logic_error );

  if( !getProblemManager().runSimulation() )
  {
    m_state = State::COMPLETED;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataRepository::Group & GeosxState::getProblemManagerAsGroup()
{ return getProblemManager(); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FieldSpecificationManager & GeosxState::getFieldSpecificationManager()
{ return getProblemManager().getFieldSpecificationManager(); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FunctionManager & GeosxState::getFunctionManager()
{ return getProblemManager().getFunctionManager(); }

} // namespace geos
