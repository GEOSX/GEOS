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
 * @file PoromechanicsInitialization.cpp
 */

#include "PoromechanicsInitialization.hpp"

#include "events/tasks/TasksManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsStatistics.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsSolver.hpp"
#include "physicsSolvers/multiphysics/LogLevelsInfo.hpp"

namespace geos
{

using namespace dataRepository;

PoromechanicsInitialization::PoromechanicsInitialization( const string & name, Group * const parent ) :
  TaskBase( name, parent ),
  m_solidMechanicsSolverName(),
  m_solidMechanicsStatistics(),
  m_solidMechanicsStateResetTask( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::solidMechanicsSolverNameString(), &m_solidMechanicsSolverName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the poromechanics solver" );

  registerWrapper( viewKeyStruct::solidMechanicsStatisticsNameString(), &m_solidMechanicsStatisticsName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "" ).
    setDescription( "Name of the solid mechanics statistics" );

  addLogLevel< logInfo::Initialization >();
}

PoromechanicsInitialization::~PoromechanicsInitialization() {}

void PoromechanicsInitialization::postInputInitialization()
{
  Group & problemManager = this->getGroupByPath( "/Problem" );
  Group & physicsSolverManager = problemManager.getGroup( "Solvers" );

  GEOS_THROW_IF( !physicsSolverManager.hasGroup( m_solidMechanicsSolverName ),
                 GEOS_FMT( "{}: solid mechanics solver named {} not found",
                           getWrapperDataContext( viewKeyStruct::solidMechanicsSolverNameString() ),
                           m_solidMechanicsSolverName ),
                 InputError );

  m_solidMechanicsSolver = &physicsSolverManager.getGroup< SolidMechanicsLagrangianFEM >( m_solidMechanicsSolverName );

  if( !m_solidMechanicsStatisticsName.empty())
  {
    TasksManager & tasksManager = problemManager.getGroup< TasksManager >( "Tasks" );

    GEOS_THROW_IF( !tasksManager.hasGroup( m_solidMechanicsStatisticsName ),
                   GEOS_FMT( "{}: {} task named {} not found",
                             getWrapperDataContext( viewKeyStruct::solidMechanicsStatisticsNameString() ),
                             SolidMechanicsStatistics::catalogName(),
                             m_solidMechanicsStatisticsName ),
                   InputError );

    m_solidMechanicsStatistics = &tasksManager.getGroup< SolidMechanicsStatistics >( m_solidMechanicsStatisticsName );
  }

  m_solidMechanicsStateResetTask.setLogLevel( getLogLevel());
  m_solidMechanicsStateResetTask.m_solidSolverName = m_solidMechanicsSolver->getName();
  m_solidMechanicsStateResetTask.postInputInitialization();
}

bool PoromechanicsInitialization::execute( real64 const time_n,
         real64 const dt,
         integer const cycleNumber,
         integer const eventCounter,
         real64 const eventProgress,
         DomainPartition & domain )
{
  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Initialization, GEOS_FMT( "Task `{}`: at time {}s, physics solver `{}` is set to perform stress initialization during the next time step(s)",
                                                                 getName(), time_n, m_solidMechanicsSolverName ) );
  m_solidMechanicsSolver->setStressInitialization( true );

  m_solidMechanicsStateResetTask.execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  m_solidMechanicsSolver->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Initialization, GEOS_FMT( "Task `{}`: at time {}s, physics solver `{}` has completed stress initialization",
                                                                 getName(), time_n + dt, m_solidMechanicsSolverName ) );
  m_solidMechanicsSolver->setStressInitialization( false );

  if( m_solidMechanicsStatistics != nullptr )
  {
    m_solidMechanicsStatistics->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }

  m_solidMechanicsStateResetTask.execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );

  // always returns false because we don't want early return (see EventManager.cpp)
  return false;
}

namespace
{
REGISTER_CATALOG_ENTRY( TaskBase, PoromechanicsInitialization, string const &, Group * const )
}

} /* namespace geos */
