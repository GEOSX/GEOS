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

#include "PhysicsSolverBase.hpp"
#include "PhysicsSolverManager.hpp"

#include "common/TimingMacros.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "mesh/DomainPartition.hpp"
#include "math/interpolation/Interpolation.hpp"
#include "common/Timer.hpp"
#include "common/Units.hpp"
#include "dataRepository/LogLevelsInfo.hpp"

#if defined(GEOS_USE_PYGEOSX)
#include "python/PySolverType.hpp"
#endif

namespace geos
{

using namespace dataRepository;

PhysicsSolverBase::PhysicsSolverBase( string const & name,
                                      Group * const parent )
  :
  ExecutableGroup( name, parent ),
  m_cflFactor(),
  m_maxStableDt{ 1e99 },
  m_nextDt( 1e99 ),
  m_numTimestepsSinceLastDtCut( -1 ),
  m_dofManager( name ),
  m_linearSolverParameters( groupKeyStruct::linearSolverParametersString(), this ),
  m_nonlinearSolverParameters( groupKeyStruct::nonlinearSolverParametersString(), this ),
  m_solverStatistics( groupKeyStruct::solverStatisticsString(), this ),
  m_systemSetupTimestamp( 0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  // This sets a flag to indicate that this object is going to select the time step size
  this->setTimesteppingBehavior( ExecutableGroup::TimesteppingBehavior::DeterminesTimeStepSize );

  registerWrapper( viewKeyStruct::cflFactorString(), &m_cflFactor ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor to apply to the `CFL condition <http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition>`_"
                    " when calculating the maximum allowable time step. Values should be in the interval (0,1] " );

  registerWrapper( viewKeyStruct::maxStableDtString(), &m_maxStableDt ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Value of the Maximum Stable Timestep for this solver." );

  this->registerWrapper( viewKeyStruct::discretizationString(), &m_discretizationName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of discretization object (defined in the :ref:`NumericalMethodsManager`) to use for this "
                    "solver. For instance, if this is a Finite Element Solver, the name of a :ref:`FiniteElement` "
                    "should be specified. If this is a Finite Volume Method, the name of a :ref:`FiniteVolume` "
                    "discretization should be specified." );

  registerWrapper( viewKeyStruct::targetRegionsString(), &m_targetRegionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Allowable regions that the solver may be applied to. Note that this does not indicate that "
                    "the solver will be applied to these regions, only that allocation will occur such that the "
                    "solver may be applied to these regions. The decision about what regions this solver will be"
                    "applied to rests in the EventManager." );

  registerWrapper( viewKeyStruct::meshTargetsString(), &m_meshTargets ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "MeshBody/Region combinations that the solver will be applied to." );

  registerWrapper( viewKeyStruct::initialDtString(), &m_nextDt ).
    setApplyDefaultValue( 1e99 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Initial time-step value required by the solver to the event manager." );

  registerWrapper( viewKeyStruct::writeLinearSystemString(), &m_writeLinearSystem ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Write matrix, rhs, solution to screen ( = 1) or file ( = 2)." );

  addLogLevel< logInfo::Fields >();
  addLogLevel< logInfo::LineSearch >();
  addLogLevel< logInfo::Solution >();
  addLogLevel< logInfo::Convergence >();
  addLogLevel< logInfo::TimeStep >();
  addLogLevel< logInfo::LinearSolver >();
  addLogLevel< logInfo::NonlinearSolver >();
  addLogLevel< logInfo::Timers >();

  registerGroup( groupKeyStruct::linearSolverParametersString(), &m_linearSolverParameters );
  registerGroup( groupKeyStruct::nonlinearSolverParametersString(), &m_nonlinearSolverParameters );
  registerGroup( groupKeyStruct::solverStatisticsString(), &m_solverStatistics );

  m_localMatrix.setName( this->getName() + "/localMatrix" );
  m_matrix.setDofManager( &m_dofManager );
}

PhysicsSolverBase::~PhysicsSolverBase() = default;

void PhysicsSolverBase::initialize_postMeshGeneration()
{
  ExecutableGroup::initialize_postMeshGeneration();
  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  generateMeshTargetsFromTargetRegions( domain.getMeshBodies());
}

void PhysicsSolverBase::generateMeshTargetsFromTargetRegions( Group const & meshBodies )
{
  for( auto const & target : m_targetRegionNames )
  {

    std::vector< string > targetTokens = stringutilities::tokenize( target, "/" );

    if( targetTokens.size()==1 ) // no MeshBody or MeshLevel specified
    {
      GEOS_ERROR_IF( meshBodies.numSubGroups() != 1,
                     getDataContext() << ": No MeshBody information is specified in" <<
                     " PhysicsSolverBase::meshTargets, but there are multiple MeshBody objects" );
      MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( 0 );
      string const meshBodyName = meshBody.getName();

      string const meshLevelName = m_discretizationName;

      string const regionName = target;
      auto const key = std::make_pair( meshBodyName, meshLevelName );
      m_meshTargets[key].emplace_back( regionName );
    }
    else if( targetTokens.size()==2 )
    {
      string const meshBodyName = targetTokens[0];
      GEOS_ERROR_IF( !meshBodies.hasGroup( meshBodyName ),
                     getWrapperDataContext( viewKeyStruct::targetRegionsString() ) << ": MeshBody (" <<
                     meshBodyName << ") is specified in targetRegions, but does not exist." );

      string const meshLevelName = m_discretizationName;

      string const regionName = targetTokens[1];


      auto const key = std::make_pair( meshBodyName, meshLevelName );
      m_meshTargets[key].emplace_back( regionName );
    }
    else
    {
      GEOS_ERROR( getDataContext() << ": Invalid specification of targetRegions" );
    }
  }
}


void PhysicsSolverBase::registerDataOnMesh( Group & meshBodies )
{
  ExecutableGroup::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      setConstitutiveNamesCallSuper( subRegion );
      setConstitutiveNames( subRegion );
    } );

  } );

}



Group * PhysicsSolverBase::createChild( string const & GEOS_UNUSED_PARAM( childKey ), string const & GEOS_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

PhysicsSolverBase::CatalogInterface::CatalogType & PhysicsSolverBase::getCatalog()
{
  static PhysicsSolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

localIndex PhysicsSolverBase::targetRegionIndex( string const & regionName ) const
{
  auto const pos = std::find( m_targetRegionNames.begin(), m_targetRegionNames.end(), regionName );
  GEOS_ERROR_IF( pos == m_targetRegionNames.end(),
                 GEOS_FMT( "{}: Region {} is not a target of the solver.",
                           getDataContext(), regionName ) );
  return std::distance( m_targetRegionNames.begin(), pos );
}

bool PhysicsSolverBase::registerCallback( void * func, const std::type_info & funcType )
{
  if( std::type_index( funcType ) == std::type_index( typeid( std::function< void( CRSMatrix< real64, globalIndex >, array1d< real64 > ) > ) ) )
  {
    m_assemblyCallback = *reinterpret_cast< std::function< void( CRSMatrix< real64, globalIndex >, array1d< real64 > ) > * >( func );
    return true;
  }

  return false;
}

real64 PhysicsSolverBase::solverStep( real64 const & time_n,
                                      real64 const & dt,
                                      const integer cycleNumber,
                                      DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // Only build the sparsity pattern if the mesh has changed
  Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

  if( meshModificationTimestamp > getSystemSetupTimestamp() )
  {
    setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );
    setSystemSetupTimestamp( meshModificationTimestamp );

    std::ostringstream oss;
    m_dofManager.printFieldInfo( oss );
    GEOS_LOG_LEVEL_INFO( logInfo::Fields, oss.str())
  }

  implicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  real64 const dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

bool PhysicsSolverBase::execute( real64 const time_n,
                                 real64 const dt,
                                 integer const cycleNumber,
                                 integer const GEOS_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  /*
   * Reset counter indicating the number of cycles since the last timestep cut
   * when the new timestep. "-1" means that no time-step cut has ocurred.
   * */
  if( dt < m_nextDt )
  {
    m_numTimestepsSinceLastDtCut = -1;
  }

  real64 dtRemaining = dt;
  real64 nextDt = dt;

  integer const maxSubSteps = m_nonlinearSolverParameters.m_maxSubSteps;

  // Keep track of substeps. It is useful to output these.
  std::vector< real64 > subStepDt( maxSubSteps, 0.0 );
  integer numOfSubSteps = 0;

  for( integer subStep = 0; subStep < maxSubSteps && dtRemaining > 0.0; ++subStep )
  {
    // reset number of nonlinear and linear iterations
    m_solverStatistics.initializeTimeStepStatistics();

    real64 const dtAccepted = solverStep( time_n + (dt - dtRemaining),
                                          nextDt,
                                          cycleNumber,
                                          domain );
    numOfSubSteps++;
    subStepDt[subStep] = dtAccepted;

    // increment the cumulative number of nonlinear and linear iterations
    m_solverStatistics.saveTimeStepStatistics();

    /*
     * Let us check convergence history of previous solve:
     * - number of nonlinear iter.
     * - if the time-step was chopped. Then we can add some heuristics to choose next dt.
     * */
    dtRemaining -= dtAccepted;

    if( dtRemaining > 0.0 )
    {
      nextDt = setNextDt( time_n + (dt - dtRemaining), dtAccepted, domain );
      if( nextDt < dtRemaining )
      {
        // better to do two equal steps than one big and one small (even potentially tiny)
        if( nextDt * 2 > dtRemaining )
        {
          nextDt = dtRemaining / 2;
          if( m_nonlinearSolverParameters.getLogLevel() > 0 )
            GEOS_LOG_RANK_0( GEOS_FMT( "{}: shortening time step to {} to cover remaining time {} in two steps", getName(), nextDt, dtRemaining ));
        }
      }
      else
      {
        nextDt = dtRemaining;
        if( m_nonlinearSolverParameters.getLogLevel() > 0 )
          GEOS_LOG_RANK_0( GEOS_FMT( "{}: shortening time step to {} to match remaining time", getName(), nextDt ));
      }
    }

    if( dtRemaining > 0.0 )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep,
                                  GEOS_FMT( "{}: sub-step = {}, accepted dt = {}, next dt = {}, remaining dt = {}", getName(), subStep, dtAccepted, nextDt, dtRemaining ) );
    }
  }

  GEOS_ERROR_IF( dtRemaining > 0.0, getDataContext() << ": Maximum allowed number of sub-steps"
                                                        " reached. Consider increasing maxSubSteps." );

  // Decide what to do with the next Dt for the event running the solver.
  m_nextDt = setNextDt( time_n + dt, nextDt, domain );

  // Increase counter to indicate how many cycles since the last timestep cut
  if( m_numTimestepsSinceLastDtCut >= 0 )
  {
    m_numTimestepsSinceLastDtCut++;
  }

  logEndOfCycleInformation( cycleNumber, numOfSubSteps, subStepDt );

  return false;
}

void PhysicsSolverBase::logEndOfCycleInformation( integer const cycleNumber,
                                                  integer const numOfSubSteps,
                                                  std::vector< real64 > const & subStepDt ) const
{
  // The formating here is a work in progress.
  GEOS_LOG_RANK_0( "\n------------------------- TIMESTEP END -------------------------" );
  GEOS_LOG_RANK_0( GEOS_FMT( "    - Cycle:      {}", cycleNumber ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "    - N substeps: {}", numOfSubSteps ) );
  std::string logMessage = "    - dt:";
  for( integer i = 0; i < numOfSubSteps; ++i )
  {
    logMessage += "  " + units::TimeFormatInfo::fromSeconds( subStepDt[i] ).toString();
  }
  // Log the complete message once
  GEOS_LOG_RANK_0( logMessage );
  GEOS_LOG_RANK_0( "------------------------------------------------------------------\n" );
}

real64 PhysicsSolverBase::setNextDt( real64 const & GEOS_UNUSED_PARAM( currentTime ),
                                     real64 const & lastDt,
                                     DomainPartition & domain )
{
  integer const minTimeStepIncreaseInterval = m_nonlinearSolverParameters.minTimeStepIncreaseInterval();
  real64 const nextDtNewton = setNextDtBasedOnNewtonIter( lastDt );
  if( m_nonlinearSolverParameters.getLogLevel() > 0 )
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on Newton iterations = {}", getName(), nextDtNewton ));
  real64 const nextDtStateChange = setNextDtBasedOnStateChange( lastDt, domain );
  if( m_nonlinearSolverParameters.getLogLevel() > 0 )
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on state change = {}", getName(), nextDtStateChange ));

  if( ( m_numTimestepsSinceLastDtCut >= 0 ) && ( m_numTimestepsSinceLastDtCut < minTimeStepIncreaseInterval ) )
  {
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: time-step size will be kept the same since it's been {} cycles since last cut.",
                                                             getName(), m_numTimestepsSinceLastDtCut ) );
    return lastDt;
  }

  if( nextDtNewton < nextDtStateChange )      // time step size decided based on convergence
  {
    if( nextDtNewton > lastDt )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: time-step required will be increased based on number of iterations.",
                                                               getName() ) );
    }
    else if( nextDtNewton < lastDt )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: time-step required will be decreased based on number of iterations.",
                                                               getName() ) );
    }
    else
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: time-step required will be kept the same based on number of iterations.",
                                                               getName() ) );
    }
  }
  else         // time step size decided based on state change
  {
    if( nextDtStateChange > lastDt )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: time-step required will be increased based on state change.",
                                                               getName()));
    }
    else if( nextDtStateChange < lastDt )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: time-step required will be decreased based on state change.",
                                                               getName()));
    }
    else
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: time-step required will be kept the same based on state change.",
                                                               getName()));
    }
  }

  return std::min( nextDtNewton, nextDtStateChange );
}

real64 PhysicsSolverBase::setNextDtBasedOnStateChange( real64 const & lastDt,
                                                       DomainPartition & domain )
{
  GEOS_UNUSED_VAR( lastDt, domain );
  return LvArray::NumericLimits< real64 >::max; // i.e., not implemented
}

real64 PhysicsSolverBase::setNextDtBasedOnNewtonIter( real64 const & lastDt )
{
  integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
  integer const iterDecreaseLimit = m_nonlinearSolverParameters.timeStepDecreaseIterLimit();
  integer const iterIncreaseLimit = m_nonlinearSolverParameters.timeStepIncreaseIterLimit();

  real64 nextDt = 0;
  if( newtonIter < iterIncreaseLimit )
  {
    // Easy convergence, let's increase the time-step.
    nextDt = lastDt * m_nonlinearSolverParameters.timeStepIncreaseFactor();
    if( m_nonlinearSolverParameters.getLogLevel() > 0 )
      GEOS_LOG_RANK_0( GEOS_FMT( "{}: number of iterations = {} is less than {}, next time step = {} (increase)", getName(), newtonIter, iterIncreaseLimit, nextDt ));
  }
  else if( newtonIter > iterDecreaseLimit )
  {
    // Tough convergence let us make the time-step smaller!
    nextDt = lastDt * m_nonlinearSolverParameters.timeStepDecreaseFactor();
    if( m_nonlinearSolverParameters.getLogLevel() > 0 )
      GEOS_LOG_RANK_0( GEOS_FMT( "{}: number of iterations = {} is more than {}, next time step = {} (decrease)", getName(), newtonIter, iterDecreaseLimit, nextDt ));
  }
  else
  {
    nextDt = lastDt;
    if( m_nonlinearSolverParameters.getLogLevel() > 0 )
      GEOS_LOG_RANK_0( GEOS_FMT( "{}: number of iterations = {} is between {} and {}, next time step = {} (no change)", getName(), newtonIter, iterIncreaseLimit, iterDecreaseLimit, nextDt ));
  }
  return nextDt;
}

real64 PhysicsSolverBase::linearImplicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                              DomainPartition & domain )
{
  // call setup for physics solver. Pre step allocations etc.
  // TODO: Nonlinear step does not call its own setup, need to decide on consistent behavior
  implicitStepSetup( time_n, dt, domain );

  {
    Timer timer( m_timers["assemble"] );

    // zero out matrix/rhs before assembly
    m_localMatrix.zero();
    m_rhs.zero();

    arrayView1d< real64 > const localRhs = m_rhs.open();

    // call assemble to fill the matrix and the rhs
    assembleSystem( time_n,
                    dt,
                    domain,
                    m_dofManager,
                    m_localMatrix.toViewConstSizes(),
                    localRhs );

    // apply boundary conditions to system
    applyBoundaryConditions( time_n,
                             dt,
                             domain,
                             m_dofManager,
                             m_localMatrix.toViewConstSizes(),
                             localRhs );

    m_rhs.close();

    if( m_assemblyCallback )
    {
      // Make a copy of LA objects and ship off to the callback
      array1d< real64 > localRhsCopy( m_rhs.localSize() );
      localRhsCopy.setValues< parallelDevicePolicy<> >( m_rhs.values() );
      m_assemblyCallback( m_localMatrix, std::move( localRhsCopy ) );
    }
  }

  {
    Timer timer( m_timers["linear solver total"] );

    // TODO: Trilinos currently requires this, re-evaluate after moving to Tpetra-based solvers
    if( m_precond )
    {
      m_precond->clear();
    }

    {
      Timer timer_create( m_timers["linear solver create"] );

      // Compose parallel LA matrix out of local matrix
      m_matrix.create( m_localMatrix.toViewConst(), m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
    }

    // Output the linear system matrix/rhs for debugging purposes
    debugOutputSystem( 0.0, 0, 0, m_matrix, m_rhs );

    // Solve the linear system
    solveLinearSystem( m_dofManager, m_matrix, m_rhs, m_solution );
  }

  // Increment the solver statistics for reporting purposes
  m_solverStatistics.logNonlinearIteration( m_linearSolverResult.numIterations );

  // Output the linear system solution for debugging purposes
  debugOutputSolution( 0.0, 0, 0, m_solution );

  {
    Timer timer( m_timers["apply solution"] );

    // apply the system solution to the fields/variables
    applySystemSolution( m_dofManager, m_solution.values(), 1.0, dt, domain );
  }

  {
    Timer timer( m_timers["update state"] );

    // update non-primary variables (constitutive models)
    updateState( domain );
  }

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dt, domain );

  // return the achieved timestep
  return dt;
}


bool PhysicsSolverBase::lineSearch( real64 const & time_n,
                                    real64 const & dt,
                                    integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                    DomainPartition & domain,
                                    DofManager const & dofManager,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    ParallelVector & rhs,
                                    ParallelVector & solution,
                                    real64 const scaleFactor,
                                    real64 & lastResidual )
{
  Timer timer( m_timers["line search"] );

  integer const maxNumberLineSearchCuts = m_nonlinearSolverParameters.m_lineSearchMaxCuts;
  real64 const lineSearchCutFactor = m_nonlinearSolverParameters.m_lineSearchCutFactor;

  // flag to determine if we should solve the system and apply the solution. If the line
  // search fails we just bail.
  bool lineSearchSuccess = false;

  real64 residualNorm = lastResidual;

  // scale factor is value applied to the previous solution. In this case we want to
  // subtract a portion of the previous solution.
  real64 localScaleFactor = -scaleFactor;
  real64 cumulativeScale = scaleFactor;

  // main loop for the line search.
  for( integer lineSearchIteration = 0; lineSearchIteration < maxNumberLineSearchCuts; ++lineSearchIteration )
  {
    // cut the scale factor by half. This means that the scale factors will
    // have values of -0.5, -0.25, -0.125, ...
    localScaleFactor *= lineSearchCutFactor;
    cumulativeScale += localScaleFactor;

    if( !checkSystemSolution( domain, dofManager, solution.values(), localScaleFactor ) )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, GEOS_FMT( "        Line search {}, solution check failed", lineSearchIteration ) );
      continue;
    }

    applySystemSolution( dofManager, solution.values(), localScaleFactor, dt, domain );

    // update non-primary variables (constitutive models)
    updateState( domain );

    // re-assemble system
    localMatrix.zero();
    rhs.zero();

    arrayView1d< real64 > const localRhs = rhs.open();
    assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    rhs.close();

    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, GEOS_FMT( "        Line search @ {:0.3f}:      ", cumulativeScale ));

    // get residual norm
    residualNorm = calculateResidualNorm( time_n, dt, domain, dofManager, rhs.values() );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );

    // if the residual norm is less than the last residual, we can proceed to the
    // solution step
    if( residualNorm < lastResidual )
    {
      lineSearchSuccess = true;
      break;
    }
  }

  lastResidual = residualNorm;
  return lineSearchSuccess;
}

bool PhysicsSolverBase::lineSearchWithParabolicInterpolation( real64 const & time_n,
                                                              real64 const & dt,
                                                              integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                                              DomainPartition & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              ParallelVector & rhs,
                                                              ParallelVector & solution,
                                                              real64 const scaleFactor,
                                                              real64 & lastResidual,
                                                              real64 & residualNormT )
{
  Timer timer( m_timers["line search"] );

  bool lineSearchSuccess = true;

  integer const maxNumberLineSearchCuts = m_nonlinearSolverParameters.m_lineSearchMaxCuts;

  real64 const sigma1 = 0.5;
  real64 const alpha = 1.e-4;

  real64 localScaleFactor = scaleFactor;
  real64 lamm = scaleFactor;
  real64 lamc = localScaleFactor;
  integer lineSearchIteration = 0;

  real64 residualNorm0 = lastResidual;

  real64 ff0 = residualNorm0*residualNorm0;
  real64 ffT = residualNormT*residualNormT;
  real64 ffm = ffT;
  real64 cumulativeScale = scaleFactor;

  while( residualNormT >= (1.0 - alpha*localScaleFactor)*residualNorm0 )
  {
    real64 const previousLocalScaleFactor = localScaleFactor;
    // Apply the three point parabolic model
    if( lineSearchIteration == 0 )
    {
      localScaleFactor *= sigma1;
    }
    else
    {
      localScaleFactor = interpolation::parabolicInterpolationThreePoints( lamc, lamm, ff0, ffT, ffm );
    }

    // Update x; keep the books on lambda
    real64 const deltaLocalScaleFactor = ( localScaleFactor - previousLocalScaleFactor );
    cumulativeScale += deltaLocalScaleFactor;

    if( !checkSystemSolution( domain, dofManager, solution.values(), deltaLocalScaleFactor ) )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, GEOS_FMT( "        Line search {}, solution check failed", lineSearchIteration ) );
      continue;
    }

    applySystemSolution( dofManager, solution.values(), deltaLocalScaleFactor, dt, domain );

    updateState( domain );

    lamm = lamc;
    lamc = localScaleFactor;

    // Keep the books on the function norms

    // re-assemble system
    // TODO: add a flag to avoid a completely useless Jacobian computation: rhs is enough
    localMatrix.zero();
    rhs.zero();

    arrayView1d< real64 > const localRhs = rhs.open();
    assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    rhs.close();

    if( logger::internal::rank==0 )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, GEOS_FMT( "        Line search @ {:0.3f}:      ", cumulativeScale ) );
    }

    // get residual norm
    residualNormT = calculateResidualNorm( time_n, dt, domain, dofManager, rhs.values() );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNormT ) );


    ffm = ffT;
    ffT = residualNormT*residualNormT;
    lineSearchIteration += 1;

    if( lineSearchIteration > maxNumberLineSearchCuts )
    {
      lineSearchSuccess = false;
      break;
    }
  }

  lastResidual = residualNormT;

  return lineSearchSuccess;
}


real64 PhysicsSolverBase::eisenstatWalker( real64 const newNewtonNorm,
                                           real64 const oldNewtonNorm,
                                           LinearSolverParameters::Krylov const & krylovParams )
{
  real64 normRatio = std::min( newNewtonNorm / oldNewtonNorm, 1.0 );
  real64 newKrylovTol = krylovParams.adaptiveGamma * std::pow( normRatio, krylovParams.adaptiveExponent );
  // the following is a safeguard to avoid too sharp tolerance reduction
  // the bound is the quadratic reduction wrt previous value
  real64 altKrylovTol = std::pow( krylovParams.relTolerance, 2.0 );

  real64 krylovTol = std::max( newKrylovTol, altKrylovTol );
  krylovTol = std::min( krylovTol, krylovParams.weakestTol );
  krylovTol = std::max( krylovTol, krylovParams.strongestTol );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LinearSolver,
                              GEOS_FMT( "        Adaptive linear tolerance = {:4.2e} (norm ratio = {:4.2e}, old tolerance = {:4.2e}, new tolerance = {:4.2e}, safeguard = {:4.2e})",
                                        krylovTol, normRatio, krylovParams.relTolerance, newKrylovTol, altKrylovTol ) );

  return krylovTol;
}

real64 PhysicsSolverBase::nonlinearImplicitStep( real64 const & time_n,
                                                 real64 const & dt,
                                                 integer const cycleNumber,
                                                 DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  // dt may be cut during the course of this step, so we are keeping a local
  // value to track the achieved dt for this step.
  real64 stepDt = dt;

  integer const maxNumberDtCuts = m_nonlinearSolverParameters.m_maxTimeStepCuts;
  real64 const dtCutFactor = m_nonlinearSolverParameters.m_timeStepCutFactor;

  bool const allowNonConverged = m_nonlinearSolverParameters.m_allowNonConverged > 0;

  integer & dtAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;

  integer const & maxConfigurationIter = m_nonlinearSolverParameters.m_maxNumConfigurationAttempts;

  integer & configurationLoopIter = m_nonlinearSolverParameters.m_numConfigurationAttempts;

  bool isConfigurationLoopConverged = false;

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      resetStateToBeginningOfStep( domain );
      resetConfigurationToBeginningOfStep( domain );
    }

    // it's the simplest configuration that can be attempted whenever Newton's fails as a last resource.
    bool attemptedSimplestConfiguration = false;

    // Configuration loop
    for( configurationLoopIter = 0; configurationLoopIter < maxConfigurationIter; ++configurationLoopIter )
    {

      outputConfigurationStatistics( domain );

      bool const isNewtonConverged = solveNonlinearSystem( time_n,
                                                           stepDt,
                                                           cycleNumber,
                                                           domain );

      if( isNewtonConverged )
      {
        isConfigurationLoopConverged = updateConfiguration( domain );

        if( isConfigurationLoopConverged )
        {
          break; // get out of configuration loop coz everything converged.
        }
        else
        {
          // increment the solver statistics for reporting purposes
          m_solverStatistics.logOuterLoopIteration();
          GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::NonlinearSolver, "---------- Configuration did not converge. Testing new configuration. ----------" );
        }
      }
      else if( !attemptedSimplestConfiguration )
      {
        resetStateToBeginningOfStep( domain );
        bool const breakLoop = resetConfigurationToDefault( domain );
        attemptedSimplestConfiguration = true;
        if( breakLoop )
        {
          break;
        }
      }
      else
      {
        break; // get out of configuration loop and cut the time-step if you can.
      }

    }  // end of configuration loop

    if( isConfigurationLoopConverged )
    {
      // get out of outer loop
      break;
    }
    else
    {
      // cut timestep, go back to beginning of step and restart the Newton loop
      stepDt *= dtCutFactor;
      m_numTimestepsSinceLastDtCut = 0;
      GEOS_LOG_LEVEL_INFO_RANK_0 ( logInfo::TimeStep, GEOS_FMT( "New dt = {}", stepDt ) );

      // notify the solver statistics counter that this is a time step cut
      m_solverStatistics.logTimeStepCut();
    }
  } // end of outer loop (dt chopping strategy)

  if( !isConfigurationLoopConverged )
  {
    GEOS_LOG_RANK_0( "Convergence not achieved." );

    if( allowNonConverged )
    {
      GEOS_LOG_RANK_0( "The accepted solution may be inaccurate." );
    }
    else
    {
      GEOS_ERROR( "Nonconverged solutions not allowed. Terminating..." );
    }
  }

  // return the achieved timestep
  return stepDt;
}

bool PhysicsSolverBase::solveNonlinearSystem( real64 const & time_n,
                                              real64 const & stepDt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{
  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer & dtAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;
  integer & configurationLoopIter = m_nonlinearSolverParameters.m_numConfigurationAttempts;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

  // keep residual from previous iteration in case we need to do a line search
  real64 lastResidual = 1e99;
  integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
  real64 scaleFactor = 1.0;

  bool isNewtonConverged = false;

  for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
  {

    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::NonlinearSolver,
                                GEOS_FMT( "    Attempt: {:2}, ConfigurationIter: {:2}, NewtonIter: {:2}", dtAttempt, configurationLoopIter, newtonIter ) );

    {
      Timer timer( m_timers["assemble"] );

      // We sync the nonlinear convergence history. The coupled solver parameters are the one being
      // used. We want to propagate the info to subsolvers. It can be important for solvers that
      // have special treatment for specific iterations.
      synchronizeNonlinearSolverParameters();

      // zero out matrix/rhs before assembly
      m_localMatrix.zero();
      m_rhs.zero();

      arrayView1d< real64 > const localRhs = m_rhs.open();

      // call assemble to fill the matrix and the rhs
      assembleSystem( time_n,
                      stepDt,
                      domain,
                      m_dofManager,
                      m_localMatrix.toViewConstSizes(),
                      localRhs );

      // apply boundary conditions to system
      applyBoundaryConditions( time_n,
                               stepDt,
                               domain,
                               m_dofManager,
                               m_localMatrix.toViewConstSizes(),
                               localRhs );

      m_rhs.close();

      if( m_assemblyCallback )
      {
        // Make a copy of LA objects and ship off to the callback
        array1d< real64 > localRhsCopy( m_rhs.localSize() );
        localRhsCopy.setValues< parallelDevicePolicy<> >( m_rhs.values() );
        m_assemblyCallback( m_localMatrix, std::move( localRhsCopy ) );
      }
    }

    real64 residualNorm = 0;
    {
      Timer timer( m_timers["convergence check"] );

      // get residual norm
      residualNorm = calculateResidualNorm( time_n, stepDt, domain, m_dofManager, m_rhs.values() );
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Convergence, GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );
    }

    // if the residual norm is less than the Newton tolerance we denote that we have
    // converged and break from the Newton loop immediately.
    if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
    {
      isNewtonConverged = true;
      break;
    }

    // if the residual norm is above the max allowed residual norm, we break from
    // the Newton loop to avoid crashes due to Newton divergence
    if( residualNorm > m_nonlinearSolverParameters.m_maxAllowedResidualNorm )
    {
      string const maxAllowedResidualNormString = NonlinearSolverParameters::viewKeysStruct::maxAllowedResidualNormString();
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Convergence, GEOS_FMT( "    The residual norm is above the {} of {}. Newton loop terminated.",
                                                                  maxAllowedResidualNormString,
                                                                  m_nonlinearSolverParameters.m_maxAllowedResidualNorm ) );
      isNewtonConverged = false;
      break;
    }

    // do line search in case residual has increased
    if( m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None
        && residualNorm > lastResidual * m_nonlinearSolverParameters.m_lineSearchResidualFactor
        && newtonIter >= m_nonlinearSolverParameters.m_lineSearchStartingIteration )
    {
      bool lineSearchSuccess = false;
      if( m_nonlinearSolverParameters.m_lineSearchInterpType == NonlinearSolverParameters::LineSearchInterpolationType::Linear )
      {
        residualNorm = lastResidual;
        lineSearchSuccess = lineSearch( time_n,
                                        stepDt,
                                        cycleNumber,
                                        domain,
                                        m_dofManager,
                                        m_localMatrix.toViewConstSizes(),
                                        m_rhs,
                                        m_solution,
                                        scaleFactor,
                                        residualNorm );
      }
      else
      {
        lineSearchSuccess = lineSearchWithParabolicInterpolation( time_n,
                                                                  stepDt,
                                                                  cycleNumber,
                                                                  domain,
                                                                  m_dofManager,
                                                                  m_localMatrix.toViewConstSizes(),
                                                                  m_rhs,
                                                                  m_solution,
                                                                  scaleFactor,
                                                                  lastResidual,
                                                                  residualNorm );
      }

      if( !lineSearchSuccess )
      {
        if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Attempt )
        {
          GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, "        Line search failed to produce reduced residual. Accepting iteration." );
        }
        else if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Require )
        {
          // if line search failed, then break out of the main Newton loop. Timestep will be cut.
          GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LineSearch, "        Line search failed to produce reduced residual. Exiting Newton Loop." );
          break;
        }
      }
    }

    {
      Timer timer( m_timers["linear solver total"] );

      // if using adaptive Krylov tolerance scheme, update tolerance.
      LinearSolverParameters::Krylov & krylovParams = m_linearSolverParameters.get().krylov;
      if( krylovParams.useAdaptiveTol )
      {
        krylovParams.relTolerance = newtonIter > 0 ? eisenstatWalker( residualNorm, lastResidual, krylovParams ) : krylovParams.weakestTol;
      }

      // TODO: Trilinos currently requires this, re-evaluate after moving to Tpetra-based solvers
      if( m_precond )
      {
        m_precond->clear();
      }

      {
        Timer timer_setup( m_timers["linear solver create"] );

        // Compose parallel LA matrix/rhs out of local LA matrix/rhs
        //
        m_matrix.create( m_localMatrix.toViewConst(), m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
      }

      // Output the linear system matrix/rhs for debugging purposes
      debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs );

      // Solve the linear system
      solveLinearSystem( m_dofManager, m_matrix, m_rhs, m_solution );

      // Increment the solver statistics for reporting purposes
      m_solverStatistics.logNonlinearIteration( m_linearSolverResult.numIterations );

      // Output the linear system solution for debugging purposes
      debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );
    }

    {
      Timer timer( m_timers["apply solution"] );

      // Compute the scaling factor for the Newton update
      scaleFactor = scalingForSystemSolution( domain, m_dofManager, m_solution.values() );

      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Global solution scaling factor = {}", getName(), scaleFactor ) );

      if( !checkSystemSolution( domain, m_dofManager, m_solution.values(), scaleFactor ) )
      {
        // TODO try chopping (similar to line search)
        GEOS_LOG_RANK_0( GEOS_FMT( "    {}: Solution check failed. Newton loop terminated.", getName()) );
        break;
      }

      // apply the system solution to the fields/variables
      applySystemSolution( m_dofManager, m_solution.values(), scaleFactor, stepDt, domain );
    }

    {
      Timer timer( m_timers["update state"] );

      // update non-primary variables (constitutive models)
      updateState( domain );
    }

    lastResidual = residualNorm;
  }

  return isNewtonConverged;
}

real64 PhysicsSolverBase::explicitStep( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                        real64 const & GEOS_UNUSED_PARAM( dt ),
                                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                        DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_THROW( "PhysicsSolverBase::ExplicitStep called!. Should be overridden.", std::runtime_error );
  return 0;
}

void PhysicsSolverBase::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                           real64 const & GEOS_UNUSED_PARAM( dt ),
                                           DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_THROW( "PhysicsSolverBase::ImplicitStepSetup called!. Should be overridden.", std::runtime_error );
}

void PhysicsSolverBase::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                   DofManager & GEOS_UNUSED_PARAM( dofManager ) ) const
{
  GEOS_ERROR( "PhysicsSolverBase::setupDofs called!. Should be overridden." );
}

void PhysicsSolverBase::setupSystem( DomainPartition & domain,
                                     DofManager & dofManager,
                                     CRSMatrix< real64, globalIndex > & localMatrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution,
                                     bool const setSparsity )
{
  GEOS_MARK_FUNCTION;

  dofManager.setDomain( domain );

  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  if( setSparsity )
  {
    SparsityPattern< globalIndex > pattern;
    dofManager.setSparsityPattern( pattern );
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  }
  localMatrix.setName( this->getName() + "/matrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );
}

void PhysicsSolverBase::assembleSystem( real64 const GEOS_UNUSED_PARAM( time ),
                                        real64 const GEOS_UNUSED_PARAM( dt ),
                                        DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                        DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                        CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                        arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_ERROR( "PhysicsSolverBase::Assemble called!. Should be overridden." );
}

void PhysicsSolverBase::applyBoundaryConditions( real64 const GEOS_UNUSED_PARAM( time ),
                                                 real64 const GEOS_UNUSED_PARAM( dt ),
                                                 DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                 DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                 CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                 arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_ERROR( "PhysicsSolverBase::applyBoundaryConditions called!. Should be overridden." );
}

namespace
{

/**
 * @brief Helper for debug output of linear algebra objects (matrices and vectors)
 * @tparam T type of LA object (must have stream insertion and .write() implemented)
 * @param obj                the object to output
 * @param cycleNumber        event cycle number
 * @param nonlinearIteration nonlinear iteration number
 * @param filePrefix          short filename prefix (e.g. "mat")
 * @param screenName           long name for screen output (e.g. "System matrix")
 * @param toScreen           whether to print on screen
 * @param toFile             whether to write to file
 */
template< typename T >
void debugOutputLAObject( T const & obj,
                          real64 const & GEOS_UNUSED_PARAM( time ),
                          integer const cycleNumber,
                          integer const nonlinearIteration,
                          string const & filePrefix,
                          string const & screenName,
                          bool const toScreen,
                          bool const toFile )
{
  if( toScreen )
  {
    string const frame( screenName.size() + 1, '=' );
    GEOS_LOG_RANK_0( frame << "\n" << screenName << ":\n" << frame );
    GEOS_LOG( obj );
  }

  if( toFile )
  {
    string const filename = GEOS_FMT( "{}_{:06}_{:02}.mtx", filePrefix.c_str(), cycleNumber, nonlinearIteration );
    obj.write( filename, LAIOutputFormat::NATIVE_ASCII );
    GEOS_LOG_RANK_0( screenName << " written to " << filename );
  }
}

}

void PhysicsSolverBase::debugOutputSystem( real64 const & time,
                                           integer const cycleNumber,
                                           integer const nonlinearIteration,
                                           ParallelMatrix const & matrix,
                                           ParallelVector const & rhs ) const
{
  // special case when flag value > 2
  if( m_writeLinearSystem > 2 && cycleNumber < m_writeLinearSystem )
    return;

  debugOutputLAObject( matrix,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_mat",
                       "System matrix",
                       m_writeLinearSystem == 1,
                       m_writeLinearSystem >= 2 );

  debugOutputLAObject( rhs,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_rhs",
                       "System right-hand side",
                       m_writeLinearSystem == 1,
                       m_writeLinearSystem >= 2 );
}

void PhysicsSolverBase::debugOutputSolution( real64 const & time,
                                             integer const cycleNumber,
                                             integer const nonlinearIteration,
                                             ParallelVector const & solution ) const
{
  // special case when flag value > 2
  if( m_writeLinearSystem > 2 && cycleNumber < m_writeLinearSystem )
    return;

  debugOutputLAObject( solution,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_sol",
                       "System solution",
                       m_writeLinearSystem == 1,
                       m_writeLinearSystem >= 2 );
}

real64
PhysicsSolverBase::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time ),
                                          real64 const & GEOS_UNUSED_PARAM( dt ),
                                          DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                          DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                          arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_ERROR( "PhysicsSolverBase::calculateResidualNorm called!. Should be overridden." );
  return 0;
}

void PhysicsSolverBase::solveLinearSystem( DofManager const & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  GEOS_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  LinearSolverParameters const & params = m_linearSolverParameters.get();
  matrix.setDofManager( &dofManager );

  if( params.solverType == LinearSolverParameters::SolverType::direct || !m_precond )
  {
    std::unique_ptr< LinearSolverBase< LAInterface > > solver = LAInterface::createSolver( params );
    {
      Timer timer_setup( m_timers["linear solver setup"] );
      solver->setup( matrix );
    }
    {
      Timer timer_setup( m_timers["linear solver solve"] );
      solver->solve( rhs, solution );
    }
    m_linearSolverResult = solver->result();
  }
  else
  {
    {
      Timer timer_setup( m_timers["linear solver setup"] );
      m_precond->setup( matrix );
    }
    std::unique_ptr< KrylovSolver< ParallelVector > > solver = KrylovSolver< ParallelVector >::create( params, matrix, *m_precond );
    {
      Timer timer_setup( m_timers["linear solver solve"] );
      solver->solve( rhs, solution );
    }
    m_linearSolverResult = solver->result();
  }

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::LinearSolver, GEOS_FMT( "        Last LinSolve(iter,res) = ( {:3}, {:4.2e} )",
                                                               m_linearSolverResult.numIterations,
                                                               m_linearSolverResult.residualReduction ) );

  if( params.stopIfError )
  {
    GEOS_ERROR_IF( m_linearSolverResult.breakdown(), getDataContext() << ": Linear solution breakdown -> simulation STOP" );
  }
  else
  {
    GEOS_WARNING_IF( !m_linearSolverResult.success(), getDataContext() << ": Linear solution failed" );
  }
}

bool PhysicsSolverBase::checkSystemSolution( DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                             DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                             arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ),
                                             real64 const GEOS_UNUSED_PARAM( scalingFactor ) )
{
  return true;
}

real64 PhysicsSolverBase::scalingForSystemSolution( DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                    DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                    arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ) )
{
  return 1.0;
}

void PhysicsSolverBase::applySystemSolution( DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                             arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( localSolution ),
                                             real64 const GEOS_UNUSED_PARAM( scalingFactor ),
                                             real64 const GEOS_UNUSED_PARAM( dt ),
                                             DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_ERROR( "PhysicsSolverBase::applySystemSolution called!. Should be overridden." );
}

void PhysicsSolverBase::updateState( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_ERROR( "PhysicsSolverBase::updateState called!. Should be overridden." );
}

bool PhysicsSolverBase::updateConfiguration( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  return true;
}

void PhysicsSolverBase::outputConfigurationStatistics( DomainPartition const & GEOS_UNUSED_PARAM( domain ) ) const
{
  // For most solvers there is nothing to do.
}

void PhysicsSolverBase::resetConfigurationToBeginningOfStep( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // For most solvers there is nothing to do.
}

void PhysicsSolverBase::resetStateToBeginningOfStep( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_ERROR( "PhysicsSolverBase::ResetStateToBeginningOfStep called!. Should be overridden." );
}

bool PhysicsSolverBase::resetConfigurationToDefault( DomainPartition & GEOS_UNUSED_PARAM( domain ) ) const
{
  // for most solvers it just breaks the loop.
  return true;
}

void PhysicsSolverBase::implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time ),
                                              real64 const & GEOS_UNUSED_PARAM( dt ),
                                              DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_ERROR( "PhysicsSolverBase::ImplicitStepComplete called!. Should be overridden." );
}

void PhysicsSolverBase::cleanup( real64 const GEOS_UNUSED_PARAM( time_n ),
                                 integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                 integer const GEOS_UNUSED_PARAM( eventCounter ),
                                 real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                 DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  m_solverStatistics.outputStatistics();

  for( auto & timer : m_timers )
  {
    real64 const time = std::chrono::duration< double >( timer.second ).count();
    real64 const minTime = MpiWrapper::min( time );
    real64 const maxTime = MpiWrapper::max( time );
    if( maxTime > 0 )
    {
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Timers, GEOS_FMT( "{}: {} time = {} s (min), {} s (max)", getName(), timer.first, minTime, maxTime ) );
    }
  }

}

Timestamp PhysicsSolverBase::getMeshModificationTimestamp( DomainPartition & domain ) const
{
  Timestamp meshModificationTimestamp = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    if( meshModificationTimestamp < mesh.getModificationTimestamp() )
    {
      meshModificationTimestamp = mesh.getModificationTimestamp();
    }
  } );
  return meshModificationTimestamp;
}

R1Tensor const PhysicsSolverBase::gravityVector() const
{
  R1Tensor rval;
  if( dynamicCast< PhysicsSolverManager const * >( &getParent() ) != nullptr )
  {
    rval = getParent().getReference< R1Tensor >( PhysicsSolverManager::viewKeyStruct::gravityVectorString() );
  }
  else
  {
    rval = {0.0, 0.0, -9.81};
  }
  return rval;
}

bool PhysicsSolverBase::checkSequentialSolutionIncrements( DomainPartition & GEOS_UNUSED_PARAM( domain ) ) const
{
  // default behavior - assume converged
  return true;
}

void PhysicsSolverBase::saveSequentialIterationState( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // up to specific solver to save what is needed
  GEOS_ERROR( "Call to PhysicsSolverBase::saveSequentialIterationState. Method should be overloaded by the solver" );
}

#if defined(GEOS_USE_PYGEOSX)
PyTypeObject * PhysicsSolverBase::getPythonType() const
{ return python::getPySolverType(); }
#endif

} // namespace geos
