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
 * @file CoupledSolver.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

#include <tuple>

namespace geos
{

template< typename ... SOLVERS >
class CoupledSolver : public SolverBase
{

public:

  /**
   * @brief main constructor for CoupledSolver Objects
   * @param name the name of this instantiation of CoupledSolver in the repository
   * @param parent the parent group of this instantiation of CoupledSolver
   */
  CoupledSolver( const string & name,
                 Group * const parent )
    : SolverBase( name, parent )
  {
    forEachArgInTuple( m_solvers, [&]( auto solver, auto idx )
    {
      using SolverType = TYPEOFPTR( solver );
      string const key = SolverType::coupledSolverAttributePrefix() + "SolverName";
      registerWrapper( key, &m_names[idx()] ).
        setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
        setInputFlag( dataRepository::InputFlags::REQUIRED ).
        setDescription( "Name of the " + SolverType::coupledSolverAttributePrefix() + " solver used by the coupled solver" );
    } );

    this->getWrapper< string >( SolverBase::viewKeyStruct::discretizationString() ).
      setInputFlag( dataRepository::InputFlags::FALSE );
  }

  /// deleted copy constructor
  CoupledSolver( CoupledSolver const & ) = delete;

  /// default move constructor
  CoupledSolver( CoupledSolver && ) = default;

  /// deleted assignment operator
  CoupledSolver & operator=( CoupledSolver const & ) = delete;

  /// deleted move operator
  CoupledSolver & operator=( CoupledSolver && ) = delete;


  /**
   * @brief Utility function to set the subsolvers pointers using the names provided by the user
   */
  void
  setSubSolvers()
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
    {
      using SolverPtr = TYPEOFREF( solver );
      using SolverType = TYPEOFPTR( SolverPtr {} );
      auto const & solverName = m_names[idx()];
      auto const & solverType = LvArray::system::demangleType< SolverType >();
      solver = this->getParent().template getGroupPointer< SolverType >( solverName );
      GEOS_THROW_IF( solver == nullptr,
                     GEOS_FMT( "{}: Could not find solver '{}' of type {}",
                               getDataContext(),
                               solverName, solverType ),
                     InputError );
      GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: found {} solver named {}", getName(), solver->getCatalogName(), solverName ) );
    } );
  }


  /**
   * @brief Utility function to set the coupling between degrees of freedom
   * @param[in] domain the domain partition
   * @param[in] dofManager the dof manager
   */
  virtual void
  setupCoupling( DomainPartition const & domain,
                 DofManager & dofManager ) const
  { GEOS_UNUSED_VAR( domain, dofManager ); }

  /**
   * @brief Utility function to compute coupling terms
   * @param[in] time_n the time at the beginning of the time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   * @paran[in] dofManager the degree of freedom manager
   * @param[in] localMatrix the local matrix
   * @param[in] localRhs the local rhs
   */
  virtual void
  assembleCouplingTerms( real64 const time_n,
                         real64 const dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs )
  { GEOS_UNUSED_VAR( time_n, dt, domain, dofManager, localMatrix, localRhs ); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->setupDofs( domain, dofManager );
    } );

    setupCoupling( domain, dofManager );
  }

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->implicitStepSetup( time_n, dt, domain );
    } );
  }

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->implicitStepComplete( time_n, dt, domain );
    } );
  }

  // general version of assembleSystem function, keep in mind many solvers will override it
  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override
  {
    /// Fully-coupled assembly.

    // 1. Assemble matrix blocks of each individual solver
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    } );

    // 2. Assemble coupling blocks
    assembleCouplingTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );
    } );
  }

  virtual void
  updateState( DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->updateState( domain );
    } );
  }

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->resetStateToBeginningOfStep( domain );
    } );
  }

  /// This method is meant to be kept final. Derived CoupledSolvers are expected, if needed,
  /// to override fullyCoupledSolverStep and/or sequentiallyCoupledSolverStep.
  real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override final
  {
    GEOS_MARK_FUNCTION;

    if( getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::FullyImplicit )
    {
      return fullyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }
    else if( getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential )
    {
      return sequentiallyCoupledSolverStep( time_n, dt, cycleNumber, domain );
    }
    else
    {
      GEOS_ERROR( getDataContext() << ": Invalid coupling type option." );
      return 0;
    }

  }


  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override
  {
    real64 norm = 0.0;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsNorm = solver->calculateResidualNorm( time_n, dt, domain, dofManager, localRhs );
      norm += singlePhysicsNorm * singlePhysicsNorm;
    } );

    return sqrt( norm );
  }

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    } );
  }

  virtual bool
  checkSystemSolution( DomainPartition & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override
  {
    bool validSolution = true;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      bool const validSinglePhysicsSolution = solver->checkSystemSolution( domain, dofManager, localSolution, scalingFactor );
      validSolution = validSolution && validSinglePhysicsSolution;
    } );
    return validSolution;
  }

  virtual real64
  scalingForSystemSolution( DomainPartition & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override
  {
    real64 scalingFactor = 1e9;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsScalingFactor = solver->scalingForSystemSolution( domain, dofManager, localSolution );
      scalingFactor = LvArray::math::min( scalingFactor, singlePhysicsScalingFactor );
    } );
    return scalingFactor;
  }

  virtual real64
  setNextDtBasedOnStateChange( real64 const & currentDt,
                               DomainPartition & domain ) override
  {
    real64 nextDt = LvArray::NumericLimits< real64 >::max;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsNextDt =
        solver->setNextDtBasedOnStateChange( currentDt, domain );
      nextDt = LvArray::math::min( singlePhysicsNextDt, nextDt );
    } );
    return nextDt;
  }

  virtual real64 setNextDtBasedOnNewtonIter( real64 const & currentDt ) override
  {
    real64 nextDt = SolverBase::setNextDtBasedOnNewtonIter( currentDt );
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      real64 const singlePhysicsNextDt =
        solver->setNextDtBasedOnNewtonIter( currentDt );
      nextDt = LvArray::math::min( singlePhysicsNextDt, nextDt );
    } );
    return nextDt;
  }

  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
    } );
    SolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
  }

  /**@}*/

  virtual bool checkSequentialSolutionIncrements( DomainPartition & domain ) const override
  {
    bool isConverged = true;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      isConverged &= solver->checkSequentialSolutionIncrements( domain );
    } );
    return isConverged;
  }

  virtual bool updateConfiguration( DomainPartition & domain ) override
  {
    bool result = true;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      result &= solver->updateConfiguration( domain );
    } );
    return result;
  }

  virtual void outputConfigurationStatistics( DomainPartition const & domain ) const override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->outputConfigurationStatistics( domain );
    } );
  }

  virtual void resetConfigurationToBeginningOfStep( DomainPartition & domain ) override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->resetConfigurationToBeginningOfStep( domain );
    } );
  }

  virtual bool resetConfigurationToDefault( DomainPartition & domain ) const override
  {
    bool result = true;
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      result &= solver->resetConfigurationToDefault( domain );
    } );
    return result;
  }

protected:

  /**
   * @brief Fully coupled solution approach solution step.
   *
   * @param time_n the current time
   * @param dt timestep size
   * @param cycleNumber
   * @param domain the domain partition
   * @return real64 size of the accepted timestep
   */
  virtual real64 fullyCoupledSolverStep( real64 const & time_n,
                                         real64 const & dt,
                                         int const cycleNumber,
                                         DomainPartition & domain )
  {
    return SolverBase::solverStep( time_n, dt, cycleNumber, domain );
  }

  /**
   * @brief Sequentially coupled solver step. It solves a nonlinear system of
   * equations using a sequential approach.
   *
   * @param time_n the current time
   * @param dt timestep size
   * @param cycleNumber
   * @param domain the domain partition
   * @return real64 size of the accepted timestep
   */
  virtual real64 sequentiallyCoupledSolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain )
  {
    GEOS_MARK_FUNCTION;

    // Only build the sparsity pattern if the mesh has changed
    Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      if( meshModificationTimestamp > solver->getSystemSetupTimestamp() )
      {
        solver->setupSystem( domain,
                             solver->getDofManager(),
                             solver->getLocalMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );
        solver->setSystemSetupTimestamp( meshModificationTimestamp );
      }
    } );

    implicitStepSetup( time_n, dt, domain );

    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
    integer const maxNumberDtCuts = solverParams.m_maxTimeStepCuts;
    real64 const dtCutFactor = solverParams.m_timeStepCutFactor;
    integer & dtAttempt = solverParams.m_numTimeStepAttempts;

    bool isConverged = false;
    // dt may be cut during the course of this step, so we are keeping a local
    // value to track the achieved dt for this step.
    real64 stepDt = dt;

    // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
    // required.
    for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
    {
      // TODO configuration loop

      // Reset the states of all solvers if any of them had to restart
      forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
      {
        solver->resetStateToBeginningOfStep( domain );
        solver->getSolverStatistics().initializeTimeStepStatistics();     // initialize counters for subsolvers
      } );
      resetStateToBeginningOfStep( domain );

      integer & iter = solverParams.m_numNewtonIterations;

      /// Sequential coupling loop
      for( iter = 0; iter < solverParams.m_maxIterNewton; iter++ )
      {
        // Increment the solver statistics for reporting purposes
        // Pass a "0" as argument (0 linear iteration) to skip the output of linear iteration stats at the end
        m_solverStatistics.logNonlinearIteration( 0 );

        startSequentialIteration( iter, domain );

        // Solve the subproblems nonlinearly
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
        {
          GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "  Iteration {:2}: {}", iter + 1, solver->getName() ) );
          real64 solverDt = solver->nonlinearImplicitStep( time_n,
                                                           stepDt,
                                                           cycleNumber,
                                                           domain );

          // save fields (e.g. pressure and temperature) after inner solve
          if( solver->getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential )
          {
            solver->saveSequentialIterationState( domain );
          }

          mapSolutionBetweenSolvers( domain, idx() );

          if( solverDt < stepDt ) // subsolver had to cut the time step
          {
            iter = 0; // restart outer loop
            stepDt = solverDt; // sync time step
          }
        } );

        // Check convergence of the outer loop
        isConverged = checkSequentialConvergence( iter,
                                                  time_n,
                                                  stepDt,
                                                  domain );

        if( isConverged )
        {
          // we still want to count current iteration
          ++iter;
          // exit outer loop
          break;
        }
        else
        {
          finishSequentialIteration( iter, domain );
        }
      }

      if( isConverged )
      {
        // Save time step statistics for the subsolvers
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {
          solver->getSolverStatistics().saveTimeStepStatistics();
        } );
        // get out of the time loop
        break;
      }
      else
      {
        // cut timestep, go back to beginning of step and restart the Newton loop
        stepDt *= dtCutFactor;
        GEOS_LOG_LEVEL_RANK_0 ( 1, GEOS_FMT( "New dt = {}", stepDt ) );

        // notify the solver statistics counter that this is a time step cut
        m_solverStatistics.logTimeStepCut();
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {
          solver->getSolverStatistics().logTimeStepCut();
        } );
      }
    }

    if( !isConverged )
    {
      GEOS_LOG_RANK_0( "Convergence not achieved." );

      if( m_nonlinearSolverParameters.m_allowNonConverged > 0 )
      {
        GEOS_LOG_RANK_0( "The accepted solution may be inaccurate." );
      }
      else
      {
        GEOS_ERROR( "Nonconverged solutions not allowed. Terminating..." );
      }
    }

    implicitStepComplete( time_n, stepDt, domain );

    return stepDt;
  }

  /**
   * @brief Maps the solution obtained from one solver to the fields used by the other solver(s)
   *
   * @param domain the domain partition
   * @param solverType the index of the solver withing this coupled solver.
   */
  virtual void mapSolutionBetweenSolvers( DomainPartition & domain,
                                          integer const solverType )
  {
    GEOS_UNUSED_VAR( domain, solverType );
  }

  virtual bool checkSequentialConvergence( int const & iter,
                                           real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition & domain )
  {
    NonlinearSolverParameters const & params = getNonlinearSolverParameters();
    bool isConverged = true;

    if( params.m_subcyclingOption == 0 )
    {
      GEOS_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****" );
    }
    else
    {
      GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "  Iteration {:2}: outer-loop convergence check", iter + 1 ) );

      if( params.sequentialConvergenceCriterion() == NonlinearSolverParameters::SequentialConvergenceCriterion::ResidualNorm )
      {
        real64 residualNorm = 0;

        // loop over all the single-physics solvers
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {

          solver->getLocalMatrix().toViewConstSizes().zero();
          solver->getSystemRhs().zero();
          arrayView1d< real64 > const localRhs = solver->getSystemRhs().open();

          // for each solver, we have to recompute the residual (and Jacobian, although not necessary)
          solver->assembleSystem( time_n,
                                  dt,
                                  domain,
                                  solver->getDofManager(),
                                  solver->getLocalMatrix().toViewConstSizes(),
                                  localRhs );
          solver->applyBoundaryConditions( time_n,
                                           dt,
                                           domain,
                                           solver->getDofManager(),
                                           solver->getLocalMatrix().toViewConstSizes(),
                                           localRhs );
          solver->getSystemRhs().close();

          // once this is done, we recompute the single-physics residual
          real64 const singlePhysicsNorm =
            solver->calculateResidualNorm( time_n,
                                           dt,
                                           domain,
                                           solver->getDofManager(),
                                           solver->getSystemRhs().values() );
          residualNorm += singlePhysicsNorm * singlePhysicsNorm;
        } );

        // finally, we perform the convergence check on the multiphysics residual
        residualNorm = sqrt( residualNorm );
        GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );
        isConverged = ( residualNorm < params.m_newtonTol );

      }
      else if( params.sequentialConvergenceCriterion() == NonlinearSolverParameters::SequentialConvergenceCriterion::NumberOfNonlinearIterations )
      {
        // TODO also make recursive?
        forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
        {
          NonlinearSolverParameters const & singlePhysicsParams = solver->getNonlinearSolverParameters();
          if( singlePhysicsParams.m_numNewtonIterations > singlePhysicsParams.m_minIterNewton )
          {
            isConverged = false;
          }
        } );
      }
      else if( params.sequentialConvergenceCriterion() == NonlinearSolverParameters::SequentialConvergenceCriterion::SolutionIncrements )
      {
        isConverged = checkSequentialSolutionIncrements( domain );
      }
      else
      {
        GEOS_ERROR( getDataContext() << ": Invalid sequential convergence criterion." );
      }

      if( isConverged )
      {
        GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "***** The iterative coupling has converged in {} iteration(s) *****", iter + 1 ) );
      }
    }
    return isConverged;
  }

  virtual void
  postInputInitialization() override
  {
    setSubSolvers();

    bool const isSequential = getNonlinearSolverParameters().couplingType() == NonlinearSolverParameters::CouplingType::Sequential;
    bool const usesLineSearch = getNonlinearSolverParameters().m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None;
    GEOS_THROW_IF( isSequential && usesLineSearch,
                   GEOS_FMT( "{}: line search is not supported by the coupled solver when {} is set to `{}`. Please set {} to `{}` to remove this error",
                             getWrapperDataContext( NonlinearSolverParameters::viewKeysStruct::couplingTypeString() ),
                             NonlinearSolverParameters::viewKeysStruct::couplingTypeString(),
                             EnumStrings< NonlinearSolverParameters::CouplingType >::toString( NonlinearSolverParameters::CouplingType::Sequential ),
                             NonlinearSolverParameters::viewKeysStruct::lineSearchActionString(),
                             EnumStrings< NonlinearSolverParameters::LineSearchAction >::toString( NonlinearSolverParameters::LineSearchAction::None ) ),
                   InputError );

    if( !isSequential )
    {
      synchronizeNonlinearSolverParameters();
    }

    if( m_nonlinearSolverParameters.m_nonlinearAccelerationType != NonlinearSolverParameters::NonlinearAccelerationType::None )
      validateNonlinearAcceleration();
  }

  virtual void validateNonlinearAcceleration()
  {
    GEOS_THROW ( GEOS_FMT( "{}: Nonlinear acceleration {} is not supported by {} solver '{}'",
                           getWrapperDataContext( NonlinearSolverParameters::viewKeysStruct::nonlinearAccelerationTypeString() ),
                           EnumStrings< NonlinearSolverParameters::NonlinearAccelerationType >::toString( m_nonlinearSolverParameters.m_nonlinearAccelerationType ),
                           getCatalogName(), getName()),
                 InputError );
  }

  virtual void
  synchronizeNonlinearSolverParameters() override
  {
    forEachArgInTuple( m_solvers, [&]( auto & solver, auto )
    {
      solver->getNonlinearSolverParameters() = getNonlinearSolverParameters();
    } );
  }

  virtual void startSequentialIteration( integer const & iter,
                                         DomainPartition & domain )
  {
    GEOS_UNUSED_VAR( iter, domain );
  }

  virtual void finishSequentialIteration( integer const & iter,
                                          DomainPartition & domain )
  {
    GEOS_UNUSED_VAR( iter, domain );
  }

protected:

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDSOLVER_HPP_ */
