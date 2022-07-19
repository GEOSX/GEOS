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
 * @file FlowProppantTransportSolver.cpp
 *
 */


#include "FlowProppantTransportSolver.hpp"

#include "mesh/DomainPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransport.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FlowProppantTransportSolver::FlowProppantTransportSolver( const string & name, Group * const parent ):
  SolverBase( name, parent ),
  m_proppantSolverName(),
  m_flowSolverName(),
  m_flowSolver{},
  m_proppantSolver{}
{
  registerWrapper( viewKeyStruct::proppantSolverNameString(), &m_proppantSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the proppant transport solver to use in the flowProppantTransport solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the flow solver to use in the flowProppantTransport solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );
}

void FlowProppantTransportSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_proppantSolver = &this->getParent().getGroup< ProppantTransport >( m_proppantSolverName );
  m_flowSolver = &this->getParent().getGroup< FlowSolverBase >( m_flowSolverName );
}

FlowProppantTransportSolver::~FlowProppantTransportSolver() = default;

void FlowProppantTransportSolver::preStepUpdate( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  if( time_n <= 0.0 )
  {
    // We need resize composition array in fractures after they are generated by SurfaceGenerator solver
    forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                 MeshLevel & mesh,
                                                 arrayView1d< string const > const & regionNames )
    {
      m_proppantSolver->resizeFractureFields( mesh, regionNames );
      // We need re-apply initial conditions to fractures after they are generated
      FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::getInstance();
      boundaryConditionManager.applyInitialConditions( mesh );
    } );
  }

  m_flowSolver->setupSystem( domain,
                             m_flowSolver->getDofManager(),
                             m_flowSolver->getLocalMatrix(),
                             m_flowSolver->getSystemRhs(),
                             m_flowSolver->getSystemSolution() );


  m_flowSolver->implicitStepSetup( time_n, dt, domain );

  m_proppantSolver->setupSystem( domain,
                                 m_proppantSolver->getDofManager(),
                                 m_proppantSolver->getLocalMatrix(),
                                 m_proppantSolver->getSystemRhs(),
                                 m_proppantSolver->getSystemSolution() );


  m_proppantSolver->implicitStepSetup( time_n, dt, domain );

  m_proppantSolver->preStepUpdate( time_n, dt, domain );
}

void FlowProppantTransportSolver::postStepUpdate( real64 const & time_n,
                                                  real64 const & dt,
                                                  DomainPartition & domain )
{
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
  m_proppantSolver->implicitStepComplete( time_n, dt, domain );
  m_proppantSolver->postStepUpdate( time_n, dt, domain );
}

void FlowProppantTransportSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{

  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_proppantSolver->resetStateToBeginningOfStep( domain );

}


real64 FlowProppantTransportSolver::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  preStepUpdate( time_n, dt, domain );

  // reset number of nonlinear iterations
  m_solverStatistics.initializeTimeStepStatistics();

  int iter = 0;
  while( iter < this->m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all sub-solvers if any of them has been reset
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );

    dtReturnTemporary = m_flowSolver->nonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    NonlinearSolverParameters const & fluidNonLinearParams = m_flowSolver->getNonlinearSolverParameters();
    if( fluidNonLinearParams.m_numNewtonIterations <= this->m_nonlinearSolverParameters.m_minIterNewton && iter > 0 )
    {
      m_solverStatistics.logNonlinearIteration();
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", Proppant Solver: " );

    dtReturnTemporary = m_proppantSolver->nonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }

  // increment the cumulative number of nonlinear iterations
  m_solverStatistics.saveTimeStepStatistics();

  postStepUpdate( time_n, dtReturn, domain );

  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, FlowProppantTransportSolver, string const &, Group * const )

} /* namespace geosx */
