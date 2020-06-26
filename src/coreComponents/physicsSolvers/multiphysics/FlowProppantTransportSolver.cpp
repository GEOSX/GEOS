/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FlowProppantTransportSolver.cpp
 *
 */


#include "FlowProppantTransportSolver.hpp"

#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/ProppantTransport.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FlowProppantTransportSolver::FlowProppantTransportSolver( const std::string & name, Group * const parent ):
  SolverBase( name, parent ),
  m_proppantSolverName(),
  m_flowSolverName(),
  m_proppantSolver{},
  m_flowSolver{}
{
  registerWrapper( viewKeyStruct::proppantSolverNameString, &m_proppantSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the proppant transport solver to use in the flowProppantTransport solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString, &m_flowSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the flow solver to use in the flowProppantTransport solver" );
}

void FlowProppantTransportSolver::RegisterDataOnMesh( dataRepository::Group * const )
{}

void FlowProppantTransportSolver::ImplicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  m_flowSolver->ImplicitStepSetup( time_n, dt, domain );
  m_proppantSolver->ImplicitStepSetup( time_n, dt, domain );
}

void FlowProppantTransportSolver::ImplicitStepComplete( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_proppantSolver->ImplicitStepComplete( time_n, dt, domain );
}

void FlowProppantTransportSolver::PostProcessInput()
{
  m_proppantSolver = this->getParent()->GetGroup< ProppantTransport >( m_proppantSolverName );
  GEOSX_ERROR_IF( m_proppantSolver == nullptr, "Invalid proppant solver name: " << m_proppantSolverName );

  m_flowSolver = this->getParent()->GetGroup< FlowSolverBase >( m_flowSolverName );
  GEOSX_ERROR_IF( m_flowSolver == nullptr, "Invalid flow solver name: " << m_flowSolverName );
}

FlowProppantTransportSolver::~FlowProppantTransportSolver() = default;

void FlowProppantTransportSolver::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  m_proppantSolver->ResetStateToBeginningOfStep( domain );
}

real64 FlowProppantTransportSolver::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  m_proppantSolver->ResizeFractureFields( *domain.getMeshBody( 0 )->getMeshLevel( 0 ) );

  if( cycleNumber == 0 )
  {
    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::get();
    boundaryConditionManager.ApplyInitialConditions( &domain );
  }

  m_flowSolver->SetupSystem( domain,
                             m_flowSolver->getDofManager(),
                             m_flowSolver->getLocalMatrix(),
                             m_flowSolver->getLocalRhs(),
                             m_flowSolver->getLocalSolution() );

  m_proppantSolver->SetupSystem( domain,
                                 m_proppantSolver->getDofManager(),
                                 m_proppantSolver->getLocalMatrix(),
                                 m_proppantSolver->getLocalRhs(),
                                 m_proppantSolver->getLocalSolution() );

  ImplicitStepSetup( time_n, dt, domain );
  m_proppantSolver->PreStepUpdate( time_n, dt, domain );

  int iter = 0;
  while( iter < this->m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      ResetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );

    dtReturnTemporary = m_flowSolver->NonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    NonlinearSolverParameters const & fluidNonLinearParams = m_flowSolver->getNonlinearSolverParameters();
    if( fluidNonLinearParams.m_numNewtonIterations <= this->m_nonlinearSolverParameters.m_minIterNewton && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", Proppant Solver: " );

    dtReturnTemporary = m_proppantSolver->NonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }

  ImplicitStepComplete( time_n, dtReturn, domain );
  m_proppantSolver->PostStepUpdate( time_n, dt, domain );
  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, FlowProppantTransportSolver, std::string const &, Group * const )

} /* namespace geosx */
