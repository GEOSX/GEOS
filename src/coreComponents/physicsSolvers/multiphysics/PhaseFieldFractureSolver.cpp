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
 * @file PhaseFieldFractureSolver.cpp
 *
 */

#include "PhaseFieldFractureSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PhaseFieldFractureSolver::PhaseFieldFractureSolver( const string & name,
                                                    Group * const parent ):
  Base( name, parent )
{
  m_couplingType = CouplingType::Sequential;  
}

void PhaseFieldFractureSolver::registerDataOnMesh( Group & meshBodies )
{
  Base::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion,
                                      FaceElementSubRegion >( [&] ( auto & elementSubRegion )
    {
      elementSubRegion.template registerWrapper< array1d< real64 > >( viewKeyStruct::totalMeanStressString() ).
        setDescription( "Total Mean Stress" );
      elementSubRegion.template registerWrapper< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() ).
        setDescription( "Total Mean Stress" );
    } );
  } );
}

void PhaseFieldFractureSolver::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );

    //***** loop over all elements and initialize the derivative arrays *****
    forAllElemsInMesh( mesh, [ &]( localIndex const er, localIndex const esr, localIndex const k )
    {
      oldTotalMeanStress[er][esr][k] = totalMeanStress[er][esr][k];
    } );
  } );
}

PhaseFieldFractureSolver::~PhaseFieldFractureSolver()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldFractureSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );

    //***** loop over all elements and initialize the derivative arrays *****
    forAllElemsInMesh( mesh, [ &]( localIndex const er,
                                   localIndex const esr,
                                   localIndex const k )
    {
      totalMeanStress[er][esr][k] = oldTotalMeanStress[er][esr][k];
    } );
  } );
}

real64 PhaseFieldFractureSolver::solverStep( real64 const & time_n,
                                             real64 const & dt,
                                             int const cycleNumber,
                                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  if( m_couplingTypeOption == CouplingTypeOption::FixedStress )
  {
    dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::TightlyCoupled )
  {
    GEOSX_ERROR( "CouplingTypeOption::FullyImplicit not yet implemented" );
  }
  return dtReturn;
}

real64 PhaseFieldFractureSolver::splitOperatorStep( real64 const & time_n,
                                                    real64 const & dt,
                                                    integer const cycleNumber,
                                                    DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsLagrangianFEM &
  solidSolver = this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  PhaseFieldDamageFEM &
  damageSolver = this->getParent().getGroup< PhaseFieldDamageFEM >( m_damageSolverName );

  damageSolver.setupSystem( domain,
                            damageSolver.getDofManager(),
                            damageSolver.getLocalMatrix(),
                            damageSolver.getSystemRhs(),
                            damageSolver.getSystemSolution(),
                            true );

  solidSolver.setupSystem( domain,
                           solidSolver.getDofManager(),
                           solidSolver.getLocalMatrix(),
                           solidSolver.getSystemRhs(),
                           solidSolver.getSystemSolution() );

  damageSolver.implicitStepSetup( time_n, dt, domain );

  solidSolver.implicitStepSetup( time_n, dt, domain );

  this->implicitStepSetup( time_n, dt, domain );

  NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
  integer & iter = solverParams.m_numNewtonIterations;
  iter = 0;
  bool isConverged = false;
  while( iter < solverParams.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      damageSolver.resetStateToBeginningOfStep( domain );
      solidSolver.resetStateToBeginningOfStep( domain );
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", MechanicsSolver: " );

    dtReturnTemporary = solidSolver.nonlinearImplicitStep( time_n,
                                                           dtReturn,
                                                           cycleNumber,
                                                           domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( solidSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }
    else if( m_subcyclingOption == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
      isConverged = true;
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", DamageSolver: " );

    dtReturnTemporary = damageSolver.nonlinearImplicitStep( time_n,
                                                            dtReturn,
                                                            cycleNumber,
                                                            domain );

    mapDamageToQuadrature( domain );

    //std::cout << "Here: " << dtReturnTemporary << std::endl;

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( m_subcyclingOption == 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
      isConverged = true;
      break;
    }

    ++iter;
  }

  GEOSX_ERROR_IF( !isConverged, "PhaseFieldFractureSolver::SplitOperatorStep() did not converge" );

  damageSolver.implicitStepComplete( time_n, dt, domain );
  solidSolver.implicitStepComplete( time_n, dt, domain );
  this->implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

void PhaseFieldFractureSolver::mapDamageToQuadrature( DomainPartition & domain )
{

  GEOSX_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    PhaseFieldDamageFEM const &
    damageSolver = this->getParent().getGroup< PhaseFieldDamageFEM >( m_damageSolverName );

    string const & damageFieldName = damageSolver.getFieldName();

    //should get reference to damage field here.
    arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

    ElementRegionManager & elemManager = mesh.getElemManager();

    // begin region loop
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [this, nodalDamage]
                                                                ( localIndex const,
                                                                CellElementSubRegion & elementSubRegion )
    {
      string const & solidModelName = elementSubRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
      constitutive::SolidBase &
      solidModel = elementSubRegion.getConstitutiveModel< constitutive::SolidBase >( solidModelName );

      ConstitutivePassThru< DamageBase >::execute( solidModel, [this, &elementSubRegion, nodalDamage]( auto & damageModel )
      {
        using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
        typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

        arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_damage;
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemNodes = elementSubRegion.nodeList();

        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( m_discretizationName );

        finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [nodalDamage, &elementSubRegion, damageFieldOnMaterial, elemNodes]( auto & finiteElement )
        {
          using FE_TYPE = TYPEOFREF( finiteElement );
          constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
          constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

          forAll< serialPolicy >( elementSubRegion.size(), [nodalDamage, damageFieldOnMaterial, elemNodes] ( localIndex const k )
          {
            for( localIndex q = 0; q < n_q_points; ++q )
            {
              real64 N[ numNodesPerElement ];
              FE_TYPE::calcN( q, N );

              damageFieldOnMaterial( k, q ) = 0;
              for( localIndex a = 0; a < numNodesPerElement; ++a )
              {
                damageFieldOnMaterial( k, q ) += N[a] * nodalDamage[elemNodes( k, a )];
                //solution is probably not going to work because the solution of the coupled solver
                //has both damage and displacements. Using the damageResult field from the Damage solver
                //is probably better
                //            std::cout<<"q, N, Dnode = "<<q<<", "<<feDiscretization->m_finiteElement->value(a, q)<<",
                // "<<nodalDamage[elemNodes(k, a)]<<std::endl;
              }
              //          std::cout<<"damage("<<k<<","<<q<<") = "<<damageFieldOnMaterial(k,q)<<std::endl;
            }
          } );
        } );
      } );
    } );
  } );

}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldFractureSolver, string const &, Group * const )

} /* namespace geosx */
