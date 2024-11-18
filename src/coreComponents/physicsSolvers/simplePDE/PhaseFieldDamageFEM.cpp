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
 * @file PhaseFieldDamageFEM.hpp
 */

#include "PhaseFieldDamageFEM.hpp"
#include "PhaseFieldDamageFEMKernels.hpp"
#include <math.h>
#include <vector>

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/Damage.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "discretizationMethods/NumericalMethodsManager.hpp"

#include "mesh/DomainPartition.hpp"

namespace geos
{

namespace dataRepository
{
namespace keys
{}
} // namespace dataRepository

using namespace dataRepository;
using namespace constitutive;

PhaseFieldDamageFEM::PhaseFieldDamageFEM( const string & name,
                                          Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_fieldName( "primaryField" )
{

  registerWrapper< TimeIntegrationOption >( PhaseFieldDamageFEMViewKeys.timeIntegrationOption.key(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "option for default time integration method" );

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.fieldVarName.key(), &m_fieldName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "name of field variable" );

  registerWrapper( viewKeyStruct::localDissipationOptionString(), &m_localDissipationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Type of local dissipation function. Can be Linear or Quadratic" );

  registerWrapper( viewKeyStruct::irreversibilityFlagString(), &m_irreversibilityFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The flag to indicate whether to apply the irreversibility constraint" );

  registerWrapper( viewKeyStruct::damageUpperBoundString(), &m_damageUpperBound ).
    setApplyDefaultValue( 1.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The upper bound of the damage" );
}

PhaseFieldDamageFEM::~PhaseFieldDamageFEM()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldDamageFEM::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodes = mesh.getNodeManager();

    nodes.registerWrapper< real64_array >( m_fieldName )
      .setApplyDefaultValue( 0.0 )
      .setPlotLevel( PlotLevel::LEVEL_0 )
      .setDescription( "Primary field variable" );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [ &]( localIndex const, CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper( viewKeyStruct::coeffNameString(), &m_coeff ).
        setApplyDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setDescription( "field variable representing the diffusion coefficient" );


      subRegion.registerWrapper< string >( viewKeyStruct::solidModelNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidModelNamesString() );
      solidMaterialName = PhysicsSolverBase::getConstitutiveName< SolidBase >( subRegion );
      GEOS_ERROR_IF( solidMaterialName.empty(), GEOS_FMT( "{}: SolidBase model not found on subregion {}",
                                                          getDataContext(), subRegion.getName() ) );

    } );
  } );
}

void PhaseFieldDamageFEM::postInputInitialization()
{
  PhysicsSolverBase::postInputInitialization();

  // Set basic parameters for solver
  // m_linearSolverParameters.logLevel = 0;
  // m_linearSolverParameters.solverType = "gmres";
  // m_linearSolverParameters.krylov.tolerance = 1e-8;
  // m_linearSolverParameters.krylov.maxIterations = 250;
  // m_linearSolverParameters.krylov.maxRestart = 250;
  // m_linearSolverParameters.preconditionerType = "amg";
  // m_linearSolverParameters.amg.smootherType = "gaussSeidel";
  // m_linearSolverParameters.amg.coarseType = "direct";
}

real64 PhaseFieldDamageFEM::solverStep( real64 const & time_n,
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitTransient )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == TimeIntegrationOption::SteadyState )
  {
    this->setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution, false );

    dtReturn = this->nonlinearImplicitStep( time_n,
                                            dt,
                                            cycleNumber,
                                            domain );
  }
  return dtReturn;
}

real64 PhaseFieldDamageFEM::explicitStep( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                          real64 const & dt,
                                          const int GEOS_UNUSED_PARAM( cycleNumber ),
                                          DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  return dt;
}

void PhaseFieldDamageFEM::implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                real64 const & GEOS_UNUSED_PARAM( dt ),
                                                DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{}

void PhaseFieldDamageFEM::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                     DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
  dofManager.addField( m_fieldName,
                       FieldLocation::Node,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );

}

void PhaseFieldDamageFEM::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > const & dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

    // Initialize all entries to zero
#if 1 // Andre...this is the new code
    localMatrix.zero();
    localRhs.zero();

    auto const localDissipation = m_localDissipationOption == LocalDissipation::Linear ?
                                  PhaseFieldDamageKernelLocalDissipation::Linear :
                                  PhaseFieldDamageKernelLocalDissipation::Quadratic;

    PhaseFieldDamageKernelFactory kernelFactory( dofIndex,
                                                 dofManager.rankOffset(),
                                                 localMatrix,
                                                 localRhs,
                                                 dt,
                                                 m_fieldName,
                                                 localDissipation );

    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy<>,
                                    constitutive::DamageBase,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            this->getDiscretizationName(),
                                                            viewKeyStruct::solidModelNamesString(),
                                                            kernelFactory );
#else // this has your changes to the old base code
    matrix.zero();
    rhs.zero();

    matrix.open();
    rhs.open();

    // begin region loop
    for( localIndex er = 0; er < elemManager->numRegions(); ++er )
    {
      ElementRegionBase & elementRegion = elemManager->GetRegion( er );

      elementRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const GEOS_UNUSED_PARAM( esr ),
                                                                            CellElementSubRegion & elementSubRegion )
      {

        constitutive::ConstitutiveBase &
        solidModel = elementSubRegion.getConstitutiveModel< constitutive::ConstitutiveBase >( m_solidModelName );

        constitutive::ConstitutivePassThru< constitutive::DamageBase >::execute( solidModel,
                                                                                 [&]( auto & damageModel )
        {
          using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
          typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

          arrayView4d< real64 const > const &
          dNdX = elementSubRegion.dNdX();

          arrayView2d< real64 const > const &
          detJ = elementSubRegion.detJ();

          localIndex const numNodesPerElement =  elementSubRegion.numNodesPerElement();
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemNodes = elementSubRegion.nodeList();

          // arrayView1d<real64 const> const &
          // coeff = elementSubRegion.getReference<array1d<real64> >(viewKeyStruct::coeffNameString());

          globalIndex_array elemDofIndex( numNodesPerElement );
          real64_array element_rhs( numNodesPerElement );
          real64_array2d element_matrix( numNodesPerElement, numNodesPerElement );

          arrayView1d< integer const > const & elemGhostRank = elementSubRegion.ghostRank();
          std::unique_ptr< FiniteElementBase > finiteElement = feDiscretization->getFiniteElement( elementSubRegion.getElementType() );
          localIndex const n_q_points = finiteElement->n_quadrature_points();

          //real64 ell = m_lengthScale;                       //phase-field length scale
          real64 ell = constitutiveUpdate.getRegularizationLength();
          //real64 Gc = m_criticalFractureEnergy;             //energy release rate
          real64 Gc = constitutiveUpdate.getCriticalFractureEnergy();

          real64 threshold = constitutiveUpdate.getEnergyThreshold();//elastic energy threshold - use when Local Dissipation is linear

          arrayView1d< real64 > const & nodalDamage = nodeManager.getReference< array1d< real64 > >( m_fieldName );
          //real64 diffusion = 1.0;
          // begin element loop, skipping ghost elements
          for( localIndex k = 0; k < elementSubRegion.size(); ++k )
          {
            if( elemGhostRank[k] < 0 )
            {
              element_rhs = 0.0;
              element_matrix = 0.0;
              for( localIndex q = 0; q < n_q_points; ++q )
              {
                real64 const strainEnergyDensity = constitutiveUpdate.calculateStrainEnergyDensity( k, q );
                real64 D = 0;                                                                 //max between threshold and
                                                                                              // Elastic energy
                if( m_localDissipationOption == LocalDissipation::Linear )
                {
                  D = std::max( threshold, strainEnergyDensity );
                  //D = max(strainEnergy(k,q), strainEnergy(k,q));//debbuging line - remove after testing
                }
                //Interpolate d and grad_d

                real64 qp_damage = 0.0;
                R1Tensor qp_grad_damage;
                R1Tensor temp;
                for( localIndex a = 0; a < numNodesPerElement; ++a )
                {
                  qp_damage += finiteElement->value( a, q ) * nodalDamage[elemNodes( k, a )];
                  temp = dNdX[k][q][a];
                  temp *= nodalDamage[elemNodes( k, a )];
                  qp_grad_damage += temp;

                }
                //std::cout << "Damage: " << qp_damage <<std::endl;
                //std::cout << "GradDamage: " << qp_grad_damage <<std::endl;
                for( localIndex a = 0; a < numNodesPerElement; ++a )
                {
                  elemDofIndex[a] = dofIndex[elemNodes( k, a )];
                  //real64 diffusion = 1.0;
                  real64 Na = finiteElement->value( a, q );
                  //element_rhs(a) += detJ[k][q] * Na * myFunc(Xq, Yq, Zq); //older reaction diffusion solver
                  if( m_localDissipationOption == LocalDissipation::Linear )
                  {
                    // element_rhs( a ) += detJ[k][q] * (Na * (ell * D - 3 * Gc / 16 )/ Gc -
                    //                                   0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( qp_grad_damage, dNdX[k][q][a] )
                    // -
                    //                                   (ell * D/Gc) * Na * qp_damage);

                    element_rhs( a ) += detJ[k][q] * ( -3 * Na / 16  -
                                                       0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( qp_grad_damage, dNdX[k][q][a] ) -
                                                       (0.5 * ell * D/Gc) * Na * constitutiveUpdate.GetDegradationDerivative( qp_damage ));

                  }
                  else
                  {
                    // element_rhs( a ) += detJ[k][q] * (Na * (2 * ell) * strainEnergyDensity / Gc -
                    //                                   (pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( qp_grad_damage, dNdX[k][q][a] ) +
                    //                                    Na * qp_damage * (1 + 2 * ell*strainEnergyDensity/Gc)) );


                    element_rhs( a ) -= detJ[k][q] * (Na * qp_damage +
                                                      (pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( qp_grad_damage, dNdX[k][q][a] ) +
                                                       Na * constitutiveUpdate.GetDegradationDerivative( qp_damage ) * (ell*strainEnergyDensity/Gc)) );
                  }

                  for( localIndex b = 0; b < numNodesPerElement; ++b )
                  {
                    real64 Nb = finiteElement->value( b, q );
                    if( m_localDissipationOption == LocalDissipation::Linear )
                    {
                      // element_matrix( a, b ) -= detJ[k][q] *
                      //                           (0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( dNdX[k][q][a], dNdX[k][q][b] ) +
                      //                            (ell * D/Gc) * Na * Nb);
                      //
                      element_matrix( a, b ) -= detJ[k][q] *
                                                (0.375*pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( dNdX[k][q][a], dNdX[k][q][b] ) +
                                                 (0.5 * ell * D/Gc) * constitutiveUpdate.GetDegradationSecondDerivative( qp_damage ) * Na * Nb);

                    }
                    else
                    {
                      // element_matrix( a, b ) -= detJ[k][q] *
                      //                           ( pow( ell, 2 ) * LvArray::tensorOps::AiBi<3>( dNdX[k][q][a], dNdX[k][q][b] ) +
                      //                               Na * Nb * (1 + 2 * ell*strainEnergyDensity/Gc )
                      //                           );

                      element_matrix( a, b ) -= detJ[k][q] *
                                                ( pow( ell, 2 ) * LvArray::tensorOps::AiBi< 3 >( dNdX[k][q][a], dNdX[k][q][b] ) +
                                                  Na * Nb * (1 + constitutiveUpdate.GetDegradationSecondDerivative( qp_damage ) * ell*strainEnergyDensity/Gc )
                                                );
                    }
                  }
                }
              }
              matrix.add( elemDofIndex, elemDofIndex, element_matrix );
              rhs.add( elemDofIndex, element_rhs );
            }
          }
        } );
      } );
    }
    matrix.close();
    rhs.close();

    if( getLogLevel() == 2 )
    {
      GEOS_LOG_RANK_0( "After PhaseFieldDamageFEM::AssembleSystem" );
      GEOS_LOG_RANK_0( "\nJacobian:\n" );
      std::cout << matrix;
      GEOS_LOG_RANK_0( "\nResidual:\n" );
      std::cout << rhs;
    }

    if( getLogLevel() >= 3 )
    {
      NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
      integer newtonIter = solverParams.m_numNewtonIterations;
#endif
  } );
}

void PhaseFieldDamageFEM::applySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               real64 const dt,
                                               DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  GEOS_MARK_FUNCTION;

  dofManager.addVectorToField( localSolution, m_fieldName, m_fieldName, scalingFactor );

  // Syncronize ghost nodes
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { m_fieldName } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         false );
  } );

}

void PhaseFieldDamageFEM::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
}

void PhaseFieldDamageFEM::applyBoundaryConditions(
  real64 const time_n,
  real64 const dt, DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  applyDirichletBCImplicit( time_n + dt, dofManager, domain, localMatrix, localRhs );

  // Apply the crack irreversibility constraint
  if( m_irreversibilityFlag )
  {
    applyIrreversibilityConstraint( dofManager, domain, localMatrix, localRhs );
  }

  if( getLogLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After PhaseFieldDamageFEM::applyBoundaryConditions" );
    GEOS_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << localMatrix.toViewConst();
    GEOS_LOG_RANK_0( "\nResidual:\n" );
    std::cout << localRhs;
  }
//
//  if( getLogLevel() >= 3 )
//  {
//    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
//    integer newtonIter = solverParams.m_numNewtonIterations;
//
//    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" +
//                          std::to_string( newtonIter ) + ".mtx";
//    matrix.write( filename_mat );
//
//    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" +
//                          std::to_string( newtonIter ) + ".mtx";
//    rhs.write( filename_rhs );
//
//    GEOS_LOG_RANK_0( "After PhaseFieldDamageFEM::applyBoundaryConditions" );
//    GEOS_LOG_RANK_0( "Jacobian: written to " << filename_mat );
//    GEOS_LOG_RANK_0( "Residual: written to " << filename_rhs );
//  }
}

real64
PhaseFieldDamageFEM::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            DomainPartition const & domain,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    const NodeManager & nodeManager = mesh.getNodeManager();
    const arrayView1d< const integer > & ghostRank = nodeManager.ghostRank();

    const arrayView1d< const globalIndex > &
    dofNumber = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );
    const globalIndex rankOffset = dofManager.rankOffset();


    forAll< parallelDevicePolicy<> >( nodeManager.size(),
                                      [localRhs, localSum, dofNumber, rankOffset, ghostRank] GEOS_HOST_DEVICE ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
        localSum += localRhs[localRow] * localRhs[localRow];
      }
    } );
  } );
  const real64 localResidualNorm[2] = { localSum.get(), 1};

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  const int rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  const int size = MpiWrapper::commSize( MPI_COMM_GEOS );
  array1d< real64 > globalValues( size * 2 );

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      2,
                      globalValues.data(),
                      2,
                      0,
                      MPI_COMM_GEOS );

  if( rank==0 )
  {
    for( int r=0; r<size; ++r )
    {
      // sum/max across all ranks
      globalResidualNorm[0] += globalValues[r*2];
      globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
    }
  }

  MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOS );


  const real64 residual = sqrt( globalResidualNorm[0] ) / ( globalResidualNorm[1] );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residual );
  }

  return residual;
}

void PhaseFieldDamageFEM::applyDirichletBCImplicit( real64 const time,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )

{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    fsManager.template apply< NodeManager >( time,
                                             mesh,
                                             m_fieldName,
                                             [&]( FieldSpecificationBase const & bc,
                                                  string const &,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  NodeManager & targetGroup,
                                                  string const GEOS_UNUSED_PARAM( fieldName ) ) -> void
    {
      bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                         parallelDevicePolicy< > >( targetSet,
                                                                    time,
                                                                    targetGroup,
                                                                    m_fieldName,
                                                                    dofManager.getKey( m_fieldName ),
                                                                    dofManager.rankOffset(),
                                                                    localMatrix,
                                                                    localRhs );
    } );

    fsManager.applyFieldValue< serialPolicy >( time, mesh, viewKeyStruct::coeffNameString() );
  } );
}

void PhaseFieldDamageFEM::applyIrreversibilityConstraint( DofManager const & dofManager,
                                                          DomainPartition & domain,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > const dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

    arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( m_fieldName );

    globalIndex const rankOffSet = dofManager.rankOffset();

    real64 const damangeUpperBound = m_damageUpperBound;

    forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const nodeIndex )
    {
      localIndex const dof = dofIndex[nodeIndex];

      if( dof > -1 )
      {
        real64 const damageAtNode = nodalDamage[nodeIndex];

        if( damageAtNode >= damangeUpperBound )
        {

          // Specify the contribution to rhs
          real64 rhsContribution;

          FieldSpecificationEqual::SpecifyFieldValue( dof,
                                                      rankOffSet,
                                                      localMatrix,
                                                      rhsContribution,
                                                      damangeUpperBound,
                                                      damageAtNode );

          globalIndex const localRow = dof - rankOffSet;

          if( localRow >= 0 && localRow < localRhs.size() )
          {
            localRhs[ localRow ] = rhsContribution;
          }
        }
      }
    } );

  } );
}

void PhaseFieldDamageFEM::saveSequentialIterationState( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // nothing to save yet
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, PhaseFieldDamageFEM, string const &,
                        Group * const )
} // namespace geos
