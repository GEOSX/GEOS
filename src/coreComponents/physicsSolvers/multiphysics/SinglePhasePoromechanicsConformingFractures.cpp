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
 * @file SinglePhasePoromechanicsConformingFractures.cpp
 */

#include "SinglePhasePoromechanicsConformingFractures.hpp"

#include "dataRepository/LogLevelsInfo.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/contact/HydraulicApertureRelationSelector.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

template< typename FLOW_SOLVER >
SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::SinglePhasePoromechanicsConformingFractures( const string & name,
                                                                                                         Group * const parent )
  : Base( name, parent )
{
  LinearSolverParameters & params = this->m_linearSolverParameters.get();
  params.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFractures;
  params.mgr.separateComponents = true;
  params.dofsPerNode = 3;
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::setupCoupling( DomainPartition const & domain,
                                                                                DofManager & dofManager ) const
{
  /// We need to add 2 coupling terms:
  // 1. Poromechanical coupling in the bulk
  Base::setupCoupling( domain, dofManager );

  // 2. Traction - pressure coupling in the fracture
  dofManager.addCoupling( SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          fields::contact::traction::key(),
                          DofManager::Connector::Elem );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::setupSystem( DomainPartition & domain,
                                                                              DofManager & dofManager,
                                                                              CRSMatrix< real64, globalIndex > & localMatrix,
                                                                              ParallelVector & rhs,
                                                                              ParallelVector & solution,
                                                                              bool const setSparsity )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( setSparsity );

  /// 1. Add all coupling terms handled directly by the DofManager
  dofManager.setDomain( domain );
  this->setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  /// 2. Add coupling terms not added by the DofManager.
  localIndex const numLocalRows = dofManager.numLocalDofs();

  SparsityPattern< globalIndex > patternOriginal;
  dofManager.setSparsityPattern( patternOriginal );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternOriginal.numRows() );
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  addTransmissibilityCouplingNNZ( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                         patternOriginal.numColumns(),
                                                         rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    globalIndex const * cols = patternOriginal.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternOriginal.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addTransmissibilityCouplingPattern( domain, dofManager, pattern.toView() );

  localMatrix.setName( this->getName() + "/matrix" );
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( numLocalRows, MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( numLocalRows, MPI_COMM_GEOS );

  setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );

  // if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  // {
  //   createPreconditioner( domain );
  // }
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::assembleSystem( real64 const time_n,
                                                                                 real64 const dt,
                                                                                 DomainPartition & domain,
                                                                                 DofManager const & dofManager,
                                                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                 arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  this->solidMechanicsSolver()->synchronizeFractureState( domain );

  assembleElementBasedContributions( time_n,
                                     dt,
                                     domain,
                                     dofManager,
                                     localMatrix,
                                     localRhs );

  // Assemble fluxes 3D/2D and get dFluidResidualDAperture
  this->flowSolver()->assembleHydrofracFluxTerms( time_n,
                                                  dt,
                                                  domain,
                                                  dofManager,
                                                  localMatrix,
                                                  localRhs,
                                                  getDerivativeFluxResidual_dNormalJump() );

  // This step must occur after the fluxes are assembled because that's when DerivativeFluxResidual_dAperture is filled.
  assembleCouplingTerms( time_n,
                         dt,
                         domain,
                         dofManager,
                         localMatrix,
                         localRhs );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::assembleElementBasedContributions( real64 const time_n,
                                                                                                    real64 const dt,
                                                                                                    DomainPartition & domain,
                                                                                                    DofManager const & dofManager,
                                                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                                    arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time_n, dt );

  /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement

  Base::assembleElementBasedTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // Flow accumulation for fractures
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion const & subRegion )
    {
      this->flowSolver()->accumulationAssemblyLaunch( dofManager, subRegion, localMatrix, localRhs );
    } );
  } );

  this->solidMechanicsSolver()->assembleContact( domain, dofManager, localMatrix, localRhs );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::assembleCouplingTerms( real64 const time_n,
                                                                                        real64 const dt,
                                                                                        DomainPartition const & domain,
                                                                                        DofManager const & dofManager,
                                                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                        arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time_n, dt );
  // These 2 steps need to occur after the fluxes are assembled because that's when DerivativeFluxResidual_dAperture is filled.
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    /// 3. assemble Force Residual w.r.t. pressure and Flow mass residual w.r.t. displacement
    assembleForceResidualDerivativeWrtPressure( mesh, regionNames, dofManager, localMatrix, localRhs );
    assembleFluidMassResidualDerivativeWrtDisplacement( mesh, regionNames, dofManager, localMatrix, localRhs );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::
setUpDflux_dApertureMatrix( DomainPartition & domain,
                            DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                            CRSMatrix< real64, globalIndex > & localMatrix )
{
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    std::unique_ptr< CRSMatrix< real64, localIndex > > & derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();

    {
      localIndex numRows = 0;
      mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                          [&]( localIndex const, FaceElementSubRegion const & subRegion )
      {
        numRows += subRegion.size();
      } );

      derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numRows );
      derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );

      derivativeFluxResidual_dAperture->reserveNonZeros( localMatrix.numNonZeros() );
      localIndex maxRowSize = -1;
      for( localIndex row = 0; row < localMatrix.numRows(); ++row )
      {
        localIndex const rowSize = localMatrix.numNonZeros( row );
        maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;
      }
      // TODO This is way too much. The With the full system rowSize is not a good estimate for this.
      for( localIndex row = 0; row < numRows; ++row )
      {
        derivativeFluxResidual_dAperture->reserveNonZeros( row, maxRowSize );
      }
    }

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
          {
            derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0], sei[iconn][k1], 0.0 );
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::
addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                DofManager const & dofManager,
                                arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &, //  meshBodyName,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & ) // regionNames
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const presDofKey = dofManager.getKey( m_pressureKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( this->solidMechanicsSolver()->getStabilizationName() );

    stabilizationMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        FaceElementSubRegion const & elementSubRegion =
          elemManager.getRegion( seri[iconn][0] ).getSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

        ArrayOfArraysView< localIndex const > const elemsToNodes = elementSubRegion.nodeList().toViewConst();

        arrayView1d< globalIndex const > const faceElementDofNumber =
          elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

        for( localIndex k0=0; k0<numFluxElems; ++k0 )
        {
          globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];
          globalIndex const rowNumber = activeFlowDOF - rankOffset;

          if( rowNumber >= 0 && rowNumber < rowLengths.size() )
          {
            for( localIndex k1=0; k1<numFluxElems; ++k1 )
            {
              // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
              // so we only add the coupling with the nodal displacements of the neighbors.
              if( k1 != k0 )
              {
                localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
                rowLengths[rowNumber] += 3*numNodesPerElement;
              }
            }
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::
addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      arrayView1d< string const > const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const presDofKey = dofManager.getKey( m_pressureKey );

    arrayView1d< globalIndex const > const &
    dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    // Get the finite volume method used to compute the stabilization
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fvDiscretization = fvManager.getFluxApproximation( this->flowSolver()->getDiscretizationName() );

    SurfaceElementRegion const & fractureRegion =
      elemManager.getRegion< SurfaceElementRegion >( this->solidMechanicsSolver()->getUniqueFractureRegionName() );
    FaceElementSubRegion const & fractureSubRegion =
      fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

    GEOS_ERROR_IF( !fractureSubRegion.hasWrapper( flow::pressure::key() ),
                   this->getDataContext() << ": The fracture subregion must contain pressure field." );

    arrayView2d< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();

    arrayView1d< globalIndex const > const &
    presDofNumber = fractureSubRegion.getReference< globalIndex_array >( presDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    fvDiscretization.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      forAll< serialPolicy >( stencil.size(), [=] ( localIndex const iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );

        // A fracture connector has to be an edge shared by two faces
        if( numFluxElems == 2 )
        {
          typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

          // First index: face element. Second index: node
          for( localIndex kf = 0; kf < 2; ++kf )
          {
            // Set row DOF index
            // Note that the 1-kf index is intentional, as this is coupling the pressure of one face cell
            // to the nodes of the adjacent cell
            localIndex const rowIndex = presDofNumber[sei[iconn][1-kf]] - rankOffset;

            if( rowIndex >= 0 && rowIndex < pattern.numRows() )
            {

              // Get fracture, face and region/subregion/element indices (for elements on both sides)
              localIndex const fractureIndex = sei[iconn][kf];

              // Get the number of nodes
              localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elem2dToFaces[fractureIndex][0] );

              // Loop over the two sides of each fracture element
              // GEOS_ERROR_IF( elem2dToFaces.sizeOfArray( fractureIndex ) != 2,
              //                "Fracture face " << fractureIndex << " has to be shared by two cells." );
              for( localIndex kf1 = 0; kf1 < 2; ++kf1 )
              {
                localIndex const faceIndex = elem2dToFaces[fractureIndex][kf1];

                // Save the list of DOF associated with nodes
                for( localIndex a=0; a<numNodesPerFace; ++a )
                {
                  for( localIndex i = 0; i < 3; ++i )
                  {
                    globalIndex const colIndex = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
                    pattern.insertNonZero( rowIndex, colIndex );
                  }
                }
              }
            }
          }
        }
      } );
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::
assembleForceResidualDerivativeWrtPressure( MeshLevel const & mesh,
                                            arrayView1d< string const > const & regionNames,
                                            DofManager const & dofManager,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  arrayView1d< real64 const > faceAreas = faceManager.faceArea();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    arrayView1d< globalIndex const > const &
    presDofNumber = subRegion.getReference< globalIndex_array >( presDofKey );
    arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( flow::pressure::key() );
    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

      real64 Nbar[3];
      Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
      Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
      Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
      LvArray::tensorOps::normalize< 3 >( Nbar );
      globalIndex rowDOF[3 * m_maxFaceNodes]; // this needs to be changed when dealing with arbitrary element types
      real64 nodeRHS[3 * m_maxFaceNodes];
      stackArray1d< real64, 3 * m_maxFaceNodes > dRdP( 3*m_maxFaceNodes );
      globalIndex colDOF[1];
      colDOF[0] = presDofNumber[kfe];

      for( localIndex kf=0; kf<2; ++kf )
      {
        localIndex const faceIndex = elemsToFaces[kfe][kf];

        // Compute local area contribution for each node
        stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
        this->solidMechanicsSolver()->computeFaceNodalArea( elemsToFaces[kfe][kf],
                                                            nodePosition,
                                                            faceToNodeMap,
                                                            faceToEdgeMap,
                                                            edgeToNodeMap,
                                                            faceCenters,
                                                            faceNormal,
                                                            faceAreas,
                                                            nodalArea );
        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          real64 const nodalForceMag = -( pressure[kfe] ) * nodalArea[a];
          real64 globalNodalForce[ 3 ];
          LvArray::tensorOps::scaledCopy< 3 >( globalNodalForce, Nbar, nodalForceMag );

          for( localIndex i=0; i<3; ++i )
          {
            rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
            // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
            nodeRHS[3*a+i] = +globalNodalForce[i] * pow( -1, kf );

            // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
            dRdP( 3*a+i ) = -nodalArea[a] * Nbar[i] * pow( -1, kf );
          }
        }

        for( localIndex idof = 0; idof < numNodesPerFace * 3; ++idof )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                        colDOF,
                                                        &dRdP[idof],
                                                        1 );
            RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], nodeRHS[idof] );
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                    arrayView1d< string const > const & regionNames,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{
  GEOS_MARK_FUNCTION;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  arrayView1d< real64 const > faceAreas = faceManager.faceArea();

  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dNormalJump = getDerivativeFluxResidual_dNormalJump().toViewConst();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
  string const & presDofKey = dofManager.getKey( m_pressureKey );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );

    SingleFluidBase const & fluid = this->template getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );
    arrayView2d< real64 const > const & density = fluid.density();

    arrayView1d< globalIndex const > const & presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

    arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();
    arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

    arrayView1d< integer const > const & fractureState = subRegion.getField< fields::contact::fractureState >();

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      globalIndex nodeDOF[2*3*m_maxFaceNodes];
      globalIndex elemDOF[1];
      elemDOF[0] = presDofNumber[kfe];

      real64 Nbar[3];
      Nbar[ 0 ] = faceNormal[kf0][0] - faceNormal[kf1][0];
      Nbar[ 1 ] = faceNormal[kf0][1] - faceNormal[kf1][1];
      Nbar[ 2 ] = faceNormal[kf0][2] - faceNormal[kf1][2];
      LvArray::tensorOps::normalize< 3 >( Nbar );

      stackArray1d< real64, 2*3*m_maxFaceNodes > dRdU( 2*3*m_maxFaceNodes );

      bool const isFractureOpen = ( fractureState[kfe] == fields::contact::FractureState::Open );

      // Accumulation derivative
      if( isFractureOpen )
      {
        for( localIndex kf=0; kf<2; ++kf )
        {
          // Compute local area contribution for each node
          stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
          this->solidMechanicsSolver()->computeFaceNodalArea( elemsToFaces[kfe][kf],
                                                              nodePosition,
                                                              faceToNodeMap,
                                                              faceToEdgeMap,
                                                              edgeToNodeMap,
                                                              faceCenters,
                                                              faceNormal,
                                                              faceAreas,
                                                              nodalArea );

          // TODO: move to something like this plus a static method.
          // localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][kf] );
          // stackArray1d<real64, 4> nodalArea( numNodesPerFace );

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            real64 const dAccumulationResidualdAperture = density[kfe][0] * nodalArea[a];
            for( localIndex i=0; i<3; ++i )
            {
              nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )]
                                                        + LvArray::integerConversion< globalIndex >( i );
              real64 const dAper_dU = -pow( -1, kf ) * Nbar[i];
              dRdU( kf*3*numNodesPerFace + 3*a+i ) = dAccumulationResidualdAperture * dAper_dU;
            }
          }
        }

        localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

        if( localRow >= 0 && localRow < localMatrix.numRows() )
        {

          localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                    nodeDOF,
                                                                    dRdU.data(),
                                                                    2 * 3 * numNodesPerFace );
        }
      }

      // flux derivative
      bool skipAssembly = true;
      localIndex const numColumns = dFluxResidual_dNormalJump.numNonZeros( kfe );
      arraySlice1d< localIndex const > const & columns = dFluxResidual_dNormalJump.getColumns( kfe );
      arraySlice1d< real64 const > const & values = dFluxResidual_dNormalJump.getEntries( kfe );

      skipAssembly &= !isFractureOpen;

      for( localIndex kfe1=0; kfe1<numColumns; ++kfe1 )
      {
        real64 const dR_dAper = values[kfe1];
        localIndex const kfe2 = columns[kfe1];

        bool const isOpen = ( fractureState[kfe2] == fields::contact::FractureState::Open );
        skipAssembly &= !isOpen;

        for( localIndex kf=0; kf<2; ++kf )
        {
          //TODO: We should avoid allocating LvArrays inside kernel
          stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
          this->solidMechanicsSolver()->computeFaceNodalArea( elemsToFaces[kfe2][kf],
                                                              nodePosition,
                                                              faceToNodeMap,
                                                              faceToEdgeMap,
                                                              edgeToNodeMap,
                                                              faceCenters,
                                                              faceNormal,
                                                              faceAreas,
                                                              nodalArea );

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              nodeDOF[ kf*3*numNodesPerFace + 3*a+i ] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe2][kf], a )]
                                                        + LvArray::integerConversion< globalIndex >( i );
              real64 const dAper_dU = -pow( -1, kf ) * Nbar[i] * ( nodalArea[a] / area[kfe2] );
              dRdU( kf*3*numNodesPerFace + 3*a+i ) = dR_dAper * dAper_dU;
            }
          }
        }

        if( !skipAssembly )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                      nodeDOF,
                                                                      dRdU.data(),
                                                                      2 * 3 * numNodesPerFace );
          }
        }
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // call base poromechanics update
  Base::updateState( domain );
  // need to call solid mechanics update separately to compute face displacement jump
  this->solidMechanicsSolver()->updateState( domain );

  // remove the contribution of the hydraulic aperture from the stencil weights
  this->flowSolver()->prepareStencilWeights( domain );

  updateHydraulicApertureAndFracturePermeability( domain );

  // update the stencil weights using the updated hydraulic aperture
  this->flowSolver()->updateStencilWeights( domain );
}

template< typename FLOW_SOLVER >
void SinglePhasePoromechanicsConformingFractures< FLOW_SOLVER >::updateHydraulicApertureAndFracturePermeability( DomainPartition & domain )
{
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 const > const dispJump           = subRegion.getField< contact::dispJump >();
      arrayView1d< real64 const > const area               = subRegion.getElementArea();
      arrayView1d< real64 const > const volume             = subRegion.getElementVolume();
      arrayView2d< real64 const > const fractureTraction   = subRegion.getField< fields::contact::traction >();
      arrayView1d< real64 const > const pressure           = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< fields::flow::aperture0 >();

      arrayView1d< real64 > const aperture                 = subRegion.getElementAperture();
      arrayView1d< real64 > const hydraulicAperture        = subRegion.getField< flow::hydraulicAperture >();
      arrayView1d< real64 > const deltaVolume              = subRegion.getField< flow::deltaVolume >();

      string const porousSolidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::solidNamesString() );
      CoupledSolidBase & porousSolid = subRegion.getConstitutiveModel< CoupledSolidBase >( porousSolidName );

      string const & hydraulicApertureRelationName = subRegion.template getReference< string >( viewKeyStruct::hydraulicApertureRelationNameString()  );
      HydraulicApertureBase const & hydraulicApertureModel = this->template getConstitutiveModel< HydraulicApertureBase >( subRegion, hydraulicApertureRelationName );

      constitutiveUpdatePassThru( hydraulicApertureModel, [&] ( auto & castedHydraulicAperture )
      {
        using HydraulicApertureType = TYPEOFREF( castedHydraulicAperture );
        typename HydraulicApertureType::KernelWrapper hydraulicApertureWrapper = castedHydraulicAperture.createKernelWrapper();

        constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
        {
          typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousMaterialWrapper = castedPorousSolid.createKernelUpdates();

          poromechanicsFracturesKernels::StateUpdateKernel::
            launch< parallelDevicePolicy<> >( subRegion.size(),
                                              porousMaterialWrapper,
                                              hydraulicApertureWrapper,
                                              dispJump,
                                              pressure,
                                              area,
                                              volume,
                                              deltaVolume,
                                              aperture,
                                              oldHydraulicAperture,
                                              hydraulicAperture,
                                              fractureTraction );

        } );
      } );
    } );
  } );
}

template class SinglePhasePoromechanicsConformingFractures<>;
template class SinglePhasePoromechanicsConformingFractures< SinglePhaseReservoirAndWells<> >;

namespace
{
typedef SinglePhasePoromechanicsConformingFractures< SinglePhaseReservoirAndWells<> > SinglePhaseReservoirPoromechanicsConformingFractures;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReservoirPoromechanicsConformingFractures, string const &, Group * const )
typedef SinglePhasePoromechanicsConformingFractures<> SinglePhasePoromechanicsConformingFractures;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhasePoromechanicsConformingFractures, string const &, Group * const )
}

} /* namespace geos */
