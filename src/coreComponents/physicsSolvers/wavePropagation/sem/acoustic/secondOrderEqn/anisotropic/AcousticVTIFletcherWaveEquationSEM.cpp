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
 * @file AcousticVTIFletcherWaveEquationSEM.cpp
 */

#include "AcousticVTIFletcherWaveEquationSEM.hpp"
#include "AcousticVTIFletcherWaveEquationSEMKernel.hpp"
#include "AcousticVTIFletcherAdjointWaveEquationSEMKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticTimeSchemeSEMKernel.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp"
#include "events/EventManager.hpp"
// #include "AcousticPMLSEMKernel.hpp" // Not working with VTI
#include "physicsSolvers/wavePropagation/shared/PrecomputeSourcesAndReceiversKernel.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

AcousticVTIFletcherWaveEquationSEM::AcousticVTIFletcherWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

}

AcousticVTIFletcherWaveEquationSEM::~AcousticVTIFletcherWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void AcousticVTIFletcherWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void AcousticVTIFletcherWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< acousticfields::ForcingRHS,
                               acousticfields::AcousticMassVector,
                               acousticfields::AcousticFreeSurfaceNodeIndicator,
                               acousticvtifields::Pressure_p_nm1,
                               acousticvtifields::Pressure_p_n,
                               acousticvtifields::Pressure_p_np1,
                               acousticvtifields::Pressure_q_nm1,
                               acousticvtifields::Pressure_q_n,
                               acousticvtifields::Pressure_q_np1,
                               acousticvtifields::DampingVector_pp,
                               acousticvtifields::DampingVector_pq,
                               acousticvtifields::DampingVector_qq,
                               acousticvtifields::DampingVector_qp,
                               acousticvtifields::StiffnessVector_p,
                               acousticvtifields::StiffnessVector_q,
                               acousticvtifields::AcousticLateralSurfaceNodeIndicator,
                               acousticvtifields::AcousticBottomSurfaceNodeIndicator >( getName() );

    /// PML is not supported
    if( m_usePML )
    {
      GEOS_ERROR( "This option (PML) is not supported" );
    }

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< acousticfields::AcousticFreeSurfaceFaceIndicator >( getName() );
    faceManager.registerField< acousticvtifields::AcousticLateralSurfaceFaceIndicator >( getName() );
    faceManager.registerField< acousticvtifields::AcousticBottomSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< acousticfields::AcousticVelocity >( getName() );
      subRegion.registerField< acousticfields::AcousticDensity >( getName() );
      subRegion.registerField< acousticfields::PartialGradient >( getName() );

      subRegion.registerField< acousticvtifields::AcousticDelta >( getName() );
      subRegion.registerField< acousticvtifields::AcousticEpsilon >( getName() );
    } );

  } );
}


void AcousticVTIFletcherWaveEquationSEM::postInputInitialization()
{
  WaveSolverBase::postInputInitialization();

  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, m_receiverCoordinates.size( 0 ) + 1 );
}

void AcousticVTIFletcherWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< globalIndex const > const nodeLocalToGlobalMap = baseMesh.getNodeManager().localToGlobalMap().toViewConst();
  ArrayOfArraysView< localIndex const > const nodesToElements = baseMesh.getNodeManager().elementList().toViewConst();
  ArrayOfArraysView< localIndex const > const facesToNodes = baseMesh.getFaceManager().nodeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodeCoords = baseMesh.getNodeManager().referencePosition();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< EXEC_POLICY >( -1 );
  receiverConstants.setValues< EXEC_POLICY >( -1 );
  receiverIsLocal.zero();

  mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                localIndex const er,
                                                                                                localIndex const esr,
                                                                                                ElementRegionBase &,
                                                                                                CellElementSubRegion & elementSubRegion )
  {
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                   getDataContext() << ": Invalid type of element, the acoustic VTI Fletcher solver is designed for hexahedral meshes only (C3D8), using the SEM formulation",
                   InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes = baseMesh.getElemManager().getRegion( er ).getSubRegion< CellElementSubRegion >( esr ).nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();
    arrayView1d< globalIndex const > const elemLocalToGlobalMap = elementSubRegion.localToGlobalMap().toViewConst();

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      {
        GEOS_MARK_SCOPE( acousticVTIFletcherWaveEquationSEMKernels::PrecomputeSourceAndReceiverKernel );
        PreComputeSourcesAndReceivers::
          Compute1DSourceAndReceiverConstants
        < EXEC_POLICY, FE_TYPE >
          ( elementSubRegion.size(),
          facesToNodes,
          nodeCoords,
          nodeLocalToGlobalMap,
          elemLocalToGlobalMap,
          nodesToElements,
          baseElemsToNodes,
          elemGhostRank,
          elemsToNodes,
          elemsToFaces,
          elemCenter,
          sourceCoordinates,
          sourceIsAccessible,
          sourceNodeIds,
          sourceConstants,
          receiverCoordinates,
          receiverIsLocal,
          receiverNodeIds,
          receiverConstants );
      }
    } );
    elementSubRegion.faceList().freeOnDevice();
    baseMesh.getElemManager().getRegion( er ).getSubRegion< CellElementSubRegion >( esr ).nodeList().freeOnDevice();
    elementSubRegion.getElementCenter().freeOnDevice();
    elementSubRegion.ghostRank().freeOnDevice();
    elementSubRegion.localToGlobalMap().freeOnDevice();
  } );
  baseMesh.getNodeManager().localToGlobalMap().freeOnDevice();
  baseMesh.getNodeManager().elementList().toView().freeOnDevice();
  baseMesh.getFaceManager().nodeList().toView().freeOnDevice();
  baseMesh.getNodeManager().referencePosition().freeOnDevice();
  facesToNodes.freeOnDevice();
  nodesToElements.freeOnDevice();
}

void AcousticVTIFletcherWaveEquationSEM::addSourceToRightHandSide( real64 const & time_n, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  real32 const timeSourceFrequency = m_timeSourceFrequency;
  real32 const timeSourceDelay = m_timeSourceDelay;
  localIndex const rickerOrder = m_rickerOrder;
  bool useSourceWaveletTables = m_useSourceWaveletTables;
  arrayView1d< TableFunction::KernelWrapper const > const sourceWaveletTableWrappers = m_sourceWaveletTableWrappers.toViewConst();
  forAll< EXEC_POLICY >( sourceConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const isrc )
  {
    if( sourceIsAccessible[isrc] == 1 )
    {
      real64 const srcValue =
        useSourceWaveletTables ? sourceWaveletTableWrappers[ isrc ].compute( &time_n ) : WaveSolverUtils::evaluateRicker( time_n, timeSourceFrequency, timeSourceDelay, rickerOrder );
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        real32 const localIncrement = sourceConstants[isrc][inode] * srcValue;
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhs[sourceNodeIds[isrc][inode]], localIncrement );
      }
    }
  } );
}

void AcousticVTIFletcherWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;
  {
    GEOS_MARK_SCOPE( WaveSolverBase::initializePostInitialConditionsPreSubGroups );
    WaveSolverBase::initializePostInitialConditionsPreSubGroups();
  }
  if( m_usePML )
  {
    AcousticVTIFletcherWaveEquationSEM::initializePML();
  }

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  applyFreeSurfaceBC( 0.0, domain );
  precomputeSurfaceFieldIndicator( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    MeshLevel & baseMesh = domain.getMeshBodies().getGroup< MeshBody >( meshBodyName ).getBaseDiscretization();
    precomputeSourceAndReceiverTerm( baseMesh, mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< acousticfields::AcousticMassVector >();
    mass.zero();

    /// damping matrices to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping_pp = nodeManager.getField< acousticvtifields::DampingVector_pp >();
    arrayView1d< real32 > const damping_pq = nodeManager.getField< acousticvtifields::DampingVector_pq >();
    arrayView1d< real32 > const damping_qp = nodeManager.getField< acousticvtifields::DampingVector_qp >();
    arrayView1d< real32 > const damping_qq = nodeManager.getField< acousticvtifields::DampingVector_qq >();
    damping_pp.zero();
    damping_pq.zero();
    damping_qp.zero();
    damping_qq.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< acousticfields::AcousticFreeSurfaceFaceIndicator >();
    arrayView1d< localIndex const > const lateralSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticLateralSurfaceFaceIndicator >();
    arrayView1d< localIndex const > const bottomSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticBottomSurfaceFaceIndicator >();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      computeTargetNodeSet( elemsToNodes, elementSubRegion.size(), fe.getNumQuadraturePoints() );

      arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfields::AcousticVelocity >();
      arrayView1d< real32 const > const density = elementSubRegion.getField< acousticfields::AcousticDensity >();
      arrayView1d< real32 const > const vti_epsilon = elementSubRegion.getField< acousticvtifields::AcousticEpsilon >();
      arrayView1d< real32 const > const vti_delta = elementSubRegion.getField< acousticvtifields::AcousticDelta >();
      arrayView1d< real32 const > const vti_sigma = elementSubRegion.getField< acousticvtifields::AcousticSigma >();

      /// Partial gradient if gradient as to be computed
      arrayView1d< real32 > grad = elementSubRegion.getField< acousticfields::PartialGradient >();

      grad.zero();

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        AcousticMatricesSEM::MassMatrix< FE_TYPE > kernelM( finiteElement );
        kernelM.template computeMassMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                          nodeCoords,
                                                                          elemsToNodes,
                                                                          velocity,
                                                                          density,
                                                                          mass );

        AcousticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );
        kernelD.template computeVTIFletcherDampingMatrices< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                                       nodeCoords,
                                                                                       elemsToFaces,
                                                                                       facesToNodes,
                                                                                       facesDomainBoundaryIndicator,
                                                                                       freeSurfaceFaceIndicator,
                                                                                       lateralSurfaceFaceIndicator,
                                                                                       bottomSurfaceFaceIndicator,
                                                                                       velocity,
                                                                                       density,
                                                                                       vti_epsilon,
                                                                                       vti_delta,
                                                                                       vti_sigma,
                                                                                       damping_pp,
                                                                                       damping_pq,
                                                                                       damping_qp,
                                                                                       damping_qq );


      } );
    } );
  } );

  if( m_timestepStabilityLimit==1 )
  {
  // TODO: adapt to VTI
    GEOS_ERROR( "This option (Time Step computation) is not supported" );
/*    real64 dtOut = 0.0;
    computeTimeStep( dtOut );
    m_timestepStabilityLimit = 0;
    m_timeStep=dtOut;*/
    
  }

  WaveSolverUtils::initTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ), m_receiverIsLocal );
  m_seismoCoeff.resize( m_receiverIsLocal.size());
  m_seismoCoeff.setValues< EXEC_POLICY >( 0.5 );
}

//This function is only to give an easy accesss to the computation of the time-step for Pygeosx interface and avoid to exit the code when
// using Pygeosx

real64 AcousticVTIFletcherWaveEquationSEM::computeTimeStep( real64 & dtOut )
{
  // TODO: adapt to VTI
  GEOS_ERROR( "This option (Time Step computation) is not supported" );
  dtOut = 0.;

  return m_timeStep * m_cflFactor;

}


void AcousticVTIFletcherWaveEquationSEM::precomputeSurfaceFieldIndicator( DomainPartition & domain )
{
  real64 const time = 0.0;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on lateral surface; 0 otherwise
  arrayView1d< localIndex > const lateralSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticLateralSurfaceFaceIndicator >();
  /// array of indicators: 1 if a node is on on lateral surface; 0 otherwise
  arrayView1d< localIndex > const lateralSurfaceNodeIndicator = nodeManager.getField< acousticvtifields::AcousticLateralSurfaceNodeIndicator >();

  /// array of indicators: 1 if a face is on on bottom surface; 0 otherwise
  arrayView1d< localIndex > const bottomSurfaceFaceIndicator = faceManager.getField< acousticvtifields::AcousticBottomSurfaceFaceIndicator >();
  /// array of indicators: 1 if a node is on on bottom surface; 0 otherwise
  arrayView1d< localIndex > const bottomSurfaceNodeIndicator = nodeManager.getField< acousticvtifields::AcousticBottomSurfaceNodeIndicator >();

  // Lateral surfaces
  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "LateralSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        lateralSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          lateralSurfaceNodeIndicator[dof] = 1;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );

  // For the Bottom surface
  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "BottomSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        bottomSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          bottomSurfaceNodeIndicator[dof] = 1;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

void AcousticVTIFletcherWaveEquationSEM::applyFreeSurfaceBC( real64 time, DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticvtifields::Pressure_p_nm1 >();
  arrayView1d< real32 > const p_n = nodeManager.getField< acousticvtifields::Pressure_p_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticvtifields::Pressure_p_np1 >();

  arrayView1d< real32 > const q_nm1 = nodeManager.getField< acousticvtifields::Pressure_q_nm1 >();
  arrayView1d< real32 > const q_n = nodeManager.getField< acousticvtifields::Pressure_q_n >();
  arrayView1d< real32 > const q_np1 = nodeManager.getField< acousticvtifields::Pressure_q_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< acousticfields::AcousticFreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< acousticfields::AcousticFreeSurfaceNodeIndicator >();

  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "FreeSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      real64 const value = bc.getScale();

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        freeSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          freeSurfaceNodeIndicator[dof] = 1;

          p_np1[dof] = value;
          p_n[dof]   = value;
          p_nm1[dof] = value;

          q_np1[dof] = value;
          q_n[dof]   = value;
          q_nm1[dof] = value;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

void AcousticVTIFletcherWaveEquationSEM::initializePML()
{  GEOS_ERROR( "This option (PML) is not supported" );
  return;
}



void AcousticVTIFletcherWaveEquationSEM::applyPML( real64 const time, DomainPartition & domain )
{
  GEOS_UNUSED_VAR(time, domain);
  GEOS_MARK_FUNCTION;
#if 0
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  parametersPML const & param = getReference< parametersPML >( viewKeyStruct::parametersPMLString() );

  /// Loop over the different mesh bodies; for wave propagation, there is only one mesh body
  /// which is the whole mesh
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    NodeManager & nodeManager = mesh.getNodeManager();

    /// Array views of the pressure p, PML auxiliary variables, and node coordinates
    arrayView1d< real32 const > const p_n = nodeManager.getField< acousticfields::Pressure_n >();
    arrayView2d< real32 const > const v_n = nodeManager.getField< acousticfields::AuxiliaryVar1PML >();
    arrayView2d< real32 > const grad_n = nodeManager.getField< acousticfields::AuxiliaryVar2PML >();
    arrayView1d< real32 > const divV_n = nodeManager.getField< acousticfields::AuxiliaryVar3PML >();
    arrayView1d< real32 const > const u_n = nodeManager.getField< acousticfields::AuxiliaryVar4PML >();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// Select the subregions concerned by the PML (specified in the xml by the Field Specification)
    /// 'targetSet' contains the indices of the elements in a given subregion
    fsManager.apply< ElementSubRegionBase,
                     PerfectlyMatchedLayer >( time,
                                              mesh,
                                              PerfectlyMatchedLayer::catalogName(),
                                              [&]( PerfectlyMatchedLayer const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )

    {

      /// Get the element to nodes mapping in the subregion
      CellElementSubRegion::NodeMapType const & elemToNodes =
        subRegion.getReference< CellElementSubRegion::NodeMapType >( CellElementSubRegion::viewKeyStruct::nodeListString() );

      /// Get a const ArrayView of the mapping above
      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst = elemToNodes.toViewConst();

      /// Array view of the wave speed
      arrayView1d< real32 const > const vel = subRegion.getReference< array1d< real32 > >( acousticfields::AcousticVelocity::key());

      /// Get the object needed to determine the type of the element in the subregion
      finiteElement::FiniteElementBase const &
      fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      real32 xMin[3];
      real32 xMax[3];
      real32 dMin[3];
      real32 dMax[3];
      real32 cMin[3];
      real32 cMax[3];
      for( integer i=0; i<3; ++i )
      {
        xMin[i] = param.xMinPML[i];
        xMax[i] = param.xMaxPML[i];
        dMin[i] = param.thicknessMinXYZPML[i];
        dMax[i] = param.thicknessMaxXYZPML[i];
        cMin[i] = param.waveSpeedMinXYZPML[i];
        cMax[i] = param.waveSpeedMaxXYZPML[i];
      }
      real32 const r = param.reflectivityPML;

      /// Get the type of the elements in the subregion
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        /// apply the PML kernel

        AcousticPMLSEM::
          PMLKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( targetSet,
          nodeCoords32,
          elemToNodesViewConst,
          vel,
          p_n,
          v_n,
          u_n,
          xMin,
          xMax,
          dMin,
          dMax,
          cMin,
          cMax,
          r,
          grad_n,
          divV_n );
        
      } );
    } );
  } );
#endif
}

real64 AcousticVTIFletcherWaveEquationSEM::explicitStepForward( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer cycleNumber,
                                                     DomainPartition & domain,
                                                     bool computeGradient )
{
  real64 dtCompute = explicitStepInternal( time_n, dt, domain, true );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & GEOS_UNUSED_PARAM ( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticvtifields::Pressure_p_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< acousticvtifields::Pressure_p_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticvtifields::Pressure_p_np1 >();

    arrayView1d< real32 > const q_nm1 = nodeManager.getField< acousticvtifields::Pressure_q_nm1 >();
    arrayView1d< real32 > const q_n = nodeManager.getField< acousticvtifields::Pressure_q_n >();
    arrayView1d< real32 > const q_np1 = nodeManager.getField< acousticvtifields::Pressure_q_np1 >();

    if( computeGradient && cycleNumber >= 0 )
    {
      GEOS_ERROR( "This option is not supported yet" );
    }

    prepareNextTimestep( mesh );
  } );

  return dtCompute;
}


real64 AcousticVTIFletcherWaveEquationSEM::explicitStepBackward( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer cycleNumber,
                                                      DomainPartition & domain,
                                                      bool computeGradient )
{
  real64 dtCompute = explicitStepInternal( time_n, dt, domain, false );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & GEOS_UNUSED_PARAM( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getField< acousticfields::AcousticMassVector >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticvtifields::Pressure_p_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< acousticvtifields::Pressure_p_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticvtifields::Pressure_p_np1 >();

    arrayView1d< real32 > const q_nm1 = nodeManager.getField< acousticvtifields::Pressure_q_nm1 >();
    arrayView1d< real32 > const q_n = nodeManager.getField< acousticvtifields::Pressure_q_n >();
    arrayView1d< real32 > const q_np1 = nodeManager.getField< acousticvtifields::Pressure_q_np1 >();


    EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
    real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
    int const maxCycle = int(round( maxTime / dt ));
    GEOS_UNUSED_VAR(maxCycle);

    if( computeGradient && cycleNumber >= 0 )
    {
      GEOS_ERROR( "This option is not supported yet" );
    }

    prepareNextTimestep( mesh );
  } );

  return dtCompute;
}

void AcousticVTIFletcherWaveEquationSEM::prepareNextTimestep( MeshLevel & mesh )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticvtifields::Pressure_p_nm1 >();
  arrayView1d< real32 > const p_n   = nodeManager.getField< acousticvtifields::Pressure_p_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticvtifields::Pressure_p_np1 >();

  arrayView1d< real32 > const q_nm1 = nodeManager.getField< acousticvtifields::Pressure_q_nm1 >();
  arrayView1d< real32 > const q_n   = nodeManager.getField< acousticvtifields::Pressure_q_n >();
  arrayView1d< real32 > const q_np1 = nodeManager.getField< acousticvtifields::Pressure_q_np1 >();

  arrayView1d< real32 > const stiffnessVector_p = nodeManager.getField< acousticvtifields::StiffnessVector_p >();
  arrayView1d< real32 > const stiffnessVector_q = nodeManager.getField< acousticvtifields::StiffnessVector_q >();
  arrayView1d< real32 > const rhs = nodeManager.getField< acousticfields::ForcingRHS >();

  SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();

  forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
  {
    localIndex const a = solverTargetNodesSet[n];

    p_nm1[a] = p_n[a];
    p_n[a]   = p_np1[a];
    q_nm1[a] = q_n[a];
    q_n[a]   = q_np1[a];

    rhs[a] = 0.0;
    stiffnessVector_p[a] = 0.0;
    stiffnessVector_q[a] = 0.0;
  } );
}

void AcousticVTIFletcherWaveEquationSEM::computeUnknowns( real64 const & time_n,
                                               real64 const & dt,
                                               DomainPartition & GEOS_UNUSED_PARAM( domain),
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames,
                                               bool const isForward )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 const > const mass = nodeManager.getField< acousticfields::AcousticMassVector >();
  arrayView1d< real32 const > const damping_pp = nodeManager.getField< acousticvtifields::DampingVector_pp >();
  arrayView1d< real32 const > const damping_pq = nodeManager.getField< acousticvtifields::DampingVector_pq >();
  arrayView1d< real32 const > const damping_qp = nodeManager.getField< acousticvtifields::DampingVector_qp >();
  arrayView1d< real32 const > const damping_qq = nodeManager.getField< acousticvtifields::DampingVector_qq >();

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticvtifields::Pressure_p_nm1 >();
  arrayView1d< real32 > const p_n = nodeManager.getField< acousticvtifields::Pressure_p_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticvtifields::Pressure_p_np1 >();

  arrayView1d< real32 > const q_nm1 = nodeManager.getField< acousticvtifields::Pressure_q_nm1 >();
  arrayView1d< real32 > const q_n = nodeManager.getField< acousticvtifields::Pressure_q_n >();
  arrayView1d< real32 > const q_np1 = nodeManager.getField< acousticvtifields::Pressure_q_np1 >();

  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< acousticfields::AcousticFreeSurfaceNodeIndicator >();
  arrayView1d< localIndex const > const lateralSurfaceNodeIndicator = nodeManager.getField< acousticvtifields::AcousticLateralSurfaceNodeIndicator >();
  arrayView1d< localIndex const > const bottomSurfaceNodeIndicator = nodeManager.getField< acousticvtifields::AcousticBottomSurfaceNodeIndicator >();
  arrayView1d< real32 > const stiffnessVector_p = nodeManager.getField< acousticvtifields::StiffnessVector_p >();
  arrayView1d< real32 > const stiffnessVector_q = nodeManager.getField< acousticvtifields::StiffnessVector_q >();
  arrayView1d< real32 > const rhs = nodeManager.getField< acousticfields::ForcingRHS >();
if( isForward )
  {
    auto kernelFactory = acousticVTIFletcherWaveEquationSEMKernels::ExplicitAcousticVTIFletcherSEMFactory( dt );

    finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            getDiscretizationName(),
                                                            "",
                                                            kernelFactory );
  }
  else
  {
    auto kernelFactory = acousticVTIFletcherAdjointWaveEquationSEMKernels::ExplicitAcousticVTIFletcherAdjointSEMFactory( dt );

    finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            getDiscretizationName(),
                                                            "",
                                                            kernelFactory );

  }
  //Modification of cycleNember useful when minTime < 0
  addSourceToRightHandSide( time_n, rhs );


  SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();
  if( !m_usePML )
  {
    GEOS_MARK_SCOPE ( updateP );
    AcousticTimeSchemeSEM::LeapFrogWithoutPML( dt, p_np1, p_n, p_nm1, mass, stiffnessVector_p, damping_pp,
                                               rhs, freeSurfaceNodeIndicator, solverTargetNodesSet );
  }
  else
  {
    GEOS_ERROR( "This option is not supported yet" );
  }
}

void AcousticVTIFletcherWaveEquationSEM::synchronizeUnknowns( real64 const & time_n,
                                                   real64 const & dt,
                                                   DomainPartition & domain,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const p_n = nodeManager.getField< acousticvtifields::Pressure_p_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticvtifields::Pressure_p_np1 >();
  arrayView1d< real32 > const q_n = nodeManager.getField< acousticvtifields::Pressure_q_n >();
  arrayView1d< real32 > const q_np1 = nodeManager.getField< acousticvtifields::Pressure_q_np1 >();

  arrayView1d< real32 > const stiffnessVector_p = nodeManager.getField< acousticvtifields::StiffnessVector_p >();
  arrayView1d< real32 > const stiffnessVector_q = nodeManager.getField< acousticvtifields::StiffnessVector_q >();
  arrayView1d< real32 > const rhs = nodeManager.getField< acousticfields::ForcingRHS >();

  /// synchronize pressure fields
  FieldIdentifiers fieldsToBeSync;
  fieldsToBeSync.addFields( FieldLocation::Node, { acousticvtifields::Pressure_p_np1::key() } );
  fieldsToBeSync.addFields( FieldLocation::Node, { acousticvtifields::Pressure_q_np1::key() } );

  if( m_usePML )
  {
    GEOS_ERROR( "This option is not supported yet" );
  }


  CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldsToBeSync,
                                mesh,
                                domain.getNeighbors(),
                                true );
  /// compute the seismic traces since last step.
  arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();

  computeAllSeismoTraces( time_n, dt, p_np1, p_n, pReceivers, m_seismoCoeff.toView(), false );
  computeAllSeismoTraces( time_n, dt, q_np1, q_n, pReceivers, m_seismoCoeff.toView(), true );
  incrementIndexSeismoTrace( time_n );

  if( m_usePML )
  {
    GEOS_ERROR( "This option is not supported yet" );
  }
}

real64 AcousticVTIFletcherWaveEquationSEM::explicitStepInternal( real64 const & time_n,
                                                      real64 const & dt,
                                                      DomainPartition & domain,
                                                      bool const isForward )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  real64 dtCompute;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    localIndex nSubSteps = (int) ceil( dt/m_timeStep );
    dtCompute = dt/nSubSteps;
    computeUnknowns( time_n, dtCompute, domain, mesh, regionNames, isForward );
    synchronizeUnknowns( time_n, dtCompute, domain, mesh, regionNames );
  } );

  return dtCompute;
}

void AcousticVTIFletcherWaveEquationSEM::cleanup( real64 const time_n,
                                       integer const cycleNumber,
                                       integer const eventCounter,
                                       real64 const eventProgress,
                                       DomainPartition & domain )
{
  // call the base class cleanup (for reporting purposes)
  PhysicsSolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    arrayView1d< real32 const > const p_n = nodeManager.getField< acousticvtifields::Pressure_p_n >();
    arrayView1d< real32 const > const p_np1 = nodeManager.getField< acousticvtifields::Pressure_p_np1 >();
    arrayView1d< real32 const > const q_n = nodeManager.getField< acousticvtifields::Pressure_q_n >();
    arrayView1d< real32 const > const q_np1 = nodeManager.getField< acousticvtifields::Pressure_q_np1 >();
    arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0.0, p_np1, p_n, pReceivers, m_seismoCoeff.toView(), false );
    computeAllSeismoTraces( time_n, 0.0, q_np1, q_n, pReceivers, m_seismoCoeff.toView(), true );

    WaveSolverUtils::writeSeismoTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                       m_receiverIsLocal, m_nsamplesSeismoTrace, pReceivers );
  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticVTIFletcherWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
