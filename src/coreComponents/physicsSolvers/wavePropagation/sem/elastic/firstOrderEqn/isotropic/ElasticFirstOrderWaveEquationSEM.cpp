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
 * @file ElasticFirstOrderWaveEquationSEM.cpp
 */

#include "ElasticFirstOrderWaveEquationSEM.hpp"
#include "ElasticFirstOrderWaveEquationSEMKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/wavePropagation/sem/elastic/shared/ElasticMatricesSEMKernel.hpp"
#include "physicsSolvers/wavePropagation/shared/PrecomputeSourcesAndReceiversKernel.hpp"
#include "events/EventManager.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

ElasticFirstOrderWaveEquationSEM::ElasticFirstOrderWaveEquationSEM( const std::string & name,
                                                                    Group * const parent ):

  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::displacementxNp1AtReceiversString(), &m_displacementxNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (x-components)" );

  registerWrapper( viewKeyStruct::displacementyNp1AtReceiversString(), &m_displacementyNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (y-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_displacementzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::sigmaxxNp1AtReceiversString(), &m_sigmaxxNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::sigmayyNp1AtReceiversString(), &m_sigmayyNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::sigmazzNp1AtReceiversString(), &m_sigmazzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::sigmaxyNp1AtReceiversString(), &m_sigmaxyNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::sigmaxzNp1AtReceiversString(), &m_sigmaxzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::sigmayzNp1AtReceiversString(), &m_sigmayzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::sourceElemString(), &m_sourceElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the sources" );

  registerWrapper( viewKeyStruct::sourceRegionString(), &m_sourceRegion ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Region containing the sources" );
}

ElasticFirstOrderWaveEquationSEM::~ElasticFirstOrderWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void ElasticFirstOrderWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void ElasticFirstOrderWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )

  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< elasticfields::Displacementx_np1,
                               elasticfields::Displacementy_np1,
                               elasticfields::Displacementz_np1,
                               elasticfields::ForcingRHS,
                               elasticfields::ElasticMassVector,
                               elasticfields::DampingVectorx,
                               elasticfields::DampingVectory,
                               elasticfields::DampingVectorz,
                               elasticfields::ElasticFreeSurfaceNodeIndicator >( getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< elasticfields::ElasticFreeSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< elasticfields::ElasticVelocityVp >( getName() );
      subRegion.registerField< elasticfields::ElasticVelocityVs >( getName() );
      subRegion.registerField< elasticfields::ElasticDensity >( getName() );
      subRegion.registerField< elasticfields::Lambda >( getName() );
      subRegion.registerField< elasticfields::Mu >( getName() );

      subRegion.registerField< elasticfields::Stresstensorxx >( getName());
      subRegion.registerField< elasticfields::Stresstensoryy >( getName());
      subRegion.registerField< elasticfields::Stresstensorzz >( getName());
      subRegion.registerField< elasticfields::Stresstensorxy >( getName());
      subRegion.registerField< elasticfields::Stresstensorxz >( getName());
      subRegion.registerField< elasticfields::Stresstensoryz >( getName());

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {

        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

        subRegion.getField< elasticfields::Stresstensorxx >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< elasticfields::Stresstensoryy >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< elasticfields::Stresstensorzz >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< elasticfields::Stresstensorxy >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< elasticfields::Stresstensorxz >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< elasticfields::Stresstensoryz >().resizeDimension< 1 >( numNodesPerElem );
      } );


    } );

  } );
}



void ElasticFirstOrderWaveEquationSEM::postInputInitialization()
{

  WaveSolverBase::postInputInitialization();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceElem.resize( numSourcesGlobal );
  m_sourceRegion.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverElem.resize( numReceiversGlobal );
  m_receiverRegion.resize( numReceiversGlobal );

  m_displacementxNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_displacementyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_displacementzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );

  m_sigmaxxNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_sigmayyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_sigmazzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_sigmaxyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_sigmaxzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_sigmayzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
}


void ElasticFirstOrderWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< globalIndex const > const nodeLocalToGlobal = baseMesh.getNodeManager().localToGlobalMap().toViewConst();
  ArrayOfArraysView< localIndex const > const nodesToElements = baseMesh.getNodeManager().elementList().toViewConst();
  ArrayOfArraysView< localIndex const > const facesToNodes = baseMesh.getFaceManager().nodeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodeCoords = baseMesh.getNodeManager().referencePosition();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex > const sourceElem = m_sourceElem.toView();
  arrayView1d< localIndex > const sourceRegion = m_sourceRegion.toView();
  sourceNodeIds.setValues< serialPolicy >( -1 );
  sourceConstants.setValues< serialPolicy >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  arrayView1d< localIndex > const receiverElem = m_receiverElem.toView();
  arrayView1d< localIndex > const receiverRegion = m_receiverRegion.toView();

  receiverNodeIds.setValues< serialPolicy >( -1 );
  receiverConstants.setValues< serialPolicy >( -1 );
  receiverIsLocal.zero();

  mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                localIndex const regionIndex,
                                                                                                localIndex const esr,
                                                                                                ElementRegionBase &,
                                                                                                CellElementSubRegion & elementSubRegion )
  {

    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                   getDataContext() << ": Invalid type of element, the elastic solver is designed for hexahedral meshes only (C3D8) ",
                   InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes = baseMesh.getElemManager().getRegion( regionIndex ).getSubRegion< CellElementSubRegion >( esr ).nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();
    arrayView1d< globalIndex const > const elemLocalToGlobal = elementSubRegion.localToGlobalMap().toViewConst();

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      PreComputeSourcesAndReceivers::
        Compute1DSourceAndReceiverConstantsWithElementsAndRegionStorage
      < EXEC_POLICY, FE_TYPE >
        ( elementSubRegion.size(),
        regionIndex,
        facesToNodes,
        nodeCoords,
        nodeLocalToGlobal,
        elemLocalToGlobal,
        nodesToElements,
        baseElemsToNodes,
        elemGhostRank,
        elemsToNodes,
        elemsToFaces,
        elemCenter,
        sourceCoordinates,
        sourceIsAccessible,
        sourceElem,
        sourceNodeIds,
        sourceConstants,
        sourceRegion,
        receiverCoordinates,
        receiverIsLocal,
        receiverElem,
        receiverNodeIds,
        receiverConstants,
        receiverRegion );
    } );
    elementSubRegion.faceList().freeOnDevice();
    baseMesh.getElemManager().getRegion( regionIndex ).getSubRegion< CellElementSubRegion >( esr ).nodeList().freeOnDevice();
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

void ElasticFirstOrderWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{

  WaveSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    MeshLevel & baseMesh = domain.getMeshBodies().getGroup< MeshBody >( meshBodyName ).getBaseDiscretization();
    precomputeSourceAndReceiverTerm( baseMesh, mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();
    arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< elasticfields::ElasticMassVector >();
    mass.zero();
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const dampingx = nodeManager.getField< elasticfields::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getField< elasticfields::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getField< elasticfields::DampingVectorz >();
    dampingx.zero();
    dampingy.zero();
    dampingz.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< elasticfields::ElasticFreeSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
      arrayView1d< real32 const > const density = elementSubRegion.getField< elasticfields::ElasticDensity >();
      arrayView1d< real32 const > const velocityVp = elementSubRegion.getField< elasticfields::ElasticVelocityVp >();
      arrayView1d< real32 const > const velocityVs = elementSubRegion.getField< elasticfields::ElasticVelocityVs >();
      arrayView1d< real32 > const lambda = elementSubRegion.getField< elasticfields::Lambda >();
      arrayView1d< real32 > const mu = elementSubRegion.getField< elasticfields::Mu >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe,
                                                                               [&]
                                                                                 ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        ElasticMatricesSEM::MassMatrix< FE_TYPE > kernelM( finiteElement );

        kernelM.template computeMassMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                          nodeCoords,
                                                                          elemsToNodes,
                                                                          density,
                                                                          mass );


        ElasticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );

        kernelD.template computeDampingMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                             nodeCoords,
                                                                             elemsToFaces,
                                                                             facesToNodes,
                                                                             facesDomainBoundaryIndicator,
                                                                             freeSurfaceFaceIndicator,
                                                                             faceNormal,
                                                                             density,
                                                                             velocityVp,
                                                                             velocityVs,
                                                                             dampingx,
                                                                             dampingy,
                                                                             dampingz );
      } );
    } );
  } );
}

real64 ElasticFirstOrderWaveEquationSEM::computeTimeStep( real64 & dtOut )
{
  GEOS_ERROR( getDataContext() << ":  Time-Step computation for the first order elastic wave propagator not yet implemented" );
  return dtOut;
}

void ElasticFirstOrderWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< elasticfields::ElasticFreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< elasticfields::ElasticFreeSurfaceNodeIndicator >();

  freeSurfaceFaceIndicator.zero();
  freeSurfaceNodeIndicator.zero();

  fsManager.apply( time,
                   domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                   WaveSolverBase::viewKeyStruct::freeSurfaceString(),
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group &,
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

          ux_np1[dof] = value;
          uy_np1[dof] = value;
          uz_np1[dof] = value;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

real64 ElasticFirstOrderWaveEquationSEM::explicitStepForward( real64 const & time_n,
                                                              real64 const & dt,
                                                              integer cycleNumber,
                                                              DomainPartition & domain,
                                                              bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}



real64 ElasticFirstOrderWaveEquationSEM::explicitStepBackward( real64 const & time_n,
                                                               real64 const & dt,
                                                               integer cycleNumber,
                                                               DomainPartition & domain,
                                                               bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  GEOS_ERROR( getDataContext() << ": Backward propagation for the first order elastic wave propagator not yet implemented" );
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}

real64 ElasticFirstOrderWaveEquationSEM::explicitStepInternal( real64 const & time_n,
                                                               real64 const & dt,
                                                               integer const cycleNumber,
                                                               DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, cycleNumber );


  arrayView2d< real64 const > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex const > const sourceElem = m_sourceElem.toView();
  arrayView1d< localIndex const > const sourceRegion = m_sourceRegion.toView();


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

  {

    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    arrayView1d< real32 const > const mass = nodeManager.getField< elasticfields::ElasticMassVector >();
    arrayView1d< real32 > const dampingx = nodeManager.getField< elasticfields::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getField< elasticfields::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getField< elasticfields::DampingVectorz >();


    arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      arrayView1d< real32 const > const velocityVp = elementSubRegion.getField< elasticfields::ElasticVelocityVp >();
      arrayView1d< real32 const > const velocityVs = elementSubRegion.getField< elasticfields::ElasticVelocityVs >();
      arrayView1d< real32 const > const density = elementSubRegion.getField< elasticfields::ElasticDensity >();

      arrayView1d< real32 > const lambda = elementSubRegion.getField< elasticfields::Lambda >();
      arrayView1d< real32 > const mu = elementSubRegion.getField< elasticfields::Mu >();

      arrayView2d< real32 > const stressxx = elementSubRegion.getField< elasticfields::Stresstensorxx >();
      arrayView2d< real32 > const stressyy = elementSubRegion.getField< elasticfields::Stresstensoryy >();
      arrayView2d< real32 > const stresszz = elementSubRegion.getField< elasticfields::Stresstensorzz >();
      arrayView2d< real32 > const stressxy = elementSubRegion.getField< elasticfields::Stresstensorxy >();
      arrayView2d< real32 > const stressxz = elementSubRegion.getField< elasticfields::Stresstensorxz >();
      arrayView2d< real32 > const stressyz = elementSubRegion.getField< elasticfields::Stresstensoryz >();


      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        //Modification of cycleNember useful when minTime < 0

        elasticFirstOrderWaveEquationSEMKernels::
          StressComputation< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          regionIndex,
          nodeCoords,
          elemsToNodes,
          ux_np1,
          uy_np1,
          uz_np1,
          density,
          velocityVp,
          velocityVs,
          lambda,
          mu,
          sourceConstants,
          sourceIsAccessible,
          sourceElem,
          sourceRegion,
          dt,
          time_n,
          m_timeSourceFrequency,
          m_timeSourceDelay,
          m_rickerOrder,
          m_useSourceWaveletTables,
          m_sourceWaveletTableWrappers,
          stressxx,
          stressyy,
          stresszz,
          stressxy,
          stressxz,
          stressyz );

        elasticFirstOrderWaveEquationSEMKernels::
          VelocityComputation< FE_TYPE > kernel2( finiteElement );
        kernel2.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          nodeManager.size(),
          nodeCoords,
          elemsToNodes,
          stressxx,
          stressyy,
          stresszz,
          stressxy,
          stressxz,
          stressyz,
          mass,
          dampingx,
          dampingy,
          dampingz,
          dt,
          ux_np1,
          uy_np1,
          uz_np1 );


      } );

      arrayView2d< real32 > const sigmaxxReceivers   = m_sigmaxxNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayyReceivers   = m_sigmayyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmazzReceivers   = m_sigmazzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxyReceivers   = m_sigmaxyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxzReceivers   = m_sigmaxzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayzReceivers   = m_sigmayzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressxx, stressxx, sigmaxxReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressyy, stressyy, sigmayyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stresszz, stresszz, sigmazzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressxy, stressxy, sigmaxyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressxz, stressxz, sigmaxzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressyz, stressyz, sigmayzReceivers );
    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { elasticfields::Displacementx_np1::key(), elasticfields::Displacementy_np1::key(), elasticfields::Displacementz_np1::key()} );
    fieldsToBeSync.addElementFields( {elasticfields::Stresstensorxx::key(), elasticfields::Stresstensoryy::key(), elasticfields::Stresstensorzz::key(),
                                      elasticfields::Stresstensorxy::key(),
                                      elasticfields::Stresstensorxz::key(), elasticfields::Stresstensoryz::key()}, regionNames );


    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const uxReceivers   = m_displacementxNp1AtReceivers.toView();
    arrayView2d< real32 > const uyReceivers   = m_displacementyNp1AtReceivers.toView();
    arrayView2d< real32 > const uzReceivers   = m_displacementzNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, dt, ux_np1, ux_np1, uxReceivers );
    computeAllSeismoTraces( time_n, dt, uy_np1, uy_np1, uyReceivers );
    computeAllSeismoTraces( time_n, dt, uz_np1, uz_np1, uzReceivers );

    incrementIndexSeismoTrace( time_n );
  } );

  return dt;

}

void ElasticFirstOrderWaveEquationSEM::cleanup( real64 const time_n,
                                                integer const cycleNumber,
                                                integer const eventCounter,
                                                real64 const eventProgress,
                                                DomainPartition & domain )
{
  PhysicsSolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< real32 const > const stressxx = elementSubRegion.getField< elasticfields::Stresstensorxx >();
      arrayView2d< real32 const > const stressyy = elementSubRegion.getField< elasticfields::Stresstensoryy >();
      arrayView2d< real32 const > const stresszz = elementSubRegion.getField< elasticfields::Stresstensorzz >();
      arrayView2d< real32 const > const stressxy = elementSubRegion.getField< elasticfields::Stresstensorxy >();
      arrayView2d< real32 const > const stressxz = elementSubRegion.getField< elasticfields::Stresstensorxz >();
      arrayView2d< real32 const > const stressyz = elementSubRegion.getField< elasticfields::Stresstensoryz >();

      arrayView2d< real32 > const sigmaxxReceivers   = m_sigmaxxNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayyReceivers   = m_sigmayyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmazzReceivers   = m_sigmazzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxyReceivers   = m_sigmaxyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxzReceivers   = m_sigmaxzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayzReceivers   = m_sigmayzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0.0, stressxx, stressxx, sigmaxxReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0.0, stressyy, stressyy, sigmayyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0.0, stresszz, stresszz, sigmazzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0.0, stressxy, stressxy, sigmaxyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0.0, stressxz, stressxz, sigmaxzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0.0, stressyz, stressyz, sigmayzReceivers );

      WaveSolverUtils::writeSeismoTraceVector( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                               m_receiverIsLocal, m_nsamplesSeismoTrace, sigmaxxReceivers, sigmayyReceivers, sigmazzReceivers );
      WaveSolverUtils::writeSeismoTraceVector( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                               m_receiverIsLocal, m_nsamplesSeismoTrace, sigmaxyReceivers, sigmaxzReceivers, sigmayzReceivers );

    } );
    arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

    // compute the seismic traces since last step.
    arrayView2d< real32 > const uxReceivers = m_displacementxNp1AtReceivers.toView();
    arrayView2d< real32 > const uyReceivers = m_displacementyNp1AtReceivers.toView();
    arrayView2d< real32 > const uzReceivers = m_displacementzNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, 0.0, ux_np1, ux_np1, uxReceivers );
    computeAllSeismoTraces( time_n, 0.0, uy_np1, uy_np1, uyReceivers );
    computeAllSeismoTraces( time_n, 0.0, uz_np1, uz_np1, uzReceivers );

    WaveSolverUtils::writeSeismoTraceVector( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                             m_receiverIsLocal, m_nsamplesSeismoTrace, uxReceivers, uyReceivers, uzReceivers );
  } );
}

void ElasticFirstOrderWaveEquationSEM::initializePML()
{
  GEOS_ERROR( getDataContext() << ": PML for the first order elastic wave propagator not yet implemented" );
}

void ElasticFirstOrderWaveEquationSEM::applyPML( real64 const, DomainPartition & )
{
  GEOS_ERROR( getDataContext() << ": PML for the first order elastic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ElasticFirstOrderWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
