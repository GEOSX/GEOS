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
 * @file EmbeddedSurfaceGenerator.cpp
 */

#include "EmbeddedSurfaceGenerator.hpp"
#include "EmbeddedSurfacesParallelSynchronization.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshFields.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/utilities/CIcomputationKernel.hpp"
#include "physicsSolvers/solidMechanics/kernels/SolidMechanicsLagrangianFEMKernels.hpp"
#include "mesh/simpleGeometricObjects/GeometricObjectManager.hpp"
#include "mesh/simpleGeometricObjects/Rectangle.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/surfaceGeneration/LogLevelsInfo.hpp"


namespace geos
{
using namespace dataRepository;
using namespace constitutive;

void NewObjectLists::insert( NewObjectLists const & newObjects )
{
  newNodes.insert( newObjects.newNodes.begin(),
                   newObjects.newNodes.end() );

  newEdges.insert( newObjects.newEdges.begin(),
                   newObjects.newEdges.end() );

  for( auto & iter : newObjects.newElements )
  {
    std::pair< localIndex, localIndex > const & key = iter.first;
    std::set< localIndex > const & values = iter.second;
    newElements[key].insert( values.begin(), values.end() );
  }

}


EmbeddedSurfaceGenerator::EmbeddedSurfaceGenerator( const string & name,
                                                    Group * const parent ):
  SolverBase( name, parent ),
  m_fractureRegionName(),
  m_mpiCommOrder( 0 )
{
  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( "FractureRegion" );

  registerWrapper( viewKeyStruct::targetObjectsNameString(), &m_targetObjectsName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "List of geometric objects that will be used to initialized the embedded surfaces/fractures." );

  registerWrapper( viewKeyStruct::mpiCommOrderString(), &m_mpiCommOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to enable MPI consistent communication ordering" );
}

EmbeddedSurfaceGenerator::~EmbeddedSurfaceGenerator()
{}

void EmbeddedSurfaceGenerator::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getBaseDiscretization();

    EmbeddedSurfaceNodeManager & nodeManager = meshLevel.getEmbSurfNodeManager();

    nodeManager.registerField< fields::parentEdgeIndex >( this->getName() );
  } );
}


void EmbeddedSurfaceGenerator::initializePostSubGroups()
{
  /*
   * Here we generate embedded elements for fractures (or faults) that already exist in the domain and
   * were specified in the input file.
   */

  // Get domain
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // Get geometric object manager
  GeometricObjectManager & geometricObjManager = GeometricObjectManager::getInstance();

  // Get meshLevel
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();

  // Get managers
  ElementRegionManager & elemManager = meshLevel.getElemManager();
  NodeManager & nodeManager = meshLevel.getNodeManager();
  EmbeddedSurfaceNodeManager & embSurfNodeManager = meshLevel.getEmbSurfNodeManager();
  EdgeManager & edgeManager = meshLevel.getEdgeManager();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager.referencePosition();

  // Get EmbeddedSurfaceSubRegions
  SurfaceElementRegion & embeddedSurfaceRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );
  EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion = embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  localIndex localNumberOfSurfaceElems         = 0;

  NewObjectLists newObjects;

  // Loop over all the fracture planes
  geometricObjManager.forGeometricObject< PlanarGeometricObject >( m_targetObjectsName, [&]( localIndex const,
                                                                                             PlanarGeometricObject & fracture )
  {
    /* 1. Find out if an element is cut by the fracture or not.
     * Loop over all the elements and for each one of them loop over the nodes and compute the
     * dot product between the distance between the plane center and the node and the normal
     * vector defining the plane. If two scalar products have different signs the plane cuts the
     * cell. If a nodes gives a 0 dot product it has to be neglected or the method won't work.
     */
    real64 const planeCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracture.getCenter() );
    real64 const normalVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracture.getNormal() );
    // Initialize variables
    globalIndex nodeIndex;
    integer isPositive, isNegative;
    real64 distVec[ 3 ];

    elemManager.forElementSubRegionsComplete< CellElementSubRegion >(
      [&]( localIndex const er, localIndex const esr, ElementRegionBase &, CellElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const cellToNodes = subRegion.nodeList();
      FixedOneToManyRelation const & cellToEdges = subRegion.edgeList();

      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      forAll< serialPolicy >( subRegion.size(), [ &, ghostRank,
                                                  cellToNodes,
                                                  nodesCoord ] ( localIndex const cellIndex )
      {
        if( ghostRank[cellIndex] < 0 )
        {
          isPositive = 0;
          isNegative = 0;
          for( localIndex kn = 0; kn < subRegion.numNodesPerElement(); kn++ )
          {
            nodeIndex = cellToNodes[cellIndex][kn];
            LvArray::tensorOps::copy< 3 >( distVec, nodesCoord[nodeIndex] );
            LvArray::tensorOps::subtract< 3 >( distVec, planeCenter );
            // check if the dot product is zero
            if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) > 0 )
            {
              isPositive = 1;
            }
            else if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) < 0 )
            {
              isNegative = 1;
            }
          } // end loop over nodes
          if( isPositive * isNegative == 1 )
          {
            bool added = embeddedSurfaceSubRegion.addNewEmbeddedSurface( cellIndex,
                                                                         er,
                                                                         esr,
                                                                         nodeManager,
                                                                         embSurfNodeManager,
                                                                         edgeManager,
                                                                         cellToEdges,
                                                                         &fracture );

            if( added )
            {
              GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::SurfaceGenerator, GEOS_FMT( "{}: element {} is fractured", subRegion.getName(), cellIndex ) );

              // Add the information to the CellElementSubRegion
              subRegion.addFracturedElement( cellIndex, localNumberOfSurfaceElems );

              newObjects.newElements[ {embeddedSurfaceRegion.getIndexInParent(), embeddedSurfaceSubRegion.getIndexInParent()} ].insert( localNumberOfSurfaceElems );

              localNumberOfSurfaceElems++;
            }
          }
        }
      } );// end loop over cells
    } );// end loop over subregions
  } );// end loop over planes

  // Launch kernel to compute connectivity index of each fractured element.
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >(
    [&]( localIndex const, localIndex const, ElementRegionBase &, CellElementSubRegion & subRegion )
  {
    finiteElement::FiniteElementBase & subRegionFE = subRegion.template getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

    finiteElement::dispatchlowOrder3D( subRegionFE, [&] ( auto const finiteElement )
    {
      using FE_TYPE = decltype( finiteElement );

      auto kernel = CIcomputationKernel< FE_TYPE >( finiteElement,
                                                    nodeManager,
                                                    subRegion,
                                                    embeddedSurfaceSubRegion );

      using KERNEL_TYPE = decltype( kernel );

      KERNEL_TYPE::template launchCIComputationKernel< parallelDevicePolicy< >, KERNEL_TYPE >( kernel );
    } );
  } );

  // add all new nodes to newObject list
  for( localIndex ni = 0; ni < embSurfNodeManager.size(); ni++ )
  {
    newObjects.newNodes.insert( ni );
  }

  // Set the ghostRank form the parent cell
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const & cellElemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  embeddedSurfaceSubRegion.inheritGhostRank( cellElemGhostRank );

  setGlobalIndices( elemManager, embSurfNodeManager, embeddedSurfaceSubRegion );

  embeddedSurfacesParallelSynchronization::sychronizeTopology( meshLevel,
                                                               domain.getNeighbors(),
                                                               newObjects,
                                                               m_mpiCommOrder,
                                                               this->m_fractureRegionName );

  addEmbeddedElementsToSets( elemManager, embeddedSurfaceSubRegion );

  EmbeddedSurfaceSubRegion::NodeMapType & embSurfToNodeMap = embeddedSurfaceSubRegion.nodeList();

  // Populate EdgeManager for embedded surfaces.
  EdgeManager & embSurfEdgeManager = meshLevel.getEmbSurfEdgeManager();

  EmbeddedSurfaceSubRegion::EdgeMapType & embSurfToEdgeMap = embeddedSurfaceSubRegion.edgeList();

  localIndex numOfPoints = embSurfNodeManager.size();

  // Create the edges
  embSurfEdgeManager.buildEdges( numOfPoints, embSurfToNodeMap.toViewConst(), embSurfToEdgeMap );
  // Node to cell map
  embSurfNodeManager.setElementMaps( elemManager );
  // Node to edge map
  embSurfNodeManager.setEdgeMaps( embSurfEdgeManager );
  embSurfNodeManager.compressRelationMaps();
}

void EmbeddedSurfaceGenerator::initializePostInitialConditionsPreSubGroups()
{}

real64 EmbeddedSurfaceGenerator::solverStep( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                             real64 const & GEOS_UNUSED_PARAM( dt ),
                                             const int GEOS_UNUSED_PARAM( cycleNumber ),
                                             DomainPartition & domain )
{
  real64 rval = 0;
  /*
   * This should be the method that generates new fracture elements based on the propagation criterion of choice.
   */
  // Add the embedded elements to the fracture stencil.
  addToFractureStencil( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );

    EmbeddedSurfaceSubRegion & fractureSubRegion = fractureRegion.getUniqueSubRegion< EmbeddedSurfaceSubRegion >();

    // Compute gravity coefficient for new elements so that gravity term is correctly computed
    real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    if( fractureSubRegion.hasField< fields::flow::gravityCoefficient >() )
    {
      arrayView2d< real64 const > const elemCenter = fractureSubRegion.getElementCenter();

      arrayView1d< real64 > const gravityCoef = fractureSubRegion.getField< fields::flow::gravityCoefficient >();

      forAll< parallelHostPolicy >( fractureSubRegion.size(), [=] ( localIndex const ei )
      {
        gravityCoef[ ei ] = LvArray::tensorOps::AiBi< 3 >( elemCenter[ ei ], gravVector );
      } );
    }
  } );

  return rval;
}

void EmbeddedSurfaceGenerator::addToFractureStencil( DomainPartition & domain )
{
  // Add embedded elements to the fracture Stencil
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getBaseDiscretization();

    for( localIndex a=0; a<fvManager.numSubGroups(); ++a )
    {
      FluxApproximationBase * const fluxApprox = fvManager.getGroupPointer< FluxApproximationBase >( a );
      if( fluxApprox!=nullptr )
      {
        fluxApprox->addEmbeddedFracturesToStencils( meshLevel, this->m_fractureRegionName );
      }
    }
  }

}

void EmbeddedSurfaceGenerator::setGlobalIndices( ElementRegionManager & elemManager,
                                                 EmbeddedSurfaceNodeManager & embSurfNodeManager,
                                                 EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )
{
  // Add new globalIndices
  int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOS );
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOS );

  localIndex_array numberOfSurfaceElemsPerRank( commSize );
  localIndex_array globalIndexOffset( commSize );
  MpiWrapper::allGather( embeddedSurfaceSubRegion.size(), numberOfSurfaceElemsPerRank );

  globalIndexOffset[0] = 0; // offSet for the globalIndex
  localIndex totalNumberOfSurfaceElements = numberOfSurfaceElemsPerRank[ 0 ];  // Sum across all ranks
  for( int rank = 1; rank < commSize; ++rank )
  {
    globalIndexOffset[rank] = globalIndexOffset[rank - 1] + numberOfSurfaceElemsPerRank[rank - 1];
    totalNumberOfSurfaceElements += numberOfSurfaceElemsPerRank[rank];
  }

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::SurfaceGenerator, "Number of embedded surface elements: " << totalNumberOfSurfaceElements );

  arrayView1d< globalIndex > const & elemLocalToGlobal = embeddedSurfaceSubRegion.localToGlobalMap();

  forAll< serialPolicy >( embeddedSurfaceSubRegion.size(), [&, elemLocalToGlobal] ( localIndex const ei )
  {
    elemLocalToGlobal( ei ) = ei + globalIndexOffset[ thisRank ] + elemManager.maxGlobalIndex() + 1;
    embeddedSurfaceSubRegion.updateGlobalToLocalMap( ei );
  } );

  embeddedSurfaceSubRegion.setMaxGlobalIndex();

  elemManager.setMaxGlobalIndex();

  // Nodes global indices
  localIndex_array numberOfNodesPerRank( commSize );
  MpiWrapper::allGather( embSurfNodeManager.size(), numberOfNodesPerRank );

  globalIndexOffset[0] = 0; // offSet for the globalIndex
  for( int rank = 1; rank < commSize; ++rank )
  {
    globalIndexOffset[rank] = globalIndexOffset[rank - 1] + numberOfNodesPerRank[rank - 1];
  }

  arrayView1d< globalIndex > const & nodesLocalToGlobal = embSurfNodeManager.localToGlobalMap();
  forAll< serialPolicy >( embSurfNodeManager.size(), [=, &embSurfNodeManager, &globalIndexOffset] ( localIndex const ni )
  {
    nodesLocalToGlobal( ni ) = ni + globalIndexOffset[ thisRank ];
    embSurfNodeManager.updateGlobalToLocalMap( ni );
  } );
}

void EmbeddedSurfaceGenerator::addEmbeddedElementsToSets( ElementRegionManager const & elemManager,
                                                          EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )
{
  // We want to create the sets for the embeddedSurfaceSubRegion which was empty when they
  // were created for the other subRegions. This way, for example, if the parent cell belongs to the set "source"
  // the embedded element will belong to the same set.

  dataRepository::Group & setGroupEmbSurf =
    embeddedSurfaceSubRegion.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    dataRepository::Group const & setGroupCell = subRegion.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );

    SortedArrayView< localIndex const > const fracturedElements = subRegion.fracturedElementsList();

    ArrayOfArraysView< localIndex const > const cellToEmbSurf = subRegion.embeddedSurfacesList().toViewConst();

    forAll< serialPolicy >( fracturedElements.size(), [&,
                                                       fracturedElements,
                                                       cellToEmbSurf ] ( localIndex const ei )
    {
      localIndex const cellIndex    = fracturedElements[ei];
      localIndex const embSurfIndex = cellToEmbSurf[cellIndex][0];
      setGroupCell.forWrappers< SortedArray< localIndex > >( [&]( auto const & wrapper )
      {
        SortedArrayView< const localIndex > const & targetSetCell = wrapper.reference();
        targetSetCell.move( hostMemorySpace );

        SortedArray< localIndex > & targetSetEmbSurf =
          setGroupEmbSurf.getWrapper< SortedArray< localIndex > >( wrapper.getName() ).reference();
        if( targetSetCell.contains( cellIndex ) )
        {
          targetSetEmbSurf.insert( embSurfIndex );
        }
      } );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        EmbeddedSurfaceGenerator,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
