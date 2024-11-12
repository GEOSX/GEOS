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
 * @file MeshLevel.cpp
 */

#include "MeshLevel.hpp"
#include "LvArray/src/memcpy.hpp"

#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"

namespace geos
{
using namespace dataRepository;

MeshLevel::MeshLevel( string const & name,
                      Group * const parent ):
  Group( name, parent ),
  m_nodeManager( new NodeManager( groupStructKeys::nodeManagerString(), this ) ),
  m_particleManager( new ParticleManager( groupStructKeys::particleManagerString(), this ) ),
  m_edgeManager( new EdgeManager( groupStructKeys::edgeManagerString(), this ) ),
  m_faceManager( new FaceManager( groupStructKeys::faceManagerString(), this ) ),
  m_elementManager( new ElementRegionManager( groupStructKeys::elemManagerString(), this ) ),
  m_embSurfNodeManager( new EmbeddedSurfaceNodeManager( groupStructKeys::embSurfNodeManagerString, this ) ),
  m_embSurfEdgeManager( new EdgeManager( groupStructKeys::embSurfEdgeManagerString, this ) ),
  m_modificationTimestamp( 0 ),
  m_isShallowCopy( false ),
  m_shallowParent( nullptr )
{

  registerGroup( groupStructKeys::nodeManagerString(), m_nodeManager );

  registerGroup( groupStructKeys::particleManagerString(), m_particleManager );

  registerGroup( groupStructKeys::edgeManagerString(), m_edgeManager );


  registerGroup< FaceManager >( groupStructKeys::faceManagerString(), m_faceManager );
  m_faceManager->nodeList().setRelatedObject( *m_nodeManager );


  registerGroup< ElementRegionManager >( groupStructKeys::elemManagerString(), m_elementManager );

  registerGroup< EdgeManager >( groupStructKeys::embSurfEdgeManagerString, m_embSurfEdgeManager );

  registerGroup< EmbeddedSurfaceNodeManager >( groupStructKeys::embSurfNodeManagerString, m_embSurfNodeManager );

  registerWrapper< integer >( viewKeys.meshLevel );

  // increment the modification timestamp at mesh level creation
  // this is to make sure that the actions that depend on this timestamp (such as system setup) are performed at the beginning of the
  // simulations
  modified();
}


MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel & source ):
  Group( name, parent ),
  m_nodeManager( source.m_nodeManager ),
  m_particleManager( source.m_particleManager ),
  m_edgeManager( source.m_edgeManager ),
  m_faceManager( source.m_faceManager ),
  m_elementManager( source.m_elementManager ),
  m_embSurfNodeManager( source.m_embSurfNodeManager ),
  m_embSurfEdgeManager( source.m_embSurfEdgeManager ),
  m_modificationTimestamp( 0 ),
  m_isShallowCopy( true ),
  m_shallowParent( &source )
{
  this->setRestartFlags( RestartFlags::NO_WRITE );

  registerGroup( groupStructKeys::nodeManagerString(), m_nodeManager );

  registerGroup( groupStructKeys::edgeManagerString(), m_edgeManager );


  registerGroup< FaceManager >( groupStructKeys::faceManagerString(), m_faceManager );
  m_faceManager->nodeList().setRelatedObject( *m_nodeManager );


  registerGroup< ElementRegionManager >( groupStructKeys::elemManagerString(), m_elementManager );

  registerGroup< EdgeManager >( groupStructKeys::embSurfEdgeManagerString, m_embSurfEdgeManager );

  registerGroup< EmbeddedSurfaceNodeManager >( groupStructKeys::embSurfNodeManagerString, m_embSurfNodeManager );

  registerWrapper< integer >( viewKeys.meshLevel );

  // increment the modification timestamp at mesh level creation
  // this is to make sure that the actions that depend on this timestamp (such as system setup) are performed at the beginning of the
  // simulations
  modified();
}


MeshLevel::MeshLevel( string const & name,
                      Group * const meshBody,
                      MeshLevel const & source,
                      int const order ):
  MeshLevel( name, meshBody )
{
  GEOS_MARK_FUNCTION;

  // constants for hex mesh
  localIndex const numNodesPerEdge = ( order+1 );
  localIndex const numNodesPerCell = ( order+1 )*( order+1 )*( order+1 );

  localIndex const numInternalNodesPerEdge = ( order-1 );
  localIndex const numInternalNodesPerFace = ( order-1 )*( order-1 );
  localIndex const numInternalNodesPerCell = ( order-1 )*( order-1 )*( order-1 );

  localIndex const numLocalVertices = source.m_nodeManager->size();
  localIndex const numLocalEdges = source.m_edgeManager->size();
  localIndex const numLocalFaces = source.m_faceManager->size();

  globalIndex const maxVertexGlobalID = source.getNodeManager().maxGlobalIndex() + 1;
  globalIndex const maxEdgeGlobalID = source.getEdgeManager().maxGlobalIndex() + 1;
  globalIndex const maxFaceGlobalID = source.getFaceManager().maxGlobalIndex() + 1;

  localIndex numLocalCells = 0;
  source.m_elementManager->forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
  {
    numLocalCells+= sourceSubRegion.size();
  } );

  ////////////////////////////////
  // Get the new number of nodes
  ////////////////////////////////
  localIndex numLocalNodes = numLocalVertices
                             + numLocalEdges * numInternalNodesPerEdge
                             + numLocalFaces * numInternalNodesPerFace
                             + numLocalCells * numInternalNodesPerCell;

  /////////////////////////
  // Nodes
  //////////////////////////

  m_nodeManager->resize( numLocalNodes );

  /////////////////////////
  // Edges
  //////////////////////////

  // the total number of nodes: to add the number of non-vertex edge nodes
  m_edgeManager->resize( numLocalEdges );
  m_edgeManager->getDomainBoundaryIndicator() = source.m_edgeManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex > edgeLocalToGlobal = m_edgeManager->localToGlobalMap();
  arrayView1d< globalIndex const > sourceEdgeLocalToGlobal = source.m_edgeManager->localToGlobalMap();
  LvArray::memcpy( edgeLocalToGlobal.toSlice(), sourceEdgeLocalToGlobal.toSlice() );

  m_edgeManager->constructGlobalToLocalMap();
  m_edgeManager->nodeList().resize( numLocalEdges, numNodesPerEdge );


  /////////////////////////
  // Faces
  //////////////////////////

  m_faceManager->resize( numLocalFaces );
  m_faceManager->faceCenter() = source.m_faceManager->faceCenter();
  m_faceManager->faceNormal() = source.m_faceManager->faceNormal();

  // copy the faces-to-edgs map from source
  m_faceManager->edgeList() = source.m_faceManager->edgeList();
  // copy the faces-to-elements map from source
  m_faceManager->elementList() = source.m_faceManager->elementList();
  m_faceManager->elementRegionList() = source.m_faceManager->elementRegionList();
  m_faceManager->elementSubRegionList() = source.m_faceManager->elementSubRegionList();
  // copy the faces-boundaryIndicator from source
  m_faceManager->getDomainBoundaryIndicator() = source.m_faceManager->getDomainBoundaryIndicator();

  arrayView1d< globalIndex > faceLocalToGlobal = m_faceManager->localToGlobalMap();
  arrayView1d< globalIndex const > sourceFaceLocalToGlobal = source.m_faceManager->localToGlobalMap();
  LvArray::memcpy( faceLocalToGlobal.toSlice(), sourceFaceLocalToGlobal.toSlice() );

  m_faceManager->constructGlobalToLocalMap();

  /////////////////////////
  // Elements
  //////////////////////////

  // check that all elements are hexahedra
  source.m_elementManager->forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    // create element region with the same name as source element region "Region"
    CellElementRegion & region = *(dynamic_cast< CellElementRegion * >( m_elementManager->createChild( sourceRegion.getCatalogName(),
                                                                                                       sourceRegion.getName() ) ) );
    // add cell block to the new element region with the same name as cell block name from source element region
    region.addCellBlockNames( sourceRegion.getCellBlockNames() );

    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {
      if( sourceSubRegion.getElementType() != ElementType::Hexahedron )
      {
        GEOS_ERROR( "Current order number "<<order<<" is higher than one are only available for hexahedral meshes." );
      }

      // create element sub region with the same name as source element sub region "cb"
      CellElementSubRegion & newSubRegion = region.getSubRegions().registerGroup< CellElementSubRegion >( sourceSubRegion.getName() );
      newSubRegion.setElementType( sourceSubRegion.getElementType() );

      // resize per elements value for the new sub region with the new number of nodes per element
      newSubRegion.resizePerElementValues( numNodesPerCell,
                                           sourceSubRegion.numEdgesPerElement(),
                                           sourceSubRegion.numFacesPerElement() );

      newSubRegion.resize( sourceSubRegion.size() );

      // copy new elemCenter map from source
      newSubRegion.getElementCenter() = sourceSubRegion.getElementCenter();

      // copy the elements-to-faces map from source
      newSubRegion.faceList() = sourceSubRegion.faceList();

      // copy the elements-to-edges map from source
      newSubRegion.edgeList() = sourceSubRegion.edgeList();

      arrayView1d< globalIndex > newSubRegionLocalToGlobal = newSubRegion.localToGlobalMap();
      arrayView1d< globalIndex const > sourceSubRegionLocalToGlobal = sourceSubRegion.localToGlobalMap();
      LvArray::memcpy( newSubRegionLocalToGlobal.toSlice(), sourceSubRegionLocalToGlobal.toSlice() );

      newSubRegion.constructGlobalToLocalMap();

    } );

  } );

  CellBlockManagerABC & cellBlockManager = meshBody->getGroup< CellBlockManagerABC >( keys::cellManager );

  cellBlockManager.generateHighOrderMaps( order,
                                          maxVertexGlobalID,
                                          maxEdgeGlobalID,
                                          maxFaceGlobalID,
                                          edgeLocalToGlobal,
                                          faceLocalToGlobal );


  /////////////////////////
  // NodesSets
  //////////////////////////

  this->generateSets();
}


MeshLevel::~MeshLevel()
{
  if( !m_isShallowCopy )
  {
    delete m_nodeManager;
    delete m_edgeManager;
    delete m_faceManager;
    delete m_elementManager;
    delete m_embSurfNodeManager;
    delete m_embSurfEdgeManager;
  }

}

void MeshLevel::initializePostInitialConditionsPostSubGroups()
{
  m_elementManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    subRegion.calculateElementGeometricQuantities( *m_nodeManager, *m_faceManager );
  } );
}

void MeshLevel::generateAdjacencyLists( arrayView1d< localIndex const > const & seedNodeList,
                                        localIndex_array & nodeAdjacencyList,
                                        localIndex_array & edgeAdjacencyList,
                                        localIndex_array & faceAdjacencyList,
                                        ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > & elementAdjacencyList,
                                        integer const depth )
{
  NodeManager const & nodeManager = getNodeManager();

  ArrayOfArraysView< localIndex const > const & nodeToElementRegionList = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementSubRegionList = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementList = nodeManager.elementList().toViewConst();

  FaceManager const & faceManager = this->getFaceManager();
  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  ElementRegionManager const & elemManager = this->getElemManager();

  std::set< localIndex > nodeAdjacencySet;
  std::set< localIndex > edgeAdjacencySet;
  std::set< localIndex > faceAdjacencySet;
  std::vector< std::vector< std::set< localIndex > > > elementAdjacencySet( elemManager.numRegions() );

  // Add the nodes, edges, and faces connected to the volumic element.
  auto const addVolumicSupport = [&]( localIndex const er,
                                      localIndex const esr,
                                      CellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList();
    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();

    for( localIndex const ei: elementAdjacencySet[er][esr] )
    {
      for( localIndex a = 0; a < elemsToNodes.size( 1 ); ++a )  // Range-based for loop not supported for GPU layout.
      {
        nodeAdjacencySet.insert( elemsToNodes( ei, a ) );
      }
      for( localIndex const fi: elemsToFaces[ei] )
      {
        faceAdjacencySet.insert( fi );
        for( auto const & edi: faceToEdges[fi] )
        {
          edgeAdjacencySet.insert( edi );
        }
      }
    }
  };

  // Add the nodes, edges, and faces connected to the fracture element.
  auto const addFractureSupport = [&]( localIndex const er,
                                       localIndex const esr,
                                       FaceElementSubRegion const & subRegion )
  {
    ArrayOfArraysView< localIndex const > const elems2dToNodes = subRegion.nodeList().toViewConst();
    ArrayOfArraysView< localIndex const > const elem2dToEdges = subRegion.edgeList().toViewConst();
    ArrayOfArraysView< localIndex const > const elems2dToFaces = subRegion.faceList().toViewConst();

    for( localIndex const ei: elementAdjacencySet[er][esr] )
    {
      for( localIndex const & ni: elems2dToNodes[ei] )
      {
        nodeAdjacencySet.insert( ni );
      }
      for( localIndex const & edi: elem2dToEdges[ei] )
      {
        edgeAdjacencySet.insert( edi );
      }
      for( localIndex const & fi: elems2dToFaces[ei] )
      {
        faceAdjacencySet.insert( fi );
      }
    }
  };

  // Add all the nodes connected/related to the well elements.
  auto const addWellSupport = [&]( localIndex const er,
                                   localIndex const esr,
                                   WellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList();

    for( localIndex const ei: elementAdjacencySet[er][esr] )
    {
      for( localIndex a = 0; a < elemsToNodes.size( 1 ); ++a )  // Range-based for loop not supported for GPU layout.
      {
        nodeAdjacencySet.insert( elemsToNodes( ei, a ) );
      }
    }
  };

  // Add all the collocated nodes of the fracture element.
  auto const addCollocatedFractureNodes = [&]( FaceElementSubRegion const & subRegion )
  {
    auto const & l2g = nodeManager.localToGlobalMap();
    auto const & g2l = nodeManager.globalToLocalMap();

    std::set< localIndex > newNodes;
    for( localIndex const & ln: nodeAdjacencySet )
    {
      globalIndex const & gn = l2g[ln];
      for( std::set< globalIndex > const & bucket: subRegion.getCollocatedNodes() )
      {
        if( bucket.find( gn ) != bucket.cend() )
        {
          for( globalIndex const & n: bucket )
          {
            auto it = g2l.find( n );
            if( it != g2l.cend() )
            {
              newNodes.insert( it->second );
            }
          }
        }
      }
    }
    nodeAdjacencySet.insert( newNodes.cbegin(), newNodes.cend() );
  };

  for( localIndex a = 0; a < elemManager.numRegions(); ++a )
  {
    elementAdjacencySet[a].resize( elemManager.getRegion( a ).numSubRegions() );
  }

  nodeAdjacencySet.insert( seedNodeList.begin(), seedNodeList.end() );

  for( integer d = 0; d < depth; ++d )
  {
    for( localIndex const nodeIndex: nodeAdjacencySet )
    {
      for( localIndex b = 0; b < nodeToElementRegionList.sizeOfArray( nodeIndex ); ++b )
      {
        localIndex const er = nodeToElementRegionList[nodeIndex][b];
        localIndex const esr = nodeToElementSubRegionList[nodeIndex][b];
        localIndex const ei = nodeToElementList[nodeIndex][b];
        elementAdjacencySet[er][esr].insert( ei );
      }
    }

    for( localIndex er = 0; er < elemManager.numRegions(); ++er )
    {
      ElementRegionBase const & elemRegion = elemManager.getRegion( er );

      elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                         CellElementSubRegion const & subRegion )
      {
        addVolumicSupport( er, esr, subRegion );
      } );
      elemRegion.forElementSubRegionsIndex< FaceElementSubRegion >( [&]( localIndex const esr,
                                                                         FaceElementSubRegion const & subRegion )
      {
        addFractureSupport( er, esr, subRegion );
        addCollocatedFractureNodes( subRegion );
      } );
      elemRegion.forElementSubRegionsIndex< WellElementSubRegion >( [&]( localIndex const esr,
                                                                         WellElementSubRegion const & subRegion )
      {
        addWellSupport( er, esr, subRegion );
      } );
    }
  }

  // Convert the `std::set` containers to `LvArray` containers.
  nodeAdjacencyList.resize( LvArray::integerConversion< localIndex >( nodeAdjacencySet.size()));
  std::copy( nodeAdjacencySet.begin(), nodeAdjacencySet.end(), nodeAdjacencyList.begin() );

  edgeAdjacencyList.resize( LvArray::integerConversion< localIndex >( edgeAdjacencySet.size()));
  std::copy( edgeAdjacencySet.begin(), edgeAdjacencySet.end(), edgeAdjacencyList.begin() );

  faceAdjacencyList.resize( LvArray::integerConversion< localIndex >( faceAdjacencySet.size()));
  std::copy( faceAdjacencySet.begin(), faceAdjacencySet.end(), faceAdjacencyList.begin() );

  for( localIndex er=0; er < elemManager.numRegions(); ++er )
  {
    ElementRegionBase const & elemRegion = elemManager.getRegion( er );

    for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
    {
      elementAdjacencyList[er][esr].get().clear();
      elementAdjacencyList[er][esr].get().resize( LvArray::integerConversion< localIndex >( elementAdjacencySet[er][esr].size()) );
      std::copy( elementAdjacencySet[er][esr].begin(),
                 elementAdjacencySet[er][esr].end(),
                 elementAdjacencyList[er][esr].get().begin() );

    }
  }
}

void MeshLevel::generateSets()
{
  GEOS_MARK_FUNCTION;

  NodeManager const & nodeManager = *m_nodeManager;

  dataRepository::Group const & nodeSets = nodeManager.sets();

  map< string, array1d< bool > > nodeInSet; // map to contain indicator of whether a node is in a set.
  string_array setNames; // just a holder for the names of the sets

  // loop over all wrappers and fill the nodeIndSet arrays for each set
  for( auto & wrapper : nodeSets.wrappers() )
  {
    string const & name = wrapper.second->getName();
    nodeInSet[name].resize( nodeManager.size() );
    nodeInSet[name].setValues< serialPolicy >( false );

    if( nodeSets.hasWrapper( name ) )
    {
      setNames.emplace_back( name );
      SortedArrayView< localIndex const > const & set = nodeSets.getReference< SortedArray< localIndex > >( name );
      for( localIndex const a : set )
      {
        nodeInSet[name][a] = true;
      }
    }
  }


  ElementRegionManager & elementRegionManager = *m_elementManager;
  elementRegionManager.forElementSubRegions( [&]( auto & subRegion )
  {
    dataRepository::Group & elementSets = subRegion.sets();

    auto const & elemToNodeMap = subRegion.nodeList();

    for( string const & setName : setNames )
    {
      arrayView1d< bool const > const nodeInCurSet = nodeInSet[setName];

      SortedArray< localIndex > & targetSet = elementSets.registerWrapper< SortedArray< localIndex > >( setName ).reference();
      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        localIndex const numNodes = subRegion.numNodesPerElement( k );

        localIndex elementInSet = true;
        for( localIndex i = 0; i < numNodes; ++i )
        {
          if( !nodeInCurSet( elemToNodeMap[ k ][ i ] ) )
          {
            elementInSet = false;
            break;
          }
        }

        if( elementInSet )
        {
          targetSet.insert( k );
        }
      }
    }
  } );
}

MeshLevel const & MeshLevel::getShallowParent() const
{
  return ( m_shallowParent!=nullptr ? *m_shallowParent : *this);
}

MeshLevel & MeshLevel::getShallowParent()
{
  return ( m_shallowParent!=nullptr ? *m_shallowParent : *this);
}


bool MeshLevel::isShallowCopyOf( MeshLevel const & comparisonLevel ) const
{
  return ( m_nodeManager == comparisonLevel.m_nodeManager ) &&
         ( m_edgeManager == comparisonLevel.m_edgeManager ) &&
         ( m_faceManager == comparisonLevel.m_faceManager ) &&
         ( m_elementManager == comparisonLevel.m_elementManager ) &&
         ( m_embSurfNodeManager == comparisonLevel.m_embSurfNodeManager) &&
         ( m_embSurfEdgeManager == comparisonLevel.m_embSurfEdgeManager ) &&
         isShallowCopy();
}


} /* namespace geos */
