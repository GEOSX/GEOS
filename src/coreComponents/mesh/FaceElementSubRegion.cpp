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

#include "FaceElementSubRegion.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{
using namespace dataRepository;


FaceElementSubRegion::FaceElementSubRegion( string const & name,
                                            dataRepository::Group * const parent ):
  SurfaceElementSubRegion( name, parent ),
  m_unmappedGlobalIndicesInToEdges(),
  m_unmappedGlobalIndicesInToFaces(),
  m_newFaceElements(),
  m_toFacesRelation()
{
  m_elementType = ElementType::Hexahedron;

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();

  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation ).
    setDescription( "Map to the faces attached to each FaceElement." ).
    reference().resize( 0, 2 );

  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString(), &m_edgesTo2dFaces ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of edge local indices to the fracture connector local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString(), &m_2dFaceToEdge ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices to edge local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString(), &m_2dFaceTo2dElems ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices face element local indices" ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::elem2dToCollocatedNodesBucketsString(), &m_2dElemToCollocatedNodesBuckets ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map eventually containing all the collocated nodes." ).
    setSizedFromParent( 1 );

#ifdef GEOS_USE_SEPARATION_COEFFICIENT
  registerWrapper( viewKeyStruct::separationCoeffString(), &m_separationCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
    setDescription( "Scalar indicator of level of separation for a fracturing face." );
#endif

  excludeWrappersFromPacking( { viewKeyStruct::faceListString() } );

  m_2dElemToElems.resize( 0, 2 );

  m_numNodesPerElement = 8;
}

void FaceElementSubRegion::copyFromCellBlock( FaceBlockABC const & faceBlock )
{
  localIndex const num2dElements = faceBlock.num2dElements();
  resize( faceBlock.num2dElements() );

  m_toNodesRelation.base() = faceBlock.get2dElemToNodes();
  m_toEdgesRelation.base() = faceBlock.get2dElemToEdges();

  // `FaceBlockABC` is designed to be heterogeneous.
  // `FaceElementSubRegion` inherits from `ElementSubRegionBase` which is meant to be homogeneous.
  // But `FaceElementSubRegion` is sometimes used as an heterogeneous sub region,
  // which emphasizes the need of a refactoring.
  // In the meantime, we try to fill the face block into the sub region and hope for the best...
  {
    auto const deduce3dElemType = [&]( integer maxSize ) -> ElementType
    {
      switch( maxSize )
      {
        case 3:
          return ElementType::Wedge;
        case 4:
          return ElementType::Hexahedron;
        case 5:
          return ElementType::Prism5;
        case 6:
          return ElementType::Prism6;
        case 7:
          return ElementType::Prism7;
        case 8:
          return ElementType::Prism8;
        case 9:
          return ElementType::Prism9;
        case 10:
          return ElementType::Prism10;
        case 11:
          return ElementType::Prism11;
        case 0:
          // In the case the fracture is empty (on this rank), then we default to hexahedron. Otherwise, there's something wrong
          GEOS_ERROR_IF_NE_MSG( num2dElements, 0, "Could not determine the element type of the fracture \"" << getName() << "\"." );
          return ElementType::Hexahedron;
        default:
          GEOS_ERROR( "Unsupported type of elements during the face element sub region creation." );
          return {};
      }
    };

    m_2dElemToCollocatedNodesBuckets = faceBlock.get2dElemsToCollocatedNodesBuckets();
    // Checking if all the 2d elements are homogeneous.
    // We rely on the number of nodes for each element to find out.
    std::vector< integer > numNodesPerElement( num2dElements );
    for( int i = 0; i < num2dElements; ++i )
    {
      numNodesPerElement[i] = m_2dElemToCollocatedNodesBuckets[i].size();
    }
    std::set< integer > const sizes( numNodesPerElement.cbegin(), numNodesPerElement.cend() );

    if( sizes.size() > 1 )
    {
      // If we have found that the input face block contains 2d elements of different types,
      // we inform the used that the situation may be at risk.
      // (We're storing the face block in a homogeneous container while it's actually heterogeneous).
      GEOS_WARNING( "Heterogeneous face element sub region found and stored as homogeneous. Use at your own risk." );
    }

    auto const it = std::max_element( sizes.cbegin(), sizes.cend() );
    integer const maxSize = it != sizes.cend() ? *it : 0;
    m_elementType = deduce3dElemType( maxSize );
    m_numNodesPerElement = maxSize;
  }

  // The `m_2dElemToElems` mappings involves element, sub regions and regions indices.
  // We store the element indices that are correct.
  // But we only have access to the cell block indices, not the sub regions indices.
  // Temporarily, and also because they share the same dimensions,
  // we store the cell block mapping at the sub region mapping location.
  // It will later be transformed into a sub regions mapping.
  // Last, we fill the regions mapping with dummy -1 values that should all be replaced eventually.
  auto const elem2dToElems = faceBlock.get2dElemToElems();
  m_2dElemToElems.resize( num2dElements, 2 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( localIndex const & j: elem2dToElems.toCellIndex[i] )
    {
      m_2dElemToElems.m_toElementIndex.emplaceBack( i, j );
    }
    for( localIndex const & j: elem2dToElems.toBlockIndex[i] )
    {
      m_2dElemToElems.m_toElementSubRegion.emplaceBack( i, j );
    }
  }

  m_toFacesRelation.base() = faceBlock.get2dElemToFaces();

  m_2dFaceToEdge = faceBlock.get2dFaceToEdge();
  m_2dFaceTo2dElems = faceBlock.get2dFaceTo2dElems();

  m_localToGlobalMap = faceBlock.localToGlobalMap();
  this->constructGlobalToLocalMap();

  for( int i = 0; i < faceBlock.num2dFaces(); ++i )
  {
    m_recalculateConnectionsFor2dFaces.insert( i );
  }

  for( localIndex i = 0; i < faceBlock.num2dElements(); ++i )
  {
    m_newFaceElements.insert( i );
  }

  // TODO We still need to be able to import fields on the FaceElementSubRegion.
}

void FaceElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

void FaceElementSubRegion::calculateSingleElementGeometricQuantities( localIndex const k,
                                                                      arrayView1d< real64 const > const & faceArea )
{
  m_elementArea[k] = faceArea[ m_toFacesRelation[k][0] ];
  m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
}

void FaceElementSubRegion::calculateElementGeometricQuantities( NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                                                                FaceManager const & faceManager )
{
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();

  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateSingleElementGeometricQuantities( k, faceArea );
  } );
}

ElementType FaceElementSubRegion::getElementType( localIndex ei ) const
{
  // We try to get the element type information from the bucket.
  // If this information does not appear to be reliable,
  // let's fall back on the 2d elem to nodes.
  auto const sizeFromBuckets = m_2dElemToCollocatedNodesBuckets[ei].size();
  auto const size = sizeFromBuckets > 0 ? sizeFromBuckets : m_toNodesRelation[ei].size() / 2;
  switch( size )
  {
    case 3:
      return ElementType::Triangle;
    case 4:
      return ElementType::Quadrilateral;
    default:
      return ElementType::Polygon;
  }
}


localIndex FaceElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex FaceElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}


template< bool DO_PACKING >
localIndex FaceElementSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > const nodeLocalToGlobal = m_toNodesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const edgeLocalToGlobal = m_toEdgesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const faceLocalToGlobal = m_toFacesRelation.relatedObjectLocalToGlobal();

  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toNodesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToNodes,
                                               packList,
                                               localToGlobal,
                                               nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toEdgesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToEdges,
                                               packList,
                                               localToGlobal,
                                               edgeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toFacesRelation.toViewConst(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::surfaceElementsToCellRegionsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_2dElemToElems,
                                               packList,
                                               m_2dElemToElems.getElementRegionManager() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elem2dToCollocatedNodesBucketsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_2dElemToCollocatedNodesBuckets.toViewConst(),
                                               packList,
                                               localToGlobal );

  return packedSize;
}


localIndex FaceElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const overwriteUpMaps,
                                                   bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::surfaceElementsToCellRegionsString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToElems,
                                     packList.toViewConst(),
                                     m_2dElemToElems.getElementRegionManager(),
                                     overwriteUpMaps );

  string elem2dToCollocatedNodesBucketsString;
  unPackedSize += bufferOps::Unpack( buffer, elem2dToCollocatedNodesBucketsString );
  GEOS_ERROR_IF_NE( elem2dToCollocatedNodesBucketsString, viewKeyStruct::elem2dToCollocatedNodesBucketsString() );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToCollocatedNodesBuckets,
                                     packList,
                                     this->globalToLocalMap() );


  GEOS_ERROR_IF_NE( m_unmappedGlobalIndicesInToNodes.size(), 0 );
  GEOS_ERROR_IF_NE( m_unmappedGlobalIndicesInToEdges.size(), 0 );
  GEOS_ERROR_IF_NE( m_unmappedGlobalIndicesInToFaces.size(), 0 );

  return unPackedSize;
}

localIndex FaceElementSubRegion::packToFaceRelationSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packToFaceRelationImpl< false >( junk, packList );
}

localIndex FaceElementSubRegion::packToFaceRelation( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  return packToFaceRelationImpl< true >( buffer, packList );
}


template< bool DO_PACKING >
localIndex FaceElementSubRegion::packToFaceRelationImpl( buffer_unit_type * & buffer,
                                                         arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > const faceLocalToGlobal = m_toFacesRelation.relatedObjectLocalToGlobal();

  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toFacesRelation.toViewConst(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );
  return packedSize;
}


localIndex FaceElementSubRegion::unpackToFaceRelation( buffer_unit_type const * & buffer,
                                                       localIndex_array & packList,
                                                       bool const GEOS_UNUSED_PARAM( overwriteUpMaps ),
                                                       bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal() );

  return unPackedSize;
}


/**
 * @brief The two mappings @p elem2dToElems3d and @p elem2dToFaces have to be consistent
 * in the sense that each @p face for a given index must "belong" to @p 3d @p element at the same index.
 * @param[in] fractureName The name of the fracture for which we're checking the mapping consistency.
 * @param[in] elem2dToElems3d A mapping.
 * @param[in,out] elem2dToFaces This mapping will be corrected if needed to match @p elem2dToElems3d.
 */
void fixNeighborMappingsInconsistency( string const & fractureName,
                                       OrderedVariableToManyElementRelation const & elem2dToElems3d,
                                       FaceElementSubRegion::FaceMapType & elem2dToFaces )
{
  {
    localIndex const num2dElems = elem2dToFaces.size();
    for( int e2d = 0; e2d < num2dElems; ++e2d )
    {
      std::set< localIndex > const sizes{
        elem2dToFaces[e2d].size(),
        elem2dToElems3d.m_toElementRegion[e2d].size(),
        elem2dToElems3d.m_toElementSubRegion[e2d].size(),
        elem2dToElems3d.m_toElementIndex[e2d].size()
      };

      if( sizes.size() != 1 || sizes.find( 2 ) == sizes.cend() )
      {
        continue;
      }

      localIndex const f0 = elem2dToFaces[e2d][0];
      localIndex const er0 = elem2dToElems3d.m_toElementRegion[e2d][0];
      localIndex const esr0 = elem2dToElems3d.m_toElementSubRegion[e2d][0];
      localIndex const ei0 = elem2dToElems3d.m_toElementIndex[e2d][0];
      auto const & faces0 = elem2dToElems3d.getElementRegionManager()->getRegion( er0 ).getSubRegion< CellElementSubRegion >( esr0 ).faceList()[ei0];

      localIndex const f1 = elem2dToFaces[e2d][1];
      localIndex const er1 = elem2dToElems3d.m_toElementRegion[e2d][1];
      localIndex const esr1 = elem2dToElems3d.m_toElementSubRegion[e2d][1];
      localIndex const ei1 = elem2dToElems3d.m_toElementIndex[e2d][1];
      auto const & faces1 = elem2dToElems3d.getElementRegionManager()->getRegion( er1 ).getSubRegion< CellElementSubRegion >( esr1 ).faceList()[ei1];

      bool const match00 = std::find( faces0.begin(), faces0.end(), f0 ) != faces0.end();
      bool const match11 = std::find( faces1.begin(), faces1.end(), f1 ) != faces1.end();
      bool const match01 = std::find( faces0.begin(), faces0.end(), f1 ) != faces0.end();
      bool const match10 = std::find( faces1.begin(), faces1.end(), f0 ) != faces1.end();

      bool const matchCrossed = !match00 && !match11 && match01 && match10;
      bool const matchStraight = match00 && match11 && !match01 && !match10;

      if( matchCrossed )
      {
        std::swap( elem2dToFaces[e2d][0], elem2dToFaces[e2d][1] );
      }
      else if( !matchStraight )
      {
        GEOS_ERROR( "Mapping neighbor inconsistency detected for fracture " << fractureName );
      }
    }
  }
}

void FaceElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );

//  fixNeighborMappingsInconsistency( getName(), m_2dElemToElems, m_toFacesRelation );
}

/**
 * @brief Returns a mapping that links any collocated node to the collocated node with the lowest index.
 * @param elem2dToCollocatedNodesBuckets All the collocated nodes buckets.
 * @return The computed mapping.
 * @details For each collocated node, the returned mapping will provide
 * the lowest id among all the collocated nodes sharing the same position.
 * That way, it's possible to know if two nodes are collocated of each other by checking if they share the same lowest id.
 */
std::map< globalIndex, globalIndex > buildReferenceCollocatedNodes( ArrayOfArrays< array1d< globalIndex > > const & elem2dToCollocatedNodesBuckets )
{
  std::map< globalIndex, globalIndex > referenceCollocatedNodes;  // Will be returned.

  // Since some 2d elem may share some nodes, some of the collocated nodes buckets will be duplicated.
  // We want to remove this.
  std::set< std::set< globalIndex > > uniqueCollocatedNodesBuckets;
  for( localIndex e2d = 0; e2d < elem2dToCollocatedNodesBuckets.size(); ++e2d )
  {
    for( integer ni = 0; ni < elem2dToCollocatedNodesBuckets[e2d].size(); ++ni )
    {
      array1d< globalIndex > const & tmp = elem2dToCollocatedNodesBuckets( e2d, ni );
      std::set< globalIndex > const tmp2( tmp.begin(), tmp.end() );
      uniqueCollocatedNodesBuckets.insert( tmp2 );
    }
  }

  // Now we can perform simple loops to compute the mapping.
  for( std::set< globalIndex > const & bucket: uniqueCollocatedNodesBuckets )
  {
    globalIndex const & refNode = *std::min_element( bucket.cbegin(), bucket.cend() );
    for( globalIndex const & n: bucket )
    {
      referenceCollocatedNodes[n] = refNode;
    }
  }

  return referenceCollocatedNodes;
}


/**
 * @brief Computes the mapping which links a pair of collocated nodes
 * to all the edges having their nodes collocated to this key pair of nodes.
 * @param referenceCollocatedNodes Mapping that link all the collocated nodes
 * to the collocated node with the lowest index (considered as a reference).
 * The reference node can then be used to make connection between geometrical objects
 * relying on collocated nodes (in the case of this function, edges).
 * @param nl2g The node local to global mapping.
 * @param edgeToNodes The mapping from local edges to local nodes.
 * @return The computed mapping.
 * @details For each edge based on collocated nodes (i.e. each edge on the fracture),
 * we consider each node and compute the reference collocated node.
 * There will be multiple edges with the same pair of reference collocated nodes.
 * This information is contained in the returned mapping.
 */
std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > >
buildCollocatedEdgeBuckets( std::map< globalIndex, globalIndex > const & referenceCollocatedNodes,
                            arrayView1d< globalIndex const > const nl2g,
                            arrayView2d< localIndex const > const edgeToNodes )
{
  GEOS_ASSERT_EQ( edgeToNodes.size( 1 ), 2 );

  // Checks if the node `gni` is handled as a collocated node on the curren rank.
  auto hasCollocatedNode = [&]( globalIndex const gni ) -> bool
  {
    return referenceCollocatedNodes.find( gni ) != referenceCollocatedNodes.cend();
  };

  // `edgeIds` is a temporary container that maps the two nodes of any edge of all the face elements,
  // to the local index of the edge itself. We fill this container using the `edgeToNodes` mapping
  // because we want to be sure to test all the combinations of nodes,
  // and therefore not to forget any possible edge.
  // It's important to note that the key of `edgeIds` are the global indices of the nodes, in no particular order.
  std::map< std::pair< globalIndex, globalIndex >, localIndex > edgesIds;
  for( localIndex lei = 0; lei < edgeToNodes.size( 0 ); ++lei )
  {
    auto const & nodes = edgeToNodes[lei];
    globalIndex const & gni0 = nl2g[nodes[0]];
    globalIndex const & gni1 = nl2g[nodes[1]];
    if( hasCollocatedNode( gni0 ) && hasCollocatedNode( gni1 ) )
    {
      edgesIds[{ gni0, gni1 }] = lei;
    }
  }

  // The key of the `collocatedEdgeBuckets` map (i.e. `std::pair< globalIndex, globalIndex >`) represents the global indices of two nodes.
  // Those two nodes are the lowest index of collocated nodes. As such, those two nodes may not form an existing edge.
  // But this trick lets us define some kind of _hash_ that allows to compare the location of the edges:
  // edges sharing the same hash lie in the same position.
  std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > collocatedEdgeBuckets;
  for( auto const & p: edgesIds )
  {
    static constexpr std::string_view nodeNotFound = "Internal error when trying to access the reference collocated node for global node {}.";

    std::pair< globalIndex, globalIndex > const & nodes = p.first;
    localIndex const & edge = p.second;

    auto it0 = referenceCollocatedNodes.find( nodes.first );
    GEOS_ERROR_IF( it0 == referenceCollocatedNodes.cend(), GEOS_FMT( nodeNotFound, nodes.first ) );
    globalIndex const n0 = it0->second;

    auto it1 = referenceCollocatedNodes.find( nodes.second );
    GEOS_ERROR_IF( it1 == referenceCollocatedNodes.cend(), GEOS_FMT( nodeNotFound, nodes.second ) );
    globalIndex const n1 = it1->second;

    std::pair< globalIndex, globalIndex > const edgeHash = std::minmax( n0, n1 );
    collocatedEdgeBuckets[edgeHash].insert( edge );
  }

  return collocatedEdgeBuckets;
}


/**
 * @brief Returns a mapping that links any collocated edge to the collocated edge with the lowest index.
 * @param collocatedEdgeBuckets Links the pairs of reference collocated nodes to the edges overlapping those nodes (even with other nodes).
 * @param edgeGhostRanks The ghost rank of the edges.
 * @return The computed mapping.
 */
std::map< localIndex, localIndex > buildReferenceCollocatedEdges( std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > const & collocatedEdgeBuckets,
                                                                  arrayView1d< integer const > const edgeGhostRanks )
{
  std::map< localIndex, localIndex > referenceCollocatedEdges;

  // We want to consider in priority the edges that are owned by the rank.
  // So this comparator takes the ghost rank into account and favors lower values of the ghost rank.
  auto comp = [=]( int e,
                   int f ) -> bool
  {
    int const re = edgeGhostRanks[e] < 0 ? 0 : 1;
    int const rf = edgeGhostRanks[f] < 0 ? 0 : 1;
    return std::tie( re, e ) < std::tie( rf, f );
  };

  for( auto const & p: collocatedEdgeBuckets )
  {
    std::set< localIndex > const & collocatedEdges = p.second;
    localIndex const refEdge = *std::min_element( collocatedEdges.cbegin(), collocatedEdges.cend(), comp );
    for( localIndex const & collocatedEdge: collocatedEdges )
    {
      referenceCollocatedEdges[collocatedEdge] = refEdge;
    }
  }

  return referenceCollocatedEdges;
}


/**
 * @brief Returns a SortedArray with sequentially increasing values, starting with value and repetitively evaluating ++value
 * @param newSize The size of the returned array.
 * @param value The starting value.
 * @return The filled array.
 */
SortedArray< localIndex > makeSortedArrayIota( localIndex newSize, localIndex value = 0 )
{
  SortedArray< localIndex > result;
  result.reserve( newSize );
  for( localIndex i = value; i < value + newSize; ++i )
  {
    result.insert( i );
  }
  return result;
}


/**
 * @brief Builds the 2d face to 2d elements mapping.
 * @param elem2dToEdges The 2d element (geometrical faces in 3d) to edges mapping.
 * @param edgesTo2dFaces The edges to 2d faces (geometrical edges in 3d) mapping.
 * @param referenceCollocatedEdges The mapping that, for a given edge index, returns the collocated edge with lowest index.
 * @return The computed mapping.
 */
ArrayOfArrays< geos::localIndex > build2dFaceTo2dElems( ArrayOfArraysView< localIndex const > const elem2dToEdges,
                                                        map< localIndex, localIndex > const & edgesTo2dFaces,
                                                        std::map< geos::localIndex, geos::localIndex > const & referenceCollocatedEdges )
{
  ArrayOfArrays< localIndex > face2dTo2dElems;

  auto const num2dElems = elem2dToEdges.size();
  auto const num2dFaces = edgesTo2dFaces.size();

  // `tmp` contains the 2d face to 2d elements mappings as a `std` container.
  // Eventually, it's copied into an `LvArray` container.
  std::vector< std::vector< localIndex > > tmp( num2dFaces );
  for( auto i = 0; i < num2dElems; ++i )
  {
    for( auto const & e: elem2dToEdges[i] )
    {
      if( e < 0 )
      {
        continue;
      }
      tmp[edgesTo2dFaces.at( referenceCollocatedEdges.at( e ) )].push_back( i );
    }
  }
  std::vector< localIndex > sizes;
  sizes.reserve( tmp.size() );
  for( std::vector< localIndex > const & t: tmp )
  {
    sizes.push_back( t.size() );
  }

  face2dTo2dElems.resizeFromCapacities< serialPolicy >( sizes.size(), sizes.data() );

  for( std::size_t i = 0; i < tmp.size(); ++i )
  {
    for( std::size_t j = 0; j < tmp[i].size(); ++j )
    {
      face2dTo2dElems.emplaceBack( i, tmp[i][j] );
    }
  }

  return face2dTo2dElems;
}


/**
 * @brief Adds some missing connections to the @p elem2dToNodes mapping.
 * @param[in] elem2dToCollocatedNodesBuckets The bucket of all the collocated nodes nodes, for each 2d element.
 * @param[in] ng2l The global to local node mapping.
 * @param[in,out] elem2dToNodes The 2d element to nodes mapping that will receive new connections
 * @details Due to the specific way 2d elements are managing nodes, some connections may be missing before the ghosting process.
 * After the ghosting has been performed, those connections can be added using this function.
 * @note This functions is meant to be called after the ghosting has occurred.
 */
void fillMissing2dElemToNodes( ArrayOfArrays< array1d< globalIndex > > const & elem2dToCollocatedNodesBuckets,
                               unordered_map< globalIndex, localIndex > const & ng2l,
                               ArrayOfArrays< localIndex > & elem2dToNodes )
{
  auto const num2dElems = elem2dToNodes.size();

  // We loop over all the collocated nodes attached to the 2d elements.
  // If a node is on the rank while not attached to the 2d element, then we add a connection.
  for( int e2d = 0; e2d < num2dElems; ++e2d )
  {
    auto bucket = elem2dToCollocatedNodesBuckets[e2d];
    for( array1d< globalIndex > const & collocatedNodes: bucket )
    {
      for( globalIndex const & collocatedNode: collocatedNodes )
      {
        auto g2l = ng2l.find( collocatedNode );
        if( g2l != ng2l.cend() )
        {
          localIndex const lni = g2l->second;
          auto nodes = elem2dToNodes[e2d];
          if( std::find( nodes.begin(), nodes.end(), lni ) == nodes.end() )
          {
            elem2dToNodes.emplaceBack( e2d, lni );
          }
        }
      }
    }
  }
}


/**
 * @brief Adds some missing connections to the @p elem2dToEdges mapping.
 * @param[in] elem2dToNodes The 2d elements to nodes mappings.
 * @param[in] nodesToEdges The nodes to edges mapping.
 * @param[in] nl2g The global to local mapping for nodes.
 * @param[in] referenceCollocatedNodes A mapping that link all each collocated node to the node that must be used as a reference.
 * @param[in] collocatedEdgeBuckets A specific mapping that links all pair of collocated nodes
 * to the reference edge built on top of the two locations.
 * @param[in,out] elem2dToEdges The 2d element to edges that will be completed/corrected.
 * @details The @p elem2dToEdges is (for the moment) used to build connection between the 2d elements.
 * Due to the specific way 2d elements are managing nodes, before the ghosting process:
 *   - some edges may be missing,
 *   - some edges may be "wrong" (in the sense where there may be multiple collocated edges,
 *     and we must be sure that all the elements use the same reference edges).
 * After the ghosting's been performed, those connections can be reset.
 */
void fillMissing2dElemToEdges( ArrayOfArraysView< localIndex const > const elem2dToNodes,
                               ArrayOfSetsView< localIndex const > const nodesToEdges,
                               arrayView1d< globalIndex const > const nl2g,
                               std::map< globalIndex, globalIndex > const & referenceCollocatedNodes,
                               std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > const & collocatedEdgeBuckets,
                               ArrayOfArrays< localIndex > & elem2dToEdges )
{
  localIndex const num2dElems = elem2dToNodes.size();
  for( localIndex e2d = 0; e2d < num2dElems; ++e2d )
  {
    auto const numNodes = elem2dToNodes.sizeOfArray( e2d );
    auto const numEdges = elem2dToEdges.sizeOfArray( e2d );
    if( 2 * numEdges == numNodes )
    {
      // If the previous conditions is true, all the information could be constructed and therefore should be OK.
      // There's no need to proceed.
      continue;
    }

    // `nodesOfEdgesTouching2dElem` deals with the edges that have at least one point touching the 2d element.
    // While `nodesOfEdgesOf2dElem` deals with the edges for which all two nodes are on the 2d element.
    // For both mappings, the key is the edge index and the values are the local nodes indices of the concerned edges.
    std::map< localIndex, std::vector< localIndex > > nodesOfEdgesTouching2dElem, nodesOfEdgesOf2dElem;
    for( localIndex const & n: elem2dToNodes[e2d] )
    {
      for( localIndex const & e: nodesToEdges[n] )
      {
        nodesOfEdgesTouching2dElem[e].push_back( n );
      }
    }
    for( auto const & ens: nodesOfEdgesTouching2dElem )
    {
      if( ens.second.size() == 2 )
      {
        nodesOfEdgesOf2dElem.insert( ens );
      }
    }
    // We are now recomputing all the edges of the 2d elements.
    // Even if some edges were already present in the initial `elem2dToEdges` mapping,
    // we'll ditch them and start with a new collection.
    std::set< localIndex > allEdgesOf2dElem;
    for( auto const & ens: nodesOfEdgesOf2dElem )
    {
      std::vector< localIndex > const & nodesOfEdge = ens.second;
      globalIndex const & gn0 = referenceCollocatedNodes.at( nl2g[ nodesOfEdge[0] ] );
      globalIndex const & gn1 = referenceCollocatedNodes.at( nl2g[ nodesOfEdge[1] ] );
      std::set< localIndex > candidateEdges = collocatedEdgeBuckets.at( std::minmax( { gn0, gn1 } ) );
      auto const min = std::min_element( candidateEdges.cbegin(), candidateEdges.cend() );
      allEdgesOf2dElem.insert( *min );
    }
    elem2dToEdges.clearArray( e2d );
    for( localIndex const & e: allEdgesOf2dElem )
    {
      elem2dToEdges.emplaceBack( e2d, e );
    }
  }
}


/**
 * @brief Builds a 2d face to edges mapping
 * @param referenceCollocatedEdges Maps all the collocated edges to a single reference collocated edge.
 * @return A 1d array that maps local index of a 2d face to the equivalent (3d) reference edge.
 * @details The @p referenceCollocatedEdges input basically contains all the edges of the fracture.
 * This function only selects the reference edges and builds a 2d face for each.
 * The 2d face ordering is more or less random (actually it's sorted on the index of the reference edge).
 * But the ordering is not critical as long as it's consistent withing the fracture.
 */
array1d< localIndex > build2dFaceToEdge( std::map< localIndex, localIndex > const & referenceCollocatedEdges )
{
  std::set< localIndex > const referenceEdges = mapValues< std::set >( referenceCollocatedEdges );

  localIndex const num2dFaces = LvArray::integerConversion< localIndex >( referenceEdges.size() );

  // For the `m_2dFaceToEdge`, we can select any values that we want.
  // But then we need to be consistent...
  array1d< localIndex > face2dToEdge;
  face2dToEdge.reserve( num2dFaces );
  for( localIndex const & refEdge: referenceEdges )
  {
    face2dToEdge.emplace_back( refEdge );
  }
  GEOS_ASSERT_EQ( num2dFaces, face2dToEdge.size() );

  return face2dToEdge;
}


/**
 * @brief Builds the edges to 2d faces mappings by inverting the 2d faces to edges mappings.
 * @param face2dToEdges The mappings to be inverted.
 * @return The mapping
 */
map< localIndex, localIndex > buildEdgesToFace2d( arrayView1d< localIndex const > const face2dToEdges )
{
  map< localIndex, localIndex > edgesToFace2d;
  localIndex const num2dFaces = face2dToEdges.size();

  for( localIndex i = 0; i < num2dFaces; ++i )
  {
    edgesToFace2d[face2dToEdges[i]] = i;
  }

  GEOS_ASSERT_EQ_MSG( LvArray::integerConversion< localIndex >( edgesToFace2d.size() ), num2dFaces, "Internal error. The mappings `edgesToFace2d` and `face2dToEdges` should have the same size" );

  return edgesToFace2d;
}


/**
 * @brief Uses the two input mappings to reorder the @p elem2dToNodes mapping.
 * @param[in] elem2dToFaces The 2d element to faces mapping.
 * @param[in] facesToNodes The face to nodes mapping.
 * @param[out] elem2dToNodes The 2d element to nodes that will be overwritten.
 * @details The @p elem2dToNodes a priori does not respect any specific nodes order.
 * But the vtk output expects the first half of the nodes to make a first face (in the correct order),
 * and the second to make a second face. But nothing imposes this in our design.
 * Even if we should have a more explicit design,
 * the current function resets this implicit information in our mappings.
 */
void fixNodesOrder( ArrayOfArraysView< localIndex const > const elem2dToFaces,
                    ArrayOfArraysView< localIndex const > const facesToNodes,
                    ArrayOfArrays< localIndex > & elem2dToNodes )
{
  localIndex const num2dElems = elem2dToNodes.size();
  for( localIndex e2d = 0; e2d < num2dElems; ++e2d )
  {
    std::vector< localIndex > nodesOfFace;
    for( localIndex fi: elem2dToFaces[e2d] )
    {
      for( localIndex ni: facesToNodes[fi] )
      {
        nodesOfFace.push_back( ni );
      }
    }
    elem2dToNodes.clearArray( e2d );
    elem2dToNodes.appendToArray( e2d, nodesOfFace.cbegin(), nodesOfFace.cend() );
  }
}


void FaceElementSubRegion::fixSecondaryMappings( NodeManager const & nodeManager,
                                                 EdgeManager const & edgeManager,
                                                 FaceManager const & faceManager,
                                                 ElementRegionManager const & elemManager )
{
  arrayView1d< globalIndex const > const nl2g = nodeManager.localToGlobalMap();
  ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();

  // First let's create the reference mappings for both nodes and edges.
  std::map< globalIndex, globalIndex > const referenceCollocatedNodes = buildReferenceCollocatedNodes( m_2dElemToCollocatedNodesBuckets );

  fillMissing2dElemToNodes( m_2dElemToCollocatedNodesBuckets, nodeManager.globalToLocalMap(), m_toNodesRelation );

  std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > const collocatedEdgeBuckets = buildCollocatedEdgeBuckets( referenceCollocatedNodes, nl2g, edgeManager.nodeList() );
  std::map< localIndex, localIndex > const referenceCollocatedEdges = buildReferenceCollocatedEdges( collocatedEdgeBuckets, edgeManager.ghostRank() );

  m_2dFaceToEdge = build2dFaceToEdge( referenceCollocatedEdges );
  m_edgesTo2dFaces = buildEdgesToFace2d( m_2dFaceToEdge.toViewConst() );

  localIndex const num2dElems = this->size();
  localIndex const num2dFaces = m_2dFaceToEdge.size();
  // Mark 2d elements and their connections as "new" so they get considered during the stencil computations.
  // This is mainly due to the dynamic process of the fracture propagation and the legacy unique way to import the fracture
  // by splitting the mesh during the first steps of the simulations.
  m_newFaceElements = makeSortedArrayIota( num2dElems );
  m_recalculateConnectionsFor2dFaces = makeSortedArrayIota( num2dFaces );

  // When a fracture element has only one (or less!) neighbor, let's try to find the other one.
  // The `ElemPath` provides all the information of a given face: obviously its face index,
  // but also the element index that touch the face, and the node indices of the face as well.
  struct ElemPath
  {
    localIndex er;
    localIndex esr;
    localIndex ei;
    localIndex face;
    std::vector< localIndex > nodes;

    bool operator<( ElemPath const & other ) const
    {
      return std::tie( er, esr, ei, face, nodes ) < std::tie( other.er, other.esr, other.ei, other.face, other.nodes );
    }
  };

  // We are building the mapping that connects all the reference (collocated) nodes of any face to the elements those nodes are touching.
  // Using this nodal information will let us reconnect the fracture 2d element to its 3d neighbor.
  std::map< std::set< globalIndex >, std::set< ElemPath > > faceRefNodesToElems;
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( [&]( localIndex const er,
                                                                         localIndex const esr,
                                                                         ElementRegionBase const & GEOS_UNUSED_PARAM( region ),
                                                                         CellElementSubRegion const & subRegion )
  {
    auto const & elemToFaces = subRegion.faceList().base();
    for( localIndex ei = 0; ei < elemToFaces.size( 0 ); ++ei )
    {
      for( auto const & face: elemToFaces[ei] )
      {
        // A set of the global indices of the nodes of the face is used as the "signature" of the face nodes.
        std::set< globalIndex > nodesOfFace;
        for( localIndex const & n: faceToNodes[face] )
        {
          auto const it = referenceCollocatedNodes.find( nl2g[n] );
          if( it != referenceCollocatedNodes.cend() )
          {
            nodesOfFace.insert( it->second );
          }
          else
          {
            // If all the nodes are not found, it makes no sense to keep on considering the current face as a candidate.
            // So we can exit the loop.
            break;
          }
        }
        // We still double-check that all the nodes of the face were found,
        // even if there's no chance for this entry to successfully match any request.
        auto const & nodes = faceToNodes[face];
        if( nodesOfFace.size() == LvArray::integerConversion< std::size_t >( nodes.size() ) )
        {
          std::vector< localIndex > const ns( nodes.begin(), nodes.end() );
          faceRefNodesToElems[nodesOfFace].insert( ElemPath{ er, esr, ei, face, ns } );
        }
      }
    }
  } );

  // Here we loop over all the elements of the fracture.
  // When there's neighbor missing, we search for a face that would lie on the collocated nodes of the fracture element.
  for( int e2d = 0; e2d < num2dElems; ++e2d )
  {
    if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) >= 2 )  // All the neighbors are known.
    {
      continue;
    }

    std::set< globalIndex > refNodes;
    if( m_toNodesRelation[e2d].size() != 0 )
    {
      for( localIndex const & n: m_toNodesRelation[e2d] )
      {
        globalIndex const & gn = nl2g[n];
        auto const it = referenceCollocatedNodes.find( gn );
        if( it != referenceCollocatedNodes.cend() )
        {
          refNodes.insert( it->second );
        }
      }
    }
    else if( m_ghostRank[e2d] < 0 )
    {
      for( integer ni = 0; ni < m_2dElemToCollocatedNodesBuckets[e2d].size(); ++ni )
      {
        array1d< globalIndex > const & bucket = m_2dElemToCollocatedNodesBuckets( e2d, ni );
        for( globalIndex const & gni: bucket )
        {
          auto const it = referenceCollocatedNodes.find( gni );
          if( it != referenceCollocatedNodes.cend() )
          {
            refNodes.insert( it->second );
          }
        }
      }
    }

    auto const match = faceRefNodesToElems.find( refNodes );
    if( match != faceRefNodesToElems.cend() )
    {
      for( ElemPath const & path: match->second )
      {
        // This `if` prevents from storing the same data twice.
        if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) == 0 || m_2dElemToElems.m_toElementIndex[e2d][0] != path.ei )
        {
          m_2dElemToElems.m_toElementRegion.emplaceBack( e2d, path.er );
          m_2dElemToElems.m_toElementSubRegion.emplaceBack( e2d, path.esr );
          m_2dElemToElems.m_toElementIndex.emplaceBack( e2d, path.ei );
          m_toFacesRelation.emplaceBack( e2d, path.face );
          for( localIndex const & n: path.nodes )
          {
            auto currentNodes = m_toNodesRelation[e2d];
            if( std::find( currentNodes.begin(), currentNodes.end(), n ) == currentNodes.end() )
            {
              m_toNodesRelation.emplaceBack( e2d, n );
            }
          }
        }
      }
    }
  }

  // Checking that each face has two neighboring elements.
  // If not, we pop up an error.
  std::vector< localIndex > isolatedFractureElements;
  for( int e2d = 0; e2d < num2dElems; ++e2d )
  {
    if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) < 2 && m_ghostRank[e2d] < 0 )
    {
      isolatedFractureElements.push_back( e2d );
    }
  }
  GEOS_ERROR_IF( !isolatedFractureElements.empty(),
                 "Fracture " << this->getName() << " has elements {" << stringutilities::join( isolatedFractureElements, ", " ) << "} with less than two neighbors." );

  fillMissing2dElemToEdges( m_toNodesRelation.toViewConst(),
                            nodeManager.edgeList().toViewConst(),
                            nl2g,
                            referenceCollocatedNodes,
                            collocatedEdgeBuckets,
                            m_toEdgesRelation );

  m_2dFaceTo2dElems = build2dFaceTo2dElems( m_toEdgesRelation.toViewConst(), m_edgesTo2dFaces, referenceCollocatedEdges );

//  fixNeighborMappingsInconsistency( getName(), m_2dElemToElems, m_toFacesRelation );

  fixNodesOrder( m_toFacesRelation.toViewConst(), faceToNodes, m_toNodesRelation );
}

void FaceElementSubRegion::inheritGhostRankFromParentFace( FaceManager const & faceManager,
                                                           std::set< localIndex > const & indices )
{
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  for( localIndex const & index: indices )
  {
    m_ghostRank[index] = faceGhostRank[m_toFacesRelation[index][0]];
  }
}

std::set< std::set< globalIndex > > FaceElementSubRegion::getCollocatedNodes() const
{
  std::set< std::set< globalIndex > > result;

  for( localIndex e2d = 0; e2d < m_2dElemToCollocatedNodesBuckets.size(); ++e2d )
  {
    for( integer ni = 0; ni < m_2dElemToCollocatedNodesBuckets[e2d].size(); ++ni )
    {
      array1d< globalIndex > const & bucket = m_2dElemToCollocatedNodesBuckets( e2d, ni );
      std::set< globalIndex > const s( bucket.begin(), bucket.end() );
      result.insert( s );
    }
  }

  return result;
}

void FaceElementSubRegion::flipFaceMap( FaceManager & faceManager,
                                        ElementRegionManager const & elemManager )
{
  ArrayOfArraysView< localIndex > const & elems2dToFaces = faceList().toView();
  arrayView2d< localIndex const > const & faceToElementRegionIndex    = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & faceToElementSubRegionIndex = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & faceToElementIndex          = faceManager.elementList();

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const cellElemGlobalIndex =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );

  forAll< parallelHostPolicy >( this->size(), [=]( localIndex const kfe )
  {
    if( elems2dToFaces.sizeOfArray( kfe ) != 2 )
    {
      return;
    }

    localIndex & f0 = elems2dToFaces[kfe][0];
    localIndex & f1 = elems2dToFaces[kfe][1];

    localIndex const er0  = faceToElementRegionIndex[f0][0];
    localIndex const esr0 = faceToElementSubRegionIndex[f0][0];
    localIndex const ek0  = faceToElementIndex[f0][0];

    localIndex const er1  = faceToElementRegionIndex[f1][0];
    localIndex const esr1 = faceToElementSubRegionIndex[f1][0];
    localIndex const ek1  = faceToElementIndex[f1][0];

    globalIndex const globalIndexElem0 = cellElemGlobalIndex[er0][esr0][ek0];
    globalIndex const globalIndexElem1 = cellElemGlobalIndex[er1][esr1][ek1];

    if( globalIndexElem0 > globalIndexElem1 )
    {
      std::swap( f0, f1 );
    }
  } );

}

void FaceElementSubRegion::fixNeighboringFacesNormals( FaceManager & faceManager,
                                                       ElementRegionManager const & elemManager )
{
  ArrayOfArraysView< localIndex > const & elems2dToFaces = faceList().toView();
  arrayView2d< localIndex const > const & faceToElementRegionIndex    = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & faceToElementSubRegionIndex = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & faceToElementIndex          = faceManager.elementList();

  arrayView2d< real64 const > const faceCenter = faceManager.faceCenter();
  FaceManager::NodeMapType & faceToNodes = faceManager.nodeList();

  auto elemCenter = elemManager.constructArrayViewAccessor< real64, 2 >( CellElementSubRegion::viewKeyStruct::elementCenterString() );

  // We need to modify the normals and the nodes ordering to be consistent.
  arrayView2d< real64 > const faceNormal = faceManager.faceNormal();
  forAll< parallelHostPolicy >( this->size(), [=, &faceToNodes]( localIndex const kfe )
  {
    if( elems2dToFaces.sizeOfArray( kfe ) != 2 )
    {
      return;
    }

    localIndex const f0 = elems2dToFaces[kfe][0];
    localIndex const f1 = elems2dToFaces[kfe][1];

    /// Note: I am assuming that the 0 element is the elementSubregion one for faces
    /// touching both a 3D and a 2D cell.
    localIndex const er0  = faceToElementRegionIndex[f0][0];
    localIndex const esr0 = faceToElementSubRegionIndex[f0][0];
    localIndex const ek0  = faceToElementIndex[f0][0];

    localIndex const er1  = faceToElementRegionIndex[f1][0];
    localIndex const esr1 = faceToElementSubRegionIndex[f1][0];
    localIndex const ek1  = faceToElementIndex[f1][0];

    real64 f0e0vector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[f0] );
    real64 f1e1vector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[f1] );

    LvArray::tensorOps::subtract< 3 >( f0e0vector, elemCenter[er0][esr0][ek0] );
    LvArray::tensorOps::subtract< 3 >( f1e1vector, elemCenter[er1][esr1][ek1] );

    // If the vector connecting the face center and the elem center is in the same
    // direction as the unit normal, we flip the normal coz it should be pointing outward
    // (i.e., towards the fracture element).
    if( LvArray::tensorOps::AiBi< 3 >( faceNormal[f0], f0e0vector ) < 0.0 )
    {
      GEOS_WARNING( GEOS_FMT( "For fracture element {}, I had to flip the normal nf0 of face {}", kfe, f0 ) );
      LvArray::tensorOps::scale< 3 >( faceNormal[f0], -1.0 );
      std::reverse( faceToNodes[f0].begin(), faceToNodes[f0].end() );
    }
    if( LvArray::tensorOps::AiBi< 3 >( faceNormal[f1], f1e1vector ) < 0.0 )
    {
      GEOS_WARNING( GEOS_FMT( "For fracture element {}, I had to flip the normal nf1 of face {}", kfe, f1 ) );
      LvArray::tensorOps::scale< 3 >( faceNormal[f1], -1.0 );
      std::reverse( faceToNodes[f1].begin(), faceToNodes[f1].end() );
    }
  } );

}

} /* namespace geos */
