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

#include "CellBlockUtilities.hpp"

#include "codingUtilities/Utilities.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/PrismUtilities.hpp"

namespace geos
{


/**
 * @brief error message to output when the number the number of node is insufficient for a given
 * face of a primitive.
 */
static const string nodeCountError = "Not enough nodes for {} element (face index = {}).\n" + string( generalMeshErrorAdvice );
/**
 * @brief error message to output when the number the number of node is insufficient for a given
 * face of a primitive.
 */
static const string faceIndexError = "Local face index out of range for {} element: face index = {}.\n" + string( generalMeshErrorAdvice );


static localIndex getFaceNodesHex( localIndex const faceNum,
                                   arraySlice1d< localIndex const, cells::NODE_MAP_USD-1 > const & elemNodes,
                                   Span< localIndex > const faceNodes )
{
  GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, "Hexahedron", faceNum ) );
  switch( faceNum )
  {
    case 0:
    {
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[1];
      faceNodes[2] = elemNodes[5];
      faceNodes[3] = elemNodes[4];
      break;
    }
    case 1:
    {
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[2];
      faceNodes[2] = elemNodes[3];
      faceNodes[3] = elemNodes[1];
      break;
    }
    case 2:
    {
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[4];
      faceNodes[2] = elemNodes[6];
      faceNodes[3] = elemNodes[2];
      break;
    }
    case 3:
    {
      faceNodes[0] = elemNodes[1];
      faceNodes[1] = elemNodes[3];
      faceNodes[2] = elemNodes[7];
      faceNodes[3] = elemNodes[5];
      break;
    }
    case 4:
    {
      faceNodes[0] = elemNodes[2];
      faceNodes[1] = elemNodes[6];
      faceNodes[2] = elemNodes[7];
      faceNodes[3] = elemNodes[3];
      break;
    }
    case 5:
    {
      faceNodes[0] = elemNodes[4];
      faceNodes[1] = elemNodes[5];
      faceNodes[2] = elemNodes[7];
      faceNodes[3] = elemNodes[6];
      break;
    }
    default:
    {
      GEOS_ERROR( GEOS_FMT( faceIndexError, "Hexahedron", faceNum ) );
    }
  }
  return 4;
}

static localIndex getFaceNodesWedge( localIndex const faceNum,
                                     arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const & elemNodes,
                                     Span< localIndex > const faceNodes )
{
  switch( faceNum )
  {
    case 0:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, "Wedge", faceNum ) );
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[1];
      faceNodes[2] = elemNodes[5];
      faceNodes[3] = elemNodes[4];
      return 4;
    }
    case 1:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, "Wedge", faceNum ) );
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[2];
      faceNodes[2] = elemNodes[3];
      faceNodes[3] = elemNodes[1];
      return 4;
    }
    case 2:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 3, GEOS_FMT( nodeCountError, "Wedge", faceNum ) );
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[4];
      faceNodes[2] = elemNodes[2];
      return 3;
    }
    case 3:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 3, GEOS_FMT( nodeCountError, "Wedge", faceNum ) );
      faceNodes[0] = elemNodes[1];
      faceNodes[1] = elemNodes[3];
      faceNodes[2] = elemNodes[5];
      return 3;
    }
    case 4:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, "Wedge", faceNum ) );
      faceNodes[0] = elemNodes[2];
      faceNodes[1] = elemNodes[4];
      faceNodes[2] = elemNodes[5];
      faceNodes[3] = elemNodes[3];
      return 4;
    }
    default:
    {
      GEOS_ERROR( GEOS_FMT( faceIndexError, "Wedge", faceNum ) );
      return 0;
    }
  }
}

static localIndex getFaceNodesTet( localIndex const faceNum,
                                   arraySlice1d< localIndex const, cells::NODE_MAP_USD-1 > const & elemNodes,
                                   Span< localIndex > const faceNodes )
{
  GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 3, GEOS_FMT( nodeCountError, "Tetrahedron", faceNum ) );
  switch( faceNum )
  {
    case 0:
    {
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[1];
      faceNodes[2] = elemNodes[3];
      break;
    }
    case 1:
    {
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[2];
      faceNodes[2] = elemNodes[1];
      break;
    }
    case 2:
    {
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[3];
      faceNodes[2] = elemNodes[2];
      break;
    }
    case 3:
    {
      faceNodes[0] = elemNodes[1];
      faceNodes[1] = elemNodes[2];
      faceNodes[2] = elemNodes[3];
      break;
    }
    default:
    {
      GEOS_ERROR( GEOS_FMT( faceIndexError, "Tetrahedron", faceNum ) );
    }
  }
  return 3;
}

static localIndex getFaceNodesPyramid( localIndex const faceNum,
                                       arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const & elemNodes,
                                       Span< localIndex > const faceNodes )
{
  switch( faceNum )
  {
    case 0:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 3, GEOS_FMT( nodeCountError, "Pyramid", faceNum ) );
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[1];
      faceNodes[2] = elemNodes[4];
      return 3;
    }
    case 1:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, "Pyramid", faceNum ) );
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[2];
      faceNodes[2] = elemNodes[3];
      faceNodes[3] = elemNodes[1];
      return 4;
    }
    case 2:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 3, GEOS_FMT( nodeCountError, "Pyramid", faceNum ) );
      faceNodes[0] = elemNodes[0];
      faceNodes[1] = elemNodes[4];
      faceNodes[2] = elemNodes[2];
      return 3;
    }
    case 3:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 3, GEOS_FMT( nodeCountError, "Pyramid", faceNum ) );
      faceNodes[0] = elemNodes[1];
      faceNodes[1] = elemNodes[3];
      faceNodes[2] = elemNodes[4];
      return 3;
    }
    case 4:
    {
      GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 3, GEOS_FMT( nodeCountError, "Pyramid", faceNum ) );
      faceNodes[0] = elemNodes[2];
      faceNodes[1] = elemNodes[4];
      faceNodes[2] = elemNodes[3];
      return 3;
    }
    default:
    {
      GEOS_ERROR( GEOS_FMT( faceIndexError, "Pyramid", faceNum ) );
      return 0;
    }
  }
}

localIndex getFaceNodes( ElementType const elementType,
                         localIndex const elemIdx,
                         localIndex const faceNumber,
                         arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elementToNodes,
                         Span< localIndex > const faceNodes )
{
  switch( elementType )
  {
    case ElementType::Hexahedron:
    {
      return getFaceNodesHex( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Wedge:
    {
      return getFaceNodesWedge( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Tetrahedron:
    {
      return getFaceNodesTet( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Pyramid:
    {
      return getFaceNodesPyramid( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Prism5:
    {
      return getFaceNodesPrism< 5 >( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Prism6:
    {
      return getFaceNodesPrism< 6 >( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Prism7:
    {
      return getFaceNodesPrism< 7 >( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Prism8:
    {
      return getFaceNodesPrism< 8 >( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Prism9:
    {
      return getFaceNodesPrism< 9 >( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Prism10:
    {
      return getFaceNodesPrism< 10 >( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    case ElementType::Prism11:
    {
      return getFaceNodesPrism< 11 >( faceNumber, elementToNodes[elemIdx], faceNodes );
    }
    default:
    {
      GEOS_ERROR( "Invalid element type " << elementType << " at face index " << faceNumber << ".\n" << generalMeshErrorAdvice );
    }
  }
  return 0;
}

/**
 * @brief Stores the data necessary to construct the various edge maps.
 *
 * These records are created and stored when visiting edges through faces.
 * Organizing and sorting them by both nodes enables identification of
 * matching edges of different faces. Face indices are preserved
 * so that face-edge maps may be constructed.
 */
struct EdgeBuilder
{
  /**
   * @brief Constructor.
   * @param [in] node1 the greater of the two node indices that comprise the edge.
   * @param [in] face the ID of the face this edge came from.
   * @param [in] edgeNum the face local index of this edge.
   */
  EdgeBuilder( localIndex const node1,
               localIndex const face,
               localIndex const edgeNum ):
    n1( node1 ),
    faceIndex( face ),
    edgeNumber( edgeNum )
  {}

  /// The larger of the two node indices that comprise the edge.
  localIndex n1;

  /// The face the edge came from.
  localIndex faceIndex;

  /// The face local index of the edge.
  localIndex edgeNumber;

private:

  /**
   * @brief Equality comparison operator.
   * @param [in] lhs left-hand side of the comparison
   * @param [in] rhs right-hand side of the comparison
   *
   * Two edges are considered equal if they share 2 nodes.
   * Because the other data structure already groups entries by lowest node,
   * here we only store and compare the second node.
   */
  friend bool operator==( EdgeBuilder const & lhs, EdgeBuilder const & rhs )
  { return lhs.n1 == rhs.n1; }
};

/**
 * @brief Populate the edgesByLowestNode map.
 * @param [in] numNodes The number of nodes.
 * @param [in] faceToNodeMap a map that associates an ordered list of nodes with each face.
 * @return a map of each node to edges that have it as their lowest index node, sorted by second node index.
 */
ArrayOfArrays< EdgeBuilder >
createEdgesByLowestNode( localIndex const numNodes,
                         ArrayOfArraysView< localIndex const > const & faceToNodeMap )
{
  array1d< localIndex > edgeCounts( numNodes );
  forAll< parallelHostPolicy >( faceToNodeMap.size(), [counts = edgeCounts.toView(),
                                                       faceToNodeMap]( localIndex const faceIndex )
  {
    // loop over all the nodes in the face. there will be an edge for each node.
    localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceIndex );
    for( localIndex a = 0; a < numNodesInFace; ++a )
    {
      // count the edge for its lowest index node
      localIndex const node0 = faceToNodeMap( faceIndex, a );
      localIndex const node1 = faceToNodeMap( faceIndex, ( a + 1 ) % numNodesInFace );
      RAJA::atomicInc< parallelHostAtomic >( &counts[ std::min( node0, node1 ) ] );
    }
  } );

  ArrayOfArrays< EdgeBuilder > edgesByLowestNode;
  edgesByLowestNode.resizeFromCapacities< parallelHostPolicy >( numNodes, edgeCounts.data() );

  forAll< parallelHostPolicy >( faceToNodeMap.size(), [&]( localIndex const faceIndex )
  {
    localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceIndex );

    // loop over all the nodes in the face. there will be an edge for each node.
    for( localIndex a = 0; a < numNodesInFace; ++a )
    {
      // sort the nodes in order of index value.
      localIndex node0 = faceToNodeMap( faceIndex, a );
      localIndex node1 = faceToNodeMap( faceIndex, ( a + 1 ) % numNodesInFace );
      if( node0 > node1 )
      {
        std::swap( node0, node1 );
      }

      // And append the edge to edgesByLowestNode.
      edgesByLowestNode.emplaceBackAtomic< parallelHostAtomic >( node0, node1, faceIndex, a );
    }
  } );

  // This comparator is not attached to EdgeBuilder because of potential inconsistencies with operator==
  auto const comp = []( EdgeBuilder const & e0, EdgeBuilder const & e1 ) -> bool
  {
    return ( e0.n1 < e1.n1 ) || ( ( e0.n1 == e1.n1 ) && ( e0.faceIndex < e1.faceIndex ) );
  };

  // Loop over all the nodes and sort the associated edges.
  forAll< parallelHostPolicy >( numNodes, [edgesByLowestNode = edgesByLowestNode.toView(),
                                           comp]( localIndex const nodeIndex )
  {
    arraySlice1d< EdgeBuilder > const edges = edgesByLowestNode[ nodeIndex ];
    std::sort( edges.begin(), edges.end(), comp );
  } );

  return edgesByLowestNode;
}


/**
 * @brief Populate the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an array containing the unique ID of the first edge associated with each node.
 * @param [in] faceToEdgeMap the map from faces to nodes.
 * @param [in,out] edgeToFaceMap the map from edges to faces.
 * @param [in,out] edgeToNodeMap the map from edges to nodes.
 */
void populateEdgeMaps( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                       arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                       ArrayOfArraysView< localIndex > const & faceToEdgeMap,
                       ArrayOfArraysView< localIndex > const & edgeToFaceMap,
                       arrayView2d< localIndex > const & edgeToNodeMap )
{
  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  GEOS_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );
  GEOS_ERROR_IF_NE( numUniqueEdges, edgeToFaceMap.size() );
  GEOS_ERROR_IF_NE( numUniqueEdges, edgeToNodeMap.size( 0 ) );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeIndex )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeIndex ];
    arraySlice1d< EdgeBuilder const > const edges = edgesByLowestNode[ nodeIndex ];
    forEqualRanges( edges.begin(), edges.end(), [&]( auto first, auto last )
    {
      // Populate the edge to node map.
      edgeToNodeMap( curEdgeID, 0 ) = nodeIndex;
      edgeToNodeMap( curEdgeID, 1 ) = first->n1;

      // Loop through all the matches and fill in the face to edge and edge to face maps.
      while( first != last )
      {
        EdgeBuilder const & edge = *first;
        faceToEdgeMap( edge.faceIndex, edge.edgeNumber ) = curEdgeID;
        edgeToFaceMap.emplaceBack( curEdgeID, edge.faceIndex );
        ++first;
      }

      ++curEdgeID;
    } );
  } );
}

/**
 * @brief Resize the edge to face and edge to node maps.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an containing the unique edge IDs for each node in edgesByLowestNode.
 * @param [out] faceToEdgeMap the map from faces to edges, resized appropriately and filled with dummy values.
 * @param [out] edgeToFaceMap the map from edges to faces, resized appropriately.
 * @param [out] edgeToNodeMap the map from edges to nodes, resized appropriately.
 */
void resizeEdgeMaps( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                     arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                     ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                     ArrayOfArrays< localIndex > & faceToEdgeMap,
                     ArrayOfArrays< localIndex > & edgeToFaceMap,
                     array2d< localIndex > & edgeToNodeMap )
{
  localIndex const numFaces = faceToNodeMap.size();
  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();

  // Count the number of faces adjacent to each edge by traversing previously built data structures
  array1d< localIndex > numFacesPerEdge( numUniqueEdges );
  forAll< parallelHostPolicy >( numNodes, [uniqueEdgeOffsets,
                                           edgesByLowestNode,
                                           numFacesPerEdge = numFacesPerEdge.toView()]( localIndex const nodeIndex )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeIndex ];
    arraySlice1d< EdgeBuilder const > const edges = edgesByLowestNode[ nodeIndex ];
    forUniqueValues( edges.begin(), edges.end(), [&]( EdgeBuilder const &, localIndex const numMatches )
    {
      numFacesPerEdge( curEdgeID++ ) = numMatches + CellBlockManagerABC::faceMapExtraSpacePerEdge();
    } );
  } );
  edgeToFaceMap.resizeFromCapacities< parallelHostPolicy >( numUniqueEdges, numFacesPerEdge.data() );

  // This relies on the fact that #nodes = #edges for each face.
  // Thus, by utilizing faceToNodeMap, we can avoid atomic counting of capacity.
  array1d< localIndex > numEdgesPerFace( numFaces );
  forAll< parallelHostPolicy >( numFaces, [faceToNodeMap,
                                           numEdgesPerFace = numEdgesPerFace.toView()] ( localIndex const faceIndex )
  {
    localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceIndex );
    numEdgesPerFace[ faceIndex ] = numNodesInFace + CellBlockManagerABC::edgeMapExtraSpacePerFace();
  } );
  faceToEdgeMap.resizeFromCapacities< parallelHostPolicy >( numFaces, numEdgesPerFace.data() );
  forAll< parallelHostPolicy >( numFaces, [numEdgesPerFace = numEdgesPerFace.toViewConst(),
                                           faceToEdgeMap = faceToEdgeMap.toView()] ( localIndex const faceIndex )
  {
    // There is no API to set size of each sub-array within its capacity in parallel
    for( localIndex i = 0; i < numEdgesPerFace[ faceIndex ] - CellBlockManagerABC::edgeMapExtraSpacePerFace(); ++i )
    {
      faceToEdgeMap.emplaceBack( faceIndex, -1 );
    }
  } );

  // Each edge is adjacent to strictly 2 nodes.
  edgeToNodeMap.resize( numUniqueEdges, 2 );
}

localIndex buildEdgeMaps( localIndex const numNodes,
                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                          ArrayOfArrays< localIndex > & faceToEdgeMap,
                          ArrayOfArrays< localIndex > & edgeToFaceMap,
                          array2d< localIndex > & edgeToNodeMap )
{
  GEOS_MARK_FUNCTION;

  ArrayOfArrays< EdgeBuilder > const edgesByLowestNode =
    createEdgesByLowestNode( numNodes, faceToNodeMap );

  array1d< localIndex > const uniqueEdgeOffsets =
    computeUniqueValueOffsets< parallelHostPolicy >( edgesByLowestNode.toViewConst() );

  resizeEdgeMaps( edgesByLowestNode.toViewConst(),
                  uniqueEdgeOffsets,
                  faceToNodeMap,
                  faceToEdgeMap,
                  edgeToFaceMap,
                  edgeToNodeMap );

  populateEdgeMaps( edgesByLowestNode.toViewConst(),
                    uniqueEdgeOffsets,
                    faceToEdgeMap.toView(),
                    edgeToFaceMap.toView(),
                    edgeToNodeMap.toView() );

  return uniqueEdgeOffsets.back();
}

}
