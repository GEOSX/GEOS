/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_MESH_GENERATORS_CELLBLOCKUTILITIES_HPP_
#define GEOS_MESH_GENERATORS_CELLBLOCKUTILITIES_HPP_

#include "mesh/ElementType.hpp"
#include "common/DataTypes.hpp"
#include "common/Span.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

/**
 * @brief Container for maps from a mesh object (node, edge or face) to cells.
 * @tparam T underlying map type
 */
template< typename T >
struct ToCellRelation
{
  T toBlockIndex; ///< Map containing a list of cell block indices for each object
  T toCellIndex;  ///< Map containing cell indices, same shape as above

  /**
   * @brief Constructor by values.
   * @param toBlockIndex_ Map containing a list of cell block indices for each object
   * @param toCellIndex_ Map containing cell indices, same shape as above
   */
  ToCellRelation( T const & toBlockIndex_,
                  T const & toCellIndex_ )
    : toBlockIndex( toBlockIndex_ ),
    toCellIndex( toCellIndex_ )
  { }

  /**
   * @brief Constructor from moved values.
   * @param toBlockIndex_ Map containing a list of cell block indices for each object
   * @param toCellIndex_ Map containing cell indices, same shape as above
   */
  ToCellRelation( T && toBlockIndex_,
                  T && toCellIndex_ )
    : toBlockIndex( toBlockIndex_ ),
    toCellIndex( toCellIndex_ )
  { }

  ToCellRelation() = default;
};

/**
 * @brief Free function that generates face to edges, edge to faces and edge to nodes mappings.
 * @param[in] numNodes The number of nodes.
 * @param[in] faceToNodeMap Face to node mappings as an input.
 * @param[out] faceToEdgeMap Face to edges will be resized and filled.
 * @param[out] edgeToFaceMap Ege to faces will be resized and filled.
 * @param[out] edgeToNodeMap Edge to nodes will be resized and filled.
 * @return The number of edges.
 */
localIndex buildEdgeMaps( localIndex numNodes,
                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                          ArrayOfArrays< localIndex > & faceToEdgeMap,
                          ArrayOfArrays< localIndex > & edgeToFaceMap,
                          array2d< localIndex > & edgeToNodeMap );

/**
 * @brief Get the local indices of the nodes in a face of the element.
 * @param[in] elementType Type of the element
 * @param[in] elemIdx the local index of the target element
 * @param[in] faceNumber the local index of the target face in the element  (this will be 0-numFacesInElement)
 * @param[in] elementToNodes Element to nodes mapping.
 * @param[out] faceNodes space to write node indices to, must be of sufficient size
 * @return number of nodes in the face
 */
localIndex getFaceNodes( ElementType const elementType,
                         localIndex const elemIdx,
                         localIndex const faceNumber,
                         arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elementToNodes,
                         Span< localIndex > const faceNodes );

/**
 * @brief Find and count ranges of repeated values in an array of sorted arrays and compute offsets.
 * @tparam POLICY execution policy
 * @tparam T value type of input arrays
 * @tparam COMP A comparator used to compare the values.
 * @param [in] sortedLists the input array of sorted arrays of values.
 * @param [in] comp the comparator used to compare the values in @p sortedLists.
 * @return an array of size @p sortedLists.size() + 1, where element i contains starting index of
 *         unique values from sub-array i in a global list of unique values, while the last value
 *         contains the total number of unique edges (i.e. an exclusive scan + total reduction).
 */
template< typename POLICY, typename T, typename COMP = std::equal_to<> >
array1d< localIndex >
computeUniqueValueOffsets( ArrayOfArraysView< T const > const & sortedLists, COMP && comp = {} )
{
  localIndex const numNodes = sortedLists.size();
  array1d< localIndex > uniqueValueOffsets( numNodes + 1 );

  // For each node, count number of unique edges that have the node as its lowest
  arrayView1d< localIndex > const numUniqueValuesView = uniqueValueOffsets.toView();
  forAll< POLICY >( numNodes, [sortedLists, numUniqueValuesView, comp=std::forward< COMP >( comp )]( localIndex const i )
  {
    arraySlice1d< T const > const list = sortedLists[ i ];
    forEqualRanges( list.begin(), list.end(), [&]( auto, auto )
    {
      ++numUniqueValuesView[ i + 1 ];
    }, comp );
  } );

  // Perform an inplace prefix-sum to get the unique edge offset.
  RAJA::inclusive_scan_inplace< POLICY >( RAJA::make_span( uniqueValueOffsets.data(), uniqueValueOffsets.size() ) );
  return uniqueValueOffsets;
}

}

#endif // GEOS_MESH_GENERATORS_CELLBLOCKUTILITIES_HPP_
