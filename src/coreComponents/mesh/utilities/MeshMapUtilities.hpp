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

/**
 * @file MeshMapUtilities.hpp
 */
#ifndef GEOS_MESH_UTILITIES_MESHMAPUTILITIES_HPP
#define GEOS_MESH_UTILITIES_MESHMAPUTILITIES_HPP

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geos
{

/**
 * @brief This namespace contains helper functions that facilitate access
 *        into the assortment of maps used by GEOSX mesh object managers
 *        (e.g. array2d/array1d(array1d)/ArrayOfArrays/ArrayOfSets, etc.)
 */
namespace meshMapUtilities
{

//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief @return the size of the map along first dimension
 * @tparam T type of map element
 * @tparam USD unit-stride dimension of the map
 * @param map reference to the map
 */
template< typename T, int USD >
GEOS_HOST_DEVICE
inline localIndex size0( arrayView2d< T, USD > const & map )
{
  return map.size( 0 );
}

/**
 * @brief @return the size of the map along first dimension
 * @tparam T type of map element
 * @param map reference to the map
 */
template< typename T >
GEOS_HOST_DEVICE
inline localIndex size0( ArrayOfArraysView< T > const & map )
{
  return map.size();
}

/**
 * @brief @return the size of the map along first dimension
 * @tparam T type of map element
 * @param map reference to the map
 */
template< typename T >
GEOS_HOST_DEVICE
inline localIndex size0( ArrayOfSetsView< T > const & map )
{
  return map.size();
}

//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Transposes an input map (array2d, ArrayOfArrays or ArrayOfSets)
 * @tparam POLICY execution policy to use
 * @tparam VIEW_TYPE type of view of the source map
 * @param srcMap the source map
 * @param dstSize number of target objects ("cols" of @p srcMap)
 * @param overAlloc overallocation (extra capacity per row) for resulting map
 * @return the transpose of @p srcMap stored as ArrayOfArrays (most general form)
 */
template< typename POLICY, typename VIEW_TYPE >
ArrayOfArrays< std::remove_const_t< typename VIEW_TYPE::ValueType > >
transposeIndexMap( VIEW_TYPE const & srcMap,
                   localIndex const dstSize,
                   localIndex const overAlloc = 0 )
{
  // Count the number of elements in each set
  array1d< localIndex > counts( dstSize );
  counts.setValues< POLICY >( overAlloc );
  forAll< POLICY >( size0( srcMap ), [srcMap, counts = counts.toView()] ( localIndex const srcIndex )
  {
    for( localIndex const dstIndex : srcMap[ srcIndex ] )
    {
      RAJA::atomicInc< AtomicPolicy< POLICY > >( &counts[ dstIndex ] );
    }
  } );

  // Allocate storage for the transpose map
  ArrayOfArrays< localIndex > dstMap;
  dstMap.resizeFromCapacities< parallelHostPolicy >( dstSize, counts.data() );

  // Fill the sub-arrays with unsorted entries
  forAll< POLICY >( size0( srcMap ), [srcMap, dstMap = dstMap.toView()] ( localIndex const srcIndex )
  {
    for( localIndex const dstIndex : srcMap[ srcIndex ] )
    {
      dstMap.emplaceBackAtomic< AtomicPolicy< POLICY > >( dstIndex, srcIndex );
    }
  } );

  return dstMap;
}

/**
 * @brief Convert ToCellRelation into ToElementRelation.
 * @tparam POLICY execution policy
 * @param blockToSubRegion a map from cell blocks to region/subregion pairs
 * @param srcMap source map (object-to-cells)
 * @param dstMap target map (object-to-elements)
 */
template< typename POLICY >
void transformCellBlockToRegionMap( arrayView2d< localIndex const > const & blockToSubRegion,
                                    ToCellRelation< ArrayOfArrays< localIndex > > const & srcMap,
                                    ToElementRelation< ArrayOfArrays< localIndex > > & dstMap )
{
  GEOS_ASSERT_EQ( blockToSubRegion.size( 1 ), 2 );
  localIndex const numObjects = srcMap.toCellIndex.size();

  localIndex const * offsets = srcMap.toCellIndex.toViewConst().getOffsets();
  dstMap.m_toElementRegion.resizeFromOffsets( numObjects, offsets );
  dstMap.m_toElementSubRegion.resizeFromOffsets( numObjects, offsets );
  dstMap.m_toElementIndex.resizeFromOffsets( numObjects, offsets );

  forAll< POLICY >( numObjects, [toCell = srcMap.toCellIndex.toViewConst(),
                                 toBlock = srcMap.toBlockIndex.toViewConst(),
                                 toRegion = dstMap.m_toElementRegion.toView(),
                                 toSubRegion = dstMap.m_toElementSubRegion.toView(),
                                 toElement = dstMap.m_toElementIndex.toView(),
                                 blockToSubRegion]( localIndex const objIndex )
  {
    arraySlice1d< localIndex const > const cells = toCell[objIndex];
    for( localIndex i = 0; i < cells.size(); ++i )
    {
      localIndex const blockIndex = toBlock( objIndex, i );
      localIndex const er = blockToSubRegion( blockIndex, 0 );
      localIndex const esr = blockToSubRegion( blockIndex, 1 );

      // Check needed because some blocks may remain unused
      if( er >= 0 && esr >= 0 )
      {
        toRegion.emplaceBack( objIndex, er );
        toSubRegion.emplaceBack( objIndex, esr );
        toElement.emplaceBack( objIndex, cells[i] );
      }
    }
  } );
}

/**
 * @brief Convert ToCellRelation into ToElementRelation.
 * @tparam POLICY execution policy
 * @param blockToSubRegion a map from cell blocks to region/subregion pairs
 * @param srcMap source map (object-to-cells)
 * @param dstMap target map (object-to-elements)
 */
template< typename POLICY, typename PERM1, typename PERM2 >
void transformCellBlockToRegionMap( arrayView2d< localIndex const > const & blockToSubRegion,
                                    ToCellRelation< array2d< localIndex, PERM1 > > const & srcMap,
                                    ToElementRelation< array2d< localIndex, PERM2 > > & dstMap )
{
  GEOS_ASSERT_EQ( blockToSubRegion.size( 1 ), 2 );
  localIndex const numObjects = srcMap.toCellIndex.size( 0 );
  localIndex const maxCellsPerObject = srcMap.toCellIndex.size( 1 );

  dstMap.m_toElementRegion.resize( numObjects, maxCellsPerObject );
  dstMap.m_toElementSubRegion.resize( numObjects, maxCellsPerObject );
  dstMap.m_toElementIndex.resize( numObjects, maxCellsPerObject );

  forAll< POLICY >( numObjects, [toCell = srcMap.toCellIndex.toViewConst(),
                                 toBlock = srcMap.toBlockIndex.toViewConst(),
                                 toRegion = dstMap.m_toElementRegion.toView(),
                                 toSubRegion = dstMap.m_toElementSubRegion.toView(),
                                 toElement = dstMap.m_toElementIndex.toView(),
                                 blockToSubRegion,
                                 maxCellsPerObject]( localIndex const objIndex )
  {
    arraySlice1d< localIndex const > const cells = toCell[objIndex];
    localIndex cellCount = 0;
    for( localIndex i = 0; i < maxCellsPerObject && cells[i] >= 0; ++i )
    {
      localIndex const blockIndex = toBlock( objIndex, i );
      localIndex const er = blockToSubRegion( blockIndex, 0 );
      localIndex const esr = blockToSubRegion( blockIndex, 1 );

      // Check needed because some blocks may remain unused
      if( er >= 0 && esr >= 0 )
      {
        toRegion( objIndex, cellCount ) = er;
        toSubRegion( objIndex, cellCount ) = esr;
        toElement( objIndex, cellCount ) = cells[i];
        ++cellCount;
      }
    }
  } );
}

} // namespace meshMapUtilities

/**
 * @brief Strucure used to hash interpolation arrays representing high-order nodes.
 * @tparam T type of node index, usually a local or global index
 */
template< typename T >
struct NodeKeyHasher
{
  /**
   * @brief @return the hash of an interpolation array representing a high-order node.
   * @param arr the array corresponding to the node to be hashed
   */
  std::size_t operator()( const std::array< T, 6 > & arr ) const
  {
    std::size_t hash = 0;
    // use a boost-style hash function
    for( auto v : arr )
    {
      hash ^= std::hash< T >{} ( v )  + 0x9e3779b9 + ( hash << 6 ) + ( hash >> 2 );
    }
    return hash;
  }
};

/**
 * @brief @return a unique interpolation array representing a high-order node coincident with a mesh vertex.
 * @tparam T type of node index, usually a local or global index
 * @param v the mesh vertex on which the high-order node lies.
 */
template< typename T >
static std::array< T, 6 > createNodeKey( T v )
{
  return std::array< T, 6 > { v, -1, -1, -1, -1, -1 };
}

/**
 * @brief @return a unique interpolation array representing a high-order node on an edge.
 * @tparam T type of node index, usually a local or global index
 * @param v1 the first mesh vertex defining the edge.
 * @param v2 the second mesh vertex defining the edge.
 * @param a the interpolation parameter, meaning that the node is 'a' steps away from v1 towards v2
 * @param order the order of the discretization
 */
template< typename T >
static std::array< T, 6 > createNodeKey( T v1, T v2, int a, int order )
{
  if( a == 0 )
    return createNodeKey( v1 );
  if( a == order )
    return createNodeKey( v2 );
  if( v1 < v2 )
  {
    return std::array< T, 6 > { v1, v2, -1, -1, a, -1 };
  }
  else
  {
    return std::array< T, 6 > { v2, v1, -1, -1, order - a, -1 };
  }
}

/**
 * @brief @return a unique interpolation array representing a high-order node on an face.
 * @tparam T type of node index, usually a local or global index
 * @param v1 the first mesh vertex defining the face.
 * @param v2 the second mesh vertex defining the face.
 * @param v3 the third mesh vertex defining the face.
 * @param v4 the fourth mesh vertex defining the face.
 * @param a the first interpolation parameter, meaning that the node is 'a' steps away from v1 towards v2 (and v3 towards v4)
 * @param b the second interpolation parameter, meaning that the node is 'b' steps away from v1 towards v2 (and v3 towards v4)
 * @param order the order of the discretization
 */
template< typename T >
static std::array< T, 6 > createNodeKey( T v1, T v2, T v3, T v4, int a, int b, int order )
{
  if( a == 0 )
    return createNodeKey( v1, v3, b, order );
  if( a == order )
    return createNodeKey( v2, v4, b, order );
  if( b == 0 )
    return createNodeKey( v1, v2, a, order );
  if( b == order )
    return createNodeKey( v3, v4, a, order );
  // arrange the vertices of the face such that v1 is the lowest value, and v2 is lower than v3
  // this ensures a coherent orientation of all face nodes
  while( v1 > v2 || v1 > v3 || v1 > v4 || v2 > v3 )
  {
    if( v1 > v2 )
    {
      std::swap( v1, v2 );
      std::swap( v3, v4 );
      a = order - a;
    }
    if( v1 > v3 )
    {
      std::swap( v1, v3 );
      std::swap( v2, v4 );
      b = order - b;
    }
    if( v1 > v4 )
    {
      std::swap( v1, v4 );
      std::swap( a, b );
      a = order - a;
      b = order - b;
    }
    if( v2 > v3 )
    {
      std::swap( v2, v3 );
      std::swap( a, b );
    }
  }
  return std::array< T, 6 > { v1, v2, v3, v4, a, b };
}

/**
 * @brief @return the unique interpolation array representing a Gauss-Lobatto node inside an element.
 * @tparam T type of node index, usually a local or global index
 * @param elemNodes indices of the nodes defining the element
 * @param q1 first interpolation parameter
 * @param q2 second interpolation parameter
 * @param q3 third interpolation parameter
 * @param order the order of the discretization
 */
template< typename T >
static std::array< T, 6 > createNodeKey( T const (&elemNodes)[ 8 ], int q1, int q2, int q3, int order )
{
  bool extremal1 = q1 == 0 || q1 == order;
  bool extremal2 = q2 == 0 || q2 == order;
  bool extremal3 = q3 == 0 || q3 == order;
  int v1 = q1/order;
  int v2 = q2/order;
  int v3 = q3/order;
  if( extremal1 && extremal2 && extremal3 )
  {
    // vertex node
    return createNodeKey( elemNodes[ v1 + 2*v2 + 4*v3 ] );
  }
  else if( extremal1 && extremal2 )
  {
    // edge node on v1, v2
    return createNodeKey( elemNodes[ v1 + 2*v2 ], elemNodes[ v1 + 2*v2 + 4 ], q3, order );
  }
  else if( extremal1 && extremal3 )
  {
    // edge node on v1, v3
    return createNodeKey( elemNodes[ v1 + 4*v3 ], elemNodes[ v1 + 2 + 4*v3 ], q2, order );
  }
  else if( extremal2 && extremal3 )
  {
    // edge node on v2, v3
    return createNodeKey( elemNodes[ 2*v2 + 4*v3 ], elemNodes[ 1 + 2*v2 + 4*v3 ], q1, order );
  }
  else if( extremal1 )
  {
    // face node on the face of type 1
    return createNodeKey( elemNodes[ v1 ], elemNodes[ v1 + 2 ], elemNodes[ v1 + 4 ], elemNodes[ v1 + 2 + 4 ], q2, q3, order );
  }
  else if( extremal2 )
  {
    // face node on the face of type 2
    return createNodeKey( elemNodes[ 2*v2 ], elemNodes[ 1 + 2*v2 ], elemNodes[ 2*v2 + 4 ], elemNodes[ 1 + 2*v2 + 4 ], q1, q3, order );
  }
  else if( extremal3 )
  {
    // face node on the face of type 3
    return createNodeKey( elemNodes[ 4*v3 ], elemNodes[ 1 + 4*v3 ], elemNodes[ 2 + 4*v3 ], elemNodes[ 1 + 2 + 4*v3 ], q1, q2, order );
  }
  else
  {
    // node internal to the cell -- no need for key, it will be created
    return createNodeKey( -1 );
  }
}

} // namespace geos

#endif //GEOS_MESH_UTILITIES_MESHMAPUTILITIES_HPP
