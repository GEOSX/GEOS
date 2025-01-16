/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */



#ifndef GEOS_DATAREPOSITORY_BUFFEROPS_HPP_
#define GEOS_DATAREPOSITORY_BUFFEROPS_HPP_

#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/traits.hpp"
#include "LvArray/src/limits.hpp"

#include <type_traits>

namespace geos
{

template< typename T >
class InterObjectRelation;

namespace bufferOps
{

template< typename T >
constexpr bool is_host_packable_scalar_v = std::is_trivial< T >::value || std::is_arithmetic< T >::value;

/// Whether an object of type T is itself packable
template< typename T >
constexpr bool is_host_packable_object_v = is_host_packable_scalar_v< T > ||
                                           traits::is_tensorT< T > ||
                                           traits::is_string< T >;

template< typename T >
constexpr bool is_container_v = !is_host_packable_object_v< T >;

namespace internal
{

// Forward decl so we can use this for contained types
template< typename T >
struct is_host_packable_helper;

}

/// Whether an object is an lvarray array/arrayview/arrayslice/arrayofarrays which ultimately contains packable objects when fully indexed
template< typename >
constexpr bool is_host_packable_array_v = false;

template< typename T, int NDIM, typename PERMUTATION >
constexpr bool is_host_packable_array_v< Array< T, NDIM, PERMUTATION > > = internal::is_host_packable_helper< T >::value;

template< typename T, int NDIM, int USD >
constexpr bool is_host_packable_array_v< ArrayView< T, NDIM, USD > > = internal::is_host_packable_helper< T >::value;

template< typename T, int NDIM, int USD >
constexpr bool is_host_packable_array_v< ArraySlice< T, NDIM, USD > > = internal::is_host_packable_helper< T >::value;

template< typename T >
constexpr bool is_host_packable_array_v< ArrayOfArrays< T > > = internal::is_host_packable_helper< T >::value;


/// Whether an object is an lvarray sortedarray which contains packable objects
template< typename >
constexpr bool is_host_packable_set_v = false;

template< typename T >
constexpr bool is_host_packable_set_v< SortedArray< T > > = internal::is_host_packable_helper< T >::value;


/// Whether an object is a map for which the keys and values are packable objects
template< typename >
constexpr bool is_host_packable_map_v = false;

template< typename T_KEY, typename T_VAL, typename SORTED >
constexpr bool is_host_packable_map_v< mapBase< T_KEY, T_VAL, SORTED > > = internal::is_host_packable_helper< T_KEY >::value &&
                                                                           internal::is_host_packable_helper< T_VAL >::value;

namespace internal
{

template< typename T >
struct is_host_packable_helper
{
  static constexpr bool value = is_host_packable_object_v< T > ||
                                is_host_packable_array_v< T > ||
                                is_host_packable_map_v< T > ||
                                is_host_packable_set_v< T >;
};

}

/// Whether the object is itself host packable
template< typename T >
constexpr bool is_host_packable_v = internal::is_host_packable_helper< T >::value;


/// Whether the object can be indexed to pack subsets of the object
template< typename T >
constexpr bool is_host_packable_by_index_v = is_host_packable_array_v< T >;


/// Whether the object is a map for which the keys are directly packable and the values can be packed by index
template< typename >
constexpr bool is_host_packable_map_by_index_v = false;

template< typename T_KEY, typename T_VAL, typename SORTED >
constexpr bool is_host_packable_map_by_index_v< mapBase< T_KEY, T_VAL, SORTED > > = is_host_packable_v< T_KEY > &&
                                                                                    is_host_packable_by_index_v< T_VAL >;


//------------------------------------------------------------------------------
// Pack(buffer,var)
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< is_host_packable_scalar_v< T >, localIndex >::type
PackData( buffer_unit_type * & buffer,
          T const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< is_host_packable_scalar_v< T >, localIndex >::type
Pack( buffer_unit_type * & buffer,
      T const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
PackData( buffer_unit_type * & buffer,
          const string & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type * & buffer,
      const string & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
PackData( buffer_unit_type * & buffer,
          SortedArray< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type * & buffer,
      SortedArray< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< traits::is_tensorT< T >, localIndex >::type
PackData( buffer_unit_type * & buffer,
          T const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< traits::is_tensorT< T >, localIndex >::type
Pack( buffer_unit_type * & buffer,
      T const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
PackData( buffer_unit_type * & buffer,
          ArrayView< T, NDIM, USD > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
Pack( buffer_unit_type * & buffer,
      ArrayView< T, NDIM, USD > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
PackData( buffer_unit_type * & buffer,
          ArrayOfArrays< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArrays< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
PackData( buffer_unit_type * & buffer,
          ArrayOfSets< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfSets< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
PackData( buffer_unit_type * & buffer,
          MAP_TYPE const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type * & buffer,
      MAP_TYPE const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
PackData( buffer_unit_type * & buffer,
          std::pair< T_FIRST, T_SECOND > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( buffer_unit_type * & buffer,
      std::pair< T_FIRST, T_SECOND > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
PackData( buffer_unit_type * & buffer,
          InterObjectRelation< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type * & buffer,
      InterObjectRelation< T > const & var );

//------------------------------------------------------------------------------
// fallthrough-implementation
template< bool DO_PACKING, typename T >
typename std::enable_if< !is_host_packable_v< T >, localIndex >::type
PackData( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
          T const & GEOS_UNUSED_PARAM( var ) )
{
  GEOS_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
// fallthrough-implementation
template< bool DO_PACKING, typename T >
typename std::enable_if< !is_host_packable_v< T >, localIndex >::type
Pack( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
      T const & GEOS_UNUSED_PARAM( var ) )
{
  GEOS_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
// PackArray(buffer,var,length)
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackPointerData( buffer_unit_type * & buffer,
                 T const * const GEOS_RESTRICT var,
                 INDEX_TYPE const length );

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackPointer( buffer_unit_type * & buffer,
             T const * const GEOS_RESTRICT var,
             INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackPointerData( buffer_unit_type * & buffer,
                 T const * const GEOS_RESTRICT var,
                 INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackPointer( buffer_unit_type * & buffer,
             T const * const GEOS_RESTRICT var,
             INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackArray( buffer_unit_type * & buffer,
           arraySlice1d< T, USD > const & var,
           INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackArray( buffer_unit_type * & buffer,
           arraySlice1d< T, USD > const & var,
           INDEX_TYPE const length );

//------------------------------------------------------------------------------
// PackByIndex(buffer,var,indices)
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_indices >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
PackDataByIndex( buffer_unit_type * & buffer,
                 ArrayView< T, NDIM, USD > const & var,
                 const T_indices & indices );

template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_indices >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
PackByIndex( buffer_unit_type * & buffer,
             ArrayView< T, NDIM, USD > const & var,
             const T_indices & indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_indices >
localIndex PackDataByIndex( buffer_unit_type * & buffer,
                            ArrayOfArrays< T > const & var,
                            T_indices const & indices );

template< bool DO_PACKING, typename T, typename T_indices >
localIndex PackByIndex( buffer_unit_type * & buffer,
                        ArrayOfArrays< T > const & var,
                        T_indices const & indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
PackDataByIndex( buffer_unit_type * & buffer,
                 MAP_TYPE const & var,
                 T_INDICES const & indices );

template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
PackByIndex( buffer_unit_type * & buffer,
             MAP_TYPE const & var,
             T_INDICES const & indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
typename std::enable_if< !is_host_packable_by_index_v< T > && !is_host_packable_map_by_index_v< T >, localIndex >::type
PackDataByIndex( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
                 T const & GEOS_UNUSED_PARAM( var ),
                 T_INDICES const & GEOS_UNUSED_PARAM( indices ) )
{
  GEOS_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") but type is not packable by index." );
  return 0;
}


template< bool DO_PACKING, typename T, typename T_INDICES >
typename std::enable_if< !is_host_packable_by_index_v< T > && !is_host_packable_map_by_index_v< T >, localIndex >::type
PackByIndex( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
             T const & GEOS_UNUSED_PARAM( var ),
             T_INDICES const & GEOS_UNUSED_PARAM( indices ) )
{
  GEOS_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
// Unpack(buffer,var)
//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< is_host_packable_scalar_v< T >, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        T & var,
        MPI_Op op );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        string & var,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< traits::is_tensorT< T >, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        T & var,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const * & buffer,
        SortedArray< T > & var,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T, int NDIM, typename PERMUTATION >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        Array< T, NDIM, PERMUTATION > & var,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
localIndex Unpack( buffer_unit_type const * & buffer,
                   ArrayOfArrays< T > & var,
                   MPI_Op op );

//------------------------------------------------------------------------------
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< array1d< globalIndex > > & var,
        localIndex const subArrayIndex,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
localIndex Unpack( buffer_unit_type const * & buffer,
                   ArrayOfSets< T > & var,
                   MPI_Op op );

//------------------------------------------------------------------------------
template< typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        MAP_TYPE & map,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T_FIRST, typename T_SECOND >
localIndex
Unpack( buffer_unit_type const * & buffer,
        std::pair< T_FIRST, T_SECOND > & var,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const * & buffer,
        InterObjectRelation< T > & var,
        MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< !is_host_packable_v< T >, localIndex >::type
Unpack( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
        T & GEOS_UNUSED_PARAM( var ),
        MPI_Op GEOS_UNUSED_PARAM( op ) )
{
  GEOS_ERROR( "Trying to unpack data type ("<<typeid(T).name()<<") but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
// UnpackArray(buffer,var,expectedLength)
//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
UnpackPointer( buffer_unit_type const * & buffer,
               T * const GEOS_RESTRICT var,
               INDEX_TYPE const expectedLength,
               MPI_Op op );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
UnpackPointer( buffer_unit_type const * & buffer,
               T * const GEOS_RESTRICT var,
               INDEX_TYPE const expectedLength,
               MPI_Op op );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
UnpackArray( buffer_unit_type const * & buffer,
             arraySlice1d< T, USD > const & var,
             INDEX_TYPE const length,
             MPI_Op op );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
UnpackArray( buffer_unit_type const * & buffer,
             arraySlice1d< T, USD > const & var,
             INDEX_TYPE const length,
             MPI_Op op );

//------------------------------------------------------------------------------
// UnpackByIndex(buffer,var,indices)
//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_indices >
localIndex
UnpackByIndex( buffer_unit_type const * & buffer,
               ArrayView< T, NDIM, USD > const & var,
               const T_indices & indices,
               MPI_Op op );

template< typename T, typename T_indices >
localIndex
UnpackByIndex( buffer_unit_type const * & buffer,
               ArrayOfArrays< T > & var,
               T_indices const & indices,
               MPI_Op op );

template< typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
UnpackByIndex( buffer_unit_type const * & buffer,
               MAP_TYPE & map,
               T_INDICES const & indices,
               MPI_Op op );

template< typename T, typename T_INDICES >
typename std::enable_if< !is_host_packable_by_index_v< T > && !is_host_packable_map_by_index_v< T >, localIndex >::type
UnpackByIndex( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
               T & GEOS_UNUSED_PARAM( var ),
               T_INDICES const & GEOS_UNUSED_PARAM( indices ),
               MPI_Op GEOS_UNUSED_PARAM( op ) )
{
  GEOS_ERROR( "Trying to unpack data type ("<<typeid(T).name()<<") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex Pack( buffer_unit_type * & buffer,
                 T const * const GEOS_RESTRICT var,
                 arraySlice1d< INDEX_TYPE const > const & indices,
                 INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
localIndex Unpack( buffer_unit_type const * & buffer,
                   T * const GEOS_RESTRICT var,
                   arraySlice1d< INDEX_TYPE const > const & indices,
                   INDEX_TYPE & length,
                   MPI_Op op );

//------------------------------------------------------------------------------
template< bool DO_PACKING, int USD >
localIndex Pack( buffer_unit_type * & buffer,
                 SortedArray< localIndex > const & var,
                 SortedArray< globalIndex > const & unmappedGlobalIndices,
                 arraySlice1d< globalIndex const, USD > const & localToGlobal );

//------------------------------------------------------------------------------
template< typename SORTED >
inline localIndex Unpack( buffer_unit_type const * & buffer,
                          SortedArray< localIndex > & var,
                          SortedArray< globalIndex > & unmappedGlobalIndices,
                          mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap,
                          bool const clearExistingSet,
                          MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfSets< T > const & var );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfSets< T > & var,
        MPI_Op op );


//------------------------------------------------------------------------------
template< bool DO_PACKING, int USD >
localIndex
Pack( buffer_unit_type * & buffer,
      arraySlice1d< localIndex const, USD > const & var,
      globalIndex const * const unmappedGlobalIndices,
      localIndex const length,
      arraySlice1d< globalIndex const > const & localToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED >
inline localIndex Unpack( buffer_unit_type const * & buffer,
                          localIndex_array & var,
                          array1d< globalIndex > & unmappedGlobalIndices,
                          mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap );

//------------------------------------------------------------------------------
inline localIndex UnpackSyncList( buffer_unit_type const * & buffer,
                                  localIndex_array & var,
                                  std::unordered_map< globalIndex, localIndex > const & globalToLocalMap );

//------------------------------------------------------------------------------
template< typename SORTED, int USD >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< localIndex, USD > & var,
        array1d< globalIndex > & unmappedGlobalIndices,
        localIndex const expectedLength,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< localIndex const > const & var,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArraysView< array1d< globalIndex > const > const & var,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< localIndex > const & var,
        array1d< localIndex > const & indices,
        mapBase< globalIndex, localIndex, SORTED0 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED1 > const & relatedObjectGlobalToLocalMap );


//------------------------------------------------------------------------------
template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< arrayView1d< localIndex const > const > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< localIndex_array > & var,
        array1d< localIndex > & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap );

//------------------------------------------------------------------------------
template< typename SORTED0 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< array1d< globalIndex > > & var,
        array1d< localIndex > & indices,
        mapBase< globalIndex, localIndex, SORTED0 > const & globalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< SortedArray< localIndex > const > const & var,
      mapBase< localIndex, SortedArray< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< SortedArray< localIndex > > & var,
        localIndex_array & indices,
        mapBase< localIndex, SortedArray< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap,
        bool const clearFlag );

//------------------------------------------------------------------------------
template< bool DO_PACKING, int USD0, int USD1 >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView2d< localIndex const, USD0 > const & var,
      arrayView1d< localIndex > const & indices,
      arraySlice1d< globalIndex const, USD1 > const & localToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED, int USD >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView2d< localIndex, USD > const & var,
        array1d< localIndex > & indices,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename SORTED, int USD0 >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView2d< localIndex const, USD0 > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arraySlice1d< globalIndex const > const & localToGlobalMap,
      arraySlice1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1, typename SORTED2, int USD >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView2d< localIndex, USD > const & var,
        localIndex_array & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type * & buffer, MAP_TYPE const & var );

//------------------------------------------------------------------------------
template< typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const * & buffer, MAP_TYPE & map, MPI_Op op );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type * & buffer, MAP_TYPE const & var, T_INDICES const & packIndices );

//------------------------------------------------------------------------------
template< typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const * & buffer, MAP_TYPE & map, T_INDICES const & unpackIndices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( buffer_unit_type * & buffer, std::pair< T_FIRST, T_SECOND > const & var );

//------------------------------------------------------------------------------
template< typename T_FIRST, typename T_SECOND >
localIndex
Unpack( buffer_unit_type const * & buffer, std::pair< T_FIRST, T_SECOND > & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type * & buffer,
      InterObjectRelation< T > const & var );

template< typename ... VARPACK >
localIndex
PackSize( VARPACK const && ... pack )
{
  buffer_unit_type * junk = nullptr;
  return Pack< false >( junk, pack ... );
}

template< typename ... VARPACK >
localIndex
PackSize( VARPACK && ... pack )
{
  buffer_unit_type * junk = nullptr;
  return Pack< false >( junk, pack ... );
}

#ifdef GEOS_USE_ARRAY_BOUNDS_CHECK
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
typename std::enable_if< !is_host_packable_by_index_v< T > &&
                         !is_host_packable_map_by_index_v< T >, localIndex >::type
Pack( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ), T const & GEOS_UNUSED_PARAM( var ), T_INDICES const & GEOS_UNUSED_PARAM( indices ) )
{
  GEOS_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, typename T_INDICES >
typename std::enable_if< !is_host_packable_by_index_v< T > &&
                         !is_host_packable_map_by_index_v< T >, localIndex >::type
Unpack( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ), T & GEOS_UNUSED_PARAM( var ), T_INDICES const & GEOS_UNUSED_PARAM( indices ) )
{
  GEOS_ERROR( "Trying to unpack data type ("<<typeid(T).name()<<") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( buffer_unit_type * & buffer,
      arraySlice1d< T > const & var,
      arraySlice1d< INDEX_TYPE > const & indices,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
localIndex
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< T > & var,
        arraySlice1d< INDEX_TYPE > const & indices,
        INDEX_TYPE & length );

#endif /* GEOS_USE_ARRAY_BOUNDS_CHECK */

} /* namespace bufferOps */
} /* namespace geos */

#include "BufferOps_inline.hpp"

#endif /* GEOS_DATAREPOSITORY_BUFFEROPS_HPP_ */
