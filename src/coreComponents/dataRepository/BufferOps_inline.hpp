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

#ifndef GEOS_DATAREPOSITORY_BUFFEROPS_INLINE_HPP_
#define GEOS_DATAREPOSITORY_BUFFEROPS_INLINE_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/traits.hpp"
#include "LvArray/src/limits.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <type_traits>

namespace geos
{

namespace bufferOps
{

template< bool DO_PACKING, typename T >
typename std::enable_if< is_host_packable_scalar_v< T >, localIndex >::type
PackData( buffer_unit_type * & buffer, T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  if( DO_PACKING )
  {
    memcpy( buffer, &var, sizeOfPackedChars );
    buffer += sizeOfPackedChars;
  }
  return sizeOfPackedChars;
}


template< bool DO_PACKING, typename T >
typename std::enable_if< is_host_packable_scalar_v< T >, localIndex >::type
Pack( buffer_unit_type * & buffer, T const & var )
{
  return PackData< DO_PACKING >( buffer, var );
}

template< bool DO_PACKING >
localIndex PackData( buffer_unit_type * & buffer, const string & var )
{
  localIndex sizeOfPackedChars = 0;
  const string::size_type varSize = var.size();
  if( DO_PACKING )
  {
    memcpy( buffer, var.data(), varSize );
    buffer += varSize;
  }
  sizeOfPackedChars += varSize;
  return sizeOfPackedChars;
}

template< bool DO_PACKING >
localIndex Pack( buffer_unit_type * & buffer, const string & var )
{
  const string::size_type varSize = var.size();
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, varSize );
  return sizeOfPackedChars + PackData< DO_PACKING >( buffer, var );
}

template< bool DO_PACKING, typename T >
localIndex PackData( buffer_unit_type * & buffer, SortedArray< T > const & var )
{
  localIndex sizeOfPackedChars = 0;
  for( T const & val : var )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, val );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
localIndex Pack( buffer_unit_type * & buffer, SortedArray< T > const & var )
{
  const localIndex length = LvArray::integerConversion< localIndex >( var.size() );
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  return sizeOfPackedChars + PackData< DO_PACKING >( buffer, var );
}

template< bool DO_PACKING, typename T, int SIZE >
localIndex PackData( buffer_unit_type * & buffer, Tensor< T, SIZE > const & var )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += PackPointer< DO_PACKING >( buffer, var.data, SIZE );
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, int SIZE >
localIndex Pack( buffer_unit_type * & buffer, Tensor< T, SIZE > const & var )
{
  return PackData< DO_PACKING >( buffer, var );
}

template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
PackData( buffer_unit_type * & buffer,
          ArrayView< T, NDIM, USD > const & var )
{
  localIndex sizeOfPackedChars = 0;
  const localIndex length = var.size();
  T const * const data = var.data();
  sizeOfPackedChars += PackPointerData< DO_PACKING >( buffer, data, length );
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
Pack( buffer_unit_type * & buffer,
      ArrayView< T, NDIM, USD > const & var )
{
  localIndex sizeOfPackedChars = PackPointer< DO_PACKING >( buffer, var.dims(), NDIM );
  sizeOfPackedChars += PackPointer< DO_PACKING >( buffer, var.strides(), NDIM );
  const localIndex length = var.size();
  T const * const data = var.data();
  sizeOfPackedChars += PackPointer< DO_PACKING >( buffer, data, length );
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
localIndex PackData( buffer_unit_type * & buffer,
                     ArrayOfArrays< T > const & var )
{
  localIndex sizeOfPackedChars = 0;
  for( localIndex a=0; a<var.size(); ++a )
  {
    T const * const data = var[a];
    sizeOfPackedChars += PackPointerData< DO_PACKING >( buffer, data, var.sizeOfArray( a ) );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
localIndex Pack( buffer_unit_type * & buffer,
                 ArrayOfArrays< T > const & var )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.size() );
  for( localIndex a=0; a<var.size(); ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.sizeOfArray( a ) );
    T const * const data = var[a];
    sizeOfPackedChars += PackPointer< DO_PACKING >( buffer, data, var.sizeOfArray( a ) );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
localIndex PackData( buffer_unit_type * & buffer,
                     ArrayOfSets< T > const & var )
{
  localIndex sizeOfPackedChars = 0;
  for( localIndex a=0; a<var.size(); ++a )
  {
    T const * const data = var[a];
    sizeOfPackedChars += PackPointerData< DO_PACKING >( buffer, data, var.sizeOfSet( a ) );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
localIndex Pack( buffer_unit_type * & buffer,
                 ArrayOfSets< T > const & var )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.size() );
  for( localIndex a=0; a<var.size(); ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.sizeOfSet( a ) );
    T const * const data = var[a];
    sizeOfPackedChars += PackPointer< DO_PACKING >( buffer, data, var.sizeOfSet( a ) );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
PackData( buffer_unit_type * & buffer, MAP_TYPE const & var )
{
  localIndex sizeOfPackedChars = 0;
  for( typename MAP_TYPE::const_iterator i = var.begin(); i != var.end(); ++i )
  {
    sizeOfPackedChars += PackData< DO_PACKING >( buffer, i->first );
    sizeOfPackedChars += PackData< DO_PACKING >( buffer, i->second );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type * & buffer, MAP_TYPE const & var )
{
  const typename MAP_TYPE::size_type length = var.size();
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  for( typename MAP_TYPE::const_iterator i = var.begin(); i != var.end(); ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->first );
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->second );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
PackData( buffer_unit_type * & buffer, std::pair< T_FIRST, T_SECOND > const & var )
{
  localIndex sizeOfPackedChars = PackData< DO_PACKING >( buffer, var.first );
  sizeOfPackedChars += PackData< DO_PACKING >( buffer, var.second );
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( buffer_unit_type * & buffer, std::pair< T_FIRST, T_SECOND > const & var )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, var.first );
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.second );
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
localIndex PackData( buffer_unit_type * & buffer, InterObjectRelation< T > const & var )
{
  return PackData< DO_PACKING >( buffer, static_cast< T const & >(var));
}

template< bool DO_PACKING, typename T >
localIndex Pack( buffer_unit_type * & buffer, InterObjectRelation< T > const & var )
{
  return Pack< DO_PACKING >( buffer, static_cast< T const & >(var));
}

//------------------------------------------------------------------------------
// PackArray(buffer,var,length)
//------------------------------------------------------------------------------

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackPointerData( buffer_unit_type * & buffer, T const * const GEOS_RESTRICT var, INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = length * sizeof(T);
  if( DO_PACKING )
  {
    memcpy( buffer, var, length * sizeof(T) );
    buffer += length * sizeof(T);
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackPointer( buffer_unit_type * & buffer, T const * const GEOS_RESTRICT var, INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  return sizeOfPackedChars + PackPointerData< DO_PACKING >( buffer, var, length );
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackPointerData( buffer_unit_type * & buffer,
                 T const * const GEOS_RESTRICT var,
                 INDEX_TYPE const length )

{
  localIndex sizeOfPackedChars = 0;
  for( INDEX_TYPE a = 0; a < length; ++a )
  {
    sizeOfPackedChars += PackData< DO_PACKING >( buffer, var[ a ] );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackPointer( buffer_unit_type * & buffer,
             T const * const GEOS_RESTRICT var,
             INDEX_TYPE const length )

{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  for( INDEX_TYPE a = 0; a < length; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ a ] );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackArray( buffer_unit_type * & buffer,
           arraySlice1d< T, USD > const & var,
           INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  sizeOfPackedChars += length * sizeof(T);
  if( DO_PACKING )
  {
    T * const GEOS_RESTRICT buffer_T = reinterpret_cast< T * >( buffer );
    for( INDEX_TYPE i = 0; i < length; ++i )
    {
      buffer_T[ i ] = var[ i ];
    }
    buffer += length * sizeof(T);
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackArray( buffer_unit_type * & buffer,
           arraySlice1d< T, USD > const & var,
           INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  for( INDEX_TYPE a = 0; a < length; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ a ] );
  }
  return sizeOfPackedChars;
}

//------------------------------------------------------------------------------
// PackByIndex(buffer,var,indices)
//------------------------------------------------------------------------------

template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_indices >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
PackDataByIndex( buffer_unit_type * & buffer,
                 ArrayView< T, NDIM, USD > const & var,
                 const T_indices & indices )
{
  localIndex sizeOfPackedChars = 0;
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    LvArray::forValuesInSlice( var[ indices[ a ] ],
                               [&sizeOfPackedChars, &buffer]( T const & value )
    {
      sizeOfPackedChars += PackData< DO_PACKING >( buffer, value );
    } );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_indices >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
PackByIndex( buffer_unit_type * & buffer,
             ArrayView< T, NDIM, USD > const & var,
             const T_indices & indices )
{
  localIndex sizeOfPackedChars = PackPointer< DO_PACKING >( buffer, var.strides(), NDIM );
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    LvArray::forValuesInSlice( var[ indices[ a ] ],
                               [&sizeOfPackedChars, &buffer]( T const & value )
    {
      sizeOfPackedChars += Pack< DO_PACKING >( buffer, value );
    } );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename T_indices >
localIndex PackDataByIndex( buffer_unit_type * & buffer,
                            ArrayOfArrays< T > const & var,
                            T_indices const & indices )
{
  localIndex sizeOfPackedChars = 0;
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    T const * const data = var[indices[a]];
    sizeOfPackedChars += PackPointerData< DO_PACKING >( buffer, data, var.sizeOfArray( indices[a] ) );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename T_indices >
localIndex PackByIndex( buffer_unit_type * & buffer,
                        ArrayOfArrays< T > const & var,
                        T_indices const & indices )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.sizeOfArray( indices[a] ) );
    T const * const data = var[indices[a]];
    sizeOfPackedChars += PackPointer< DO_PACKING >( buffer, data, var.sizeOfArray( indices[a] ) );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
PackDataByIndex( buffer_unit_type * & buffer,
                 MAP_TYPE const & var,
                 T_INDICES const & indices )
{
  localIndex sizeOfPackedChars = 0;
  for( typename MAP_TYPE::const_iterator i = var.begin(); i != var.end(); ++i )
  {
    sizeOfPackedChars += PackData< DO_PACKING >( buffer, i->first );
    sizeOfPackedChars += PackDataByIndex< DO_PACKING >( buffer, i->second, indices );
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
PackByIndex( buffer_unit_type * & buffer,
             MAP_TYPE const & var,
             T_INDICES const & indices )
{
  const typename MAP_TYPE::size_type length = var.size();
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  for( typename MAP_TYPE::const_iterator i = var.begin(); i != var.end(); ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->first );
    sizeOfPackedChars += PackByIndex< DO_PACKING >( buffer, i->second, indices );
  }
  return sizeOfPackedChars;
}

//------------------------------------------------------------------------------
// Unpack(buffer,var)
//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< is_host_packable_scalar_v< T >, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        T & var,
        MPI_Op op )
{
  localIndex const sizeOfUnpackedChars = sizeof(T);
  T const castValue = *reinterpret_cast< T const * >( buffer );
  if( op == MPI_REPLACE )
  {
    var = castValue;
  }
  else if( op == MPI_SUM )
  {
    if constexpr ( traits::has_plus_equal_v< T > )
    {
      var += castValue;
    }
    else
    {
      GEOS_ERROR( GEOS_FMT( "Unsupported unpack operator+= for type {}!", typeid(T).name() ) );
    }
  }
  else if( op == MPI_MAX )
  {
    if constexpr ( traits::has_less_than_v< T > ) // && traits::has_equality_v< T >
    {
      var = std::max( var, castValue );
    }
    else
    {
      GEOS_ERROR( GEOS_FMT( "Unsupported unpack operator< for type {}!", typeid(T).name() ) );
    }
  }
  else
  {
    GEOS_ERROR( "Unsupported MPI operator! MPI_SUM, MPI_REPLACE and MPI_MAX are supported." );
  }
  memcpy( &var, buffer, sizeOfUnpackedChars );
  buffer += sizeOfUnpackedChars;
  return sizeOfUnpackedChars;
}

inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        string & var,
        MPI_Op op )
{
  string::size_type stringsize = 0;
  localIndex sizeOfUnpackedChars = Unpack( buffer, stringsize, MPI_REPLACE );
  var.resize( stringsize );
  string castValue( stringsize, ' ' );
  memcpy( &castValue[0], buffer, stringsize );
  if( op == MPI_REPLACE )
  {
    var = castValue;
  }
  else if( op == MPI_SUM )
  {
    var += castValue;
  }
  else
  {
    GEOS_ERROR( "Unsupported MPI operator! MPI_SUM and MPI_REPLACE are supported." );
  }
  buffer += stringsize;
  sizeOfUnpackedChars += stringsize;
  return sizeOfUnpackedChars;
}

template< typename T, int SIZE >
localIndex
Unpack( buffer_unit_type const * & buffer,
        Tensor< T, SIZE > & var,
        MPI_Op op )
{
  localIndex sizeOfUnpackedChars = 0;
  sizeOfUnpackedChars += UnpackPointer( buffer, var.data, SIZE, op );
  return sizeOfUnpackedChars;
}

template< typename T >
localIndex
Unpack( buffer_unit_type const * & buffer,
        SortedArray< T > & var,
        MPI_Op op )
{
  var.clear();
  localIndex set_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, set_length, MPI_REPLACE );
  for( localIndex a=0; a<set_length; ++a )
  {
    T temp;
    sizeOfUnpackedChars += Unpack( buffer, temp, op );
    var.insert( temp );
  }
  return sizeOfUnpackedChars;
}

template< typename T, int NDIM, typename PERMUTATION >
typename std::enable_if< is_host_packable_v< T >, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        Array< T, NDIM, PERMUTATION > & var,
        MPI_Op op )
{
  localIndex dims[NDIM];
  localIndex sizeOfUnpackedChars = UnpackPointer( buffer, dims, NDIM, MPI_REPLACE );
  var.resize( NDIM, dims );

  localIndex strides[NDIM];
  sizeOfUnpackedChars += UnpackPointer( buffer, strides, NDIM, MPI_REPLACE );
  for( int i=0; i<NDIM; ++i )
  {
    GEOS_ASSERT_EQ( strides[i], var.strides()[i] );
  }

  sizeOfUnpackedChars += UnpackPointer( buffer, var.data(), var.size(), op );
  return sizeOfUnpackedChars;
}

template< typename T >
localIndex Unpack( buffer_unit_type const * & buffer,
                   ArrayOfArrays< T > & var,
                   MPI_Op op )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex numOfArrays;
  sizeOfUnpackedChars += Unpack( buffer, numOfArrays, op );
  var.resize( numOfArrays );
  for( localIndex a=0; a<numOfArrays; ++a )
  {
    localIndex sizeOfArray;
    sizeOfUnpackedChars += Unpack( buffer, sizeOfArray, op );
    var.resizeArray( a, sizeOfArray );
    T * data = var[a];
    sizeOfUnpackedChars += UnpackPointer( buffer, data, sizeOfArray, op );
  }
  return sizeOfUnpackedChars;
}

inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< array1d< globalIndex > > & var,
        localIndex const subArrayIndex,
        MPI_Op op )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = bufferOps::Unpack( buffer, length, op );

  var.resizeArray( subArrayIndex, length );

  for( localIndex a = 0; a < length; ++a )
  {
    array1d< globalIndex > & tmp = var( subArrayIndex, a );
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, tmp, op );
  }

  return sizeOfUnpackedChars;
}


template< typename T >
localIndex Unpack( buffer_unit_type const * & buffer,
                   ArrayOfSets< T > & var,
                   MPI_Op op )
{
  ArrayOfArrays< T > varAsArray;
  localIndex sizeOfUnpackedChars = Unpack( buffer, varAsArray, op );
  var.template assimilate< parallelHostPolicy >( std::move( varAsArray ), LvArray::sortedArrayManipulation::SORTED_UNIQUE );
  return sizeOfUnpackedChars;
}

template< typename MAP_TYPE >
typename std::enable_if< is_host_packable_map_v< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        MAP_TYPE & map,
        MPI_Op op )
{
  map.clear();
  typename MAP_TYPE::size_type map_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, map_length, MPI_REPLACE );
  for( typename MAP_TYPE::size_type a = 0; a < map_length; ++a )
  {
    typename MAP_TYPE::key_type key;
    typename MAP_TYPE::mapped_type value;
    sizeOfUnpackedChars += Unpack( buffer, key, MPI_REPLACE );
    sizeOfUnpackedChars += Unpack( buffer, value, op );
    map[key] = std::move( value );
  }
  return sizeOfUnpackedChars;
}

template< typename T_FIRST, typename T_SECOND >
localIndex Unpack( buffer_unit_type const * & buffer,
                   std::pair< T_FIRST, T_SECOND > & var,
                   MPI_Op op )
{
  localIndex sizeOfUnpackedChars = Unpack( buffer, var.first, op );
  sizeOfUnpackedChars += Unpack( buffer, var.second, op );
  return sizeOfUnpackedChars;
}

template< typename T >
localIndex Unpack( buffer_unit_type const * & buffer,
                   InterObjectRelation< T > & var,
                   MPI_Op op )
{
  return Unpack( buffer, static_cast< T & >(var), op );
}

//------------------------------------------------------------------------------
// UnpackArray(buffer,var,expectedLength)
//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
UnpackPointer( buffer_unit_type const * & buffer,
               T * const GEOS_RESTRICT var,
               INDEX_TYPE const expectedLength,
               MPI_Op op )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  GEOS_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                   expectedLength << " != " << length );
  GEOS_DEBUG_VAR( expectedLength );
  T const * castBuffer = reinterpret_cast< T const * >( buffer );
  if( op == MPI_REPLACE )
  {
    for( int ii = 0; ii < length; ++ii )
    {
      var[ ii ] = castBuffer[ ii ];
    }
  }
  else if( op == MPI_SUM )
  {
    for( int ii = 0; ii < length; ++ii )
    {
      var[ ii ] += castBuffer[ ii ];
    }
  }
  else if( op == MPI_MAX )
  {
    for( int ii = 0; ii < length; ++ii )
    {
      var[ ii ] = std::max( var[ii], castBuffer[ ii ] );
    }
  }
  else
  {
    GEOS_ERROR( "Unsupported MPI operator! MPI_SUM, MPI_REPLACE and MPI_MAX are supported." );
  }
  sizeOfUnpackedChars += length * sizeof(T);
  buffer += length * sizeof(T);
  return sizeOfUnpackedChars;
}

template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
UnpackPointer( buffer_unit_type const * & buffer,
               T * const GEOS_RESTRICT var,
               INDEX_TYPE const expectedLength,
               MPI_Op op )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  GEOS_ASSERT_EQ( length, expectedLength );
  GEOS_DEBUG_VAR( expectedLength );
  for( INDEX_TYPE a=0; a<length; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a], op );
  }

  return sizeOfUnpackedChars;
}

template< typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
UnpackArray( buffer_unit_type const * & buffer,
             arraySlice1d< T, USD > const & var,
             INDEX_TYPE const expectedLength,
             MPI_Op op )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  GEOS_DEBUG_VAR( expectedLength );
  GEOS_ASSERT_EQ( length, expectedLength );

  T const * castBuffer = reinterpret_cast< T const * >( buffer );
  if( op == MPI_REPLACE )
  {
    for( int ii = 0; ii < length; ++ii )
    {
      var[ ii ] = castBuffer[ ii ];
    }
  }
  else if( op == MPI_SUM )
  {
    for( int ii = 0; ii < length; ++ii )
    {
      var[ ii ] += castBuffer[ ii ];
    }
  }
  else if( op == MPI_MAX )
  {
    for( int ii = 0; ii < length; ++ii )
    {
      var[ ii ] = std::max( var[ii], castBuffer[ ii ] );
    }
  }
  else
  {
    GEOS_ERROR( "Unsupported MPI operator! MPI_SUM, MPI_REPLACE and MPI_MAX are supported." );
  }

  buffer += length * sizeof(T);
  sizeOfUnpackedChars += length * sizeof(T);
  return sizeOfUnpackedChars;
}

template< typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
UnpackArray( buffer_unit_type const * & buffer,
             arraySlice1d< T, USD > const & var,
             INDEX_TYPE const expectedLength,
             MPI_Op op )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  GEOS_DEBUG_VAR( expectedLength );
  GEOS_ASSERT_EQ( length, expectedLength );

  for( INDEX_TYPE a=0; a<length; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a], op );
  }

  return sizeOfUnpackedChars;
}

//------------------------------------------------------------------------------
// UnpackByIndex(buffer,var,indices)
//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_indices >
localIndex
UnpackByIndex( buffer_unit_type const * & buffer,
               ArrayView< T, NDIM, USD > const & var,
               const T_indices & indices,
               MPI_Op op )
{
  localIndex strides[NDIM];
  localIndex sizeOfUnpackedChars = UnpackPointer( buffer, strides, NDIM, MPI_REPLACE );

  for( localIndex a=0; a<indices.size(); ++a )
  {
    LvArray::forValuesInSlice( var[ indices[ a ] ],
                               [&sizeOfUnpackedChars, &buffer, &op] ( T & value )
    {
      sizeOfUnpackedChars += Unpack( buffer, value, op );
    } );
  }
  return sizeOfUnpackedChars;
}

template< typename T, typename T_indices >
localIndex
UnpackByIndex( buffer_unit_type const * & buffer,
               ArrayOfArrays< T > & var,
               T_indices const & indices,
               MPI_Op op )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex numUnpackedIndices = 0;
  sizeOfUnpackedChars += Unpack( buffer, numUnpackedIndices, MPI_REPLACE );
  GEOS_ERROR_IF( numUnpackedIndices != indices.size(), "number of unpacked indices does not equal expected number" );
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    localIndex sizeOfSubArray;
    sizeOfUnpackedChars += Unpack( buffer, sizeOfSubArray, MPI_REPLACE );
    var.resizeArray( indices[a], sizeOfSubArray );
    sizeOfUnpackedChars += UnpackArray( buffer, var[indices[a]], sizeOfSubArray, op );
  }
  return sizeOfUnpackedChars;
}

template< typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
UnpackByIndex( buffer_unit_type const * & buffer,
               MAP_TYPE & map,
               T_INDICES const & unpackIndices,
               MPI_Op op )
{
  map.clear();
  typename MAP_TYPE::size_type map_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, map_length, MPI_REPLACE );
  for( typename MAP_TYPE::size_type a = 0; a < map_length; ++a )
  {
    typename MAP_TYPE::key_type key;
    typename MAP_TYPE::mapped_type value;
    sizeOfUnpackedChars += Unpack( buffer, key, MPI_REPLACE );
    sizeOfUnpackedChars += UnpackByIndex( buffer, value, unpackIndices, op );
    map[key] = std::move( value );
  }
  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex Pack( buffer_unit_type * & buffer,
                 T const * const GEOS_RESTRICT var,
                 arraySlice1d< INDEX_TYPE const > const & indices,
                 INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( localIndex a=0; a<length; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}

template< typename T, typename INDEX_TYPE >
localIndex Unpack( buffer_unit_type const * & buffer,
                   T * const GEOS_RESTRICT var,
                   arraySlice1d< INDEX_TYPE const > const & indices,
                   INDEX_TYPE & length,
                   MPI_Op op )
{
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );

  for( INDEX_TYPE a=0; a<length; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[ indices[a] ], op );
  }

  return sizeOfUnpackedChars;
}

#ifdef GEOS_USE_ARRAY_BOUNDS_CHECK

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
Pack( buffer_unit_type * & buffer,
      arraySlice1d< T > const & var,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( INDEX_TYPE a=0; a<length; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ a ] );
  }

  return sizeOfPackedChars;
}

template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< T > & var,
        INDEX_TYPE const expectedLength,
        MPI_Op op )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  GEOS_DEBUG_VAR( expectedLength );
  GEOS_ASSERT_EQ( length, expectedLength );

  T const * const GEOS_RESTRICT castBuffer = reinterpret_cast< T const * >( buffer );
  if( op == MPI_REPLACE )
  {
    for( INDEX_TYPE ii = 0; ii < length; ++ii )
    {
      var[ ii ] = castBuffer[ ii ];
    }
  }
  else if( op == MPI_SUM )
  {
    for( INDEX_TYPE ii = 0; ii < length; ++ii )
    {
      var[ ii ] += castBuffer[ ii ];
    }
  }
  else if( op == MPI_MAX )
  {
    for( INDEX_TYPE ii = 0; ii < length; ++ii )
    {
      var[ ii ] = std::max( var[ ii ], castBuffer[ ii ] );
    }
  }
  else
  {
    GEOS_ERROR( "Unsupported MPI operator! MPI_SUM, MPI_REPLACE and MPI_MAX are supported." );
  }

  buffer += length * sizeof(T);
  sizeOfUnpackedChars += length * sizeof(T);

  return sizeOfUnpackedChars;
}

template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< T > & var,
        INDEX_TYPE const expectedLength,
        MPI_Op op )
{
  INDEX_TYPE length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  GEOS_DEBUG_VAR( expectedLength );
  GEOS_ASSERT_EQ( length, expectedLength );

  for( INDEX_TYPE a=0; a<length; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[a], op );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( buffer_unit_type * & buffer,
      arraySlice1d< T > const & var,
      arraySlice1d< INDEX_TYPE > const & indices,
      INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( INDEX_TYPE a=0; a<length; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}


template< typename T, typename INDEX_TYPE >

localIndex
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< T > & var,
        arraySlice1d< INDEX_TYPE > const & indices,
        INDEX_TYPE & length,
        MPI_Op op )
{
  localIndex sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += Unpack( buffer, length, MPI_REPLACE );

  for( INDEX_TYPE a=0; a<length; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[ indices[a] ], op );
  }

  return sizeOfUnpackedChars;
}

#endif /* GEOS_USE_ARRAY_BOUNDS_CHECK */

template< bool DO_PACKING, int USD >
localIndex Pack( buffer_unit_type * & buffer,
                 SortedArray< localIndex > const & var,
                 SortedArray< globalIndex > const & unmappedGlobalIndices,
                 arraySlice1d< globalIndex const, USD > const & localToGlobal )
{
  const localIndex length = LvArray::integerConversion< localIndex >( var.size()+unmappedGlobalIndices.size());
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( localIndex const lid : var )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobal[ lid ] );
  }

  for( globalIndex const gid : unmappedGlobalIndices )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, gid );
  }


  return sizeOfPackedChars;
}

template< typename SORTED >
inline
localIndex Unpack( buffer_unit_type const * & buffer,
                   SortedArray< localIndex > & var,
                   SortedArray< globalIndex > & unmappedGlobalIndices,
                   mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap,
                   bool const clearExistingSet )
{
  if( clearExistingSet )
  {
    var.clear();
  }
  localIndex set_length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, set_length, MPI_REPLACE );

  for( localIndex a=0; a<set_length; ++a )
  {
    globalIndex temp;
    sizeOfUnpackedChars += Unpack( buffer, temp, MPI_REPLACE );
    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator iter = globalToLocalMap.find( temp );
    if( iter==globalToLocalMap.end() )
    {
      unmappedGlobalIndices.insert( temp );
    }
    else
    {
      var.insert( iter->second );
    }
  }

  return sizeOfUnpackedChars;
}

template< bool DO_PACKING >
localIndex Pack( buffer_unit_type * & buffer,
                 SortedArrayView< localIndex const > const & var,
                 arrayView1d< localIndex const > const & packList,
                 arraySlice1d< globalIndex const > const & localToGlobal )
{
  localIndex length = 0;
  for( auto a : packList )
  {
    length += var.count( a );
  }

  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( localIndex a=0; a< packList.size(); ++a )
  {
    if( var.count( packList[ a ] ) )
    {
      sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobal[packList[a]] );
    }
  }

  return sizeOfPackedChars;
}

template< bool DO_PACKING >
localIndex Pack( buffer_unit_type * & buffer,
                 SortedArrayView< localIndex const > const & var,
                 arrayView1d< localIndex const > const & packList,
                 SortedArrayView< globalIndex const > const & unmappedGlobalIndices,
                 arraySlice1d< globalIndex const > const & localToGlobal )
{

  localIndex length = 0;
  array1d< localIndex > temp( var.size() );

  for( auto a : packList )
  {
    if( var.count( a ) )
    {
      temp[length] = a;
      ++length;
    }
  }
  temp.resize( length );
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );


  for( localIndex a=0; a<length; ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobal[temp[a]] );
  }

  for( globalIndex const gid : unmappedGlobalIndices )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, gid );
  }

  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename T_indices >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArrays< T > const & var,
      T_indices const & indices )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0; a<indices.size(); ++a )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, var.sizeOfArray( indices[a] ) );
    T const * const data = var[indices[a]];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, data, var.sizeOfArray( indices[a] ) );
  }
  return sizeOfPackedChars;
}


template< typename T, typename T_indices >
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< T > & var,
        T_indices const & indices )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex numUnpackedIndices;
  sizeOfUnpackedChars += Unpack( buffer, numUnpackedIndices, MPI_REPLACE );

  GEOS_ERROR_IF( numUnpackedIndices!=indices.size(), "number of unpacked indices does not equal expected number" );

  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex sizeOfSubArray;
    sizeOfUnpackedChars += Unpack( buffer, sizeOfSubArray, MPI_REPLACE );
    var.resizeArray( indices[a], sizeOfSubArray );
    sizeOfUnpackedChars += Unpack( buffer, var[indices[a]], sizeOfSubArray, MPI_REPLACE );
  }


  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, int USD >
localIndex Pack( buffer_unit_type * & buffer,
                 arraySlice1d< localIndex const, USD > const & var,
                 globalIndex const * const unmappedGlobalIndices,
                 localIndex const length,
                 arraySlice1d< globalIndex const > const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );
  sizeOfPackedChars += length*sizeof(globalIndex);

  if( DO_PACKING )
  {
    globalIndex * const buffer_GI = reinterpret_cast< globalIndex * >(buffer);
    for( localIndex a=0; a<length; ++a )
    {
      if( var[a] != unmappedLocalIndexValue )
      {
        buffer_GI[a] = localToGlobalMap[var[a]];
      }
      else
      {
        buffer_GI[a] = unmappedGlobalIndices[a];
      }
    }

    buffer += length * sizeof(globalIndex);
  }

  return sizeOfPackedChars;
}

template< typename SORTED >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        localIndex_array & var,
        array1d< globalIndex > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  var.resize( length );
  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices.setValues< serialPolicy >( unmappedLocalIndexValue );

  bool unpackedGlobalFlag = false;
  for( localIndex a=0; a<length; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex, MPI_REPLACE );

    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
      iter = globalToLocalMap.find( unpackedGlobalIndex );
    if( iter == globalToLocalMap.end() )
    {
      var[a] = unmappedLocalIndexValue;
      unmappedGlobalIndices[a] = unpackedGlobalIndex;
      unpackedGlobalFlag = true;
    }
    else
    {
      var[a] = iter->second;
    }
  }
  if( !unpackedGlobalFlag )
  {
    unmappedGlobalIndices.clear();
  }

  return sizeOfUnpackedChars;
}

inline
localIndex
UnpackSyncList( buffer_unit_type const * & buffer,
                localIndex_array & var,
                std::unordered_map< globalIndex, localIndex > const & globalToLocalMap )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );
  var.resize( length );

  for( localIndex a=0; a<length; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex, MPI_REPLACE );
    var[a] = globalToLocalMap.at( unpackedGlobalIndex );
  }

  return sizeOfUnpackedChars;
}

template< typename SORTED >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< localIndex > & var,
        localIndex const subArrayIndex,
        array1d< globalIndex > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = Unpack( buffer, length, MPI_REPLACE );

  var.resizeArray( subArrayIndex, length );
  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices.setValues< serialPolicy >( unmappedLocalIndexValue );

  bool unpackedGlobalFlag = false;
  for( localIndex a=0; a<length; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex, MPI_REPLACE );

    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
      iter = globalToLocalMap.find( unpackedGlobalIndex );
    if( iter == globalToLocalMap.end() )
    {
      var( subArrayIndex, a ) = unmappedLocalIndexValue;
      unmappedGlobalIndices[a] = unpackedGlobalIndex;
      unpackedGlobalFlag = true;
    }
    else
    {
      var( subArrayIndex, a ) = iter->second;
    }
  }
  if( !unpackedGlobalFlag )
  {
    unmappedGlobalIndices.clear();
  }

  return sizeOfUnpackedChars;
}

template< typename SORTED, int USD >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arraySlice1d< localIndex, USD > & var,
        array1d< globalIndex > & unmappedGlobalIndices,
        localIndex const expectedLength,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex length;
  sizeOfUnpackedChars += Unpack( buffer, length, MPI_REPLACE );

  GEOS_ASSERT_EQ( length, expectedLength );
  GEOS_DEBUG_VAR( expectedLength );

  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices.setValues< serialPolicy >( unmappedLocalIndexValue );

  bool unpackedGlobalFlag = false;
  for( localIndex a=0; a<length; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, unpackedGlobalIndex, MPI_REPLACE );

    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
      iter = globalToLocalMap.find( unpackedGlobalIndex );
    if( iter == globalToLocalMap.end() )
    {
      var[a] = unmappedLocalIndexValue;
      unmappedGlobalIndices[a] = unpackedGlobalIndex;
      unpackedGlobalFlag = true;
    }
    else
    {
      var[a] = iter->second;
    }
  }
  if( !unpackedGlobalFlag )
  {
    unmappedGlobalIndices.clear();
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< localIndex const > const & var,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );
    if( var[li] != -1 )
    {
      sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                               relatedObjectLocalToGlobalMap[var[li]] );
    }
    else
    {
      sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                               globalIndex( -1 ) );
    }
  }

  return sizeOfPackedChars;
}

template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArraysView< array1d< globalIndex > const > const & var,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    sizeOfPackedChars += bufferOps::PackArray< DO_PACKING >( buffer,
                                                             var[li],
                                                             var.sizeOfArray( li ) );
  }

  return sizeOfPackedChars;
}

template< typename SORTED0, typename SORTED1 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< localIndex > const & var,
        array1d< localIndex > const & indices,
        mapBase< globalIndex, localIndex, SORTED0 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED1 > const & relatedObjectGlobalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  for( localIndex a=0; a<indices.size(); ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi, MPI_REPLACE );
    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    globalIndex mappedGlobalIndex;
    sizeOfUnpackedChars += Unpack( buffer, mappedGlobalIndex, MPI_REPLACE );
    if( mappedGlobalIndex != -1 )
    {
      var[li] = relatedObjectGlobalToLocalMap.at( mappedGlobalIndex );
    }
    else
    {
      var[li] = -1;
    }

  }
  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< arrayView1d< localIndex const > const > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, array1d< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > junk;
    array1d< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                junk :
                                                iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li].toSliceConst(),
                                             unmappedGI.data(),
                                             var[li].size(),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< localIndex_array > & var,
        array1d< localIndex > & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0; a<indices.size(); ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi, MPI_REPLACE );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    array1d< globalIndex > unmappedIndices;
    sizeOfUnpackedChars += Unpack( buffer,
                                   var[li],
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap );

    if( unmappedIndices.size() > 0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }
  }
  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArraysView< localIndex const > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;
  array1d< globalIndex > junk;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, array1d< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                junk :
                                                iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.data(),
                                             var.sizeOfArray( li ),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      ArrayOfArraysView< localIndex const > const & var,
      mapBase< localIndex, SortedArray< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;
  SortedArray< globalIndex > junk;

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, SortedArray< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    SortedArray< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                    junk :
                                                    iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.data(),
                                             var.sizeOfArray( li ),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< localIndex > & var,
        array1d< localIndex > & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );
  array1d< globalIndex > unmappedIndices;

  for( localIndex a=0; a<indices.size(); ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi, MPI_REPLACE );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    unmappedIndices.resize( 0 );
    sizeOfUnpackedChars += Unpack( buffer,
                                   var,
                                   li,
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap );

    if( unmappedIndices.size() > 0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }
  }
  return sizeOfUnpackedChars;
}

template< typename SORTED0 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfArrays< array1d< globalIndex > > & var,
        array1d< localIndex > & indices,
        mapBase< globalIndex, localIndex, SORTED0 > const & globalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = bufferOps::Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn != 0 && numIndicesUnpacked != indices.size(),
                 "number of unpacked indices(" << numIndicesUnpacked << ") does not equal size of "
                                                                        "indices passed into Unpack function(" << sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );
  array1d< globalIndex > unmappedIndices;

  for( localIndex a=0; a<indices.size(); ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, gi, MPI_REPLACE );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li != globalToLocalMap.at( gi ),
                     "global index " << gi << " unpacked from buffer does not equal the lookup "
                                     << li << " for localIndex " << li << " on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    unmappedIndices.resize( 0 );
    sizeOfUnpackedChars += Unpack( buffer,
                                   var,
                                   li,
                                   MPI_REPLACE );
  }
  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename SORTED >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView1d< SortedArray< localIndex > const > const & var,
      mapBase< localIndex, SortedArray< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, indices.size() );

  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, SortedArray< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    SortedArray< globalIndex > junk;
    SortedArray< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                    junk :
                                                    iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI,
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView1d< SortedArray< localIndex > > & var,
        localIndex_array & indices,
        mapBase< localIndex, SortedArray< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap,
        bool const clearFlag )
{
  localIndex sizeOfUnpackedChars=0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0; a<numIndicesUnpacked; ++a )
  {

    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi, MPI_REPLACE );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    SortedArray< globalIndex > unmappedIndices;
    sizeOfUnpackedChars += Unpack( buffer,
                                   var[li],
                                   unmappedIndices,
                                   relatedObjectGlobalToLocalMap,
                                   clearFlag );

    if( unmappedIndices.size() > 0 )
    {
      unmappedGlobalIndices[li].insert( unmappedIndices.data(), unmappedIndices.size() );
    }
  }
  return sizeOfUnpackedChars;
}

template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        ArrayOfSets< localIndex > & var,
        localIndex_array & indices,
        mapBase< localIndex, SortedArray< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap,
        bool const clearFlag )
{
  localIndex sizeOfUnpackedChars=0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  // for objects related to the above local index li (e.g. up/down mappings)
  // global indices not yet known on the local rank
  std::vector< globalIndex > unmapped;

  // local indices of known global indices
  std::vector< localIndex > mapped;

  // local indices of known objects not yet present in the map
  std::vector< localIndex > mappedNew;

  // storage for new values that don't fit into existing capacity
  array1d< localIndex > indiciesToInsert;
  ArrayOfSets< localIndex > valuesToInsert;
  indiciesToInsert.reserve( numIndicesUnpacked );
  valuesToInsert.reserve( numIndicesUnpacked );
  valuesToInsert.reserveValues( numIndicesUnpacked * 12 ); // guesstimate

  for( localIndex a=0; a<numIndicesUnpacked; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi, MPI_REPLACE );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    if( clearFlag )
    {
      var.clearSet( li );
    }

    localIndex set_length;
    sizeOfUnpackedChars += Unpack( buffer, set_length, MPI_REPLACE );

    mapped.clear();
    mappedNew.clear();
    unmapped.clear();

    // again, the global indices being unpacked here are for
    //  objects related to the global index recvd and
    //  mapped to a local index above (e.g. up/down maps)
    for( localIndex b = 0; b < set_length; ++b )
    {
      globalIndex temp;
      sizeOfUnpackedChars += Unpack( buffer, temp, MPI_REPLACE );
      auto iter = relatedObjectGlobalToLocalMap.find( temp );
      // if we have no existing global-to-local information
      //  for the recv'd global index
      if( iter == relatedObjectGlobalToLocalMap.end() )
      {
        unmapped.push_back( temp );
      }
      // if we have existing global-to-local information
      //  use that mapping and store the local index
      else
      {
        mapped.push_back( iter->second );
      }
    }

    // insert known local indices into the set of indices related to the local index
    mapped.resize( LvArray::sortedArrayManipulation::makeSortedUnique( mapped.begin(), mapped.end() ) );
    std::set_difference( mapped.begin(), mapped.end(), var[li].begin(), var[li].end(), std::back_inserter( mappedNew ) );
    localIndex const numNewValues = LvArray::integerConversion< localIndex >( mappedNew.size() );

    // check if we have enough capacity to insert new indices
    if( numNewValues <= var.capacityOfSet( li ) - var.sizeOfSet( li ) )
    {
      // all good, just insert new entries
      var.insertIntoSet( li, mappedNew.begin(), mappedNew.end() );
    }
    else
    {
      // need to stash away and rebuild the map later
      localIndex const k = indiciesToInsert.size();
      indiciesToInsert.emplace_back( li );
      valuesToInsert.appendSet( numNewValues );
      valuesToInsert.insertIntoSet( k, mappedNew.begin(), mappedNew.end() );
    }

    // insert unknown global indices related to the local index into an additional mapping to resolve externally
    unmapped.resize( LvArray::sortedArrayManipulation::makeSortedUnique( unmapped.begin(), unmapped.end() ) );
    if( unmapped.size() > 0 )
    {
      unmappedGlobalIndices[li].insert( unmapped.begin(), unmapped.end() );
    }
  }

  // If there were element lists that didn't fit in the map, rebuild the whole thing
  if( !indiciesToInsert.empty() )
  {
    // Copy old capacities (no way to direct copy, must kernel launch)
    array1d< localIndex > newCapacities( var.size() );
    forAll< parallelHostPolicy >( var.size(), [var = var.toViewConst(), newCapacities = newCapacities.toView()]( localIndex const k )
    {
      newCapacities[k] = var.capacityOfSet( k );
    } );

    // Add new capacities where needed
    for( localIndex i = 0; i < indiciesToInsert.size(); ++i )
    {
      newCapacities[indiciesToInsert[i]] += valuesToInsert.sizeOfSet( i );
    }

    // Allocate new map
    ArrayOfSets< localIndex > newVar;
    newVar.resizeFromCapacities< parallelHostPolicy >( var.size(), newCapacities.data() );

    // Fill new map with old values
    forAll< parallelHostPolicy >( var.size(), [var = var.toViewConst(), newVar = newVar.toView()]( localIndex const k )
    {
      newVar.insertIntoSet( k, var[k].begin(), var[k].end() );
    } );

    // Insert new values
    for( localIndex i = 0; i < indiciesToInsert.size(); ++i )
    {
      localIndex const k = indiciesToInsert[i];
      newVar.insertIntoSet( k, valuesToInsert[i].begin(), valuesToInsert[i].end() );
    }

    // Replace the old map
    var = std::move( newVar );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, int USD0, int USD1 >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView2d< localIndex const, USD0 > const & var,
      arrayView1d< localIndex > const & indices,
      arraySlice1d< globalIndex const, USD1 > const & localToGlobalMap )
{
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    sizeOfPackedChars += PackArray< DO_PACKING >( buffer, var[li], var.size( 1 ) );
  }

  return sizeOfPackedChars;
}

template< typename SORTED, int USD >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView2d< localIndex, USD > const & var,
        array1d< localIndex > & indices,
        mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );

  for( localIndex a=0; a<numIndicesUnpacked; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi, MPI_REPLACE );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    localIndex * const varSlice = var[li];
    sizeOfUnpackedChars += UnpackPointer( buffer, varSlice, var.size( 1 ), MPI_REPLACE );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename SORTED, int USD0 >
localIndex
Pack( buffer_unit_type * & buffer,
      arrayView2d< localIndex const, USD0 > const & var,
      mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
      arrayView1d< localIndex const > const & indices,
      arraySlice1d< globalIndex const > const & localToGlobalMap,
      arraySlice1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars = 0;
  array1d< globalIndex > invalidGlobalIndices( var.size( 1 ) );
  for( localIndex a=0; a<var.size( 1 ); ++a )
  {
    invalidGlobalIndices[a] = -1;
  }

  sizeOfPackedChars += Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a=0; a<indices.size(); ++a )
  {
    localIndex li = indices[a];
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    typename mapBase< localIndex, array1d< globalIndex >, SORTED >::const_iterator
      iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > const & unmappedGI = iterUnmappedGI==unmappedGlobalIndices.end() ?
                                                invalidGlobalIndices :
                                                iterUnmappedGI->second;

    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.data(),
                                             var.size( 1 ),
                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}


template< typename SORTED0, typename SORTED1, typename SORTED2, int USD >
inline
localIndex
Unpack( buffer_unit_type const * & buffer,
        arrayView2d< localIndex, USD > const & var,
        localIndex_array & indices,
        mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
        mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
        mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Unpack( buffer, numIndicesUnpacked, MPI_REPLACE );
  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );
  array1d< globalIndex > unmappedIndices;

  for( localIndex a=0; a<numIndicesUnpacked; ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += Unpack( buffer, gi, MPI_REPLACE );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
      li = globalToLocalMap.at( gi );
    }

    arraySlice1d< localIndex, USD - 1 > varSlice = var[li];

    unmappedIndices.resize( 0 );
    sizeOfUnpackedChars += Unpack( buffer,
                                   varSlice,
                                   unmappedIndices,
                                   var.size( 1 ),
                                   relatedObjectGlobalToLocalMap );

    if( unmappedIndices.size()>0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }

  }

  return sizeOfUnpackedChars;
}

template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type * & buffer, MAP_TYPE const & var, T_INDICES const & packIndices )
{
  typename MAP_TYPE::size_type const length = var.size();
  localIndex sizeOfPackedChars = Pack< DO_PACKING >( buffer, length );

  for( typename MAP_TYPE::const_iterator i=var.begin(); i!=var.end(); ++i )
  {
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->first );
    sizeOfPackedChars += Pack< DO_PACKING >( buffer, i->second, packIndices );
  }

  return sizeOfPackedChars;
}


template< typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_host_packable_map_by_index_v< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const * & buffer, MAP_TYPE & map, T_INDICES const & unpackIndices )
{
  map.clear();
  typename MAP_TYPE::size_type map_length;

  localIndex sizeOfUnpackedChars = Unpack( buffer, map_length, MPI_REPLACE );

  for( typename MAP_TYPE::size_type a=0; a<map_length; ++a )
  {
    typename MAP_TYPE::key_type key;
    typename MAP_TYPE::mapped_type value;
    sizeOfUnpackedChars += Unpack( buffer, key, MPI_REPLACE );
    sizeOfUnpackedChars += Unpack( buffer, value, unpackIndices, MPI_REPLACE );

    map[key] = value;
  }

  return sizeOfUnpackedChars;
}

} /* namespace bufferOps */
} /* namespace geos */

#endif
