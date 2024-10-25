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
 * @file DataTypes.hpp
 *
 * This file contains various aliases and functions that provide operations regarding the
 * use of the static types.
 */

#ifndef GEOS_COMMON_DATATYPES_HPP
#define GEOS_COMMON_DATATYPES_HPP

// Source includes
#include "common/GeosxConfig.hpp"
#include "GeosxMacros.hpp"
#include "BufferAllocator.hpp"
#include "DataLayouts.hpp"
#include "Tensor.hpp"
#include "LvArray/src/Macros.hpp"
#include "LvArray/src/Array.hpp"
#include "LvArray/src/ArrayOfArrays.hpp"
#include "LvArray/src/ArrayOfSets.hpp"
#include "LvArray/src/SparsityPattern.hpp"
#include "LvArray/src/CRSMatrix.hpp"
#include "LvArray/src/SortedArray.hpp"
#include "LvArray/src/StackBuffer.hpp"
#include "LvArray/src/ChaiBuffer.hpp"

#include "Path.hpp"

// TPL includes
#include <camp/camp.hpp>

// System includes
#ifdef GEOS_USE_MPI
  #include <mpi.h>
#endif

#include <cassert>
//#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <vector>

/*
 * top level geos namespace contains all code that is specific to GEOSX
 */
namespace geos
{

/**
 * @name Basic data types used in GEOSX.
 */
///@{

/// Unsigned size type.
using size_t      = std::size_t;

/// Signed integer type.
using integer     = std::int32_t;

/// Local index type (for indexing objects within an MPI partition).
using localIndex  = GEOS_LOCALINDEX_TYPE;

/// Global index type (for indexing objects across MPI partitions).
using globalIndex = GEOS_GLOBALINDEX_TYPE;

/// String type.
using string      = std::string;

/// String type.
using string_view = std::string_view;

/// 32-bit floating point type.
using real32 = float;
/// 64-bit floating point type.
using real64 = double;

///@}

/**
 * @name Binary buffer data types.
 */
///@{

/// Type stored in communication buffers.
using buffer_unit_type = signed char;

#ifdef GEOS_USE_CHAI
/// Type of storage for communication buffers.
using buffer_type = std::vector< buffer_unit_type, BufferAllocator< buffer_unit_type > >;
#else
/// Type of storage for communication buffers.
using buffer_type = std::vector< buffer_unit_type >;
#endif

///@}

/**
 * @name Alias for timestamp in GEOSX.
 */
///@{

/// Timestamp type (used to perform actions such a sparsity pattern computation after mesh modifications)
using Timestamp = unsigned long long int;

///@}

//START_SPHINX_INCLUDE_00

/**
 * @name Aliases for LvArray::Array class family.
 */
///@{

/// Multidimensional array type. See LvArray:Array for details.
template< typename T,
          int NDIM,
          typename PERMUTATION=camp::make_idx_seq_t< NDIM > >
using Array = LvArray::Array< T, NDIM, PERMUTATION, localIndex, LvArray::ChaiBuffer >;

/// Multidimensional array view type. See LvArray:ArrayView for details.
template< typename T,
          int NDIM,
          int USD = NDIM - 1 >
using ArrayView = LvArray::ArrayView< T, NDIM, USD, localIndex, LvArray::ChaiBuffer >;

/// Multidimensional array slice type. See LvArray:ArraySlice for details.
template< typename T, int NDIM, int USD = NDIM - 1 >
using ArraySlice = LvArray::ArraySlice< T, NDIM, USD, localIndex >;

/// Multidimensional stack-based array type. See LvArray:StackArray for details.
template< typename T, int NDIM, int MAXSIZE, typename PERMUTATION=camp::make_idx_seq_t< NDIM > >
using StackArray = LvArray::StackArray< T, NDIM, PERMUTATION, localIndex, MAXSIZE >;

///@}

/**
 * @name Short-hand aliases for commonly used array types.
 */
///@{

/// Alias for a local (stack-based) rank-1 tensor type
using R1Tensor = Tensor< real64, 3 >;
/// Alias for a local (stack-based) rank-1 tensor type using 32 bits integers
using R1Tensor32 = Tensor< real32, 3 >;

/// Alias for a local (stack-based) rank-2 Voigt tensor type
using R2SymTensor = Tensor< real64, 6 >;


/// Alias for 1D array.
template< typename T >
using array1d = Array< T, 1 >;

/// Alias for 1D array view.
template< typename T >
using arrayView1d = ArrayView< T, 1 >;

/// Alias for 1D array slice.
template< typename T, int USD = 0 >
using arraySlice1d = ArraySlice< T, 1, USD >;

/// Alias for 1D stack array.
template< typename T, int MAXSIZE >
using stackArray1d = StackArray< T, 1, MAXSIZE >;

/// Alias for 2D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 2 > >
using array2d = Array< T, 2, PERMUTATION >;

/// Alias for 2D array view.
template< typename T, int USD = 1 >
using arrayView2d = ArrayView< T, 2, USD >;

/// Alias for 2D array slice.
template< typename T, int USD = 1 >
using arraySlice2d = ArraySlice< T, 2, USD >;

/// Alias for 2D stack array.
template< typename T, int MAXSIZE >
using stackArray2d = StackArray< T, 2, MAXSIZE >;

/// Alias for 3D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 3 > >
using array3d = Array< T, 3, PERMUTATION >;

/// Alias for 3D array view.
template< typename T, int USD=2 >
using arrayView3d = ArrayView< T, 3, USD >;

/// Alias for 3D array slice.
template< typename T, int USD=2 >
using arraySlice3d = ArraySlice< T, 3, USD >;

/// Alias for 3D stack array.
template< typename T, int MAXSIZE >
using stackArray3d = StackArray< T, 3, MAXSIZE >;

/// Alias for 4D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 4 > >
using array4d = Array< T, 4, PERMUTATION >;

/// Alias for 4D array view.
template< typename T, int USD=3 >
using arrayView4d = ArrayView< T, 4, USD >;

/// Alias for 4D array slice.
template< typename T, int USD=3 >
using arraySlice4d = ArraySlice< T, 4, USD >;

/// Alias for 4D stack array.
template< typename T, int MAXSIZE >
using stackArray4d = StackArray< T, 4, MAXSIZE >;

/// Alias for 5D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 5 > >
using array5d = Array< T, 5, PERMUTATION >;

/// Alias for 5D array view.
template< typename T, int USD=4 >
using arrayView5d = ArrayView< T, 5, USD >;

/// Alias for 5D array slice.
template< typename T, int USD=4 >
using arraySlice5d = ArraySlice< T, 5, 4 >;

/// Alias for 5D stack array.
template< typename T, int MAXSIZE >
using stackArray5d = StackArray< T, 5, MAXSIZE >;

///@}

/**
 * @name Aliases for sorted arrays and set types.
 */
///@{

/// A set of local indices.
template< typename T >
using set = std::set< T >;

/// A sorted array of local indices.
template< typename T >
using SortedArray = LvArray::SortedArray< T, localIndex, LvArray::ChaiBuffer >;

/// A sorted array view of local indices.
template< typename T >
using SortedArrayView = LvArray::SortedArrayView< T, localIndex, LvArray::ChaiBuffer >;

///@}

/**
 * @name Aliases for LvArray::ArrayOfArrays class family.
 */
///@{

/// Array of variable-sized arrays. See LvArray::ArrayOfArrays for details.
template< typename T, typename INDEX_TYPE=localIndex >
using ArrayOfArrays = LvArray::ArrayOfArrays< T, INDEX_TYPE, LvArray::ChaiBuffer >;

/// View of array of variable-sized arrays. See LvArray::ArrayOfArraysView for details.
template< typename T, typename INDEX_TYPE=localIndex, bool CONST_SIZES=std::is_const< T >::value >
using ArrayOfArraysView = LvArray::ArrayOfArraysView< T, INDEX_TYPE const, CONST_SIZES, LvArray::ChaiBuffer >;

/// Array of variable-sized sets. See LvArray::ArrayOfSets for details.
template< typename T, typename INDEX_TYPE=localIndex >
using ArrayOfSets = LvArray::ArrayOfSets< T, INDEX_TYPE, LvArray::ChaiBuffer >;

/// View of array of variable-sized sets. See LvArray::ArrayOfSetsView for details.
template< typename T, typename INDEX_TYPE=localIndex >
using ArrayOfSetsView = LvArray::ArrayOfSetsView< T, INDEX_TYPE const, LvArray::ChaiBuffer >;

/// Alias for Sparsity pattern class.
template< typename COL_INDEX, typename INDEX_TYPE=localIndex >
using SparsityPattern = LvArray::SparsityPattern< COL_INDEX, INDEX_TYPE, LvArray::ChaiBuffer >;

/// Alias for Sparsity pattern View.
template< typename COL_INDEX, typename INDEX_TYPE=localIndex >
using SparsityPatternView = LvArray::SparsityPatternView< COL_INDEX, INDEX_TYPE const, LvArray::ChaiBuffer >;

/// Alias for CRS Matrix class.
template< typename T, typename COL_INDEX=globalIndex >
using CRSMatrix = LvArray::CRSMatrix< T, COL_INDEX, localIndex, LvArray::ChaiBuffer >;

/// Alias for CRS Matrix View.
template< typename T, typename COL_INDEX=globalIndex >
using CRSMatrixView = LvArray::CRSMatrixView< T, COL_INDEX, localIndex const, LvArray::ChaiBuffer >;

///@}

//END_SPHINX_INCLUDE_00

/**
 * @name Ordered and unordered map types.
 */
///@{

/**
 * @brief Base template for ordered and unordered maps.
 * @tparam TKEY key type
 * @tparam TVAL value type
 * @tparam SORTED a @p std::integral_constant<bool> indicating whether map is ordered
 */
template< typename TKEY, typename TVAL, typename SORTED >
class mapBase
{};

/// @cond DO_NOT_DOCUMENT
template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, true > > : public std::map< TKEY, TVAL >
{
  using std::map< TKEY, TVAL >::map; // enable list initialization
};

template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, false > > : public std::unordered_map< TKEY, TVAL >
{
  using std::unordered_map< TKEY, TVAL >::unordered_map; // enable list initialization
};
/// @endcond

/**
 * @brief Stream output operator for map types.
 * @tparam K key type
 * @tparam V value type
 * @tparam SORTED
 * @param stream output stream
 * @param map the map to print
 * @return reference to output stream
 */
template< typename K, typename V, typename SORTED >
inline
std::ostream & operator<< ( std::ostream & stream, mapBase< K, V, SORTED > const & map )
{
  stream << "{\n";
  for( auto const & pair : map )
  {
    stream << pair.first << " : " << pair.second << "\n";
  }
  stream << "}";
  return stream;
}

/// Ordered map type.
template< typename TKEY, typename TVAL >
using map = mapBase< TKEY, TVAL, std::integral_constant< bool, true > >;

/// Unordered map type.
template< typename TKEY, typename TVAL >
using unordered_map = mapBase< TKEY, TVAL, std::integral_constant< bool, false > >;

///@}

/**
 * @name Aliases for commonly used array types.
 */
///@{

/// A 1-dimensional array of geos::integer types.
using integer_array = array1d< integer >;

/// A 1-dimensional array of geos::real32 types.
using real32_array = array1d< real32 >;

/// A 1-dimensional array of geos::real64 types.
using real64_array = array1d< real64 >;

/// A 1-dimensional array of geos::string types.
using string_array = array1d< string >;

/// A 1-dimensional array of geos::Path types
using path_array = array1d< Path >;

/// A 1-dimensional array of geos::localIndex types
using localIndex_array = array1d< localIndex >;

/// A 1-dimensional array of geos::globalIndex types
using globalIndex_array = array1d< globalIndex >;


/// A 2-dimensional array of geos::integer types.
using integer_array2d = array2d< integer >;

/// A 2-dimensional array of geos::real32 types.
using real32_array2d = array2d< real32 >;

/// A 2-dimensional array of geos::real64 types.
using real64_array2d = array2d< real64 >;

/// A 2-dimensional array of geos::localIndex types
using localIndex_array2d = array2d< localIndex >;

/// A 2-dimensional array of geos::globalIndex types
using globalIndex_array2d = array2d< globalIndex >;


/// A 3-dimensional array of geos::integer types.
using integer_array3d = array3d< integer >;

/// A 3-dimensional array of geos::real32 types.
using real32_array3d = array3d< real32 >;

/// A 3-dimensional array of geos::real64 types.
using real64_array3d = array3d< real64 >;

/// A 3-dimensional array of geos::localIndex types.
using localIndex_array3d = array3d< localIndex >;

/// A 3-dimensional array of geos::globalIndex types.
using globalIndex_array3d = array3d< globalIndex >;


/// A 4-dimensional array of geos::integer types.
using integer_array4d = array4d< integer >;

/// A 4-dimensional array of geos::real32 types.
using real32_array4d = array4d< real32 >;

/// A 4-dimensional array of geos::real64 types.
using real64_array4d = array4d< real64 >;

/// A 4-dimensional array of geos::localIndex types.
using localIndex_array4d = array4d< localIndex >;

/// A 4-dimensional array of geos::globalIndex types.
using globalIndex_array4d = array4d< globalIndex >;

///@}

/// A variable for the maximum value of a geos::globalIndex.
constexpr static auto GLOBALINDEX_MAX = std::numeric_limits< globalIndex >::max();

/// A variable for the maximum value of a geos::localIndex.
constexpr static auto LOCALINDEX_MAX = std::numeric_limits< localIndex >::max();

/// A global variable for the value of a object that has not been assigned a geos::globalIndex.
constexpr static localIndex unmappedLocalIndexValue = -1;

} /* namespace geos */

#endif /* GEOS_COMMON_DATATYPES_HPP */
