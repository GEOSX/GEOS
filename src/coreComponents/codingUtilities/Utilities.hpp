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
 * @file Utilities.hpp
 */

#ifndef GEOS_CODINGUTILITIES_UTILITIES_H_
#define GEOS_CODINGUTILITIES_UTILITIES_H_

#include "common/format/StringUtilities.hpp"
#include "common/logger/Logger.hpp"
#include "common/DataTypes.hpp"
#include "LvArray/src/limits.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

/**
 * @brief Compare two real values with a tolerance.
 * @tparam type of real value
 * @param val1 first value
 * @param val2 second value
 * @param relTol relative tolerance for comparison
 * @return @p true if @p val1 and @p val2 are within tolerance of each other
 */
template< typename T >
GEOS_FORCE_INLINE GEOS_HOST_DEVICE constexpr
bool isEqual( T const val1, T const val2, T const relTol = 0.0 )
{
  T const absTol = ( relTol > 10 * LvArray::NumericLimits< T >::epsilon ) ? relTol * ( fabs( val1 ) + fabs( val2 ) ) * 0.5 : 0.0;
  return ( val2 - absTol ) <= val1 && val1 <= ( val2 + absTol );
}

/**
 * @brief Test if a real value is (almost) zero.
 * @tparam T type of real value
 * @param val the value to test
 * @param tol absolute tolerance for comparison
 * @return @p true if @p val is within @p tol of zero
 */
template< typename T >
GEOS_FORCE_INLINE GEOS_HOST_DEVICE constexpr
bool isZero( T const val, T const tol = LvArray::NumericLimits< T >::epsilon )
{
  return -tol <= val && val <= tol;
}

template< typename T >
GEOS_FORCE_INLINE GEOS_HOST_DEVICE constexpr
bool isOdd( T x )
{
  static_assert( std::is_integral< T >::value, "Not meaningful for non-integral types" );
  return (x & 1);
}

template< typename T >
GEOS_FORCE_INLINE GEOS_HOST_DEVICE constexpr
bool isEven( T x )
{
  return !isOdd( x );
}

template< typename T1, typename T2, typename SORTED >
T2 & stlMapLookup( mapBase< T1, T2, SORTED > & Map, const T1 & key )
{
  typename mapBase< T1, T2, SORTED >::iterator MapIter = Map.find( key );
  GEOS_ERROR_IF( MapIter==Map.end(), "Key not found: " << key );
  return MapIter->second;
}


template< typename T1, typename T2, typename SORTED >
const T2 & stlMapLookup( const mapBase< T1, T2, SORTED > & Map, const T1 & key )
{
  return (stlMapLookup( const_cast< mapBase< T1, T2, SORTED > & >(Map), key ));
}

template< typename T1, typename T2, typename SORTED, typename LAMBDA >
bool executeOnMapValue( mapBase< T1, T2, SORTED > const & Map, const T1 & key, LAMBDA && lambda )
{
  bool rval = false;
  typename mapBase< T1, T2, SORTED >::const_iterator MapIter = Map.find( key );
  if( MapIter!=Map.end() )
  {
    rval = true;
    lambda( MapIter->second );
  }

  return rval;
}

template< typename T_KEY, typename T_VALUE, typename SORTED >
T_VALUE softMapLookup( mapBase< T_KEY, T_VALUE, SORTED > const & theMap,
                       T_KEY const & key,
                       T_VALUE const failValue )
{
  T_VALUE rvalue;
  typename mapBase< T_KEY, T_VALUE, SORTED >::const_iterator iter = theMap.find( key );
  if( iter==theMap.end() )
  {
    rvalue = failValue;
  }
  else
  {
    rvalue = iter->second;
  }
  return rvalue;
}

/**
 * @brief Call a user function on sub-ranges of repeating values.
 * @tparam ITER type of range iterator
 * @tparam FUNC type of user function
 * @tparam COMP type of comparison function
 * @param first iterator to start of the input range
 * @param last iterator past the end of the input range
 * @param func the function to call, must be callable with a pair of iterators
 *
 * User function will be called once per consecutive group of equal values, as defined
 * by @p comp, with a pair of iterators to the beginning and past the end of each group.
 * For example, given an input [1,1,2,2,2,3,3,1], @p func will be called with iterators
 * to sub-ranges [1,1], [2,2,2], [3,3], [1].
 *
 * @note @p comp is an equality comparison, not a less-than type predicate used in sorting.
 *       For a range sorted with some predicate P, one can use comp(x,y) = !(P(x,y) || P(y,x)),
 *       but in most cases default value (std::equal_to<>) should be fine.
 */
template< typename ITER, typename FUNC, typename COMP = std::equal_to<> >
void forEqualRanges( ITER first, ITER const last, FUNC && func, COMP && comp = {} )
{
  using R = typename std::iterator_traits< ITER >::reference;
  while( first != last )
  {
    R curr = *first;
    auto const pred = [&]( R v ) { return comp( curr, v ); };
    ITER const next = std::find_if_not( std::next( first ), last, pred );
    func( first, next );
    first = next;
  }
}

/**
 * @brief Execute a user function on unique values in a range.
 * @tparam ITER type of range iterator
 * @tparam FUNC type of user function
 * @tparam COMP type of comparison function
 * @param first iterator to start of the range
 * @param last iterator past the end of the range
 * @param func the function to call, must be callable with value and size (as int)
 *
 * User function will be called once per consecutive group of equal values, as defined
 * by @p comp, with a reference to one of the values from each group and the group size.
 * If the range is (partially) ordered w.r.t. to a predicate compatible with @p comp,
 * the function will be called once per unique value in the entire range.
 * For example, given an input [a,a,b,b,b,c,c,a], @p func will be called with
 * (a,2), (b,3), (c,2), (a,1), where a,b,c will be references to one of the
 * corresponding elements in the input (which one exactly is unspecified).
 *
 * @note @p comp is an equality comparison, not a less-than type predicate used in sorting.
 *       For a range sorted with some predicate P, one can use comp(x,y) = !(P(x,y) || P(y,x)),
 *       but in most cases default value (std::equal_to<>) should be fine.
 */
template< typename ITER, typename FUNC, typename COMP = std::equal_to<> >
void forUniqueValues( ITER first, ITER const last, FUNC && func, COMP && comp = {} )
{
  auto const f = [&func]( ITER r_first, ITER r_last )
  {
    func( *r_first, std::distance( r_first, r_last ) );
  };
  forEqualRanges( first, last, f, std::forward< COMP >( comp ) );
}

/**
 * @brief Perform lookup in a map of options and throw a user-friendly exception if not found.
 * @tparam KEY map key type
 * @tparam VAL map value type
 * @tparam SORTED whether map is ordered or unordered
 * @param map the option map
 * @param option the lookep key
 * @param optionName name of the option to use in exception error message
 * @param contextName name of the lookup context (e.g. the data repository group)
 */
template< typename KEY, typename VAL, typename SORTED >
VAL findOption( mapBase< KEY, VAL, SORTED > const & map,
                KEY const & option,
                string const & optionName,
                string const & contextName )
{
  auto const iter = map.find( option );
  GEOS_THROW_IF( iter == map.end(),
                 GEOS_FMT( "{}: unsupported option '{}' for {}.\nSupported options are: {}",
                           contextName,
                           const_cast< typename std::remove_const< KEY & >::type >( option ),
                           optionName,
                           stringutilities::join( mapKeys( map ), ", " ) ),
                 InputError );
  return iter->second;
}

namespace
{

/**
 * @brief Apply functor @p transformer onto the map elements and store the results into a container of type @p C.
 * @tparam MAP Type of the considered map.
 * @tparam C Type of the container holding the keys.
 * @tparam TRANSFORMER Type of the unary functor.
 * @param[in] map The map from which keys will be extracted.
 * @param transformer Which unary functor to apply to the @p map elements.
 * @return The container with the results for the @p transformer application.
 */
template< template< typename ... > class C, typename MAP, typename TRANSFORMER >
auto mapTransformer( MAP const & map,
                     TRANSFORMER const & transformer )
{
  using v = std::invoke_result_t< TRANSFORMER, typename MAP::const_reference >;
  C< v > result;
  auto inserter = std::inserter( result, result.end() );
  std::transform( map.begin(), map.end(), inserter, transformer );
  return result;
}

}

/**
 * @brief Extract the keys from the given map.
 * @tparam MAP Type of the considered map.
 * @tparam C Type of the container holding the keys.
 * @param[in] map The map from which keys will be extracted.
 * @return The container with the keys.
 */
template< template< typename ... > class C = std::vector, typename MAP >
C< typename MAP::key_type > mapKeys( MAP const & map )
{
  auto transformer = []( auto const & p ) -> typename MAP::key_type
  { return p.first; };
  return mapTransformer< C >( map, transformer );
}

/**
 * @brief Extract the values from the given map.
 * @tparam MAP Type of the considered map.
 * @tparam C Type of the container holding the values.
 * @param[in] map The map from which values will be extracted.
 * @return The container with the values.
 */
template< template< typename ... > class C = std::vector, typename MAP >
C< typename MAP::mapped_type > mapValues( MAP const & map )
{
  auto transformer = []( typename MAP::const_reference p ) -> typename MAP::mapped_type
  { return p.second; };
  return mapTransformer< C >( map, transformer );
}

namespace internal
{
template< class F, class ... Ts, std::size_t ... Is >
void forEachArgInTuple( std::tuple< Ts ... > const & tuple, F && func, std::index_sequence< Is ... > )
{
  using expander = int[];
  (void) expander { 0, ( (void)func( std::get< Is >( tuple ), std::integral_constant< size_t, Is >{} ), 0 )... };
}

template< class F, class ... Ts, std::size_t ... Is >
void forEachArgInTuple( std::tuple< Ts ... > & tuple, F && func, std::index_sequence< Is ... > )
{
  using expander = int[];
  (void) expander { 0, ( (void)func( std::get< Is >( tuple ), std::integral_constant< size_t, Is >{} ), 0 )... };
}

}

/**
 * @brief Visit every element in a tuple applying a function.
 * @tparam F type of function
 * @tparam Ts types of tuple elements
 * @param tuple the target tuple
 * @param func the function to apply
 *
 * The function will be called with a reference to the tuple element and
 * a compile-time (std::integral_constant) index of the tuple element.
 */
template< class F, class ... Ts >
void forEachArgInTuple( std::tuple< Ts ... > const & tuple, F && func )
{
  internal::forEachArgInTuple( tuple, std::forward< F >( func ), std::make_index_sequence< sizeof...( Ts ) >() );
}

/**
 * @brief Visit every element in a tuple applying a function.
 * @tparam F type of function
 * @tparam Ts types of tuple elements
 * @param tuple the target tuple
 * @param func the function to apply
 *
 * The function will be called with a reference to the tuple element and
 * a compile-time (std::integral_constant) index of the tuple element.
 */
template< class F, class ... Ts >
void forEachArgInTuple( std::tuple< Ts ... > & tuple, F && func )
{
  internal::forEachArgInTuple( tuple, std::forward< F >( func ), std::make_index_sequence< sizeof...( Ts ) >() );
}

/**
 * @brief Utility function to convert the value of an enumerator to its underlying type (integer).
 * @tparam ENUMERATION the type of the enumeration
 * @param[in] value the value of the enumerator
 * @return the integer conversion of @p value
 */
template< typename ENUMERATION >
constexpr std::underlying_type_t< ENUMERATION > toUnderlying( ENUMERATION const value )
{
  return static_cast< std::underlying_type_t< ENUMERATION > >( value );
}

/**
 * @brief Utility function to convert a pointer to an enumeration to a pointer to its underlying type (integer).
 * @tparam ENUMERATION the type of the enumeration
 * @param[in] enumPtr the pointer to the enumeration
 * @return the pointer to the enumeration underlying type
 */
template< typename ENUMERATION >
std::underlying_type_t< ENUMERATION > * toUnderlyingPtr( ENUMERATION * const enumPtr )
{
  return reinterpret_cast< std::underlying_type_t< ENUMERATION > * >( enumPtr );
}

// The code below should work with any subscriptable vector/matrix types

/**
 * @brief Utility function to copy a vector into another vector
 * @tparam VEC1 the type of the source vector
 * @tparam VEC2 the type of the destination vector
 * @param[in] N the number of values to copy
 * @param[in] v1 the source vector
 * @param[out] v2 the destination vector
 * @param[in] offset the first index of v2 at which we start copying values
 */
template< typename VEC1, typename VEC2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void copy( integer const N,
           VEC1 const & v1,
           VEC2 const & v2,
           integer const offset = 0 )
{
  for( integer i = 0; i < N; ++i )
  {
    v2[offset+i] = v1[i];
  }
}

/**
 * @brief Utility function to apply the chain rule to compute df_dx as a function of df_dy and dy_dx
 * @tparam MATRIX the type of the matrix of derivatives
 * @tparam VEC1 the type of the source vector
 * @tparam VEC2 the type of the destination vector
 * @param[in] N the number of derivative values
 * @param[in] dy_dx the matrix of derivatives of y(x) wrt x
 * @param[in] df_dy the derivatives of f wrt y
 * @param[out] df_dx the computed derivatives of f wrt x
 * @param[in] firstDerivativeOffset the first derivative offset in df_dy
 */
template< typename MATRIX, typename VEC1, typename VEC2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void applyChainRule( integer const N,
                     MATRIX const & dy_dx,
                     VEC1 const & df_dy,
                     VEC2 && df_dx,
                     integer const firstDerivativeOffset = 0 )
{
  // this could use some dense linear algebra
  for( integer i = 0; i < N; ++i )
  {
    df_dx[i] = 0.0;
    for( integer j = 0; j < N; ++j )
    {
      df_dx[i] += df_dy[firstDerivativeOffset+j] * dy_dx[j][i];
    }
  }
}

/**
 * @brief Utility function to apply the chain rule to compute df_dxy in place
 * @tparam MATRIX the type of the matrix of derivatives
 * @tparam VEC1 the type of the source vector
 * @tparam VEC2 the type of the utility vector
 * @param[in] N the number of derivative values
 * @param[in] dy_dx the matrix of derivatives of y(x) wrt x
 * @param[out] df_dxy the derivatives of f wrt y
 * @param[out] work utility array
 * @param[in] offset the first derivative offset in df_dy
 */
template< typename MATRIX, typename VEC1, typename VEC2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void applyChainRuleInPlace( integer const N,
                            MATRIX const & dy_dx,
                            VEC1 && df_dxy,
                            VEC2 && work,
                            integer const firstDeriv = 0 )
{
  applyChainRule( N, dy_dx, df_dxy, work, firstDeriv );
  copy( N, work, df_dxy, firstDeriv );
}

/**
 * @brief Internal struct to provide no-op defaults used in the inclusion
 *   of lambda functions into kernel component functions.
 * @struct NoOpFunc
 */
struct NoOpFunc
{
  template< typename ... Ts >
  GEOS_HOST_DEVICE
  constexpr void
  operator()( Ts && ... ) const {}
};

template< typename FLAGS_ENUM >
struct BitFlags
{
  GEOS_HOST_DEVICE
  void set( FLAGS_ENUM flag )
  {
    m_FlagValue |= (integer)flag;
  }

  GEOS_HOST_DEVICE
  bool isSet( FLAGS_ENUM flag ) const
  {
    return (m_FlagValue & (integer)flag) == (integer)flag;
  }

private:

  using FlagsType = std::underlying_type_t< FLAGS_ENUM >;
  FlagsType m_FlagValue{};

};

// Temporary functions (axpy, scale, dot) used in Aitken's acceleration. Will be removed once the nonlinear
// acceleration implementation scheme will use LAI vectors. See issue #2891
// (https://github.com/GEOS-DEV/GEOS/issues/2891)

template< typename VECTOR_TYPE, typename SCALAR_TYPE >
VECTOR_TYPE axpy( VECTOR_TYPE const & vec1,
                  VECTOR_TYPE const & vec2,
                  SCALAR_TYPE const alpha )
{
  GEOS_ASSERT( vec1.size() == vec2.size() );
  const localIndex N = vec1.size();
  VECTOR_TYPE result( N );
  RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, N ),
                                      [&] GEOS_HOST ( localIndex const i )
    {
      result[i] = vec1[i] + alpha * vec2[i];
    } );
  return result;
}

template< typename VECTOR_TYPE, typename SCALAR_TYPE >
VECTOR_TYPE scale( VECTOR_TYPE const & vec,
                   SCALAR_TYPE const scalarMult )
{
  const localIndex N = vec.size();
  VECTOR_TYPE result( N );
  RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, N ),
                                      [&] GEOS_HOST ( localIndex const i )
    {
      result[i] = scalarMult * vec[i];
    } );
  return result;
}

template< typename VECTOR_TYPE >
real64 dot( VECTOR_TYPE const & vec1,
            VECTOR_TYPE const & vec2 )
{
  GEOS_ASSERT( vec1.size() == vec2.size());
  RAJA::ReduceSum< parallelHostReduce, real64 > result( 0.0 );
  const localIndex N = vec1.size();
  RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, N ),
                                      [&] GEOS_HOST ( localIndex const i )
    {
      result += vec1[i] * vec2[i];
    } );
  return result.get();
}

} // namespace geos

#endif /* GEOS_CODINGUTILITIES_UTILITIES_H_ */
