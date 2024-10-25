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
 * @file TypeDispatch.hpp
 *
 * Collection of tools to help dispatch templated functions on types
 */

#ifndef GEOS_COMMON_TYPEDISPATCH_HPP
#define GEOS_COMMON_TYPEDISPATCH_HPP

#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"

#include <camp/camp.hpp>

#include <unordered_map>

namespace geos
{

/**
 * @brief Namespace containing type dispatching facilities.
 */
namespace types
{

namespace internal
{

template< camp::idx_t NDIM >
struct LayoutList
{
  static_assert( NDIM < 5, "Dispatching on 5D arrays and higher not implemented" );
  using types = camp::list<>;
};

template<>
struct LayoutList< 1 >
{
  using types = camp::list< RAJA::PERM_I >;
};

template<>
struct LayoutList< 2 >
{
  using types = camp::list< RAJA::PERM_IJ,
                            RAJA::PERM_JI >;
};

template<>
struct LayoutList< 3 >
{
  using types = camp::list< RAJA::PERM_IJK,
                            RAJA::PERM_IKJ,
                            RAJA::PERM_JIK,
                            RAJA::PERM_JKI,
                            RAJA::PERM_KIJ,
                            RAJA::PERM_KJI >;
};

template<>
struct LayoutList< 4 >
{
  using types = camp::list< RAJA::PERM_IJKL,
                            RAJA::PERM_IJLK,
                            RAJA::PERM_IKJL,
                            RAJA::PERM_IKLJ,
                            RAJA::PERM_ILJK,
                            RAJA::PERM_ILKJ,
                            RAJA::PERM_JIKL,
                            RAJA::PERM_JILK,
                            RAJA::PERM_JKIL,
                            RAJA::PERM_JKLI,
                            RAJA::PERM_JLIK,
                            RAJA::PERM_JLKI,
                            RAJA::PERM_KIJL,
                            RAJA::PERM_KILJ,
                            RAJA::PERM_KJIL,
                            RAJA::PERM_KJLI,
                            RAJA::PERM_KLIJ,
                            RAJA::PERM_KLJI,
                            RAJA::PERM_LIJK,
                            RAJA::PERM_LIKJ,
                            RAJA::PERM_LJIK,
                            RAJA::PERM_LJKI,
                            RAJA::PERM_LKIJ,
                            RAJA::PERM_LKJI >;
};

template< int NDIM >
using Layouts = typename LayoutList< NDIM >::types;

// Expands a list of dimensions with all corresponding permutations
template< typename T >
struct AddLayoutsImpl;

template< camp::idx_t ... NDIMS >
struct AddLayoutsImpl< camp::list< camp::num< NDIMS >... > >
{
  using types = typename camp::join< camp::cartesian_product< camp::list< camp::num< NDIMS > >, Layouts< NDIMS > > ... >::type;
};

template< typename T >
using AddLayouts = typename AddLayoutsImpl< T >::types;

// Helper to convert a constructed camp::list of type tuples into an actual Array type
template< typename T >
struct ArrayType;

template< typename T, typename NDIM, typename LAYOUT >
struct ArrayType< camp::list< T, camp::list< NDIM, LAYOUT > > >
{
  using type = ::geos::Array< T, NDIM::value, LAYOUT >;
};

// Helper to apply a template to all types in a list
template< template< typename ... > class F, typename T >
struct ApplyImpl;

template< template< typename ... > class F, typename ... Ts >
struct ApplyImpl< F, camp::list< Ts ... > >
{
  using types = camp::list< typename F< Ts >::type ... >;
};

template< template< typename ... > class F, typename T >
using Apply = typename ApplyImpl< F, T >::types;

// Helper to increment values in an compile-time integer sequence
template< typename T, int N >
struct IncrementImpl;

template< camp::idx_t ... Is, int N >
struct IncrementImpl< camp::idx_seq< Is... >, N >
{
  using type = camp::idx_seq< Is + N ... >;
};

template< typename T, int N >
using Increment = typename IncrementImpl< T, N >::type;

} // namespace internal

/**
 * @brief Generate a sequence of integers from @p Begin up to (and including) @p End.
 */
template< integer Begin, integer End >
using IntegerSequence = internal::Increment< camp::make_idx_seq_t< End - Begin + 1 >, Begin >;

/**
 * @brief Construct a list of types.
 */
template< typename ... Ts >
using TypeList = camp::list< Ts... >;

/**
 * @brief Construct a list of list type.
 */
template< typename LIST >
using ListofTypeList = internal::Apply< camp::list, LIST >;

/**
 * @brief Concatenate multiple type lists.
 */
template< typename ... Ls >
using Join = typename camp::join< Ls ... >::type;

/**
 * @brief Construct a list of GEOSX multidimensional array types (geos::Array), containing all
 *        value types in type list @p TYPES and all dimensionalities in list of integral constants @p NDIMS.
 */
template< typename TYPES, typename NDIMS >
using ArrayTypes = internal::Apply< internal::ArrayType, camp::cartesian_product< TYPES, internal::AddLayouts< NDIMS > > >;

/**
 * @brief List of major integral types used in GEOSX type system.
 */
using IntegralTypes = TypeList< integer, localIndex, globalIndex >;

/**
 * @brief List of major real-valued (floating point) types used in GEOSX type system.
 */
using RealTypes = TypeList< real32, real64 >;

/**
 * @brief Generate a list of types representing array dimensionalities from @p M up to (and including) @p N.
 */
template< int M, int N >
using DimsRange = camp::as_list< IntegerSequence< M, N > >;

/**
 * @brief Generate a list of types representing array dimensionality exactly @p N.
 */
template< int N >
using DimsSingle = DimsRange< N, N >;

/**
 * @brief Generate a list of types representing array dimensionalities up to (and including) @p N.
 */
template< int N >
using DimsUpTo = DimsRange< 1, N >;

/**
 * @brief List of real-valued array types typically used in GEOSX (dimensions up to 4).
 */
using RealArrays = ArrayTypes< RealTypes, DimsUpTo< 4 > >;

/**
 * @brief List of integer-valued array types typically used in GEOSX (dimensions 1 and 2).
 */
using IntegralArrays = ArrayTypes< IntegralTypes, DimsUpTo< 2 > >;

/**
 * @brief List of all array types (real- and integer-valued) typically used in GEOSX.
 */
using StandardArrays = Join< RealArrays, IntegralArrays >;

namespace internal
{
/**
 * @brief struct to define the hash of a tuple
 *
 */
struct tuple_hash
{
  template< class... Ts >
  std::size_t operator()( const std::tuple< Ts... > & t ) const
  {
    std::size_t hash = 0;
    std::apply( [&hash]( auto &&... args )
    {
      ((hash ^= std::hash< std::decay_t< decltype(args) > >{} (args)), ...);
    }, t );
    return hash;
  }
};

/**
 * @brief Function to create a tuple
 * @tparam Ts tuple types
 */
template< typename ... Ts >
auto createTypeIndexTuple( TypeList< Ts... > )
{
  return std::make_tuple( std::type_index( typeid(Ts))... );
}

/**
 * @brief Get the static/singleton instance of the type map
 * @tparam LIST The
 * @tparam Is integer sequence
 * @return reference to the static map
 */
template< typename LIST, std::size_t ... Is >
auto const & getTypeMap( LIST, std::integer_sequence< std::size_t, Is... > )
{
  using KeyType = decltype( createTypeIndexTuple( camp::first< LIST > {} ) );
  static std::unordered_map< KeyType, std::size_t, tuple_hash > const result = { { createTypeIndexTuple( camp::at_t< LIST, camp::num< Is > >{} ), Is } ... };
  return result;
}


template< typename T >
struct TypeIdentity
{
  using type = T;
};

/**
 * @brief Function to output string containing the types of a TypeList
 * @tparam Ts The types contained in the TypeList
 * @param pre string to prepend the printer with
 * @param post string to postpend the printer with
 * @param printer a function that prints the type Ts
 * @return A string containing the types in the TypeList
 */
template< typename ... Ts, typename P >
string listToString( TypeList< Ts... >,
                     string const & pre,
                     string const & post,
                     P printer )
{
  return ( ( pre + printer( TypeIdentity< Ts >{} ) + post ) + ... );
}

HAS_MEMBER_FUNCTION_NO_RTYPE( getTypeId, );

/**
 * @brief Function to return the typeid() of a type or a wrapped type.
 * @tparam T the type of object or the wrapper
 * @param[in] object The object or wrapped object
 * @return the result of typeid()
 */
template< typename T >
std::type_info const & typeIdWrapper( T const & object )
{
  if constexpr ( HasMemberFunction_getTypeId< T > )
    return object.getTypeId();
  else
    return typeid(object);
}

/**
 * @brief
 *
 * @tparam TypeTuples
 * @tparam LAMBDA
 * @param combinations
 * @param type_index
 * @param lambda
 * @return true
 * @return false
 */
template< typename ... TypeTuples, typename LAMBDA, typename Index_type >
bool dispatchViaTable( TypeList< TypeTuples... > const combinations,
                       LAMBDA && lambda,
                       Index_type const type_index )
{
  static_assert( sizeof...(TypeTuples) > 0, "Dispatching on empty type list not supported" );

  // Initialize a table of handlers, once per unique combination of type list and lambda
  using Handler = void (*)( LAMBDA && );
  static Handler const handlers[] = { []( LAMBDA && f ){ f( TypeTuples{} ); } ... };

  auto const & typeIndexMap = getTypeMap( combinations, std::index_sequence_for< TypeTuples... >{} );
  auto const it = typeIndexMap.find( type_index );

  if( it != typeIndexMap.end() )
  {
    handlers[ it->second ]( std::forward< LAMBDA >( lambda ) );
    return true;
  }
  return false;
}

} // namespace internal

/**
 * @brief Dispatch a generic worker function @p lambda based on runtime type.
 *
 * @tparam LIST type of the list of types
 * @tparam LAMBDA  type of user-provided function or lambda, must have one auto argument
 * @tparam Ts types of the objects to be dispatched.
 * @param combinations list of types
 * @param lambda lambda user-provided callable
 * @param objects objects to be dispatched
 * @return @p true iff type has been dispatch
 *
 * todo Do we want errorIfTypeNotFound parameter? Options:
 *       - make it a template parameter (caller always knows whether or not it wants hard errors)
 *       - make the caller process return value and raise error if needed (and however they want)
 */
template< typename LIST, typename LAMBDA, typename ... Ts >
bool dispatch( LIST const combinations,
               LAMBDA && lambda,
               Ts && ... objects )
{
  bool const success = internal::dispatchViaTable( combinations,
                                                   std::forward< LAMBDA >( lambda ),
                                                   std::make_tuple( std::type_index( internal::typeIdWrapper( objects ))... ) );

  if( !success )
  {
    auto typePrinter = []( auto t ){ return LvArray::system::demangle( typeid( typename decltype(t)::type ).name() ); };
    auto typeListPrinter = [typePrinter]( auto tlist ){ return internal::listToString( typename decltype( tlist )::type{}, "\n  ", "", typePrinter ); };

    GEOS_ERROR( "Types were not dispatched. The types of the input objects are:\n" <<
                "( "<<(  ( "\n  " + LvArray::system::demangle( internal::typeIdWrapper( objects ).name() ) ) + ... )<<" \n)\n"<<
                "and the dispatch options are:\n"<<
                internal::listToString( combinations, "\n(", "\n)", typeListPrinter ) );
  }
  return success;
}

} // namespace types

} // namespace geos

#endif //GEOS_COMMON_TYPEDISPATCH_HPP
