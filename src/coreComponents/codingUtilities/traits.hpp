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
 * @file traits.hpp
 */

#ifndef GEOS_CODINGUTILITIES_TRAITS_HPP_
#define GEOS_CODINGUTILITIES_TRAITS_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "LvArray/src/typeManipulation.hpp"
#include "LvArray/src/bufferManipulation.hpp"
#include "SFINAE_Macros.hpp"

// System includes
#include <type_traits>

namespace geos
{

namespace traits
{

/**
 * @brief Defines a static constexpr bool HasMemberFunction_data< @p CLASS >
 *        that is true iff the method @p CLASS ::data() exists and the return value is convertable to a pointer.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( data, void const *, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_begin< @p CLASS >
 *        that is true iff the method @p CLASS ::begin() exists and the return value is convertable to a pointer.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( begin, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_end< @p CLASS >
 *        that is true iff the method @p CLASS ::end() exists and the return value is convertable to a pointer.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( end, );

/// True iff @p T is a range-like type, i.e. has begin() and end() member functions
template< typename T >
constexpr bool is_range_like = HasMemberFunction_begin< T > && HasMemberFunction_end< T >;

/**
 * @brief Defines a static constexpr bool HasMemberFunction_move< @p CLASS >
 *        that is true iff the method @p CLASS ::move(LvArray::MemorySpace, bool) exists.
 * @tparam CLASS The type to test.
 */
template< typename CLASS >
static constexpr bool HasMemberFunction_move = LvArray::bufferManipulation::HasMemberFunction_move< CLASS >;

/**
 * @brief Defines a static constexpr bool HasMemberFunction_getPreviousSpace< @p CLASS >
 *        that is true iff the method @p CLASS ::getPreviousSpace() exists and the return value is convertable to a LvArray::MemorySpace.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( getPreviousSpace, LvArray::MemorySpace, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_registerTouch< @p CLASS >
 *        that is true iff the method @p CLASS ::registerTouch( LvArray::MemorySpace ) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( registerTouch, LvArray::MemorySpace() );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_checkTouch< @p CLASS >
 *        that is true iff the method @p CLASS :checkTouch( ) exists and the return value is converatble to bool.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( checkTouch, bool, );


/**
 * @brief Defines a static constexpr bool HasMemorySpaceFunctions< @p CLASS >
 *        that is true iff the class exposes the set a memory space movement functions defined above.
 * @tparam CLASS The type to test.
 */
template< typename CLASS >
static constexpr bool HasMemorySpaceFunctions = HasMemberFunction_move< CLASS > &&
                                                HasMemberFunction_registerTouch< CLASS > &&
                                                HasMemberFunction_getPreviousSpace< CLASS >;

/**
 * @brief Defines a static constexpr bool HasMemberFunction_setName< @p CLASS >
 *        that is true iff the method @p CLASS ::setName( string ) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( setName, string() );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_size< @p CLASS >
 *        that is true iff the method @p CLASS ::size() exists and the return value is convertable to a localIndex.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( size, localIndex, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_capacity< @p CLASS >
 *        that is true iff the method @p CLASS ::capacity() exists and the return value is convertable to a localIndex.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( capacity, localIndex, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_resize< @p CLASS >
 *        that is True iff the method @p CLASS ::resize( int ) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( resize, 0 );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_resizeDefault< @p CLASS >
 *        that is true iff the method @p CLASS ::resizeDefault( int, int, int) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( resizeDefault, 0, 0 );


/**
 * @brief Defines a static constexpr bool HasMemberFunction_reserve< @p CLASS >
 *        that is true iff the method @p CLASS ::reserve( localIndex ) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( reserve, localIndex( 55 ) );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_toView< @p CLASS >
 *        that is true iff the method @p CLASS ::toView() exists.
 * @tparam CLASS The type to test.
 */
template< typename CLASS >
static constexpr bool HasMemberFunction_toView = LvArray::typeManipulation::HasMemberFunction_toView< CLASS >;

/**
 * @brief Defines a static constexpr bool with two template parameter CanStreamInto
 *        that is true iff the operation std::declval< SRC & >() >> src::declval< DST & >() exists.
 * @tparam SRC The type of the source.
 * @tparam DST The type of the destination.
 */
IS_VALID_EXPRESSION_2( CanStreamInto, SRC, DST, std::declval< SRC & >() >> std::declval< DST & >() );

/**
 * @brief Defines a static constexpr bool HasMemberType_value_type< @p CLASS >
 *        that is true iff @p CLASS ::value_type is valid and not an enum.
 * @tparam CLASS The type to test.
 */
HAS_ALIAS( value_type );

namespace internal
{

/**
 * @brief Trait to check if `operator+=` is defined between two types.
 *
 * This struct is a type trait that determines if the `operator+=` is defined between two types, T and U.
 * It inherits from `std::false_type` by default, indicating that the operator is not defined.
 * A specialized version inherits from `std::true_type` if the operator is indeed defined.
 *
 * @tparam T The left-hand side type for the `operator+=`.
 * @tparam U The right-hand side type for the `operator+=`, defaults to the same as T.
 * @tparam Void Helper type for SFINAE.
 */
template< typename T, typename U, typename = void >
struct has_plus_equal : std::false_type {};

/**
 * @brief Specialization of has_plus_equal when `operator+=` is defined.
 *
 * This specialization of has_plus_equal uses SFINAE to detect if `T& operator+=(U)` exists.
 * If the operator is found, this struct inherits from `std::true_type`.
 *
 * @tparam T The left-hand side type for the `operator+=`.
 * @tparam U The right-hand side type for the `operator+=`.
 */
template< typename T, typename U >
struct has_plus_equal< T, U, std::void_t< decltype(std::declval< T & >() += std::declval< U >()) > > : std::true_type {};

} // namespace internal

/**
 * @brief Helper variable template to simplify usage of has_plus_equal trait.
 *
 * This variable template provides a more convenient way to use the has_plus_equal trait.
 * It evaluates to `true` if `T` supports `operator+=` with `U`, and `false` otherwise.
 *
 * @tparam T The left-hand side type for the `operator+=`.
 * @tparam U The right-hand side type for the `operator+=`, defaults to the same as T.
 */
template< typename T, typename U=T >
constexpr bool has_plus_equal_v = internal::has_plus_equal< T, U >::value;

namespace internal
{

/**
 * @brief Trait to check if `operator<` is defined between two types.
 *
 * This struct is a type trait that determines if the `operator<` is defined between two types, T and U.
 * It inherits from `std::false_type` by default, indicating that the operator is not defined.
 * A specialized version inherits from `std::true_type` if the operator is indeed defined.
 *
 * @tparam T The left-hand side type for the `operator<`.
 * @tparam U The right-hand side type for the `operator<`, defaults to the same as T.
 * @tparam Void Helper type for SFINAE.
 */
template< typename T, typename U, typename = void >
struct has_less_than : std::false_type {};

/**
 * @brief Specialization of has_less_than when `operator<` is defined.
 *
 * This specialization of has_less_than uses SFINAE to detect if `T& operator<(U)` exists.
 * If the operator is found, this struct inherits from `std::true_type`.
 *
 * @tparam T The left-hand side type for the `operator<`.
 * @tparam U The right-hand side type for the `operator<`.
 */
template< typename T, typename U >
struct has_less_than<T, U, std::void_t<decltype(std::declval< T & >() < std::declval< U >())>> : std::true_type {};

} // namespace internal

/**
 * @brief Helper variable template to simplify usage of has_less_than trait.
 *
 * This variable template provides a more convenient way to use the has_less_than trait.
 * It evaluates to `true` if `T` supports `operator<` with `U`, and `false` otherwise.
 *
 * @tparam T The left-hand side type for the `operator<`.
 * @tparam U The right-hand side type for the `operator<`, defaults to the same as T.
 */
template< typename T, typename U=T >
constexpr bool has_less_than_v = internal::has_less_than< T, U >::value;

namespace internal
{
template< class T,
          bool HAS_DATA_METHOD = HasMemberFunction_data< T > >
struct GetPointerType
{
  using Pointer = T *;
  using ConstPointer = T const *;
};

template< class T >
struct GetPointerType< T, true >
{
  using Pointer = decltype( std::declval< T >().data() );
  using ConstPointer = std::remove_pointer_t< Pointer > const *;
};

} // namespace internal

/// Type aliased to whatever T::data() returns or T * if that method doesn't exist.
template< typename T >
using Pointer = typename internal::GetPointerType< T >::Pointer;

/// The const version of Pointer.
template< typename T >
using ConstPointer = typename internal::GetPointerType< T >::ConstPointer;

/// Type aliased to whatever T::toView() returns or T & if that method doesn't exist.
template< typename T >
using ViewType = LvArray::typeManipulation::ViewType< T >;

/// Type aliased to whatever T::toViewConst() returns or T const & if that method doesn't exist.
template< typename T >
using ViewTypeConst = LvArray::typeManipulation::ViewTypeConst< T >;

/// True if T is or inherits from string.
template< typename T >
constexpr bool is_string = std::is_base_of< string, T >::value;

/// True if T is an instantiation of LvArray::Array.
template< typename T >
constexpr bool is_array = LvArray::isArray< T >;

/// True if T is an instantiation of LvArray::ArrayView.
template< typename T >
constexpr bool is_array_view = LvArray::isArrayView< T >;

/// True if T is an instantiation of LvArray::SortedArray.
template< typename T >
constexpr bool is_sorted_array = LvArray::isSortedArray< T >;

/// True if T is an instantiation of LvArray:SortedArrayView.
template< typename T >
constexpr bool is_sorted_array_view = LvArray::isSortedArrayView< T >;

/// True if T is an instantiation of LvArray::ArrayView or LvArray::Array
template< typename T >
constexpr bool is_array_type = traits::is_array_view< T > || traits::is_array< T >;

/// True if T is an instantiation of LvArray::SortedArrayView or LvArray::SortedArray
template< typename T >
constexpr bool is_sorted_array_type = traits::is_sorted_array_view< T > || traits::is_sorted_array< T >;

/// True if T is a Tensor class.
template< typename T >
constexpr bool is_tensorT = std::is_same< std::remove_const_t< T >, R1Tensor >::value ||
                            std::is_same< std::remove_const_t< T >, R1Tensor32 >::value ||
                            std::is_same< std::remove_const_t< T >, R2SymTensor >::value;

/// True of T has operator=() defined.
template< typename T >
struct hasCopyAssignmentOperatorImpl
{
private:
  template< typename U > static constexpr auto test( int )->decltype( std::declval< U >()=std::declval< U >(), bool () )
  { return true; }

  template< typename U > static constexpr auto test( ... )->bool
  { return false; }
public:
  static constexpr bool value = test< T >( 0 );
};
/// True if T has operator= defined, or it is arithmetic or an enum.
template< typename T >
static constexpr bool hasCopyAssignmentOp = hasCopyAssignmentOperatorImpl< T >::value ||
                                            std::is_arithmetic< T >::value ||
                                            std::is_enum< T >::value;

namespace internal
{

template< typename T, template< typename ... > class L >
constexpr std::size_t type_list_index( T, L<> ) { return 0; }

template< typename T, template< typename ... > class L, typename ... Ts >
constexpr std::size_t type_list_index( T, L< T, Ts ... > ) { return 0; }

template< typename T, template< typename ... > class L, typename U, typename ... Ts >
constexpr std::size_t type_list_index( T, L< U, Ts ... > ) { return 1 + type_list_index( T{}, L< Ts ... >{} ); }

} // namespace internal

/**
 * @brief Index of a given type in a type list.
 * @tparam T type to find
 * @tparam LIST list of types (any variadic type, such as std::tuple<>, camp::list<>, etc.)
 *
 *  If @p T is not found, returns the size of the tuple.
 *  This is consistent with @p end iterator for STL algorithms.
 */
template< typename T, typename LIST >
constexpr std::size_t type_list_index = internal::type_list_index( T{}, LIST{} );

} /* namespace traits */

template< typename T, bool COND >
struct add_const_if
{
  using type = typename std::conditional< COND, typename std::add_const< T >::type, T >::type;
};

template< typename T, bool COND >
using add_const_if_t = typename add_const_if< T, COND >::type;

} /* namespace geos */

#endif /* GEOS_CODINGUTILITIES_TRAITS_HPP_ */
