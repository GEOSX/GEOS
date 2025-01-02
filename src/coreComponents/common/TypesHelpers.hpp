/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TypesHelpers.hpp
 *
 */

#ifndef TYPES_HELPERS_HPP
#define TYPES_HELPERS_HPP

#include <type_traits>

namespace geos
{

namespace internal
{
/**
 * @brief Trait to determine if a type defines a `value_type` member.
 *
 * This primary template defaults to `std::false_type`, indicating that
 * the type `T` does not define a `value_type` member.
 *
 * @tparam T The type to check.
 * @tparam void A SFINAE parameter used to specialize the trait.
 */
template< typename T, typename = void >
struct has_value_type : std::false_type {};

/**
 * @brief Specialization of `has_value_type` for types with a `value_type` member.
 *
 * If the type `T` defines a `value_type` member, this specialization
 * is used, which inherits from `std::true_type`.
 *
 * @tparam T The type to check.
 */
template< typename T >
struct has_value_type< T, std::void_t< typename T::value_type > > : std::true_type {};

/**
 * @brief Trait to determine if a type defines a `ValueType` member.
 *
 * This primary template defaults to `std::false_type`, indicating that
 * the type `T` does not define a `ValueType` member.
 *
 * @tparam T The type to check.
 * @tparam void A SFINAE parameter used to specialize the trait.
 */
template< typename T, typename = void >
struct has_ValueType : std::false_type {};

/**
 * @brief Specialization of `has_ValueType` for types with a `ValueType` member.
 *
 * If the type `T` defines a `ValueType` member, this specialization
 * is used, which inherits from `std::true_type`.
 *
 * @tparam T The type to check.
 */
template< typename T >
struct has_ValueType< T, std::void_t< typename T::ValueType > > : std::true_type {};

}   // namespace internal

/**
 * @brief Trait to retrieve the `value_type` or `ValueType` of a type `T`.
 *
 * This primary template provides a static assertion error if `T` does not
 * define either `value_type` or `ValueType`.
 *
 * @tparam T The type from which to extract the type alias.
 * @tparam Enable A SFINAE parameter used for specialization.
 */
template< typename T, typename Enable = void >
struct get_value_type
{
  static_assert( sizeof(T) == 0, "T must define either value_type or ValueType." );
};

/**
 * @brief Specialization of `get_value_type` for types with a `value_type` member.
 *
 * If the type `T` defines a `value_type` member, this specialization
 * retrieves it as the alias `type`.
 *
 * @tparam T The type from which to extract `value_type`.
 */
template< typename T >
struct get_value_type< T, std::enable_if_t< internal::has_value_type< T >::value > >
{
  using type = typename T::value_type;
};

/**
 * @brief Specialization of `get_value_type` for types with a `ValueType` member.
 *
 * If the type `T` does not define a `value_type` but defines a `ValueType`,
 * this specialization retrieves it as the alias `type`.
 *
 * @tparam T The type from which to extract `ValueType`.
 */
template< typename T >
struct get_value_type< T, std::enable_if_t< !internal::has_value_type< T >::value && internal::has_ValueType< T >::value > >
{
  using type = typename T::ValueType;
};

} // namespace geos

#endif /* TYPES_HELPERS_HPP */
