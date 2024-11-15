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
 * @file TypeHelpers.hpp
 *
 */

#ifndef TYPES_HELPERS_HPP
#define TYPES_HELPERS_HPP

#include <type_traits>

namespace geos
{

namespace internal
{

template< typename T, typename = void >
struct has_value_type : std::false_type {};

template< typename T >
struct has_value_type< T, std::void_t< typename T::value_type > > : std::true_type {};

template< typename T, typename = void >
struct has_ValueType : std::false_type {};

template< typename T >
struct has_ValueType< T, std::void_t< typename T::ValueType > > : std::true_type {};

}   // namespace internal

template< typename T, typename Enable = void >
struct get_value_type
{
  static_assert( sizeof(T) == 0, "CONTAINER_TYPE must define either value_type or ValueType." );
};

template< typename T >
struct get_value_type< T, std::enable_if_t< internal::has_value_type< T >::value > >
{
  using type = typename T::value_type;
};

template< typename T >
struct get_value_type< T, std::enable_if_t< !internal::has_value_type< T >::value && internal::has_ValueType< T >::value > >
{
  using type = typename T::ValueType;
};

} // namespace geos

#endif /* TYPES_HELPERS_HPP */
