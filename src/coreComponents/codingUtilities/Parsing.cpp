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
 * @file Parsing.hpp
 */

#include <fast_float.h>

#include <cstdlib>

namespace geos
{

namespace
{

template< typename T, std::enable_if_t< std::is_floating_point< T >::value > * = nullptr >
char const * parseValueImpl( char const * const first,
                             char const * const last,
                             T & value )
{
  using namespace fast_float;
  from_chars_result const res = from_chars( first, last, value, chars_format::general );
  return res.ec == std::errc() && !std::isinf( value ) ? res.ptr : first;
}

template< typename T, std::enable_if_t< std::is_integral< T >::value > * = nullptr >
char const * parseValueImpl( char const * const first,
                             char const * const last,
                             T & value )
{
  if( first == last )
  {
    return first;
  }

  errno = 0;
  char * tmp{};
  long long const v = std::strtoll( first, &tmp, 0 ); // strtol is not const-correct
  char const * const ptr = tmp;

  // Error handling from strtol is a bit quirky
  if( tmp == nullptr || std::distance( ptr, last ) <= 0 || ( ( v == LLONG_MIN || v == LLONG_MAX ) && errno == ERANGE ) )
  {
    return first;
  }
  value = static_cast< T >( v );
  return ptr;
}

}

template< typename T >
char const * parseValue( char const * const first,
                         char const * const last,
                         T & value )
{
  return parseValueImpl( first, last, value );
}

#define INST_PARSEVALUE( T ) \
  template char const * parseValue< T >( char const * const first, char const * const last, T & value )

INST_PARSEVALUE( float );
INST_PARSEVALUE( double );
INST_PARSEVALUE( short );
INST_PARSEVALUE( int );
INST_PARSEVALUE( long );
INST_PARSEVALUE( long long );
// Add other types as needed

#undef INST_PARSEVALUE

} // namespace geos
