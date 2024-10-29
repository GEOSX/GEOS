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
 * @file Parsing.hpp
 */

#ifndef GEOS_CODINGUTILITIES_PARSING_HPP_
#define GEOS_CODINGUTILITIES_PARSING_HPP_

#include "common/DataTypes.hpp"

#include <fstream>

namespace geos
{

/**
 * @brief Parse an numeric value from a character sequence.
 * @tparam T type of value
 * @param first pointer to start of the character sequence
 * @param last pointer past-the-end of the character sequence
 * @param value the parsed value, or unchanged if parsing failed
 * @return pointer to the first character following the parsed value,
 *         or @p first if parsing did no succeed for any reason
 */
template< typename T >
char const * parseValue( char const * const first,
                         char const * const last,
                         T & value );

/**
 * @brief Parse a sequence of values into a container
 * @tparam CONTAINER type of container that supports emplace_back()
 * @tparam SEPFUNC type of function that determines separator chars
 * @param first pointer to start of buffer
 * @param last pointer past-the-end of buffer
 * @param target the container to fill
 * @param issep function that returns true if given character is a value separator
 * @return @p last if the entire buffer has been processed, or pointer to
 *         the start of the unprocessed part if a parsing error occurred
 */
template< typename CONTAINER, typename SEPFUNC >
char const * parseBuffer( char const * first,
                          char const * const last,
                          CONTAINER & target,
                          SEPFUNC issep = [] ( char const c ){ return std::isspace( c ); } )
{
  using T = typename CONTAINER::value_type;
  static_assert( std::is_arithmetic< T >::value && !std::is_same< T, char >::value,
                 "Only valid for arithmetic types except char" );

  while( true )
  {
    while( first != last && issep( *first ) )
    {
      ++first;
    }
    T value;
    char const * const ptr = parseValue( first, last, value );
    if( ptr == first )
    {
      break;
    }
    target.emplace_back( value );
    first = ptr;
  }
  return first;
}

/**
 * @brief Read a sequence of values from file into a container
 * @tparam CONTAINER type of container that supports emplace_back()
 * @tparam SEPFUNC type of function that determines separator chars
 * @param first pointer to start of buffer
 * @param last pointer past-the-end of buffer
 * @param target the container to fill
 * @param issep function that returns true if given character is a value separator
 * @return @p last if the entire buffer has been processed, or pointer to
 *         the start of the unprocessed part if a parsing error occurred
 * @throws std::runtime_error if file IO or parsing error occurred
 */
template< typename CONTAINER, typename SEPFUNC >
void parseFile( string const & filename,
                CONTAINER & target,
                SEPFUNC issep = [] ( char const c ){ return std::isspace( c ); } )
{
  // Read file and process in 16kb chunks
  std::size_t constexpr BUF_SIZE = 16384;
  char buf[BUF_SIZE+1];
  buf[BUF_SIZE] = '\0'; // safe padding for strtol

  std::ifstream inputStream( filename );
  while( inputStream )
  {
    inputStream.read( buf, BUF_SIZE );
    std::streamsize count = inputStream.gcount();

    // Rewind the stream a bit until we find a separator character.
    // This guarantees we pause/resume reading on a value boundary.
    if( !inputStream.eof() )
    {
      while( count > 0 && !issep( buf[count - 1] ) )
      {
        --count;
      }
      if( count > 0 )
      {
        inputStream.seekg( count - inputStream.gcount(), std::ios_base::cur );
      }
    }

    // Process the buffer
    char const * const end = buf + count;
    char const * const ptr = parseBuffer( buf, end, target, issep );

    // If buffer not exhausted, show an error with a snippet of unprocessed part
    if( ptr != end )
    {
      std::ptrdiff_t const left = std::distance( ptr, end );
      GEOS_THROW( GEOS_FMT( "Unable to parse value in file {} at position {}: {}...",
                            filename, static_cast< std::streamoff >( inputStream.tellg() ) - left,
                            string( ptr, std::min( left, std::ptrdiff_t{32} ) ) ),
                  std::runtime_error );
    }
  }

  GEOS_THROW_IF( inputStream.fail() && !inputStream.eof(),
                 GEOS_FMT( "Error while reading file {}: {}", filename, std::strerror( errno ) ),
                 std::runtime_error );
}

} // namespace geos

#endif //GEOS_CODINGUTILITIES_PARSING_HPP_
