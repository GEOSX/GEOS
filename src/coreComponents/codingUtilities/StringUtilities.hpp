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
 * @file StringUtilities.hpp
 */

#ifndef GEOS_CODINGUTILITIES_STRINGUTILITIES_HPP_
#define GEOS_CODINGUTILITIES_STRINGUTILITIES_HPP_

#include "common/DataTypes.hpp"

#include <iomanip>
#include <sstream>

namespace geos
{
namespace stringutilities
{

/**
 * @brief Return a copy of the string in lower case.
 * @param input The input string which is not modified.
 * @return A new string instance.
 */
string toLower( string const & input );

/**
 * @brief Join strings or other printable objects with a delimiter.
 * @tparam IT   type of iterator into the range of objects to join
 * @tparam S    type of delimiter, usually char, char const * or string
 * @param first iterator to start of the range
 * @param last  iterator past-the-end of the range
 * @param delim delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename IT, typename S = char >
string join( IT first, IT last, S const & delim = S() )
{
  if( first == last )
  {
    return {};
  }
  std::ostringstream oss;
  oss << *first;
  while( ++first != last )
  {
    oss << delim << *first;
  }
  return oss.str();
}

/**
 * @brief Join strings or other printable objects with a delimiter.
 * @tparam CONTAINER type of container to join
 * @tparam S    type of delimiter, usually char, char const * or string
 * @param container container to join
 * @param delim delimiter used to glue together strings
 * @return a string containing input values concatenated with a delimiter
 */
template< typename CONTAINER, typename S = char >
string join( CONTAINER const & cont, S const & delim = S() )
{
  return join( std::begin( cont ), std::end( cont ), delim );
}

/**
 * @brief Concatenate variadic arguments into a string with a delimiter.
 * @tparam S type of delimiter (printable to std::ostringstream)
 * @tparam T type of first argument (printable to std::ostringstream)
 * @tparam Ts types of remaining arguments (printable to std::ostringstream)
 * @param delim delimiter
 * @param v first value
 * @param vs remaining values
 * @return string containing concatenated printed arguments
 */
template< typename S = char, typename T, typename ... Ts >
string concat( S const & delim, T const & v, Ts const & ... vs )
{
  std::ostringstream oss;
  oss << v;
  // Use array initializer and comma trick to get "fold expression" pre C++-17
  using expander = int[];
  (void) expander{ 0, ( void ( oss << delim << vs ), 0) ... };
  return oss.str();
}

/**
 * @brief Subdivide the string in substrings by the specified delimiters.
 * @tparam CONTAINER The templated class of the results container (std::vector by default).
 * @param str The string to subdivide.
 * @param delimiters String that contains the list of possible delimiters.
 * @param treatConsecutiveDelimAsOne If enabled, consecutive delimiters will be treated as one.
 *                                   If not enabled, consecutive delimiters will result in empty entries.
 * @param preTrimStr If enabled, delimiters at the borders of the string will be ignored.
 *                   If not enabled, those delimiters will result in in empty entries.
 * @return The container of the divided substrings.
 */
template< template< class ... > class CONTAINER = std::vector >
CONTAINER< string > tokenize( string const & str,
                              string const & delimiters,
                              bool const treatConsecutiveDelimAsOne = true,
                              bool const preTrimStr = false )
{
  CONTAINER< string > tokens;
  string::size_type tokenBegin, tokenEnd, strEnd;

  if( preTrimStr )
  {
    tokenBegin = str.find_first_not_of( delimiters );
    strEnd = str.find_last_not_of( delimiters ) + 1;
  }
  else
  {
    tokenBegin = 0;
    strEnd = str.size();
  }

  while( ( ( tokenEnd = str.find_first_of( delimiters, tokenBegin ) ) < strEnd ) && tokenBegin < strEnd )
  {
    tokens.emplace_back( str.substr( tokenBegin, tokenEnd - tokenBegin ) );
    tokenBegin = !treatConsecutiveDelimAsOne ? tokenEnd + 1 : str.find_first_not_of( delimiters, tokenEnd );
  }

  if( tokenBegin < strEnd )
  {
    tokens.emplace_back( str.substr( tokenBegin, strEnd-tokenBegin ));
  }
  else if( !preTrimStr && str.find_first_of( delimiters, strEnd - 1 ) != string::npos )
  {
    tokens.emplace_back( "" );
  }

  return tokens;
}

/**
 * @brief Subdivide the string in substrings by whitespaces separators (see std::isspace()).
 *        Do not create any empty substrings.
 * @tparam CONTAINER The templated class of the results container (std::vector by default).
 * @param str The string to subdivide.
 * @return CONTAINER< string > The list of the subdivided substrings (std::vector< string > for instance).
 */
template< template< class ... > class CONTAINER = std::vector >
CONTAINER< string > tokenizeBySpaces( string const & str )
{
  return tokenize< CONTAINER >( str, " \f\n\r\t\v", true, true );
}

/**
 * @brief Trim the string
 * @param[in] str the string to trim
 * @param[in] charsToRemove the list of characters to remove
 * @return the trimmed string
 */
string_view trim( string_view str,
                  string_view charsToRemove );

/**
 * @brief Trim the string so it does not starts nor ends with any whitespaces
 * @param[in] str the string to trim
 * @return the trimmed string
 */
string_view trimSpaces( string_view str );

/**
 * @brief Search for a string in the line, and return the line truncated before the string
 * @param[in] str the line to truncate
 * @param[in] strToRemove the string to search for in the line
 * @return the new (truncated) string
 */
string removeStringAndFollowingContent( string const & str,
                                        string const & strToRemove );

/**
 * @brief Take a string, and return a array1d with the cast values
 * @tparam T the type to which the string will be cast
 * @param[in] str the string to turn into an array1d
 * @return the array1d that stores the cast values
 */
template< typename T >
array1d< T > fromStringToArray( string const & str )
{
  array1d< T > v;
  T sub;

  std::istringstream iss( str );
  while( iss >> sub )
  {
    v.emplace_back( sub );
  }
  return v;
}

/**
 * @brief Take a numerical value and convert/scale it to a string with a metric
 *  prefix. i.e. Kilo, Mega, Giga, Tera, Peta, Exa
 *
 * @tparam T Type of the value to be converted
 * @param value The value to be converted
 * @return String containging the scaled value.
 */
template< typename T >
string toMetricPrefixString( T const & value );

/**
 * @brief Compute the length of a constant string at compile-time.
 */
// TODO c++17: this function is to remove in favor of std::string_view
constexpr size_t cstrlen( char const * const str )
{
  if( str )
  {
    char const * ptr = str;
    for(; *ptr != '\0'; ++ptr )
    {}
    return ptr - str;
  }
  else
  {
    return 0;
  }
}

/**
 * @return true if the string starts with the prefix.
 * @param str The string to compare
 * @param prefix A prefix we want to know if the string starts with.
 */
constexpr bool startsWith( std::string_view str, std::string_view prefix )
{
  return str.size() >= prefix.size() &&
         str.compare( 0, prefix.size(), prefix ) == 0;
}

/**
 * @return true if the string ends with the suffix.
 * @param str The string to compare
 * @param suffix A suffix we want to know if the string ends with.
 */
constexpr bool endsWith( std::string_view str, std::string_view suffix )
{
  return str.size() >= suffix.size() &&
         str.compare( str.size()-suffix.size(), suffix.size(), suffix ) == 0;
}

/**
 * @brief Overloading operator (<<) for std::optional<T>.
 *
 * This function displays the value contained in a std::optional<T> object if one exists.
 * Otherwise, it produces no output.
 *
 * @tparam T The type of the value contained std::optional.
 * @param os An output stream (for example, std::cout).
 * @param optValue std::optional<T> value to display.
 * @return The output stream
 */
template< typename T >
std::ostream & operator<<( std::ostream & os, std::optional< T > const & optValue )
{
  if( optValue )
  {
    os << optValue.value();
  }
  return os;
}

} // namespace stringutilities
} // namespace geos

#endif /* GEOS_CODINGUTILITIES_STRINGUTILITIES_HPP_ */
