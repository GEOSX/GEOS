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
 * @file StringUtilities.cpp
 */

#include "StringUtilities.hpp"
#include "limits.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>

namespace geos
{
namespace stringutilities
{

string toLower( string const & input )
{
  string output;
  output.resize( input.size() );
  auto const toLowerCase = []( unsigned char c )
  { return std::tolower( c ); };
  std::transform( input.cbegin(), input.cend(), output.begin(), toLowerCase );
  return output;
}

string_view trim( string_view str,
                  string_view charsToRemove )
{
  std::size_t const first = str.find_first_not_of( charsToRemove );
  if( first != string::npos )
  {
    std::size_t const last = str.find_last_not_of( charsToRemove );
    return str.substr( first, ( last - first + 1 ) );
  }
  return {};
}
string_view trimSpaces( string_view str )
{
  return trim( str, " \f\n\r\t\v" );
}


string removeStringAndFollowingContent( string const & str,
                                        string const & strToRemove )
{
  string newStr = str;

  // check if the line contains the string to remove
  std::size_t const pos = newStr.find( strToRemove );

  if( pos != string::npos )
  {
    // remove the character and everything afterwards
    newStr = newStr.substr( 0, pos );
  }
  return newStr;
}

// put definition here so we can control the allowable values of T and
// modication of this function triggers a whole code recompile...which
// should be avoided.
template< typename T >
string toMetricPrefixString( T const & value )
{
  // These are the metric prefixes corrosponding to kilo, mega, giga...etc.
  char const prefixes[12] = { 'f', 'p', 'n', 'u', 'm', ' ', 'K', 'M', 'G', 'T', 'P', 'E'};
  string rval;

  int const power = floor( log10( std::abs( (double)value ) ) );
  int const a = floor( power / 3.0 );

  real64 const scaledValue = value * pow( 10.0, -a * 3 );

  // format the output of the value to 3 significant digits and append the
  // metric prefix.
  int const p = 2-std::abs( power - a * 3 );
  char temp[10];
  snprintf( temp, 8, "%5.*f %c", p, scaledValue, prefixes[a+5] );
  rval = temp;

  GEOS_ERROR_IF( rval.empty(),
                 GEOS_FMT( "The value of {} was not able to be converted with a metric prefix", value ) );


  return rval;
}
template string toMetricPrefixString( int const & );
template string toMetricPrefixString( long int const & );
template string toMetricPrefixString( long long int const & );
template string toMetricPrefixString( unsigned long int const & );
template string toMetricPrefixString( unsigned long long int const & );
template string toMetricPrefixString( float const & );
template string toMetricPrefixString( double const & );


}
}
