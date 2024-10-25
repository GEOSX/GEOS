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

// Source includes
#include "../StringUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace stringutilities;

struct TokenizeTest
{
  string const strToTest;
  string const delims;
  bool const treatConsecutiveDelimAsOne;
  bool const preTrimStr;
  std::vector< string > const expected;

  TokenizeTest( string const & strToTest_,
                string const & delims_,
                bool treatConsecutiveDelimAsOne_,
                bool preTrimStr_,
                std::vector< string > const & expected_ ):
    strToTest( strToTest_ ),
    delims( delims_ ),
    treatConsecutiveDelimAsOne( treatConsecutiveDelimAsOne_ ),
    preTrimStr( preTrimStr_ ),
    expected( expected_ )
  {}

  friend std::ostream & operator<<( std::ostream & os, TokenizeTest const & test )
  {
    os << "TokenizeTests( \"" << test.strToTest
       << "\", \"" << test.delims << '\"'
       << ", delimsAsOne=" << test.treatConsecutiveDelimAsOne
       << ", preTrim=" << test.preTrimStr << " )";
    return os;
  }
};

struct TokenizeBySpacesTest
{
  string const strToTest;
  std::vector< string > const expected;

  TokenizeBySpacesTest( string const & strToTest_,
                        std::vector< string > const & expected_ ):
    strToTest( strToTest_ ),
    expected( expected_ )
  {}

  friend std::ostream & operator<<( std::ostream & os, TokenizeBySpacesTest const & test )
  {
    os << "TokenizeBySpacesTest( \"" << test.strToTest << "\" )";
    return os;
  }
};

TEST( testStringUtilities, tokenize )
{

  // Path tokenizing test
  {
    map< string, std::pair< std::vector< string >, std::vector< string > > > const
    entries =
    {
      { "//entry0//entry1//entry2", { { "", "entry0", "entry1", "entry2" },
          { "", "", "entry0", "", "entry1", "", "entry2" } } },
      { "entry0//entry1/entry2", { { "entry0", "entry1", "entry2" },
          { "entry0", "", "entry1", "entry2" } } }
    };

    for( auto const & entry : entries )
    {
      string const & key = entry.first;
      std::vector< string > const & values0 = entry.second.first;
      std::vector< string > const & values1 = entry.second.second;

      std::vector< string > tokens0 = stringutilities::tokenize( key, "/", true );
      std::vector< string > tokens1 = stringutilities::tokenize( key, "/", false );


      EXPECT_TRUE( tokens0==values0 );
      EXPECT_TRUE( tokens1==values1 );
    }
  }

  // Various strings tokenizing test
  {
    std::vector< TokenizeTest > const tokenizeTests = {
      TokenizeTest( string( "a|b" ), "|", false, false, { "a", "b" } ),
      TokenizeTest( string( "a|b" ), "|", true, false, { "a", "b" } ),
      TokenizeTest( string( "a|b" ), "|", false, true, { "a", "b" } ),
      TokenizeTest( string( "a|b" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "|a|b|" ), "|", false, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "|a|b|" ), "|", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "|a|b|" ), "|", false, true, { "a", "b" } ),
      TokenizeTest( string( "|a|b|" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "|a||b|" ), "|", false, false, { "", "a", "", "b", "" } ),
      TokenizeTest( string( "|a||b|" ), "|", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "|a||b|" ), "|", false, true, { "a", "", "b" } ),
      TokenizeTest( string( "|a||b|" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "a|||b" ), "|", false, false, { "a", "", "", "b" } ),
      TokenizeTest( string( "a|||b" ), "|", true, false, { "a", "b" } ),
      TokenizeTest( string( "a|||b" ), "|", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( "a|||b" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "||a|||b||" ), "|", false, false, { "", "", "a", "", "", "b", "", "" } ),
      TokenizeTest( string( "||a|||b||" ), "|", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "||a|||b||" ), "|", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( "||a|||b||" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "|" ), "|", false, false, { "", "" } ),
      TokenizeTest( string( "|" ), "|", true, false, { "", "" } ),
      TokenizeTest( string( "|" ), "|", false, true, { } ),
      TokenizeTest( string( "|" ), "|", true, true, { } ),

      TokenizeTest( string( "|||" ), "|", false, false, { "", "", "", "" } ),
      TokenizeTest( string( "|||" ), "|", true, false, { "", "" } ),
      TokenizeTest( string( "|||" ), "|", false, true, { } ),
      TokenizeTest( string( "|||" ), "|", true, true, { } ),

      TokenizeTest( string( "ab" ), "|", false, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "|", true, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "|", false, true, { "ab" } ),
      TokenizeTest( string( "ab" ), "|", true, true, { "ab" } ),

      TokenizeTest( string( "" ), "|", false, false, { } ),
      TokenizeTest( string( "" ), "|", true, false, { } ),
      TokenizeTest( string( "" ), "|", false, true, { } ),
      TokenizeTest( string( "" ), "|", true, true, { } ),

      TokenizeTest( string( "" ), "", false, false, { } ),
      TokenizeTest( string( "" ), "", true, false, { } ),
      TokenizeTest( string( "" ), "", false, true, { } ),
      TokenizeTest( string( "" ), "", true, true, { } ),

      TokenizeTest( string( "ab" ), "", false, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "", true, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "", false, true, { "ab" } ),
      TokenizeTest( string( "ab" ), "", true, true, { "ab" } ),

      TokenizeTest( string( "| a|  b |" ), "| ", false, false, { "", "", "a", "", "", "b", "", "" } ),
      TokenizeTest( string( "| a|  b |" ), "| ", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "| a|  b |" ), "| ", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( "| a|  b |" ), "| ", true, true, { "a", "b" } ),

      TokenizeTest( string( " |a ||b| " ), "| ", false, false, { "", "", "a", "", "", "b", "", "" } ),
      TokenizeTest( string( " |a ||b| " ), "| ", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( " |a ||b| " ), "| ", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( " |a ||b| " ), "| ", true, true, { "a", "b" } ),
    };

    for( TokenizeTest test : tokenizeTests )
    {
      std::vector< string > const r = stringutilities::tokenize( test.strToTest,
                                                                 test.delims,
                                                                 test.treatConsecutiveDelimAsOne,
                                                                 test.preTrimStr );

      if( r != test.expected )
      {
        string const result = "{ '" + stringutilities::join( r, "' , '" ) + "' }";
        string const expected = "{ '" + stringutilities::join( test.expected, "' , '" ) + "' }";
        FAIL() << test << " failed: " << result <<" results instead of " << expected << ".";
      }
    }

  }

  // Spaces tokenizing test
  {
    std::vector< TokenizeBySpacesTest > const tokenizeBSTests = {
      TokenizeBySpacesTest( string( "a b" ), { "a", "b" } ),
      TokenizeBySpacesTest( string( " a b " ), { "a", "b" } ),
      TokenizeBySpacesTest( string( "  a  b  " ), { "a", "b" } ),
      TokenizeBySpacesTest( string( "\ta \tb\n\r" ), { "a", "b" } ),
      TokenizeBySpacesTest( string( " " ), { } ),
      TokenizeBySpacesTest( string( "   " ), { } ),
      TokenizeBySpacesTest( string( "ab" ), { "ab" } ),
      TokenizeBySpacesTest( string( "\t ab \n\r" ), { "ab" } ),
    };

    for( TokenizeBySpacesTest test : tokenizeBSTests )
    {
      std::vector< string > const r = stringutilities::tokenizeBySpaces( test.strToTest );

      if( r != test.expected )
      {
        string const result = "{ '" + stringutilities::concat( r, "' , '" ) + "' }";
        string const expected = "{ '" + stringutilities::concat( test.expected, "' , '" ) + "' }";
        FAIL() << test << " failed: " << result << " results instead of " << expected << ".";
      }
    }
  }

}


TEST( testStringUtilities, toMetricPrefixString )
{
  double const values[36] = { 1.234567890e-15,
                              1.234567890e-14,
                              1.234567890e-13,
                              1.234567890e-12,
                              1.234567890e-11,
                              1.234567890e-10,
                              1.234567890e-9,
                              1.234567890e-8,
                              1.234567890e-7,
                              1.234567890e-6,
                              1.234567890e-5,
                              1.234567890e-4,
                              1.234567890e-3,
                              1.234567890e-2,
                              1.234567890e-1,
                              1.234567890e0,
                              1.234567890e1,
                              1.234567890e2,
                              1.234567890e3,
                              1.234567890e4,
                              1.234567890e5,
                              1.234567890e6,
                              1.234567890e7,
                              1.234567890e8,
                              1.234567890e9,
                              1.234567890e10,
                              1.234567890e11,
                              1.234567890e12,
                              1.234567890e13,
                              1.234567890e14,
                              1.234567890e15,
                              1.234567890e16,
                              1.234567890e17,
                              1.234567890e18,
                              1.234567890e19,
                              1.234567890e20 };

  string const answer[36] = { " 1.23 f",
                              " 12.3 f",
                              "  123 f",
                              " 1.23 p",
                              " 12.3 p",
                              "  123 p",
                              " 1.23 n",
                              " 12.3 n",
                              "  123 n",
                              " 1.23 u",
                              " 12.3 u",
                              "  123 u",
                              " 1.23 m",
                              " 12.3 m",
                              "  123 m",
                              " 1.23  ",
                              " 12.3  ",
                              "  123  ",
                              " 1.23 K",
                              " 12.3 K",
                              "  123 K",
                              " 1.23 M",
                              " 12.3 M",
                              "  123 M",
                              " 1.23 G",
                              " 12.3 G",
                              "  123 G",
                              " 1.23 T",
                              " 12.3 T",
                              "  123 T",
                              " 1.23 P",
                              " 12.3 P",
                              "  123 P",
                              " 1.23 E",
                              " 12.3 E",
                              "  123 E" };


  for( int a=0; a<36; ++a )
  {
    std::string const result = toMetricPrefixString( values[a] );
    std::string const negResult = toMetricPrefixString( -values[a] );

    EXPECT_STRCASEEQ( result.c_str(), answer[a].c_str() );
    std::string negAnswer = answer[a];
    int const sign = negAnswer.find_first_not_of( ' ' );
    negAnswer[sign-1] = '-';
    EXPECT_STRCASEEQ( negResult.c_str(), negAnswer.c_str() );
  }
}

TEST( testStringUtilities, testStartsAndEndsWith )
{
  // classic use cases
  EXPECT_TRUE( stringutilities::startsWith( "Hello World", "Hello" ) );
  EXPECT_TRUE( stringutilities::endsWith( "Hello World", "World" ) );

  // inverted prefix & suffix
  EXPECT_FALSE( stringutilities::endsWith( "Hello World", "Hello" ) );
  EXPECT_FALSE( stringutilities::startsWith( "Hello World", "World" ) );
  EXPECT_FALSE( stringutilities::endsWith( "Hello World", "H" ) );
  EXPECT_FALSE( stringutilities::startsWith( "Hello World", "d" ) );

  // If prefix / suffix equals input string, then it must return true
  EXPECT_TRUE( stringutilities::startsWith( "Hello World", "Hello World" ) );
  EXPECT_TRUE( stringutilities::endsWith( "Hello World", "Hello World" ) );
  EXPECT_TRUE( stringutilities::startsWith( "H", "H" ) );
  EXPECT_TRUE( stringutilities::endsWith( "d", "d" ) );
  EXPECT_TRUE( stringutilities::startsWith( "", "" ) );
  EXPECT_TRUE( stringutilities::endsWith( "", "" ) );

  // Empty prefix / suffix are expected to work
  EXPECT_TRUE( stringutilities::startsWith( "Hello World", "" ) );
  EXPECT_TRUE( stringutilities::endsWith( "Hello World", "" ) );

  // the prefix / suffix is longer than the input string: return false (inverted parameters mistake?)
  EXPECT_FALSE( stringutilities::startsWith( "Hello", "Hello World" ) );
  EXPECT_FALSE( stringutilities::endsWith( "World", "Hello World" ) );

  // the prefix / suffix is longer than the input string: return false (inverted parameters mistake?)
  EXPECT_FALSE( stringutilities::startsWith( "Hello", "Hello World" ) );
  EXPECT_FALSE( stringutilities::endsWith( "World", "Hello World" ) );
  EXPECT_FALSE( stringutilities::startsWith( "", "Hello World" ) );
  EXPECT_FALSE( stringutilities::endsWith( "", "Hello World" ) );
}
