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

#include "codingUtilities/Parsing.hpp"

#include <gtest/gtest.h>

#include <random>
#include <sstream>

template< typename T, typename ENABLE = void >
struct distribution_helper
{};

template< typename T >
struct distribution_helper< T, std::enable_if_t< std::is_floating_point< T >::value > >
{
  static auto make()
  {
    return std::normal_distribution< T >();
  }
};

template< typename T >
struct distribution_helper< T, std::enable_if_t< std::is_integral< T >::value > >
{
  static auto make()
  {
    return std::uniform_int_distribution< T >( std::numeric_limits< T >::min(), std::numeric_limits< T >::max() );
  }
};

template< typename T >
void compare( T const lhs, T const rhs )
{
  EXPECT_EQ( lhs, rhs );
}

void compare( float const lhs, float const rhs )
{
  EXPECT_FLOAT_EQ( lhs, rhs );
}

void compare( double const lhs, double const rhs )
{
  EXPECT_DOUBLE_EQ( lhs, rhs );
}

template< typename T >
class ParsingTest : public ::testing::TestWithParam< T >
{
protected:

  static constexpr std::array< char, 4 > const separators = { ' ', ',', '\n', ';' };

  static bool issep( char const c )
  {
    return std::find( separators.begin(), separators.end(), c ) != separators.end();
  };

  std::vector< T > reference;
  std::string input;

  void SetUp() override
  {
    std::mt19937_64 gen( 2022 );
    auto dist = distribution_helper< T >::make();
    for( int i = 0; i < 256; ++i )
    {
      reference.push_back( dist( gen ) );
    }

    std::ostringstream os;
    os << std::scientific << std::setprecision( 16 );
    for( std::size_t i = 0; i < reference.size(); ++i )
    {
      os << reference[i] << separators[i % separators.size()];
    }
    input = os.str();
  }

  template< typename CONTAINER >
  void compareToReference( CONTAINER const & values ) const
  {
    ASSERT_EQ( values.size(), reference.size() );
    for( std::size_t i = 0; i < reference.size(); ++i )
    {
      compare( values[i], reference[i] );
    }
  }

  void testParseBuffer() const
  {
    std::vector< T > vec;
    char const * ptr = geos::parseBuffer( input.data(), input.data() + input.size(), vec, issep );
    EXPECT_EQ( ptr, input.data() + input.size() );
    compareToReference( vec );
  }

  void testParseBufferInvalid() const
  {
    std::vector< T > vec;
    auto const issep_invalid = []( char const c ){ return c == '|'; };
    char const * ptr = geos::parseBuffer( input.data(), input.data() + input.size(), vec, issep_invalid );
    EXPECT_NE( ptr, input.data() + input.size() );
  }

  void testParseFile() const
  {
    std::string const fname = GEOS_FMT( "testParsing_{}_input", typeid(T).name() );
    std::ofstream os( fname );
    os << input;
    os.close();

    std::vector< T > vec;
    geos::parseFile( fname, vec, issep );
    compareToReference( vec );

    std::remove( fname.c_str() );
  }

  void testParseFileInvalid() const
  {
    std::string const fname = GEOS_FMT( "testParsing_{}_input_invalid", typeid(T).name() );
    std::ofstream os( fname );
    os << input;
    os.close();

    std::vector< T > vec;
    auto const issep_invalid = []( char const c ){ return c == '|'; };
    EXPECT_THROW( geos::parseFile( fname, vec, issep_invalid ), std::runtime_error );

    std::remove( fname.c_str() );
  }
};

template< typename T >
std::array< char, 4 > const ParsingTest< T >::separators;

using Types = ::testing::Types<
  float,
  double,
  int,
  long,
  long long
  >;

TYPED_TEST_SUITE( ParsingTest, Types, );

TYPED_TEST( ParsingTest, parseBuffer )
{
  this->testParseBuffer();
}

TYPED_TEST( ParsingTest, parseBufferInvalid )
{
  this->testParseBufferInvalid();
}

TYPED_TEST( ParsingTest, parseFile )
{
  this->testParseFile();
}

TYPED_TEST( ParsingTest, parseFileInvalid )
{
  this->testParseFileInvalid();
}

int main( int argc, char * argv[] )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
