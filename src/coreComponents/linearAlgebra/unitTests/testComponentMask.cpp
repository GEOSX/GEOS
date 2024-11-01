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
 * @file testComponentMask.cpp
 */

#include "linearAlgebra/utilities/ComponentMask.hpp"

#include "gtest/gtest.h"

#include <vector>
#include <algorithm>

using namespace geos;

template< int N >
void compare( ComponentMask< N > const & mask, std::vector< int > const & expected )
{
  EXPECT_EQ( mask.empty(), expected.empty() );
  EXPECT_EQ( mask.size(), expected.size() );
  EXPECT_EQ( std::vector< int >( mask.begin(), mask.end() ), expected );
}

std::vector< int > make_range( int lo, int const hi, int const step = 1 )
{
  std::vector< int > r( hi - lo );
  std::generate( r.begin(), r.end(), [&lo, step] { int const ret = lo; lo += step; return ret; } );
  return r;
}

std::vector< int > join_range( std::vector< int > lhs, std::vector< int > const & rhs )
{
  lhs.insert( lhs.end(), rhs.begin(), rhs.end() );
  return lhs;
}

template< typename CONSTANT >
class ComponentMaskTest : public ::testing::Test
{
protected:

  static constexpr int MAX_COMP = CONSTANT::value;
  using CompMask = ComponentMask< MAX_COMP >;
};


using MaskSizes = ::testing::Types<
  std::integral_constant< int, 0 >,
  std::integral_constant< int, 1 >,
  std::integral_constant< int, 4 >,
  std::integral_constant< int, 8 >,
  std::integral_constant< int, 12 >,
  std::integral_constant< int, 16 >,
  std::integral_constant< int, 24 >,
  std::integral_constant< int, 32 >,
  std::integral_constant< int, 48 >,
  std::integral_constant< int, 64 >
  >;

TYPED_TEST_SUITE( ComponentMaskTest, MaskSizes, );

TYPED_TEST( ComponentMaskTest, SimpleConstruct_ZeroComp_NoneIncluded )
{
  using Mask = typename TestFixture::CompMask;
  compare( Mask( 0, false ), {} );
}

TYPED_TEST( ComponentMaskTest, SimpleConstruct_ZeroComp_AllIncluded )
{
  using Mask = typename TestFixture::CompMask;
  compare( Mask( 0, true ), {} );
}

TYPED_TEST( ComponentMaskTest, SimpleConstruct_HalfComp_NoneIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N/2, false ), {} );
}

TYPED_TEST( ComponentMaskTest, SimpleConstruct_HalfComp_AllIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N/2, true ), make_range( 0, N/2 ) );
}

TYPED_TEST( ComponentMaskTest, SimpleConstruct_MaxComp_NoneIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N, false ), {} );
}

TYPED_TEST( ComponentMaskTest, SimpleConstruct_MaxComp_AllIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N, true ), make_range( 0, N ) );
}

TYPED_TEST( ComponentMaskTest, RangeConstruct_MaxComp_AllIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N, 0, N ), make_range( 0, N ) );
}

TYPED_TEST( ComponentMaskTest, RangeConstruct_MaxComp_FirstHalfIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N, 0, N/2 ), make_range( 0, N/2 ) );
}

TYPED_TEST( ComponentMaskTest, RangeConstruct_MaxComp_SecondHalfIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N, N/2, N ), make_range( N/2, N ) );
}

TYPED_TEST( ComponentMaskTest, RangeConstruct_MaxComp_MiddleHalfIncluded )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  compare( Mask( N, N/4, 3*N/4 ), make_range( N/4, 3*N/4 ) );
}

TYPED_TEST( ComponentMaskTest, ConversionConstruct_HalfComp_NoneIncluded )
{
  using Mask = typename TestFixture::CompMask;
  int constexpr N = TestFixture::MAX_COMP;
  ComponentMask< N/2 > src( N/2, false );
  compare( Mask( src ), {} );
}

TYPED_TEST( ComponentMaskTest, ConversionConstruct_HalfComp_AllIncluded )
{
  using Mask = typename TestFixture::CompMask;
  int constexpr N = TestFixture::MAX_COMP;
  ComponentMask< N/2 > src( N/2, true );
  compare( Mask( src ), make_range( 0, N/2 ) );
}

TYPED_TEST( ComponentMaskTest, Set_SingleComp )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  if( N > 0 )
  {
    Mask mask( N );
    mask.set( N / 3 );
    compare( mask, make_range( N / 3, N / 3 + 1 ) );
  }
}

TYPED_TEST( ComponentMaskTest, Set_EveryOtherComp )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  Mask mask( N );
  std::vector< int > expected;
  for( int i = 0; i < N; i += 2 )
  {
    mask.set( i );
    expected.push_back( i );
  }
  compare( mask, expected );
}

TYPED_TEST( ComponentMaskTest, Set_EveryComp )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  Mask mask( N );
  std::vector< int > expected;
  for( int i = 0; i < N; i += 1 )
  {
    mask.set( i );
    expected.push_back( i );
  }
  compare( mask, expected );
}

TYPED_TEST( ComponentMaskTest, Unset_SingleComp )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  if( N > 0 )
  {
    Mask mask( N, true );
    mask.unset( N / 3 );
    compare( mask, join_range( make_range( 0, N / 3 ), make_range( N / 3 + 1, N ) ) );
  }
}

TYPED_TEST( ComponentMaskTest, Unset_EveryOtherComp )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  Mask mask( N, true );
  std::vector< int > expected;
  for( int i = 0; i < N; i += 2 )
  {
    mask.unset( i );
    if( i + 1 < N )
    {
      expected.push_back( i + 1 );
    }
  }
  compare( mask, expected );
}

TYPED_TEST( ComponentMaskTest, Unset_EveryComp )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  Mask mask( N, true );
  for( int i = 0; i < N; i += 1 )
  {
    mask.unset( i );
  }
  compare( mask, {} );
}

TYPED_TEST( ComponentMaskTest, Invert_RemoveFirstHalf )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  Mask mask( N, 0, N/2 );
  mask.invert();
  compare( mask, make_range( N/2, N ) );
}

TYPED_TEST( ComponentMaskTest, Invert_RemoveMiddleHalf )
{
  using Mask = typename TestFixture::CompMask;
  constexpr int N = TestFixture::MAX_COMP;
  Mask mask( N, N/4, 3*N/4 );
  mask.invert();
  compare( mask, join_range( make_range( 0, N/4 ), make_range( 3*N/4, N ) ) );
}
