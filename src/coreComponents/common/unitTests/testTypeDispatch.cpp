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

#include "common/TypeDispatch.hpp"

#include <gtest/gtest.h>

using namespace geos;

static_assert( std::is_same< types::DimsRange< 3, 5 >, camp::as_list< camp::idx_seq< 3, 4, 5 > > >::value, "DimsRange failed" );
static_assert( std::is_same< types::DimsSingle< 4 >, camp::as_list< camp::idx_seq< 4 > > >::value, "DimsSingle failed" );
static_assert( std::is_same< types::DimsUpTo< 4 >, camp::as_list< camp::idx_seq< 1, 2, 3, 4 > > >::value, "DimsUpTo failed" );

static_assert( std::is_same< types::ArrayTypes< types::TypeList< int >, types::DimsRange< 2, 3 > >,
                             types::TypeList<
                               Array< int, 2, RAJA::PERM_IJ >,
                               Array< int, 2, RAJA::PERM_JI >,
                               Array< int, 3, RAJA::PERM_IJK >,
                               Array< int, 3, RAJA::PERM_IKJ >,
                               Array< int, 3, RAJA::PERM_JIK >,
                               Array< int, 3, RAJA::PERM_JKI >,
                               Array< int, 3, RAJA::PERM_KIJ >,
                               Array< int, 3, RAJA::PERM_KJI >
                               >
                             >::value, "ArrayTypes< <int>, <2, 3> > failed" );

static_assert( std::is_same< types::ArrayTypes< types::TypeList< int, double >, types::DimsSingle< 2 > >,
                             types::TypeList<
                               Array< int, 2, RAJA::PERM_IJ >,
                               Array< int, 2, RAJA::PERM_JI >,
                               Array< double, 2, RAJA::PERM_IJ >,
                               Array< double, 2, RAJA::PERM_JI >
                               >
                             >::value, "ArrayTypes< <int, double>, <2> > failed" );

static_assert( std::is_same< types::ArrayTypes< types::TypeList< int, double >, types::DimsUpTo< 2 > >,
                             types::TypeList<
                               Array< int, 1, RAJA::PERM_I >,
                               Array< int, 2, RAJA::PERM_IJ >,
                               Array< int, 2, RAJA::PERM_JI >,
                               Array< double, 1, RAJA::PERM_I >,
                               Array< double, 2, RAJA::PERM_IJ >,
                               Array< double, 2, RAJA::PERM_JI >
                               >
                             >::value, "ArrayTypes< <int, double>, <1, 2> > failed" );

namespace typeDispatchTest
{
class A
{
public:
  A() = default;
  virtual ~A() = default;
  virtual string getTypeName() const {return "A";}
};

class B : public A
{
public:
  virtual string getTypeName() const { return "B";}
};

class C : public A
{
public:
  virtual string getTypeName() const { return "C";}
};

class D : public A
{
public:
  virtual string getTypeName() const { return "D";}
};

}

template< typename TRUE_TYPES_LIST, size_t I = 0, typename ... Ts >
typename std::enable_if< I == sizeof...(Ts), void >::type
checkType( types::TypeList< Ts... > )
{
  return;
}

template< typename TRUE_TYPES_LIST, size_t I = 0, typename ... Ts >
typename std::enable_if< (I < sizeof...(Ts)), void>::type
checkType( types::TypeList< Ts... > list )
{
  EXPECT_TRUE( ( std::is_same< camp::at_t< decltype(list), camp::num< I > >,
                               camp::at_t< TRUE_TYPES_LIST, camp::num< I > > >::value ) ) << "Dispatch matched the wrong type";
  // check next element
  checkType< TRUE_TYPES_LIST, I + 1 >( list );
}

template< typename TRUE_TYPES_LIST, typename ... Ts, typename ... Vs >
void testDispatch( types::TypeList< Ts... > const list,
                   Vs && ... objects )
{
  bool const result = types::dispatch( list, []( auto tupleOfTypes )
  {
    checkType< TRUE_TYPES_LIST >( tupleOfTypes );
  }, std::forward< Vs >( objects ) ... );
  EXPECT_TRUE( result ) << "Dispatch failed to match the type";
}

TEST( testDispatch, DispatchSimpleTypesSingles )
{
  using Types = types::ListofTypeList< types::TypeList< int, float > >;

  int a = 0;
  testDispatch< types::TypeList< int > >( Types{}, a );

  float b = 0;
  testDispatch< types::TypeList< float > >( Types{}, b );

}


TEST( testDispatch, DispatchSimpleTypes )
{
  using Types = types::TypeList< types::TypeList< int, int >,
                                 types::TypeList< int, float > >;

  int a = 0;
  float b = 0.0;

  testDispatch< types::TypeList< int, float > >( Types{}, a, b );

}

TEST( testDispatch, DispatchVirtualTypePairs )
{

  using namespace typeDispatchTest;
  using Types = types::TypeList< types::TypeList< B, B >,
                                 types::TypeList< B, C >,
                                 types::TypeList< B, D > >;

  std::unique_ptr< A > b = std::make_unique< B >();
  std::unique_ptr< A > c = std::make_unique< C >();

  testDispatch< types::TypeList< B, C > >( Types{}, *b, *c );
}

TEST( testDispatch, DispatchVirtualTypeTriplets )
{

  using namespace typeDispatchTest;
  using Types = types::TypeList< types::TypeList< B, B, B >,
                                 types::TypeList< B, C, C >,
                                 types::TypeList< B, C, D > >;

  std::unique_ptr< A > b = std::make_unique< B >();
  std::unique_ptr< A > c = std::make_unique< C >();
  std::unique_ptr< A > d = std::make_unique< D >();

  testDispatch< types::TypeList< B, C, D > >( Types{}, *b, *c, *d );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
