/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "common/TypeDispatch.hpp"

#include <gtest/gtest.h>

using namespace geos;

namespace typeDispatchTest
{
class ConstitutiveBase
{
public:
  ConstitutiveBase() = default;
  int m_a;
};

class Model1 : public ConstitutiveBase
{
public:
  Model1():
    ConstitutiveBase()
  {m_a = 1;};
};

class Model2 : public ConstitutiveBase
{
public:
  Model2():
    ConstitutiveBase()
  {m_a = 2; };
};

class FEMBase {};

class QuadP1 : public FEMBase
{
  public:
  static constexpr int num_nodes = 4;
};

class TriangleP1 : public FEMBase
{
  public:
  static constexpr int num_nodes = 3;
};

}

template< typename ... Ts >
void testDispatchSimple( types::TypeList< Ts... > const list )
{
  bool const result = types::dispatchCombinations( list, []( auto tupleOfTypes )
  {
    using FirstType  = camp::first< decltype(tupleOfTypes) >;
    using SecondType = camp::second< decltype(tupleOfTypes) >;
    std::cout << "test_dispatch_simple:  dispatched with " << typeid(FirstType).name() << ", " << typeid(SecondType).name() << std::endl;
  }, 0.0f, 0 );
  EXPECT_TRUE( result ) << "Dispatch failed to match the type";
}

template< typename ... Ts >
void testDispatchVirtual( types::TypeList< Ts... > const list,
                          typeDispatchTest::FEMBase & fem,
                          typeDispatchTest::ConstitutiveBase & constitutiveModel )
{
  bool const result = types::dispatchCombinations( list, [&]( auto tupleOfTypes )
  {
    using FemType   = camp::first< decltype(tupleOfTypes) >;
    using ModelType = camp::second< decltype(tupleOfTypes) >;
    FemType & femCasted = static_cast< FemType & >(fem);
    ModelType & modelCasted = static_cast< ModelType & >(constitutiveModel);
    std::cout << "test_dispatch_virtual: dispatched with fem.num_nodes: " << femCasted.num_nodes << ", model.m_a = " << modelCasted.m_a << std::endl;
  }, fem, constitutiveModel );
  EXPECT_TRUE( result ) << "Dispatch failed to match the type";
}


TEST( testDispatchSimple, DispatchSimpleTypes )
{
  using Types = types::TypeList< types::TypeList< int, int >,
                                 types::TypeList< float, int > >;

  testDispatchSimple( Types{} );
}

TEST( testDispatchVirtual, DispatchVirtualTypes )
{

  using namespace typeDispatchTest;
  using Types = types::TypeList< types::TypeList< QuadP1, Model1 >,
                                 types::TypeList< TriangleP1, Model2 > >;

  std::unique_ptr< FEMBase >         fem    = std::make_unique< QuadP1 >();
  std::unique_ptr< ConstitutiveBase > model = std::make_unique< Model1 >();

  testDispatchVirtual( Types{}, *fem, *model );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
