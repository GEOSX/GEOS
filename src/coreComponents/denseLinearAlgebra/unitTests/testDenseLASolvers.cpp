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

#include "denseLinearAlgebra/denseLASolvers.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "testUtils.hpp"

// TPL includes
#include <gtest/gtest.h>

#include <random>

namespace geos
{
namespace denseLinearAlgebra
{
namespace testing
{

constexpr real64 machinePrecision = 1.0e3 * LvArray::NumericLimits< real64 >::epsilon;

template< std::ptrdiff_t N >
class LinearSystem
{
public:
  real64 matrix[N][N]{};
  real64 solution[N]{};
  real64 rhs[N]{};
};

template< std::ptrdiff_t N >
class InvertibleLinearSystem : public LinearSystem< N >
{
public:
  using LinearSystem< N >::matrix;
  using LinearSystem< N >::solution;
  using LinearSystem< N >::rhs;

  InvertibleLinearSystem( short int seed )
  {
    std::mt19937 generator( seed );
    std::uniform_real_distribution< real64 > distribution( -10.0, 10.0 );
    std::uniform_real_distribution< real64 > perturbation( -20.0, 20.0 );
    for( ptrdiff_t i=0; i<N; ++i )
    {
      solution[i] = distribution( generator );
      for( ptrdiff_t j=0; j<N; ++j )
      {
        matrix[i][j] = distribution( generator );
      }
      // add a perturbation on the diagonal to avoid singular matrices.
      matrix[i][i] += perturbation( generator );
    }

    // Compute rhs as matrix * solution
    for( std::ptrdiff_t i = 0; i < N; ++i )
    {
      rhs[i] = 0.0;
      for( std::ptrdiff_t j = 0; j < N; ++j )
      {
        rhs[i] += matrix[i][j] * solution[j];
      }
    }
  }
};

template< std::ptrdiff_t N >
class SingularLinearSystem : public LinearSystem< N >
{
public:
  using LinearSystem< N >::matrix;
  using LinearSystem< N >::solution;
  using LinearSystem< N >::rhs;

  SingularLinearSystem( short int seed )
  {
    std::mt19937 generator( seed );
    std::uniform_real_distribution< real64 > distribution( -10.0, 10.0 );
    for( ptrdiff_t i=0; i<N-1; ++i )
    {
      solution[i] = 0.0;
      for( ptrdiff_t j=0; j < N; ++j )
      {
        matrix[i][j] = distribution( generator );
      }
    }
    for( ptrdiff_t j=0; j<N; ++j )
    {
      matrix[N-1][j] = 2.0*matrix[0][j];
    }

    // Compute rhs as matrix * solution
    for( std::ptrdiff_t i = 0; i < N; ++i )
    {
      rhs[i] = 0.0;
      for( std::ptrdiff_t j = 0; j < N; ++j )
      {
        rhs[i] += matrix[i][j] * solution[j];
      }
    }
  }
};

template< typename N >
class DenseLinearSolverTest : public ::testing::Test
{
public:

  static constexpr std::ptrdiff_t size = N::value;

  DenseLinearSolverTest() = default;
  ~DenseLinearSolverTest() override = default;

  void test_solve()
  {
    InvertibleLinearSystem< size > LS( 2024 );

    // GEOS_UNUSED_VAR(LS);

    forAll< parallelDevicePolicy<> >( 1, [=] GEOS_HOST_DEVICE ( int )
    {
      real64 sol[size]{};
      real64 matrix[size][size]{};
      real64 rhs[size]{};

      LvArray::tensorOps::copy< size, size >( matrix, LS.matrix );
      LvArray::tensorOps::copy< size >( rhs, LS.rhs );
      bool const success = denseLinearAlgebra::solve< size >( matrix, rhs, sol );

      PORTABLE_EXPECT_TRUE( success );

      for( std::ptrdiff_t i = 0; i < size; ++i )
      {
        PORTABLE_EXPECT_NEAR( sol[i],
                              LS.solution[i],
                              machinePrecision );
      }
    } );
  }
  void test_singularSystem()
  {
    SingularLinearSystem< size > LS( 2024 );

    forAll< parallelDevicePolicy<> >( 1, [=] GEOS_HOST_DEVICE ( int )
    {

      real64 sol[size]{};
      real64 matrix[size][size]{};
      real64 rhs[size]{};

      LvArray::tensorOps::copy< size, size >( matrix, LS.matrix );
      LvArray::tensorOps::copy< size >( rhs, LS.rhs );
      bool const success = denseLinearAlgebra::solve< size >( matrix, rhs, sol );

      PORTABLE_EXPECT_FALSE( success );
    } );
  }

};

using Dimensions = ::testing::Types< std::integral_constant< std::ptrdiff_t, 2 >,
                                     std::integral_constant< std::ptrdiff_t, 3 >,
                                     std::integral_constant< std::ptrdiff_t, 4 >,
                                     std::integral_constant< std::ptrdiff_t, 5 >,
                                     std::integral_constant< std::ptrdiff_t, 6 >,
                                     std::integral_constant< std::ptrdiff_t, 7 >,
                                     std::integral_constant< std::ptrdiff_t, 8 >,
                                     std::integral_constant< std::ptrdiff_t, 9 > >;


class NameGenerator
{
public:
  template< typename T >
  static std::string GetName( int )
  {
    if constexpr (T::value == 2) return "TwoByTwo";
    if constexpr (T::value == 3) return "ThreeByThree";
    if constexpr (T::value == 4) return "FourByFour";
    if constexpr (T::value == 5) return "FiveByFive";
    if constexpr (T::value == 6) return "SixBySix";
    if constexpr (T::value == 7) return "SevenBySeven";
    if constexpr (T::value == 8) return "EightByEight";
    if constexpr (T::value == 9) return "NineByNine";
  }
};

TYPED_TEST_SUITE( DenseLinearSolverTest, Dimensions, NameGenerator );

TYPED_TEST( DenseLinearSolverTest, testDenseLA )
{
  this->test_solve();
  this->test_singularSystem();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}

} // testing

} // denseLinearAlgebra

} // namespace geos
