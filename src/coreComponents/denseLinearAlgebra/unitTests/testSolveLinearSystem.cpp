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

/**
 * @file testSolveLinearSystem.cpp
 */

// Source includes
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

#include "gtest/gtest.h"

using namespace geos;

constexpr int MAX_SIZE = 20;
constexpr real64 machinePrecision = 1.0e2 * LvArray::NumericLimits< real64 >::epsilon;

// Test matrices
enum TestMatrixType
{
  LAPLACE,
  SAMPLING,
  GRCAR
};

template< TestMatrixType MATRIX_TYPE >
struct TestMatrix {};

template<>
struct TestMatrix< TestMatrixType::LAPLACE >
{
  template< int USD >
  static void create( arraySlice2d< real64, USD > const & A )
  {
    int const N = LvArray::integerConversion< int >( A.size( 0 ) );
    GEOS_ASSERT( 2 <= N );
    LvArray::forValuesInSlice( A, []( real64 & a ){ a = 0.0; } );
    for( int i = 0; i < N-1; i++ )
    {
      A( i+1, i ) = -1.0;
      A( i+1, i+1 ) = 2.0;
      A( i, i+1 ) = -1.0;
    }
    A( 0, 0 ) = 2.0;
  }
};

template<>
struct TestMatrix< TestMatrixType::SAMPLING >
{
  template< int USD >
  static void create( arraySlice2d< real64, USD > const & A )
  {
    int const N = LvArray::integerConversion< int >( A.size( 0 ) );
    GEOS_ASSERT( 2 <= N );

    real64 minajj = LvArray::NumericLimits< real64 >::max;
    for( int j = 1; j <= N; j++ )
    {
      real64 ajj = 0.0;
      for( int i = 1; i <= N; i++ )
      {
        if( i != j )
        {
          A( i-1, j-1 ) = static_cast< real64 >(i)/static_cast< real64 >(i-j);
          ajj += A( i-1, j-1 );
        }
      }
      A( j-1, j-1 ) = ajj;
      minajj = LvArray::math::min( ajj, minajj );
    }
    for( int j = 0; j < N; j++ )
    {
      A( j, j ) -= minajj;
    }
  }
};

template<>
struct TestMatrix< TestMatrixType::GRCAR >
{
  template< int USD >
  static void create( arraySlice2d< real64, USD > const & A )
  {
    int const N = LvArray::integerConversion< int >( A.size( 0 ) );
    GEOS_ASSERT( 4 < N );
    LvArray::forValuesInSlice( A, []( real64 & a ){ a = 0.0; } );
    for( int i = 0; i < N-1; i++ )
    {
      A( i+1, i ) = -1.0;
    }
    for( int c = 0; c < 4; c++ )
    {
      for( int i = 0; i < N-c; i++ )
      {
        A( i, i+c ) = 1.0;
      }
    }
  }
};

// Randomly reorder a matrix
template< int USD >
void random_permutation( arraySlice2d< real64, USD > const & A )
{
  int const N = LvArray::integerConversion< int >( A.size( 0 ) );
  GEOS_ASSERT( 2 <= N );

  StackArray< real64, 1, MAX_SIZE > ordering( N );
  BlasLapackLA::vectorRand( ordering.toSlice() );

  for( int j = 0; j < N; j++ )
  {
    int const i = static_cast< int >(ordering[j]*N);
    if( i != j )
    {
      for( int k = 0; k < N; ++k )
      {
        std::swap( A( i, k ), A( j, k ));
      }
    }
  }
}

// Local naive matrix-vector multiply
template< int USD >
void matrix_vector_multiply( arraySlice2d< real64 const, USD > const & A,
                             arraySlice1d< real64 const > const & x,
                             arraySlice1d< real64 > const & b )
{
  int const N = LvArray::integerConversion< int >( A.size( 0 ) );
  for( int i = 0; i < N; ++i )
  {
    real64 bi = 0.0;
    for( int j = 0; j < N; ++j )
    {
      bi += A( i, j )*x( j );
    }
    b( i ) = bi;
  }
}

// Local naive matrix-matrix multiply
template< int USD >
void matrix_matrix_multiply( arraySlice2d< real64 const, USD > const & A,
                             arraySlice2d< real64 const, USD > const & X,
                             arraySlice2d< real64, USD > const & B )
{
  int const K = LvArray::integerConversion< int >( A.size( 1 ) );
  int const M = LvArray::integerConversion< int >( X.size( 0 ) );
  int const N = LvArray::integerConversion< int >( X.size( 1 ) );
  for( int i = 0; i < M; ++i )
  {
    for( int j = 0; j < N; ++j )
    {
      real64 bij = 0.0;
      for( int k = 0; k < K; ++k )
      {
        bij += A( i, k )*X( k, j );
      }
      B( i, j ) = bij;
    }
  }
}

template< typename MatrixType >
struct ArrayType {};

template<>
struct ArrayType< Array< real64, 2, MatrixLayout::COL_MAJOR_PERM > >
{
  using type = Array< real64, 1 >;
};
template<>
struct ArrayType< Array< real64, 2, MatrixLayout::ROW_MAJOR_PERM > >
{
  using type = Array< real64, 1 >;
};
template<>
struct ArrayType< StackArray< real64, 2, MAX_SIZE *MAX_SIZE, MatrixLayout::COL_MAJOR_PERM > >
{
  using type = StackArray< real64, 1, MAX_SIZE >;
};
template<>
struct ArrayType< StackArray< real64, 2, MAX_SIZE *MAX_SIZE, MatrixLayout::ROW_MAJOR_PERM > >
{
  using type = StackArray< real64, 1, MAX_SIZE >;
};

template< typename MatrixType >
class LinearSolveFixture : public ::testing::Test
{
public:
  using VectorType = typename ArrayType< MatrixType >::type;
public:
  LinearSolveFixture() = default;
  ~LinearSolveFixture() override = default;

  template< TestMatrixType TEST_MATRIX, int N >
  void test_matrix_vector_solve() const
  {
    static_assert( 1 < N && N <= MAX_SIZE );

    MatrixType A( N, N );
    TestMatrix< TEST_MATRIX >::create( A.toSlice() );
    random_permutation( A.toSlice() );

    VectorType b ( N );
    VectorType x ( N );
    VectorType x0 ( N );

    // Save selected values of A to check later
    real64 const a00 = A( 0, 0 );
    real64 const a01 = A( 0, 1 );
    real64 const a10 = A( N-1, N-2 );
    real64 const a11 = A( N-1, N-1 );

    // Create actual solution
    LvArray::forValuesInSlice( x0.toSlice(), []( real64 & a ){ a = 1.0; } );

    // Multiply to get rhs
    matrix_vector_multiply( A.toSliceConst(), x0.toSliceConst(), b.toSlice() );

    // Save selected values of b to check later
    real64 const b0 = b( 0 );
    real64 const b1 = b( N-1 );

    // Solve
    BlasLapackLA::solveLinearSystem( A.toSliceConst(), b.toSliceConst(), x.toSlice() );

    // Check solution
    for( int i = 0; i < N; ++i )
    {
      EXPECT_NEAR( x( i ), x0( i ), machinePrecision );
    }

    // Check that we have not destroyed A
    EXPECT_NEAR( A( 0, 0 ), a00, machinePrecision );
    EXPECT_NEAR( A( 0, 1 ), a01, machinePrecision );
    EXPECT_NEAR( A( N-1, N-2 ), a10, machinePrecision );
    EXPECT_NEAR( A( N-1, N-1 ), a11, machinePrecision );

    // Check that we have not destroyed b
    EXPECT_NEAR( b( 0 ), b0, machinePrecision );
    EXPECT_NEAR( b( N-1 ), b1, machinePrecision );
  }

  template< TestMatrixType TEST_MATRIX, int N >
  void test_matrix_vector_solve_inplace( ) const
  {
    static_assert( 1 < N && N <= MAX_SIZE );

    MatrixType A( N, N );
    TestMatrix< TEST_MATRIX >::create( A.toSlice() );
    random_permutation( A.toSlice() );

    VectorType x ( N );
    VectorType x0 ( N );

    // Create actual solution
    LvArray::forValuesInSlice( x0.toSlice(), []( real64 & a ){ a = 1.0; } );

    // Multiply to get rhs
    matrix_vector_multiply( A.toSliceConst(), x0.toSliceConst(), x.toSlice() );

    // Solve
    BlasLapackLA::solveLinearSystem( A, x.toSlice() );

    // Check in place
    for( int i = 0; i < N; ++i )
    {
      EXPECT_NEAR( x( i ), x0( i ), machinePrecision );
    }
  }

  template< TestMatrixType TEST_MATRIX, int N, int M >
  void test_matrix_matrix_solve( ) const
  {
    static_assert( 1 < N && N <= MAX_SIZE );
    static_assert( 1 <= M && M <= MAX_SIZE );

    MatrixType A ( N, N );
    TestMatrix< TEST_MATRIX >::create( A.toSlice() );
    random_permutation( A.toSlice() );

    MatrixType B ( N, M );
    MatrixType X ( N, M );
    MatrixType X0 ( N, M );

    real64 const a00 = A( 0, 0 );
    real64 const a01 = A( 0, 1 );
    real64 const a10 = A( N-1, N-2 );
    real64 const a11 = A( N-1, N-1 );

    // Create actual solution
    // Populate matrix with random coefficients
    BlasLapackLA::matrixRand( X0,
                              BlasLapackLA::RandomNumberDistribution::UNIFORM_m1p1 );

    // Multiply to get rhs
    matrix_matrix_multiply( A.toSliceConst(), X0.toSliceConst(), B.toSlice() );

    // Save selected values of B to check later
    real64 const b00 = B( 0, 0 );
    real64 const b01 = B( N-1, M-1 );

    // Solve
    BlasLapackLA::solveLinearSystem( A.toSliceConst(), B.toSliceConst(), X.toSlice() );

    // Check
    for( int i = 0; i < N; ++i )
    {
      for( int j = 0; j < M; ++j )
      {
        EXPECT_NEAR( X( i, j ), X0( i, j ), machinePrecision );
      }
    }

    // Check that we have not destroyed A
    EXPECT_NEAR( A( 0, 0 ), a00, machinePrecision );
    EXPECT_NEAR( A( 0, 1 ), a01, machinePrecision );
    EXPECT_NEAR( A( N-1, N-2 ), a10, machinePrecision );
    EXPECT_NEAR( A( N-1, N-1 ), a11, machinePrecision );

    // Check that we have not destroyed b
    EXPECT_NEAR( B( 0, 0 ), b00, machinePrecision );
    EXPECT_NEAR( B( N-1, M-1 ), b01, machinePrecision );
  }

  template< TestMatrixType TEST_MATRIX, int N, int M >
  void test_matrix_matrix_solve_inplace( ) const
  {
    static_assert( 1 < N && N <= MAX_SIZE );
    static_assert( 1 <= M && M <= MAX_SIZE );

    MatrixType A ( N, N );
    TestMatrix< TEST_MATRIX >::create( A.toSlice() );
    random_permutation( A.toSlice() );

    MatrixType X ( N, M );
    MatrixType X0 ( N, M );

    // Create actual solution
    // Populate matrix with random coefficients
    BlasLapackLA::matrixRand( X0,
                              BlasLapackLA::RandomNumberDistribution::UNIFORM_m1p1 );

    // Multiply to get rhs
    matrix_matrix_multiply( A.toSliceConst(), X0.toSliceConst(), X.toSlice() );

    // Solve
    BlasLapackLA::solveLinearSystem( A.toSlice(), X.toSlice() );

    // Check
    for( int i = 0; i < N; ++i )
    {
      for( int j = 0; j < M; ++j )
      {
        EXPECT_NEAR( X( i, j ), X0( i, j ), machinePrecision );
      }
    }
  }
};

using LinearSolveTypes = ::testing::Types<
  Array< real64, 2, MatrixLayout::COL_MAJOR_PERM >,
  Array< real64, 2, MatrixLayout::ROW_MAJOR_PERM >,
  StackArray< real64, 2, MAX_SIZE *MAX_SIZE, MatrixLayout::COL_MAJOR_PERM >,
  StackArray< real64, 2, MAX_SIZE *MAX_SIZE, MatrixLayout::ROW_MAJOR_PERM > >;

class NameGenerator
{
public:
  template< typename T >
  static std::string GetName( int )
  {
    if constexpr (std::is_same_v< T, Array< real64, 2, MatrixLayout::COL_MAJOR_PERM > >) return "col-major-heap-array";
    if constexpr (std::is_same_v< T, Array< real64, 2, MatrixLayout::ROW_MAJOR_PERM > >) return "row-major-heap-array";
    if constexpr (std::is_same_v< T, StackArray< real64, 2, MAX_SIZE *MAX_SIZE, MatrixLayout::COL_MAJOR_PERM > >) return "col-major-stack-array";
    if constexpr (std::is_same_v< T, StackArray< real64, 2, MAX_SIZE *MAX_SIZE, MatrixLayout::ROW_MAJOR_PERM > >) return "row-major-stack-array";
  }
};

TYPED_TEST_SUITE( LinearSolveFixture, LinearSolveTypes, NameGenerator );

TYPED_TEST( LinearSolveFixture, matrix_vector_solve_laplace )
{
  this->template test_matrix_vector_solve< TestMatrixType::LAPLACE, 5 >();
  this->template test_matrix_vector_solve_inplace< TestMatrixType::LAPLACE, 5 >();
  this->template test_matrix_vector_solve< TestMatrixType::LAPLACE, 12 >();
  this->template test_matrix_vector_solve_inplace< TestMatrixType::LAPLACE, 12 >();
}

TYPED_TEST( LinearSolveFixture, matrix_vector_solve_grcar )
{
  this->template test_matrix_vector_solve< TestMatrixType::GRCAR, 6 >();
  this->template test_matrix_vector_solve_inplace< TestMatrixType::GRCAR, 6 >();
  this->template test_matrix_vector_solve< TestMatrixType::GRCAR, 10 >();
  this->template test_matrix_vector_solve_inplace< TestMatrixType::GRCAR, 10 >();
}

TYPED_TEST( LinearSolveFixture, matrix_vector_solve_sampling )
{
  this->template test_matrix_vector_solve< TestMatrixType::SAMPLING, 3 >();
  this->template test_matrix_vector_solve_inplace< TestMatrixType::SAMPLING, 3 >();
  this->template test_matrix_vector_solve< TestMatrixType::SAMPLING, 20 >();
  this->template test_matrix_vector_solve_inplace< TestMatrixType::SAMPLING, 20 >();
}

TYPED_TEST( LinearSolveFixture, matrix_matrix_solve_laplace )
{
  this->template test_matrix_matrix_solve< TestMatrixType::LAPLACE, 5, 1 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::LAPLACE, 5, 1 >();
  this->template test_matrix_matrix_solve< TestMatrixType::LAPLACE, 5, 3 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::LAPLACE, 5, 3 >();
  this->template test_matrix_matrix_solve< TestMatrixType::LAPLACE, 5, 5 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::LAPLACE, 5, 5 >();
  this->template test_matrix_matrix_solve< TestMatrixType::LAPLACE, 5, 12 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::LAPLACE, 5, 12 >();
}

TYPED_TEST( LinearSolveFixture, matrix_matrix_solve_grcar )
{
  this->template test_matrix_matrix_solve< TestMatrixType::GRCAR, 5, 1 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::GRCAR, 5, 1 >();
  this->template test_matrix_matrix_solve< TestMatrixType::GRCAR, 5, 3 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::GRCAR, 5, 3 >();
  this->template test_matrix_matrix_solve< TestMatrixType::GRCAR, 5, 5 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::GRCAR, 5, 5 >();
  this->template test_matrix_matrix_solve< TestMatrixType::GRCAR, 5, 12 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::GRCAR, 5, 12 >();
}

TYPED_TEST( LinearSolveFixture, matrix_matrix_solve_sampling )
{
  this->template test_matrix_matrix_solve< TestMatrixType::SAMPLING, 5, 1 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::SAMPLING, 5, 1 >();
  this->template test_matrix_matrix_solve< TestMatrixType::SAMPLING, 12, 3 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::SAMPLING, 5, 3 >();
  this->template test_matrix_matrix_solve< TestMatrixType::SAMPLING, 5, 5 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::SAMPLING, 5, 5 >();
  this->template test_matrix_matrix_solve< TestMatrixType::SAMPLING, 5, 12 >();
  this->template test_matrix_matrix_solve_inplace< TestMatrixType::SAMPLING, 5, 12 >();
}
