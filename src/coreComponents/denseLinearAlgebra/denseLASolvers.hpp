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
 * @file denseLASolvers.hpp
 */
#ifndef GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_
#define GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_

#include "common/DataTypes.hpp"
#include "denseLinearAlgebra/common/layouts.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/logger/Logger.hpp"

#include <complex>

namespace geos
{

namespace denseLinearAlgebra
{

namespace details
{

constexpr real64 singularMatrixTolerance = 1e2*LvArray::NumericLimits< real64 >::epsilon;

/**
 * @brief Solves a 2x2 linear system A * x = b.
 *
 * This function solves a linear system of the form A * x = b, where A is a 2x2 matrix,
 * b is a 2x1 vector, and x is the solution vector. The function checks the sizes
 * of the inputs to ensure they conform to the expected dimensions. It also checks that
 * the determinant of matrix A is not near zero to avoid solving a singular system.
 *
 * @tparam MATRIX_TYPE The type of the matrix A. Must support indexing with `A[i][j]`.
 * @tparam RHS_TYPE The type of the right-hand side vector b. Must support indexing with `b[i]`.
 * @tparam SOL_TYPE The type of the solution vector x. Must support indexing with `x[i]`.
 *
 * @param[in] A The 2x2 matrix representing the system of equations. Must have size 2x2.
 * @param[in] b The 2-element vector representing the right-hand side of the equation.
 * @param[out] x The 2-element vector that will store the solution to the system.
 * @return bool that sepcifies whether the solve succeeded (1) or not (0).
 */
template< typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveTwoByTwoSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  LvArray::tensorOps::internal::checkSizes< 2, 2 >( A );
  LvArray::tensorOps::internal::checkSizes< 2 >( b );
  LvArray::tensorOps::internal::checkSizes< 2 >( x );

  real64 const detA = LvArray::tensorOps::determinant< 2 >( A );

  if( LvArray::math::abs( detA ) < singularMatrixTolerance )
    return false;

  real64 const invA = 1.0 / detA;

  x[0] = ( A[1][1] * b[0] - A[0][1] * b[1] ) * invA;
  x[1] = ( A[0][0] * b[1] - A[1][0] * b[0] ) * invA;

  return true;
}

/**
 * @brief Solves a 3x3 linear system A * x = b.
 *
 * This function solves a linear system of the form A * x = b, where A is a 3x3 matrix,
 * b is a 3x1 vector, and x is the solution vector. The function checks the sizes
 * of the inputs to ensure they conform to the expected dimensions. It also checks that
 * the determinant of matrix A is not near zero to avoid solving a singular system.
 *
 * @tparam MATRIX_TYPE The type of the matrix A. Must support indexing with `A[i][j]`.
 * @tparam RHS_TYPE The type of the right-hand side vector b. Must support indexing with `b[i]`.
 * @tparam SOL_TYPE The type of the solution vector x. Must support indexing with `x[i]`.
 *
 * @param[in] A The 3x3 matrix representing the system of equations. Must have size 3x3.
 * @param[in] b The 3-element vector representing the right-hand side of the equation.
 * @param[out] x The 3-element vector that will store the solution to the system.
 * @return bool that sepcifies whether the solve succeeded (1) or not (0).
 */
template< typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveThreeByThreeSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  LvArray::tensorOps::internal::checkSizes< 3, 3 >( A );
  LvArray::tensorOps::internal::checkSizes< 3 >( b );
  LvArray::tensorOps::internal::checkSizes< 3 >( x );

  real64 const detA = LvArray::tensorOps::determinant< 3 >( A );

  if( LvArray::math::abs( detA ) < singularMatrixTolerance )
    return false;

  real64 const invA = 1.0 / detA;

  real64 const detX0 = b[0] * ( A[1][1] * A[2][2] - A[2][1] * A[1][2] ) -
                       b[1] * ( A[0][1] * A[2][2] - A[0][2] * A[2][1] ) +
                       b[2] * ( A[0][1] * A[1][2] - A[0][2] * A[1][1] );

  real64 const detX1 = A[0][0] * ( b[1] * A[2][2] - b[2] * A[1][2] ) -
                       A[1][0] * ( b[0] * A[2][2] - b[2] * A[0][2] ) +
                       A[2][0] * ( b[0] * A[1][2] - b[1] * A[0][2] );

  real64 const detX2 = A[0][0] * ( A[1][1] * b[2] - A[2][1] * b[1] ) -
                       A[1][0] * ( A[0][1] * b[2] - A[2][1] * b[0] ) +
                       A[2][0] * ( A[0][1] * b[1] - A[1][1] * b[0] );

  x[0] = detX0 * invA;
  x[1] = detX1 * invA;
  x[2] = detX2 * invA;

  return true;
}

/**
 * @brief Solves a linear system where the matrix is upper triangular using back substitution.
 *
 * This function solves the linear system `Ax = b`, where `A` is an upper triangular matrix, using
 * back substitution. The solution `x` is computed and stored in the provided output vector.
 *
 * @tparam N The size of the square matrix `A`.
 * @tparam MATRIX_TYPE The type of the matrix `A`.
 * @tparam RHS_TYPE The type of the right-hand side vector `b`.
 * @tparam SOL_TYPE The type of the solution vector `x`.
 * @param[in] A The upper triangular matrix representing the coefficients of the system.
 * @param[in] b The right-hand side vector. It is used to compute the solution.
 * @param[out] x The solution vector. The result of solving the system `Ax = b` using back substitution.
 */
template< std::ptrdiff_t N,
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
void solveUpperTriangularSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  for( std::ptrdiff_t i = N - 1; i >= 0; --i )
  {
    real64 sum = b[i];
    for( std::ptrdiff_t j = i + 1; j < N; ++j )
    {
      sum -= A[i][j] * x[j];
    }
    x[i] = sum / A[i][i];
  }
}

/**
 * @brief Solves a linear system using Gaussian elimination.
 *
 * This function performs Gaussian elimination on the given matrix `A` and right-hand side vector `b`.
 * It transforms the matrix `A` boolo an upper triangular matrix and then solves for the solution `x`
 * using back substitution.
 *
 * @tparam N The size of the square matrix `A`.
 * @tparam MATRIX_TYPE The type of the matrix `A`.
 * @tparam RHS_TYPE The type of the right-hand side vector `b`.
 * @tparam SOL_TYPE The type of the solution vector `x`.
 * @param[in,out] A The matrix to be transformed boolo an upper triangular matrix. Modified in place.
 * @param[in,out] b The right-hand side vector. Modified in place to reflect the transformed system.
 * @param[out] x The solution vector. The result of solving the system `Ax = b`.
 * @return bool that sepcifies whether the solve succeeded (1) or not (0).
 */
template< std::ptrdiff_t N,
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveGaussianElimination( MATRIX_TYPE & A, RHS_TYPE & b, SOL_TYPE && x )
{
  static_assert( N > 0, "N must be greater than 0." );
  LvArray::tensorOps::internal::checkSizes< N, N >( A );
  LvArray::tensorOps::internal::checkSizes< N >( b );
  LvArray::tensorOps::internal::checkSizes< N >( x );

  // Step 1: Transform  boolo an upper triangular matrix

  // 1.a. Find the pivot
  for( std::ptrdiff_t i = 0; i < N; ++i )
  {
    std::ptrdiff_t max_row = i;
    for( std::ptrdiff_t k = i + 1; k < N; ++k )
    {
      if( LvArray::math::abs( A[k][i] ) > LvArray::math::abs( A[max_row][i] ))
      {
        max_row = k;
      }
    }

    // 1.b. Swap rows
    for( std::ptrdiff_t k = i; k < N; ++k )
    {
      // std::swap( A[i][k], A[max_row][k] );
      real64 const temp = A[max_row][k];
      A[max_row][k] = A[i][k];
      A[i][k] = temp;
    }
    // std::swap( b[i], b[max_row] ); cannot be done on device
    real64 const temp = b[i];
    b[i] =  b[max_row];
    b[max_row] = temp;


    if( LvArray::math::abs( A[i][i] ) < singularMatrixTolerance )
      return false;

    // 1.c Eliminate entries below the pivot
    for( std::ptrdiff_t k = i + 1; k < N; ++k )
    {
      real64 const scaling = A[k][i] / A[i][i];
      for( std::ptrdiff_t j = i; j < N; ++j )
      {
        A[k][j] -= scaling * A[i][j];
      }
      b[k] -= scaling * b[i];
    }
  }

  // Step 2: Backward substitution
  solveUpperTriangularSystem< N >( A, b, std::forward< SOL_TYPE >( x ) );

  return true;
}

} // details namespace

/**
 * @brief Solves a linear system using the most appropriate method based on the size of the system.
 *
 * This function determines the appropriate method for solving a linear system `Ax = b` based on
 * the size of the matrix `A`. For 2x2 and 3x3 systems, specialized solvers are used. For larger systems,
 * Gaussian elimination is employed. The matrix and the rhs are modified by the function.
 *
 * @tparam N The size of the square matrix `A`.
 * @tparam MATRIX_TYPE The type of the matrix `A`.
 * @tparam RHS_TYPE The type of the right-hand side vector `b`.
 * @tparam SOL_TYPE The type of the solution vector `x`.
 * @tparam MODIFY_MATRIX boolean flag indicating whether the input matrix `A` and vector `b` should be modified.
 *                       If `1`, the matrix `A` and vector `b` are modified in place. If `0`, copies of
 *                       `A` and `b` are made, and the original data is left unchanged.
 * @param[in] A The matrix representing the coefficients of the system.
 * @param[in] b The right-hand side vector.
 * @param[out] x The solution vector. The result of solving the system `Ax = b`.
 * @return bool that sepcifies whether the solve succeeded (1) or not (0).
 */
template< std::ptrdiff_t N,
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE,
          bool MODIFY_MATRIX = 1 >
GEOS_HOST_DEVICE
inline
bool solve( MATRIX_TYPE & A, RHS_TYPE & b, SOL_TYPE && x )
{
  static_assert( N > 0, "N must be greater than 0." );
  static_assert( N < 10, "N cannot be larger than 9" );
  LvArray::tensorOps::internal::checkSizes< N, N >( A );
  LvArray::tensorOps::internal::checkSizes< N >( b );
  LvArray::tensorOps::internal::checkSizes< N >( x );

  if constexpr ( N == 2 )
  {
    return details::solveTwoByTwoSystem( A, b, std::forward< SOL_TYPE >( x ) );
  }
  else if constexpr ( N == 3 )
  {
    return details::solveThreeByThreeSystem( A, b, std::forward< SOL_TYPE >( x ) );
  }
  else
  {
    if constexpr ( MODIFY_MATRIX )
    {
      return details::solveGaussianElimination< N >( A, b, std::forward< SOL_TYPE >( x ) );
    }
    else
    {
      real64 A_copy[N][N]{};
      real64 b_copy[N]{};

      for( std::ptrdiff_t i=0; i < N; ++i )
      {
        b_copy[i] = b[i];
        for( std::ptrdiff_t j=0; j < N; ++j )
        {
          A_copy[i][j] = A[i][j];
        }
      }
      return details::solveGaussianElimination< N >( A_copy, b_copy, std::forward< SOL_TYPE >( x ) );
    }
  }
}

} // denseLinearAlgebra

} // geos


#endif /*GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_*/
