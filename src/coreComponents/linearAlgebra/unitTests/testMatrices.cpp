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
 * @file testMatrices.cpp
 */

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"

#include <gtest/gtest.h>

using namespace geos;

template< typename LAI >
class MatrixTest : public ::testing::Test
{};

TYPED_TEST_SUITE_P( MatrixTest );

#if 0
TYPED_TEST_P( MatrixTest, MatrixFunctions )
{
  // Define aliases
  using Vector = typename TypeParam::ParallelVector;
  using Matrix = typename TypeParam::ParallelMatrix;

  // Get the MPI rank
  int numranks = MpiWrapper::commSize( MPI_COMM_GEOSX );
  int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  std::cout << "*** Rank: " << rank << std::endl;

  // Dummy vector and Matrix
  Matrix C;
  Matrix D;
  {
    // Test matrix-matrix product: C = A*B
    Matrix A;
    geos::testing::compute2DLaplaceOperator( MPI_COMM_GEOSX, 2 * numranks, A );
    Matrix B( A );

    A.multiply( B, C );
    A.leftMultiplyTranspose( A, D );
  }

  // Define some vectors, matrices
  Vector vec1, vec2, vec3;
  Matrix mat1, mat2, mat3, mat4;
  mat1.createWithLocalSize( 2, 2, MPI_COMM_GEOSX ); // 2*numranks x 2*numranks
  mat2.createWithGlobalSize( 2, 2, MPI_COMM_GEOSX ); // 2x2
  mat3.createWithLocalSize( 2, 3, 3, MPI_COMM_GEOSX ); // 2*numranks x 3*numranks
  mat4.createWithGlobalSize( 3, 4, 3, MPI_COMM_GEOSX ); // 3x4

  // Testing create, globalRows, globalCols
  EXPECT_EQ( mat1.numGlobalRows(), 2 * numranks );
  EXPECT_EQ( mat1.numGlobalCols(), 2 * numranks );
  EXPECT_EQ( mat2.numGlobalRows(), 2 );
  EXPECT_EQ( mat2.numGlobalCols(), 2 );
  EXPECT_EQ( mat3.numGlobalRows(), 2 * numranks );
  EXPECT_EQ( mat3.numGlobalCols(), 3 * numranks );
  EXPECT_EQ( mat4.numGlobalRows(), 3 );
  EXPECT_EQ( mat4.numGlobalCols(), 4 );

  // Testing add/set/insert element
  //  mat1.insert( 1, 0, .5 );
  //  mat1.close();
  //  mat1.set( 1, 0, 5 );
  //  mat1.close();
  mat1.open();
  mat1.add( 1, 0, 1 );
  mat1.add( 1, 0, 2 );
  mat1.close();

  // Testing add/set/insert c-style, getRowCopy
  globalIndex inds1[2] = { 0, 2 };
  globalIndex inds2[1] = { 0 };
  globalIndex inds3[3] = { 0, 1, 2 };
  real64 vals1[2] = { 5, 10 };
  real64 vals2[1] = { 1 };
  real64 vals3[3] = { .5, 1, 2 };

  globalIndex iRow = 1;

  if( ( mat4.ilower() <= iRow ) && ( iRow < mat4.iupper() ) )
  {
    mat4.insert( iRow, inds3, vals3, 3 );
    //    mat4.close();
    //    mat4.open();
    mat4.set( iRow, inds1, vals1, 2 );
    //    mat4.close();
    //    mat4.open();
    mat4.add( iRow, inds2, vals2, 1 );
    //    mat4.close();
  }
  mat4.close();

  array1d< real64 > colvals_CHECK( 3 );
  colvals_CHECK( 0 ) = 6;
  colvals_CHECK( 1 ) = 1;
  colvals_CHECK( 2 ) = 10;

  if( ( mat4.ilower() <= iRow ) && ( iRow < mat4.iupper() ) )
  {
    localIndex const rowLength = mat4.getGlobalRowLength( iRow );
    EXPECT_EQ( rowLength, colvals_CHECK.size() );
    array1d< real64 > colvals( rowLength );
    array1d< globalIndex > colinds( rowLength );
    mat4.getRowCopy( iRow, colinds, colvals );
    for( int i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( colvals( colinds[i] ), colvals_CHECK( i ) ); //HYPRE does not return sorted cols!
    }
  }
  // Testing add/set/insert array1d
  Matrix mat6;
  mat6.createWithGlobalSize( 4, 4, MPI_COMM_GEOSX );
  array1d< real64 > vals6( 3 );
  array1d< real64 > vals7( 3 );
  array1d< globalIndex > inds6( 3 );
  vals6[0] = 1;
  vals6[1] = .5;
  vals6[2] = -3;
  vals7[0] = 1;
  vals7[1] = 1;
  vals7[2] = 1;
  inds6[0] = 0;
  inds6[1] = 1;
  inds6[2] = 3;

  iRow = 0;
  if( ( mat6.ilower() <= iRow ) && ( iRow < mat6.iupper() ) )
  {
    mat6.insert( iRow, inds6, vals6 );
    //	  mat6.close();
    mat6.set( iRow, inds6, vals7 );
    //	  mat6.close();
    mat6.add( iRow, inds6, vals6 );
  }
  mat6.close();

  // Testing add/set/insert array2d
  Matrix mat7;
  mat7.createWithGlobalSize( 4, 4, MPI_COMM_GEOSX );
  array1d< globalIndex > rows( 2 );
  array1d< globalIndex > cols( 2 );
  array2d< real64 > vals8( 2, 2 );
  rows[0] = 0;
  rows[1] = 2;
  cols[0] = 1;
  cols[1] = 3;
  vals8[0][0] = 1;
  vals8[0][1] = 2;
  vals8[1][0] = 3;
  vals8[1][1] = 4;
  if( ( mat7.ilower() <= *std::min_element( rows.data(), rows.data() + rows.size() ) ) &&
      ( *std::max_element( rows.data(), rows.data() + rows.size() ) < mat7.iupper() ) )
  {
    mat7.insert( rows, cols, vals8 );
    //    mat7.close();
    mat7.add( rows, cols, vals8 );
    //    mat7.close();
  }
  mat7.close();

  // Testing set and zero
  mat7.open();
  mat7.set( 2 );
  mat7.close();

  mat7.open();
  mat1.zero();
  mat1.close();

  // Testing vector multiply, matrix multiply, MatrixMatrixMultiply
  vec1.createWithGlobalSize( 2, MPI_COMM_GEOSX );
  vec2.createWithGlobalSize( 2, MPI_COMM_GEOSX );
  vec1.set( 1 );
  vec1.close();
  globalIndex inds4[2] = { 0, 1 };
  real64 vals4[2] = { 1, 3 };
  real64 vals5[2] = { 2, 1 };

  if( ( mat2.ilower() <= 0 ) && ( 0 < mat2.iupper() ) )
  {
    mat2.insert( 0, inds4, vals4, 2 );
  }
  if( ( mat2.ilower() <= 1 ) && ( 1 < mat2.iupper() ) )
  {
    mat2.insert( 1, inds4, vals5, 2 );
  }
  mat2.close();

  mat2.multiply( vec1, vec2 );

  if( ( vec2.ilower() <= 0 ) && ( 0 < vec2.iupper() ) )
  {
    EXPECT_DOUBLE_EQ( vec2.get( 0 ), 4 );
  }
  if( ( vec2.ilower() <= 1 ) && ( 1 < vec2.iupper() ) )
  {
    EXPECT_DOUBLE_EQ( vec2.get( 1 ), 3 );
  }

  // Matrix-Matrix multiply
  Matrix mat2mat2;

  {
    Matrix mat22( mat2 );

    mat22.multiply( mat22, mat2mat2 );
  }
}
#endif

TYPED_TEST_P( MatrixTest, MatrixMatrixOperations )
{
  using Matrix = typename TypeParam::ParallelMatrix;

  globalIndex const n = 100;
  Matrix A;
  geos::testing::compute2DLaplaceOperator( MPI_COMM_GEOSX, n, A );

  Matrix A_squared;
  A.multiply( A, A_squared );

  real64 const a = A.normInf();
  real64 const b = A_squared.normInf();

  EXPECT_DOUBLE_EQ( a, 8.0 );
  EXPECT_DOUBLE_EQ( b, 64.0 );
}

TYPED_TEST_P( MatrixTest, RectangularMatrixOperations )
{
  using Matrix = typename TypeParam::ParallelMatrix;

  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );

  // Set a size that allows to run with arbitrary number of processes
  globalIndex const nRows = std::max( 100, mpiSize );
  globalIndex const nCols = 2 * nRows;

  Matrix A;
  A.createWithGlobalSize( nRows, nCols, 2, MPI_COMM_GEOSX );

  A.open();
  for( globalIndex i = A.ilower(); i < A.iupper(); ++i )
  {
    real64 const entry0 = static_cast< real64 >( i + 1 );
    real64 const entry1 = -entry0;
    A.insert( i, 2 * i, entry0 );
    A.insert( i, 2 * i + 1, entry1 );
  }
  A.close();

  // Check on sizes
  EXPECT_EQ( A.numGlobalRows(), nRows );
  EXPECT_EQ( A.numGlobalCols(), nCols );

  // Check on norms
  real64 const a = A.norm1();
  real64 const b = A.normInf();
  real64 const c = A.normFrobenius();

  EXPECT_DOUBLE_EQ( a, static_cast< real64 >( nRows ) );
  EXPECT_DOUBLE_EQ( b, static_cast< real64 >( nCols ) );
  EXPECT_DOUBLE_EQ( c, std::sqrt( static_cast< real64 >( nRows * ( nRows + 1 ) * ( 2 * nRows + 1 ) ) / 3.0 ) );
}

REGISTER_TYPED_TEST_SUITE_P( MatrixTest,
                             MatrixMatrixOperations,
                             RectangularMatrixOperations );

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, MatrixTest, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, MatrixTest, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, MatrixTest, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
