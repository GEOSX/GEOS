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
 * @file testVectors.cpp
 */

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "common/GEOS_RAJA_Interface.hpp"


#include <gtest/gtest.h>

using namespace geos;

/** ---------------------- Helpers ---------------------- **/

template< typename POLICY, typename VEC >
void createAndAssemble( localIndex const startSize, VEC & x )
{
  // Create a vector of local sizes:
  // - startSize      on rank 0;
  // - startSize + 1  on rank 1;
  // - etc.
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  localIndex const localSize = rank + startSize;
  globalIndex const rankOffset = rank * ( rank - 1 ) / 2 + rank * startSize;

  x.create( localSize, MPI_COMM_GEOSX );
  arrayView1d< real64 > const values = x.open();
  forAll< POLICY >( localSize, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    globalIndex const row = rankOffset + i;
    values[i] = std::pow( -1.0, row ) * ( 1.0 + row );
  } );
  x.close();
}

namespace ops
{

auto const identity = []( real64 const x ){ return x; };

struct multiply
{
  real64 factor;
  real64 operator()( real64 const x ) const { return factor * x; }
};

auto const reciprocal = []( real64 const x ){ return 1.0 / x; };

}

template< typename OP_X = decltype( ops::identity ),
          typename OP_Y = decltype( ops::identity ) >
void compareValues( arrayView1d< real64 const > const & x,
                    arrayView1d< real64 const > const & y,
                    bool exact = true,
                    OP_X const op_x = ops::identity,
                    OP_Y const op_y = ops::identity )
{
  EXPECT_EQ( x.size(), y.size() );
  x.move( hostMemorySpace, false );
  y.move( hostMemorySpace, false );
  for( localIndex i = 0; i < x.size(); ++i )
  {
    if( exact )
    {
      EXPECT_EQ( op_x( x[i] ), op_y( y[i] ) );
    }
    else
    {
      EXPECT_DOUBLE_EQ( op_x( x[i] ), op_y( y[i] ) );
    }
  }
}

template< typename VEC,
          typename OP_X = decltype( ops::identity ),
          typename OP_Y = decltype( ops::identity ) >
void compareVectors( VEC const & x,
                     VEC const & y,
                     bool exact = true,
                     OP_X op_x = ops::identity,
                     OP_Y op_y = ops::identity )
{
  EXPECT_EQ( y.created(), x.created() );
  if( !x.created() || !y.created() )
  {
    return;
  }

  EXPECT_TRUE( MpiWrapper::commCompare( y.comm(), x.comm() ) );
  EXPECT_EQ( y.localSize(), x.localSize() );
  EXPECT_EQ( y.globalSize(), x.globalSize() );

  EXPECT_EQ( y.closed(), x.closed() );
  if( !x.closed() || !y.closed() )
  {
    return;
  }

  compareValues( x.values(), y.values(), exact, std::forward< OP_X >( op_x ), std::forward< OP_Y >( op_y ) );
}

/** ---------------------- Tests ---------------------- **/

template< typename LAI >
class VectorTest : public ::testing::Test
{};

TYPED_TEST_SUITE_P( VectorTest );

TYPED_TEST_P( VectorTest, create )
{
  using Vector = typename TypeParam::ParallelVector;

  MPI_Comm const comm = MPI_COMM_GEOSX;
  int const rank = MpiWrapper::commRank( comm );
  int const nproc = MpiWrapper::commSize( comm );

  localIndex const localSize = rank + 1;
  globalIndex const globalSize = ( nproc + 1 ) * nproc / 2;
  globalIndex const rankOffset = ( rank + 1 ) * rank / 2;

  Vector x;
  x.create( localSize, comm );

  EXPECT_TRUE( x.ready() );
  EXPECT_TRUE( MpiWrapper::commCompare( x.comm(), comm ) );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );
  EXPECT_EQ( x.ilower(), rankOffset );
  EXPECT_EQ( x.iupper(), rankOffset + localSize );
}

TYPED_TEST_P( VectorTest, copyConstruction )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );
  Vector y( x );

  // Test that values are equal after copy construction
  compareVectors( x, y );

  // Test that vectors don't share data after copy construction
  real64 const factor = 0.5;
  x.scale( factor );

  compareVectors( x, y, true, ops::identity, ops::multiply{factor} );
}

TYPED_TEST_P( VectorTest, moveConstruction )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );
  localIndex const localSize = x.localSize();
  globalIndex const globalSize = x.globalSize();

  array1d< real64 > values( x.localSize() );
  values.template setValues< geos::parallelDevicePolicy<> >( x.values() );

  Vector y( std::move( x ) );

  EXPECT_TRUE( y.ready() );
  EXPECT_TRUE( MpiWrapper::commCompare( y.comm(), MPI_COMM_GEOSX ) );
  EXPECT_EQ( y.localSize(), localSize );
  EXPECT_EQ( y.globalSize(), globalSize );
  compareValues( y.values(), values );
}

TYPED_TEST_P( VectorTest, copy )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  Vector y;
  y.create( x.localSize(), x.comm() );
  y.copy( x );

  compareVectors( x, y );
}

TYPED_TEST_P( VectorTest, setAllValues )
{
  using Vector = typename TypeParam::ParallelVector;

  localIndex const localSize = 3;
  real64 const value = 1.23;

  Vector x;
  x.create( localSize, MPI_COMM_GEOSX );
  x.set( value );

  arrayView1d< real64 const > const values = x.values();
  values.move( hostMemorySpace, false );
  for( localIndex i = 0; i < localSize; ++i )
  {
    EXPECT_EQ( values[i], value );
  }
}

TYPED_TEST_P( VectorTest, zeroAllValues )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  x.zero();

  arrayView1d< real64 const > const values = x.values();
  values.move( hostMemorySpace, false );
  for( localIndex i = 0; i < values.size(); ++i )
  {
    EXPECT_EQ( values[i], 0.0 );
  }
}

TYPED_TEST_P( VectorTest, scaleValues )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  array1d< real64 > values( x.localSize() );
  values.template setValues< geos::parallelDevicePolicy<> >( x.values() );

  real64 const factor = 0.5;
  Vector y( x );
  y.scale( factor );

  compareVectors( x, y, true, ops::multiply{factor}, ops::identity );
}

TYPED_TEST_P( VectorTest, reciprocal )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  Vector y( x );
  y.reciprocal();

  compareVectors( x, y, true, ops::reciprocal, ops::identity );
}

TYPED_TEST_P( VectorTest, dotProduct )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  Vector y( x );
  y.reciprocal();

  real64 const dp = x.dot( y );
  EXPECT_DOUBLE_EQ( dp, x.globalSize() );
}

TYPED_TEST_P( VectorTest, axpy )
{
  using Vector = typename TypeParam::ParallelVector;

  real64 const alpha = 1.23;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  Vector y( x );
  x.axpy( alpha, y );

  compareVectors( x, y, false, ops::identity, ops::multiply{ 1.0 + alpha } );
}

TYPED_TEST_P( VectorTest, axpby )
{
  using Vector = typename TypeParam::ParallelVector;

  real64 const alpha = 2.0;
  real64 const beta = 3.0;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  Vector y( x );
  x.axpby( alpha, y, beta );

  compareVectors( x, y, false, ops::identity, ops::multiply{ alpha + beta } );
}

TYPED_TEST_P( VectorTest, norm1 )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  real64 const normTrue = ( x.globalSize() + 1 ) * x.globalSize() / 2;
  EXPECT_DOUBLE_EQ( x.norm1(), normTrue );
  x.scale( -1.0 );
  EXPECT_DOUBLE_EQ( x.norm1(), normTrue );
}

TYPED_TEST_P( VectorTest, norm2 )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  real64 const normTrue = std::sqrt( ( 2 * x.globalSize() + 1 ) * ( x.globalSize() + 1 ) * x.globalSize() / 6 );
  EXPECT_DOUBLE_EQ( x.norm2(), normTrue );
  x.scale( -1.0 );
  EXPECT_DOUBLE_EQ( x.norm2(), normTrue );
}

TYPED_TEST_P( VectorTest, normInf )
{
  using Vector = typename TypeParam::ParallelVector;

  Vector x;
  createAndAssemble< geos::parallelDevicePolicy<> >( 3, x );

  real64 const normTrue = x.globalSize();
  EXPECT_DOUBLE_EQ( x.normInf(), normTrue );
  x.scale( -1.0 );
  EXPECT_DOUBLE_EQ( x.normInf(), normTrue );
}

REGISTER_TYPED_TEST_SUITE_P( VectorTest,
                             create,
                             copyConstruction,
                             moveConstruction,
                             copy,
                             setAllValues,
                             zeroAllValues,
                             scaleValues,
                             reciprocal,
                             dotProduct,
                             axpy,
                             axpby,
                             norm1,
                             norm2,
                             normInf );

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, VectorTest, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, VectorTest, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, VectorTest, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
