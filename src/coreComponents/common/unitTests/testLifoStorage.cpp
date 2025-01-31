/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "mainInterface/initialization.hpp"
#define LIFO_DISABLE_CALIPER
#include "common/LifoStorage.hpp"
#include "LvArray/src/Array.hpp"
#include "LvArray/src/MallocBuffer.hpp"
#if defined(LVARRAY_USE_CHAI)
  #include "LvArray/src/ChaiBuffer.hpp"
#endif

#include <gtest/gtest.h>

#ifndef __CUDA_ARCH__
#define PORTABLE_EXPECT_EQ( L, R ) EXPECT_EQ( L, R )
#define PORTABLE_EXPECT_NEAR( L, R, EPSILON ) EXPECT_LE( math::abs( ( L ) -( R ) ), EPSILON ) << \
    STRINGIZE( L ) " = " << ( L ) << "\n" << STRINGIZE( R ) " = " << ( R );
#else
#define PORTABLE_EXPECT_EQ( L, R ) LVARRAY_ERROR_IF_NE( L, R )
#define PORTABLE_EXPECT_NEAR( L, R, EPSILON ) LVARRAY_ERROR_IF_GE_MSG( math::abs( ( L ) -( R ) ), EPSILON, \
                                                                       STRINGIZE( L ) " = " << ( L ) << "\n" << STRINGIZE( R ) " = " << ( R ) );
#endif

namespace geos
{

namespace local
{

template< typename >
struct RAJAHelper
{};

using serialPolicy = RAJA::seq_exec;

template<>
struct RAJAHelper< serialPolicy >
{
  using ReducePolicy = RAJA::seq_reduce;
  using AtomicPolicy = RAJA::seq_atomic;
  static constexpr LvArray::MemorySpace space = LvArray::MemorySpace::host;
};

#if defined(RAJA_ENABLE_OPENMP)

using parallelHostPolicy = RAJA::omp_parallel_for_exec;

template<>
struct RAJAHelper< parallelHostPolicy >
{
  using ReducePolicy = RAJA::omp_reduce;
  using AtomicPolicy = RAJA::omp_atomic;
  static constexpr LvArray::MemorySpace space = LvArray::MemorySpace::host;
};

#endif

#if defined(GEOS_USE_CUDA)

template< unsigned long THREADS_PER_BLOCK >
using devicePolicy = RAJA::cuda_exec< THREADS_PER_BLOCK >;

template< typename X, typename Y, typename C, size_t BLOCK_SIZE, bool ASYNC >
struct RAJAHelper< RAJA::policy::cuda::cuda_exec_explicit< X, Y, C, BLOCK_SIZE, ASYNC > >
{
  using ReducePolicy = RAJA::cuda_reduce;
  using AtomicPolicy = RAJA::cuda_atomic;
  static constexpr LvArray::MemorySpace space = LvArray::MemorySpace::cuda;
};

#endif

}

template< typename POLICY >
void testLifoStorage( int elemCnt, int numberOfElementsOnDevice, int numberOfElementsOnHost, int totalNumberOfBuffers )
{

  array1d< float > array( elemCnt );
  array.move( local::RAJAHelper< POLICY >::space );
  LifoStorage< float, localIndex > lifo( "lifo", array, numberOfElementsOnDevice, numberOfElementsOnHost, totalNumberOfBuffers );

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {

    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, j, elemCnt] GEOS_HOST_DEVICE ( int i ) { dataPointer[ i ] = j*elemCnt+i; } );
    lifo.push( array );
  }

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {
    lifo.pop( array );
    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, totalNumberOfBuffers, j, elemCnt] GEOS_HOST_DEVICE ( int i )
    {
      PORTABLE_EXPECT_EQ( dataPointer[ i ], (float)(totalNumberOfBuffers-j-1)*elemCnt+i );
    } );
  }
}

template< typename POLICY >
void testLifoStorageBig( int elemCnt, int numberOfElementsOnDevice, int numberOfElementsOnHost, int totalNumberOfBuffers )
{

  array1d< float > array( elemCnt );
  array.move( local::RAJAHelper< POLICY >::space );
  LifoStorage< float, localIndex > lifo( "lifo", array, numberOfElementsOnDevice, numberOfElementsOnHost, totalNumberOfBuffers );

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {

    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, j, elemCnt] GEOS_HOST_DEVICE ( int i ) { dataPointer[ i ] = j*elemCnt+i; } );
    lifo.push( array );
  }

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {
    lifo.pop( array );
    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, j, elemCnt] GEOS_HOST_DEVICE ( int i )
    {
      PORTABLE_EXPECT_EQ( dataPointer[ i ], (float)(10-j-1)*elemCnt+i );
    } );
  }
}

template< typename POLICY >
void testLifoStorageAsync( int elemCnt, int numberOfElementsOnDevice, int numberOfElementsOnHost, int totalNumberOfBuffers )
{
  array1d< float > array( elemCnt );
  array.move( local::RAJAHelper< POLICY >::space );
  LifoStorage< float, localIndex > lifo( "lifo", array, numberOfElementsOnDevice, numberOfElementsOnHost, totalNumberOfBuffers );

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {

    float * dataPointer = array.data();
    lifo.pushWait( );
    forAll< POLICY >( elemCnt, [dataPointer, j, elemCnt] GEOS_HOST_DEVICE ( int i ) { dataPointer[ i ] = j*elemCnt+i; } );
    lifo.pushAsync( array );
  }

  for( int j = 0; j < totalNumberOfBuffers; j++ )
  {
    lifo.popAsync( array );
    lifo.popWait( );
    float * dataPointer = array.data();
    forAll< POLICY >( elemCnt, [dataPointer, totalNumberOfBuffers, j, elemCnt] GEOS_HOST_DEVICE ( int i )
    {
      PORTABLE_EXPECT_EQ( dataPointer[ i ], (float)(totalNumberOfBuffers-j-1)*elemCnt+i );
    } );
  }
}



#ifdef GEOS_USE_CUDA
// running tests on GPUs
TEST( LifoStorageTest, LifoStorageBufferOnCUDA )
{
  testLifoStorage< local::devicePolicy< 32 > >( 10, 2, 3, 10 );
}

TEST( LifoStorageTest, LifoStorageBufferOnCUDAlarge )
{
  testLifoStorageBig< parallelDevicePolicy< > >( 1000000, 2, 3, 10 );
}

TEST( LifoStorageTest, LifoStorageBufferOnCUDAlargeAutoSizeHost )
{
  testLifoStorageBig< parallelDevicePolicy< > >( 1000000, 2, -80, 10 );
}

TEST( LifoStorageTest, LifoStorageBufferOnCUDAlargeAutoSizeDevice )
{
  testLifoStorageBig< parallelDevicePolicy< > >( 1000000, -80, 3, 10 );
}

TEST( LifoStorageTest, LifoStorageBufferOnCUDAlargeAutoSizeBoth )
{
  testLifoStorageBig< parallelDevicePolicy< > >( 1000000, -80, -80, 10 );
}

// The following test is disabled for now, as it produces a random assertion error
// that affects the current CI (see issue https://github.com/GEOS-DEV/GEOS/issues/3355).
// The error appears randomly on some configurations, and is probably related to the
// size-0 device buffer case that this test covers. This case should not be frequent
// in practise, but the issue should be resolved and the test reactivated as soon as
// a solution is found.
//TEST( LifoStorageTest, LifoStorageBufferOnCUDANoDeviceBuffer )
//{
//  testLifoStorage< local::devicePolicy< 32 > >( 10, 0, 3, 10 );
//}

TEST( LifoStorageTest, LifoStorageAsyncBufferOnCUDA )
{
  testLifoStorageAsync< local::devicePolicy< 32 > >( 10, 2, 3, 10 );
}

#else
// running tests on CPUs
TEST( LifoStorageTest, LifoStorageBufferOnHost )
{
  testLifoStorage< local::serialPolicy >( 10, 2, 3, 10 );
}

TEST( LifoStorageTest, LifoStorageBufferOnHostNoDeviceBuffer )
{
  testLifoStorage< local::serialPolicy >( 10, 0, 3, 10 );
}

TEST( LifoStorageTest, LifoStorageAsyncBufferOnHost )
{
  testLifoStorageAsync< local::serialPolicy >( 10, 2, 3, 10 );
}

#endif

}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
