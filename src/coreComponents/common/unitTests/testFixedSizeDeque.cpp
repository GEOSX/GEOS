#include "mainInterface/initialization.hpp"
#include "common/fixedSizeDeque.hpp"
#include "LvArray/src/Array.hpp"
#include "LvArray/src/MallocBuffer.hpp"
#if defined(LVARRAY_USE_CHAI)
  #include "LvArray/src/ChaiBuffer.hpp"
#endif

#include <gtest/gtest.h>

using namespace geos;

TEST( FixedSizeDequeTest, ZeroSizedDeque )
{
  int maxArray = 0;
  int elemCnt = 10;
  camp::resources::Resource stream = { camp::resources::Host{} };
  fixedSizeDeque< float, int > empty_deque( maxArray, elemCnt, LvArray::MemorySpace::host, stream );
  EXPECT_EQ( true, empty_deque.empty());
  EXPECT_EQ( true, empty_deque.full());
}

TEST( FixedSizeDequeTest, emplace_back )
{
  int maxArray = 2;
  int elemCnt = 10;
  camp::resources::Resource stream = { camp::resources::Host{} };
  fixedSizeDeque< float, int > deque( maxArray, elemCnt, LvArray::MemorySpace::host, stream );
  EXPECT_EQ( true, deque.empty());
  EXPECT_EQ( false, deque.full());

  array1d< float > array( elemCnt );
  deque.emplace_back( array.toSliceConst() );
  EXPECT_EQ( false, deque.empty());
  EXPECT_EQ( false, deque.full());

  deque.emplace_back( array.toSliceConst() );
  EXPECT_EQ( false, deque.empty());
  EXPECT_EQ( true, deque.full());
}

TEST( FixedSizeDequeTest, emplace_front_and_back )
{
  int maxArray = 2;
  int elemCnt = 10;
  camp::resources::Resource stream = { camp::resources::Host{} };
  fixedSizeDeque< float, int > deque( maxArray, elemCnt, LvArray::MemorySpace::host, stream );
  EXPECT_EQ( true, deque.empty());
  EXPECT_EQ( false, deque.full());

  array1d< float > array( elemCnt );
  deque.emplace_front( array.toSliceConst() );
  EXPECT_EQ( false, deque.empty());
  EXPECT_EQ( false, deque.full());

  deque.emplace_back( array.toSliceConst() );
  EXPECT_EQ( false, deque.empty());
  EXPECT_EQ( true, deque.full());
}

TEST( FixedSizeDequeTest, emplace_and_pop )
{
  int maxArray = 2;
  int elemCnt = 10;
  camp::resources::Resource stream = { camp::resources::Host{} };
  fixedSizeDeque< float, int > deque( maxArray, elemCnt, LvArray::MemorySpace::host, stream );
  EXPECT_EQ( true, deque.empty());
  EXPECT_EQ( false, deque.full());

  array1d< float > array( elemCnt );
  for( int i = 0; i < elemCnt; i++ )
    array[i] = i;
  deque.emplace_back( array.toSliceConst() );
  EXPECT_EQ( false, deque.empty());
  EXPECT_EQ( false, deque.full());

  for( int i = 0; i < elemCnt; i++ )
    array[i] = 2*i;
  deque.emplace_front( array.toSliceConst() );
  EXPECT_EQ( false, deque.empty());
  EXPECT_EQ( true, deque.full());

  for( int i = 0; i < elemCnt; i++ )
    EXPECT_EQ( deque.front()[i], 2*i );
  for( int i = 0; i < elemCnt; i++ )
    EXPECT_EQ( deque.back()[i], i );


  deque.pop_back();
  EXPECT_EQ( false, deque.empty());
  EXPECT_EQ( false, deque.full());

  for( int i = 0; i < elemCnt; i++ )
    EXPECT_EQ( deque.front()[i], 2*i );
  for( int i = 0; i < elemCnt; i++ )
    EXPECT_EQ( deque.back()[i], 2*i );

  deque.pop_front();
  EXPECT_EQ( true, deque.empty());
  EXPECT_EQ( false, deque.full());
}

TEST( FixedSizeDequeTest, emplace_and_pop_front )
{
  int maxArray = 10;
  int elemCnt = 10;
  camp::resources::Resource stream = { camp::resources::Host{} };
  fixedSizeDeque< float, int > deque( maxArray, elemCnt, LvArray::MemorySpace::host, stream );
  EXPECT_EQ( true, deque.empty());
  EXPECT_EQ( false, deque.full());

  array1d< float > array( elemCnt );
  for( int j = 0; j < maxArray; j++ )
  {
    for( int i = 0; i < elemCnt; i++ )
      array[i] = i+maxArray*j;
    deque.emplace_front( array.toSliceConst() );
    if( j+1 < maxArray )
    {
      EXPECT_EQ( false, deque.empty());
      EXPECT_EQ( false, deque.full());
    }
    else
    {
      EXPECT_EQ( false, deque.empty());
      EXPECT_EQ( true, deque.full());
    }
  }
  for( int j = 0; j < maxArray; j++ )
  {
    for( int i = 0; i < elemCnt; i++ )
      EXPECT_EQ( deque.front()[i], i+(maxArray-j-1)*maxArray );

    deque.pop_front();
    if( j + 1 < maxArray )
    {
      EXPECT_EQ( false, deque.empty());
      EXPECT_EQ( false, deque.full());
    }
    else
    {
      EXPECT_EQ( true, deque.empty());
      EXPECT_EQ( false, deque.full());
    }
  }
}

#ifdef GEOSX_USE_CUDA
TEST( FixedSizeDequeTest, emplace_and_pop_front_cuda )
{
  int maxArray = 10;
  int elemCnt = 10;
  camp::resources::Resource stream = { camp::resources::Cuda{} };
  fixedSizeDeque< float, int > deque( maxArray, elemCnt, LvArray::MemorySpace::cuda, stream );
  EXPECT_EQ( true, deque.empty());
  EXPECT_EQ( false, deque.full());

  array1d< float > array( elemCnt );
  for( int j = 0; j < maxArray; j++ )
  {
    for( int i = 0; i < elemCnt; i++ )
      array[i] = i + j *maxArray;
    camp::resources::Event e = deque.emplace_front( array.toSliceConst() );
    stream.wait_for( &e );
    if( j+1 < maxArray )
    {
      EXPECT_EQ( false, deque.empty());
      EXPECT_EQ( false, deque.full());
    }
    else
    {
      EXPECT_EQ( false, deque.empty());
      EXPECT_EQ( true, deque.full());
    }
  }

  for( int j = 0; j < maxArray; j++ )
  {
    LvArray::memcpy( array.toSlice(), deque.front() );
    array.move( LvArray::MemorySpace::host, false );
    for( int i = 0; i < elemCnt; i++ )
      EXPECT_EQ( array[i], i+(maxArray-j-1)*maxArray );

    deque.pop_front();
    if( j + 1 < maxArray )
    {
      EXPECT_EQ( false, deque.empty());
      EXPECT_EQ( false, deque.full());
    }
    else
    {
      EXPECT_EQ( true, deque.empty());
      EXPECT_EQ( false, deque.full());
    }
  }
}
#endif

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
