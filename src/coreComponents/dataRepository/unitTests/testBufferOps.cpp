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

#include <gtest/gtest.h>

#include "dataRepository/BufferOps.hpp"

using namespace geos;
using namespace bufferOps;

class NonTrivialType
{
public:
  NonTrivialType() { }
  ~NonTrivialType() { }
};

TEST( testGeosxTraits, test_is_host_packable_object )
{
  static_assert( is_host_packable_object_v< int >, "Should be true." );
  static_assert( is_host_packable_object_v< double >, "Should be true." );
  static_assert( is_host_packable_object_v< R1Tensor >, "Should be true." );
  static_assert( is_host_packable_object_v< string >, "Should be true." );

  static_assert( !is_host_packable_object_v< void >, "Should be false." );
  static_assert( !is_host_packable_object_v< array1d< double > >, "Should be false." );
  static_assert( !is_host_packable_object_v< SortedArray< double > >, "Should be false." );
  static_assert( !is_host_packable_object_v< map< string, int > >, "Should be false." );
  static_assert( !is_host_packable_object_v< std::pair< string, int > >, "Should be false." );
}


TEST( testGeosxTraits, test_is_array_packable )
{
  static_assert( is_host_packable_array_v< array2d< real64, RAJA::PERM_IJ > >, "Should be true." );
  static_assert( is_host_packable_array_v< array2d< real64, RAJA::PERM_JI > >, "Should be true." );

  static_assert( !is_host_packable_array_v< int >, "Should be false." );
  static_assert( !is_host_packable_array_v< double >, "Should be false." );
  static_assert( !is_host_packable_array_v< void >, "Should be false." );
}


TEST( testGeosxTraits, test_is_host_packable_map )
{
  static_assert( is_host_packable_map_v< map< string, int > >, "Should be true." );
  static_assert( is_host_packable_map_v< map< string, array1d< int > > >, "Should be true." );
  static_assert( !is_host_packable_map_v< map< string, std::pair< int, int > > >, "Should be false" );
}

TEST( testGeosxTraits, test_is_host_packable_scalar )
{
  static_assert( is_host_packable_scalar_v< int >, "Should be true." );
  static_assert( is_host_packable_scalar_v< double >, "Should be true." );
  static_assert( !is_host_packable_scalar_v< R1Tensor >, "Should be false" );
  static_assert( !is_host_packable_scalar_v< string >, "Should be false" );
  static_assert( !is_host_packable_scalar_v< array1d< double > >, "Should be false." );
}

TEST( testGeosxTraits, test_is_host_packable )
{
  static_assert( !is_host_packable_v< NonTrivialType >, "Should be false" );
}

TEST( testBufferOps, test_pack_unpack_data_array )
{
  array1d< int > val1( 1 );
  val1[0] = 8675309;
  auto view = val1.toView();
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackData< false >( buff_head, view );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackData< true >( buff_head, view );
  // There is no UnpackData as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integer itself
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 8675309 );

}

TEST( testBufferOps, test_pack_unpack_data_arrayofarrays )
{
  ArrayOfArrays< int > val1( 2, 1 );
  val1.emplaceBack( 0, 8675309 );
  val1.emplaceBack( 1, 9035768 );
  ASSERT_EQ( val1.size(), 2 );
  ASSERT_EQ( val1.sizeOfArray( 0 ), 1 );
  ASSERT_EQ( val1.sizeOfArray( 1 ), 1 );
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackData< false >( buff_head, val1 );
  ASSERT_EQ( sz, 8 );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackData< true >( buff_head, val1 );
  // There is no UnpackData as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integers themselves
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 8675309 );
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[4] ), 9035768 );
}

TEST( testBufferOps, test_pack_unpack_data_arrayofsets )
{
  ArrayOfSets< int > val1( 2, 1 );
  val1.insertIntoSet( 0, 8675309 );
  val1.insertIntoSet( 1, 9035768 );
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackData< false >( buff_head, val1 );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackData< true >( buff_head, val1 );
  // There is no UnpackData as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integers themselves
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 8675309 );
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[4] ), 9035768 );
}

TEST( testBufferOps, test_pack_unpack_data_maps )
{
  map< int, float > val1;
  val1[8675309] = 0.5f; // 0.5 doesn't round per IEEE754
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackData< false >( buff_head, val1 );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackData< true >( buff_head, val1 );
  // There is no UnpackData as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integer and float themselves
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 8675309 );
  ASSERT_EQ( *reinterpret_cast< float * >( &buff[4] ), 0.5f ); // 0.5 doesn't round per IEEE754
}

TEST( testBufferOps, test_pack_unpack_pairs )
{
  std::pair< int, float > val1;
  val1.first = 8675309;
  val1.second = 0.5f; // 0.5 doesn't round per IEEE754
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackData< false >( buff_head, val1 );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackData< true >( buff_head, val1 );
  ASSERT_EQ( sz1, sz );
  // There is no UnpackData as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integer and float themselves
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 8675309 );
  ASSERT_EQ( *reinterpret_cast< float * >( &buff[4] ), 0.5f ); // 0.5 doesn't round per IEEE754
}

TEST( testBufferOps, test_pack_unpack_nontrivial_pointer )
{
  std::vector< array1d< int > > val1( 1 );
  val1[0].resize( 1 );
  val1[0][0] = 8675309;
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackPointerData< false >( buff_head, &val1[0], 1 );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackPointerData< true >( buff_head, &val1[0], 1 );
  ASSERT_EQ( sz1, sz );
  // There is no UnpackPointerData as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integer itself
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 8675309 );
}

TEST( testBufferOps, test_pack_unpack_data_by_index_array )
{
  array1d< int > val1( 2 );
  array1d< localIndex > idx( 1 );
  val1[1] = 8675309;
  idx[0] = 1;
  auto view = val1.toView();
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackDataByIndex< false >( buff_head, view, idx );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackDataByIndex< true >( buff_head, view, idx );
  ASSERT_EQ( sz1, sz );
  // There is no UnpackDataByIndex as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integer itself
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 8675309 );
}

TEST( testBufferOps, test_pack_unpack_data_by_index_arrayofarrays )
{
  ArrayOfArrays< int > val1( 2, 1 );
  val1.emplaceBack( 0, 8675309 );
  val1.emplaceBack( 1, 9035768 );
  ASSERT_EQ( val1.size(), 2 );
  ASSERT_EQ( val1.sizeOfArray( 0 ), 1 );
  ASSERT_EQ( val1.sizeOfArray( 1 ), 1 );
  array1d< localIndex > idx( 1 );
  idx[0] = 1;
  ASSERT_EQ( idx.size(), 1 );
  ASSERT_EQ( idx[0], 1 );
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackDataByIndex< false >( buff_head, val1, idx );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackDataByIndex< true >( buff_head, val1, idx );
  // There is no UnpackDataByIndex as we only use the Pack*Data functions for history output
  // Packing raw data should result in the buffer containing only the integer itself
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( *reinterpret_cast< int * >( &buff[0] ), 9035768 );
}

TEST( testBufferOps, test_pack_unpack_by_index_arrayofarrays )
{
  ArrayOfArrays< int > val1( 2, 1 );
  val1.emplaceBack( 0, 8675309 );
  val1.emplaceBack( 1, 9035768 );
  ASSERT_EQ( val1.size(), 2 );
  ASSERT_EQ( val1.sizeOfArray( 0 ), 1 );
  ASSERT_EQ( val1.sizeOfArray( 1 ), 1 );
  array1d< localIndex > idx( 1 );
  idx[0] = 1;
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackByIndex< false >( buff_head, val1, idx );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackByIndex< true >( buff_head, val1, idx );
  ArrayOfArrays< int > val2( 2, 1 );
  const buffer_unit_type * c_buff_head = &buff[0];
  localIndex sz2 = UnpackByIndex( c_buff_head, val2, idx, MPI_REPLACE );
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( sz2, sz );
  ASSERT_EQ( val2[1][0], 9035768 );
  c_buff_head = &buff[0];
  UnpackByIndex( c_buff_head, val2, idx, MPI_SUM );
  ASSERT_EQ( val2[1][0], 9035768*2 );
  val2[1][0] = 9035768 - 1;
  c_buff_head = &buff[0];
  UnpackByIndex( c_buff_head, val2, idx, MPI_MAX );
  ASSERT_EQ( val2[1][0], 9035768 );
}

TEST( testBufferOps, test_pack_unpack_arrayslice )
{
  array1d< int > val1( 1 );
  val1[0] = 8675309;
  auto view1 = val1.toView();
  auto slice1 = view1.toSlice();
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackArray< false >( buff_head, slice1, 1 );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackArray< true >( buff_head, slice1, 1 );
  array1d< int > val2( 1 );
  auto view2 = val2.toView();
  auto slice2 = view2.toSlice();
  const buffer_unit_type * c_buff_head = &buff[0];
  localIndex sz2 = UnpackArray( c_buff_head, slice2, 1, MPI_REPLACE );
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( sz2, sz );
  ASSERT_EQ( val2[0], 8675309 );
  c_buff_head = &buff[0];
  UnpackArray( c_buff_head, slice2, 1, MPI_SUM );
  ASSERT_EQ( slice2[0], 8675309*2 );
  slice2[0] = 8675309 - 1;
  c_buff_head = &buff[0];
  UnpackArray( c_buff_head, slice2, 1, MPI_MAX );
  ASSERT_EQ( slice2[0], 8675309 );
}

TEST( testBufferOps, test_pack_unpack_arrayslice_nontrivial )
{
  array1d< std::string > val1( 1 );
  val1[0] = std::string( "8675309" );
  auto view = val1.toView();
  auto slice1 = view.toSlice();
  buffer_unit_type * buff_head = nullptr;
  localIndex sz = PackArray< false >( buff_head, slice1, 1 );
  buffer_type buff( sz );
  buff_head = &buff[0];
  localIndex sz1 = PackArray< true >( buff_head, slice1, 1 );
  array1d< std::string > val2( 1 );
  auto view2 = val2.toView();
  auto slice2 = view2.toSlice();
  const buffer_unit_type * c_buff_head = &buff[0];
  localIndex sz2 = UnpackArray( c_buff_head, slice2, 1, MPI_REPLACE );
  ASSERT_EQ( sz1, sz );
  ASSERT_EQ( sz2, sz );
  ASSERT_EQ( slice2[0], std::string( "8675309" ) );
  c_buff_head = &buff[0];
  UnpackArray( c_buff_head, slice2, 1, MPI_SUM );
  ASSERT_EQ( slice2[0], std::string( "86753098675309" ) );
  // c_buff_head = &buff[0];
  //UnpackByIndex( c_buff_head, slice2, idx, MPI_MAX );
}

TEST( testBufferOps, test_packdata_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    buffer_unit_type * buff_head = nullptr;
    PackData< false >( buff_head, type );
  }, GEOS_FMT( "Trying to pack data type \\({}\\) but type is not packable.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}

TEST( testBufferOps, test_pack_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    buffer_unit_type * buff_head = nullptr;
    Pack< false >( buff_head, type );
  }, GEOS_FMT( "Trying to pack data type \\({}\\) but type is not packable.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}

TEST( testBufferOps, test_packdatabyindex_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    buffer_unit_type * buff_head = nullptr;
    PackDataByIndex< false >( buff_head, type, (int *)nullptr );
  }, GEOS_FMT( "Trying to pack data type \\({}\\) but type is not packable by index.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}

TEST( testBufferOps, test_packbyindex_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    buffer_unit_type * buff_head = nullptr;
    PackByIndex< false >( buff_head, type, (int *)nullptr );
  }, GEOS_FMT( "Trying to pack data type \\({}\\) but type is not packable by index.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}

TEST( testBufferOps, test_unpack_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    const buffer_unit_type * buff_head = nullptr;
    Unpack( buff_head, type, MPI_REPLACE );
  }, GEOS_FMT( "Trying to unpack data type \\({}\\) but type is not packable.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}

TEST( testBufferOps, test_unpackbyindex_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    const buffer_unit_type * buff_head = nullptr;
    UnpackByIndex( buff_head, type, (int *)nullptr, MPI_REPLACE );
  }, GEOS_FMT( "Trying to unpack data type \\({}\\) but type is not packable by index.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}

#ifdef GEOSX_USE_ARRAY_BOUNDS_CHECK

TEST( testBufferOps, test_pack_map_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    buffer_unit_type * buff_head = nullptr;
    Pack< false >( buff_head, type, (int *)nullptr );
  }, GEOS_FMT( "Trying to unpack data type \\({}\\) but type is not packable by index.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}

TEST( testBufferOps, test_unpack_map_failure )
{
  // redirect cout to cerr
  std::streambuf * orig_cout_streambuf = std::cout.rdbuf();
  std::cout.rdbuf( std::cerr.rdbuf());
  NonTrivialType type;
  ASSERT_DEATH( {
    buffer_unit_type * buff_head = nullptr;
    Unpack< false >( buff_head, type, (int *)nullptr );
  }, GEOS_FMT( "Trying to unpack data type \\({}\\) but type is not packable by index.", typeid(NonTrivialType).name() ) );
  std::cout.rdbuf( orig_cout_streambuf );
}


#endif
