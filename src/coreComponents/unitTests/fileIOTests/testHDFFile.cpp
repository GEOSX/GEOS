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

#include "fileIO/timeHistory/HDFHistoryIO.hpp"
#include "fileIO/timeHistory/HDFFile.hpp"
#include "mainInterface/initialization.hpp"
#include "dataRepository/BufferOpsDevice.hpp"

#include <gtest/gtest.h>

using namespace geos;

class HDFFileIOTest : public ::testing::TestWithParam<bool> {
};

TEST_P( HDFFileIOTest, HDFFile )
{
  GEOS_MARK_FUNCTION;
  bool useMPIO = GetParam();
  {
    HDFFile file( "empty", true, useMPIO, MPI_COMM_GEOSX );
  }
  hid_t file_id = 0;
  if( useMPIO )
  {
    file_id = H5Fcreate( "empty.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  }
  else
  {
    file_id = H5Fcreate( "empty.0.hdf5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  }
  herr_t err = H5Fclose( file_id );
  GEOS_ERROR_IF( err < 0, "H5Fclose failed." );
}

TEST_P( HDFFileIOTest, SingleValueHistory )
{
  string filename( GEOS_FMT( "single_value_{}", static_cast<int>( GetParam() ) ) );
  HistoryMetadata spec( "Time History", 1, std::type_index( typeid(real64)));

  real64 time = 0.0;
  HDFHistoryIO io( filename, GetParam(), spec );
  io.init( true );
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    time += 0.333;
    buffer_unit_type * buffer = io.getBufferHead( );
    memcpy( buffer, &time, sizeof(real64));
  }
  io.write( );
  io.compressInFile();
}

TEST_P( HDFFileIOTest, ArrayHistory )
{
  srand( time( NULL ));

  {
    string filename( GEOS_FMT( "array1d_history_{}", static_cast<int>( GetParam() ) ) ); 
    Array< real64, 1 > arr( 4096 );
    real64 count = 0.0;
    forValuesInSlice( arr.toSlice(), [&count]( real64 & value )
    {
      value = count++;
    } );

    HistoryMetadata spec = getHistoryMetadata( "Array1d History", arr.toViewConst( ), 1 );
    HDFHistoryIO io( filename, GetParam(), spec );
    io.init( true );

    buffer_unit_type * buffer = io.getBufferHead( );
    parallelDeviceEvents packEvents;
    bufferOps::PackDataDevice< true >( buffer, arr.toViewConst( ), packEvents );
    waitAllDeviceEvents( packEvents );

    io.write( );

    //read and check the data using hdf api
    // remove( filename.c_str() );
  }
  {
    string filename( "array2d_history" );
    Array< real64, 2 > arr( 1024, 4 );
    real64 count = 0.0;
    forValuesInSlice( arr.toSlice(), [&count]( real64 & value )
    {
      value = count++;
    } );

    HistoryMetadata spec = getHistoryMetadata( "Array2d History", arr.toViewConst( ), 4 );
    HDFHistoryIO io( filename, GetParam(), spec );
    io.init( true );

    buffer_unit_type * buffer = io.getBufferHead( );
    parallelDeviceEvents packEvents;
    bufferOps::PackDataDevice< true >( buffer, arr.toViewConst( ), packEvents );
    waitAllDeviceEvents( packEvents );

    io.write( );
    io.compressInFile();

    //read and check the data using hdf api
    // remove( filename.c_str() );
  }
}

TEST_P( HDFFileIOTest, IdxArrayHistory )
{
  srand( time( NULL ));
  {
    string filename( GEOS_FMT( "array1d_idx_history_{}", static_cast<int>( GetParam() ) ) );
    Array< localIndex, 1 > idx( 256 );
    Array< real64, 2 > arr( 1024, 4 );
    real64 count = 0.0;
    forValuesInSlice( arr.toSlice(), [&count]( real64 & value )
    {
      value = count++;
    } );
    forValuesInSlice( idx.toSlice(), []( localIndex & value )
    {
      value = rand() % 1024;
    } );

    HistoryMetadata spec = getHistoryMetadata( "Array1d Idx History", arr.toViewConst( ), 4, idx.size( ));
    HDFHistoryIO io( filename, GetParam(), spec );
    io.init( true );

    buffer_unit_type * buffer = io.getBufferHead( );
    parallelDeviceEvents packEvents;
    bufferOps::PackDataByIndexDevice< true >( buffer, arr.toViewConst( ), idx.toViewConst( ), packEvents );
    waitAllDeviceEvents( packEvents );

    io.write( );
    io.compressInFile();
  }
}

INSTANTIATE_TEST_SUITE_P(
    HDFFileIOTests,
    HDFFileIOTest,
    ::testing::Values( true, false )
);

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geos::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
