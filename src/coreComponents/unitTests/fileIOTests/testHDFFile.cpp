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

#include "fileIO/timeHistory/HDFHistoryIO.hpp"
#include "fileIO/timeHistory/HDFFile.hpp"
#include "mainInterface/initialization.hpp"
#include "dataRepository/BufferOpsDevice.hpp"

#include <gtest/gtest.h>

using namespace geos;

TEST( testHDFIO, HDFFile )
{
  GEOS_MARK_FUNCTION;
  HDFFile file( "empty", true, false, MPI_COMM_GEOS );
  hid_t file_id = H5Fcreate( "empty", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  H5Fclose( file_id );
}

TEST( testHDFIO, SingleValueHistory )
{
  string filename( "single_value" );
  HistoryMetadata spec( "Time History", 1, std::type_index( typeid(real64)));

  real64 time = 0.0;
  HDFHistoryIO io( filename, spec );
  io.init( true );
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    time += 0.333;
    buffer_unit_type * buffer = io.getBufferHead( );
    memcpy( buffer, &time, sizeof(real64));
  }
  io.write( );
}

TEST( testHDFIO, ArrayHistory )
{
  srand( time( NULL ));

  {
    string filename( "array1d_history" );
    Array< real64, 1 > arr( 4096 );
    real64 count = 0.0;
    forValuesInSlice( arr.toSlice(), [&count]( real64 & value )
    {
      value = count++;
    } );

    HistoryMetadata spec = getHistoryMetadata( "Array1d History", arr.toViewConst( ), 1 );
    HDFHistoryIO io( filename, spec );
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
    HDFHistoryIO io( filename, spec );
    io.init( true );

    buffer_unit_type * buffer = io.getBufferHead( );
    parallelDeviceEvents packEvents;
    bufferOps::PackDataDevice< true >( buffer, arr.toViewConst( ), packEvents );
    waitAllDeviceEvents( packEvents );

    io.write( );

    //read and check the data using hdf api
    // remove( filename.c_str() );
  }
}

TEST( testHDFIO, IdxArrayHistory )
{
  srand( time( NULL ));
  {
    string filename( "array1d_idx_history" );
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
    HDFHistoryIO io( filename, spec );
    io.init( true );

    buffer_unit_type * buffer = io.getBufferHead( );
    parallelDeviceEvents packEvents;
    bufferOps::PackDataByIndexDevice< true >( buffer, arr.toViewConst( ), idx.toViewConst( ), packEvents );
    waitAllDeviceEvents( packEvents );

    io.write( );
  }
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geos::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
