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
#include "common/TimingMacros.hpp"
#include "mainInterface/initialization.hpp"

#include <gtest/gtest.h>

using namespace geos;

class HDFParallelFileIOTest : public ::testing::TestWithParam<bool> {
};

TEST_P( HDFParallelFileIOTest, SingleValueHistory )
{
  GEOS_MARK_FUNCTION;
  string filename( GEOS_FMT( "single_value_parallel_{}", static_cast<int>( GetParam() ) ) );
  HistoryMetadata spec( "Time History", 1, std::type_index( typeid(real64)));

  int rank = MpiWrapper::commRank( );

  HDFHistoryIO io( filename, GetParam(), spec );
  io.init( true );
  real64 val = 0.0;
  for( localIndex tidx = 0; tidx < 100; ++tidx )
  {
    val += 0.5 * (rank+1);
    buffer_unit_type * buffer = io.getBufferHead( );
    memcpy( buffer, &val, sizeof(real64));
  }

  io.write( );
  io.compressInFile();
}


TEST_P( HDFParallelFileIOTest, DynamicHistory )
{
  GEOS_MARK_FUNCTION;
  string filename( GEOS_FMT( "dynamic_parallel_{}", static_cast<int>( GetParam() ) ) );
  HistoryMetadata singleItemSpec( "some_real", 1, std::type_index( typeid(real64)));

  int rnk = MpiWrapper::commRank( );

  HDFHistoryIO io( filename, GetParam(), singleItemSpec );
  io.init( true );
  real64 val = 0.0;
  for( localIndex tidx = 0; tidx < 4; ++tidx )
  {
    val += 0.5 * ( rnk + 1 );
    if( ( rnk + tidx ) % 2 == 0 )
    {
      io.updateCollectingCount( 1 );
      buffer_unit_type * buffer = io.getBufferHead( );
      memcpy( buffer, &val, sizeof(real64));
    }
    else
    {
      io.updateCollectingCount( 0 );
      io.getBufferHead( );
    }
  }
  io.write( );
  val = 0.0;
  for( localIndex tidx = 0; tidx < 8; ++tidx )
  {
    val += 0.5 * ( rnk + 1 );
    if( ( rnk + tidx ) % 2 == 1 )
    {
      io.updateCollectingCount( 1 );
      buffer_unit_type * buffer = io.getBufferHead( );
      memcpy( buffer, &val, sizeof(real64) );
    }
    else
    {
      io.updateCollectingCount( 0 );
      io.getBufferHead( );
    }
  }
  io.write( );
  io.compressInFile();
}

INSTANTIATE_TEST_SUITE_P(
    HDFParallelFileIOTests,
    HDFParallelFileIOTest,
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
