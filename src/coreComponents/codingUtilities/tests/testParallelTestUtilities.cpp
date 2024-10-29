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

#include "mainInterface/initialization.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

#include <gtest/gtest.h>

TEST( UnitTestUtilities, skipInSerial )
{
  SKIP_TEST_IN_SERIAL( "This test should not be run in serial." );

  if( geos::MpiWrapper::commSize() == 1 )
  {
    GTEST_FAIL();
  }
  else
  {
    GTEST_SUCCEED();
  }
}

TEST( UnitTestUtilities, skipInParallel )
{
  SKIP_TEST_IN_PARALLEL( "This test should not be run in parallel." );

  if( geos::MpiWrapper::commSize() != 1 )
  {
    GTEST_FAIL();
  }
  else
  {
    GTEST_SUCCEED();
  }
}

TEST( UnitTestUtilities, expected )
{
  using namespace geos;
  using namespace geos::testing;

  if( MpiWrapper::commSize() == 1 )
  {
    ASSERT_EQ( 1, expected( 1, {} ) );
    ASSERT_EQ( 1, expected( 1, { 2 } ) );
    ASSERT_EQ( 1, expected( 1, { 2, 3 } ) );
  }
  else
  {
    ASSERT_EQ( MpiWrapper::commRank() == 0 ? 2 : 3, expected( 1, { 2, 3 } ) );
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
