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

#include "common/GeosxConfig.hpp"
#include <gtest/gtest.h>

TEST( Toolchain, NDEBUGfromTPls )
{
  /*
   * This test guards against spurious propagation of -DNDEBUG preprocessor flag
   * (HDF5 from the TPLs in our case), which has the bogus effect of disabling LvArray assertions:
   * we check that we are in RelWithDebInfo or Release build type when NDEBUG is defined and in Debug
   * configuration when NDEBUG is not defined: thus, LvArray assertions remain in Debug builds.
   */
  bool constexpr isDebug = std::string_view( GEOS_CMAKE_BUILD_TYPE ) == std::string_view( "Debug" );

#ifdef NDEBUG
  ASSERT_FALSE( isDebug );  // RelWithDebInfo or Release builds only
#else
  ASSERT_TRUE( isDebug );  // Debug builds only
#endif
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  int const result = RUN_ALL_TESTS();

  return result;
}
