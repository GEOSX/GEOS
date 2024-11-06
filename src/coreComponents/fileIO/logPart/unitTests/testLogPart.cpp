/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "common/DataTypes.hpp"
#include "fileIO/logPart/LogPart.hpp"
#include <gtest/gtest.h>

using namespace geos;

TEST( testSection, sectionWithTitle )
{
  std::ostringstream oss;
  LogPart logPart( "section name" );
  logPart.begin( oss );
  std::cout << oss.str() << std::endl;
  EXPECT_EQ( oss.str(),
             "\n######################################################################\n"
             "##                           section name                           ##\n"
             "######################################################################\n\n"
             );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );
  std::cout << oss.str() << std::endl;
  EXPECT_EQ( oss.str(),
             "\n##                       End of section name                        ##\n"
             "######################################################################\n\n"
             );
  oss.clear();
}

TEST( testSection, sectionWithTitleAndOneDescription )
{
  std::ostringstream oss;
  LogPart logPart( "section name" );
  logPart.addDescription( "description name" );
  logPart.begin( oss );
  std::cout << oss.str() << std::endl;
  EXPECT_EQ( oss.str(),
             "\n######################################################################\n"
             "##                           section name                           ##\n"
             "######################################################################\n"
             "##  description name                                                ##\n\n"
             );
  oss.clear();
}

TEST( testSection, sectionWithSetWidth )
{
  std::ostringstream oss;
  LogPart logPart( "section name" );
  logPart.addDescription( "description name 1" );
  logPart.addDescription( "description name 2" );
  logPart.setMinWidth( 100 );
  logPart.begin( oss );
  std::cout << oss.str() << std::endl;
  EXPECT_EQ( oss.str(),
             "\n####################################################################################################\n"
             "##                                          section name                                          ##\n"
             "####################################################################################################\n"
             "##  description name 1                                                                            ##\n"
             "##  description name 2                                                                            ##\n\n"
             );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );
  std::cout << oss.str() << std::endl;
  EXPECT_EQ( oss.str(),
             "\n##                                      End of section name                                       ##\n"
             "####################################################################################################\n\n"
             );
  oss.clear();
}

TEST( testSection, sectionMultipleDescriptions )
{
  std::ostringstream oss;
  LogPart logPart( "TIMESTEP START" );
  logPart.addDescription( "Time", "00h08m20s out of 2d, 21h26m40s (0% completed)", "500 s / 250000 s" );
  logPart.addDescription( "Delta Time", "00h16m40s (1000 s)" );
  logPart.addDescription( "- Cycle: 1" );
  logPart.setMinWidth( 70 );
  logPart.begin( oss );
  std::cout << oss.str() << std::endl;
  EXPECT_EQ ( oss.str(),
              "\n######################################################################\n"
              "##                          TIMESTEP START                          ##\n"
              "######################################################################\n"
              "##  - Time: 00h08m20s out of 2d, 21h26m40s (0% completed)           ##\n"
              "##          500 s / 250000 s                                        ##\n"
              "##  - Delta Time: 00h16m40s (1000 s)                                ##\n"
              "##  - Cycle: 1                                                      ##\n\n"
              );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );
  std::cout << oss.str() << std::endl;
  EXPECT_EQ( oss.str(),
             "\n##                      End of TIMESTEP START                       ##\n"
             "######################################################################\n\n"
             );
  oss.clear();
}

TEST( testSection, sectionEndDescription )
{
  std::ostringstream oss;
  LogPart logPart( "TIMESTEP START" );
  logPart.addDescription( "description" );
  logPart.addEndDescription( "test end description" );
  logPart.setMinWidth( 70 );
  logPart.begin( oss );
  oss.clear();
  oss.str( "" );

  logPart.end( oss );

  std::cout << oss.str() << std::endl;

  EXPECT_EQ( oss.str(),
             "\n##  test end description                                            ##\n"
             "######################################################################\n"
             "##                      End of TIMESTEP START                       ##\n"
             "######################################################################\n\n"
             );
  oss.clear();
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();;
}
