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

/**
 * @file testGeometricObjects.cpp
 */

#include "mesh/simpleGeometricObjects/Cylinder.hpp"

#include <gtest/gtest.h>

namespace geos
{
using namespace dataRepository;

void setCylinderParameters( Cylinder & cylinder,
                            real64 const (&inputPoint1)[3],
                            real64 const (&inputPoint2)[3],
                            real64 const & inputRadius,
                            real64 const & inputInnerRadius )
{
  auto & point1 = cylinder.getReference< R1Tensor >( Cylinder::viewKeyStruct::point1String() );
  point1[0] = inputPoint1[0]; point1[1] = inputPoint1[1]; point1[2] = inputPoint1[2];
  auto & point2 = cylinder.getReference< R1Tensor >( Cylinder::viewKeyStruct::point2String() );
  point2[0] = inputPoint2[0]; point2[1] = inputPoint2[1]; point2[2] = inputPoint2[2];
  auto & radius = cylinder.getReference< real64 >( Cylinder::viewKeyStruct::radiusString() );
  radius = inputRadius;
  auto & innerRadius = cylinder.getReference< real64 >( Cylinder::viewKeyStruct::innerRadiusString() );
  innerRadius = inputInnerRadius;
}

TEST( GeometricObjectTests, Cylinder )
{
  real64 testCoord[3]{};
  integer const numSamples = 10;

  conduit::Node node;
  dataRepository::Group parent( "testGroup", node );

  // Step 1: checks for Cylinder 1

  Cylinder cylinder1( "cylinder1", &parent );

  real64 const inputPoint1Cylinder1[3] = { 0, 0, -2650.1 };
  real64 const inputPoint2Cylinder1[3] = { 0, 0, -2500.0 };
  real64 const inputRadiusCylinder1 = 750.0;
  real64 const inputInnerRadiusCylinder1 = 1.0;
  setCylinderParameters( cylinder1,
                         inputPoint1Cylinder1, inputPoint2Cylinder1,
                         inputRadiusCylinder1, inputInnerRadiusCylinder1 );

  real64 startX = -100;
  real64 endX   =  100;
  real64 startY = -100;
  real64 endY   =  100;
  real64 startZ = -2649;
  real64 endZ   = -2501;

  testCoord[1] = -100;
  testCoord[2] = -2480;
  for( integer i = 0; i < numSamples; ++i )
  {
    testCoord[0] = startX + i * (endX-startX)/numSamples;
    EXPECT_FALSE( cylinder1.isCoordInObject( testCoord ) );
  }

  testCoord[0] = -100;
  testCoord[2] = -2660;
  for( integer j = 0; j < numSamples; ++j )
  {
    testCoord[1] = startY + j * (endY-startY)/numSamples;
    EXPECT_FALSE( cylinder1.isCoordInObject( testCoord ) );
  }

  testCoord[0] = 0.24;
  testCoord[1] = 0.1;
  for( integer k = 0; k < numSamples; ++k )
  {
    testCoord[2] = startZ + k * (endZ-startZ)/numSamples;
    EXPECT_FALSE( cylinder1.isCoordInObject( testCoord ) );
  }

  testCoord[0] = 10;
  testCoord[1] = 70;
  for( integer k = 0; k < numSamples; ++k )
  {
    testCoord[2] = startZ + k * (endZ-startZ)/numSamples;
    EXPECT_TRUE( cylinder1.isCoordInObject( testCoord ) );
  }

  // Step 2: checks for Cylinder 2

  Cylinder cylinder2( "cylinder2", &parent );

  real64 const inputPoint1Cylinder2[3] = {  1, 1, 0 };
  real64 const inputPoint2Cylinder2[3] = { -1, -1, 0 };
  real64 const inputRadiusCylinder2 = 1.0;
  real64 const inputInnerRadiusCylinder2 = 0.1;
  setCylinderParameters( cylinder2,
                         inputPoint1Cylinder2, inputPoint2Cylinder2,
                         inputRadiusCylinder2, inputInnerRadiusCylinder2 );

  startX = -0.9;
  endX   =  0.9;
  startY = -0.9;
  endY   =  0.9;

  testCoord[1] = -1;
  testCoord[2] = -2;
  for( integer i = 0; i < numSamples; ++i )
  {
    testCoord[0] = startX + i * (endX-startX)/numSamples;
    EXPECT_FALSE( cylinder2.isCoordInObject( testCoord ) );
  }

  testCoord[0] = -1;
  testCoord[2] = 1.5;
  for( integer j = 0; j < numSamples; ++j )
  {
    testCoord[1] = startY + j * (endY-startY)/numSamples;
    EXPECT_FALSE( cylinder2.isCoordInObject( testCoord ) );
  }

  testCoord[2] = 0;
  for( integer j = 0; j < numSamples; ++j )
  {
    testCoord[1] = startY + j * (endY-startY)/numSamples;
    testCoord[0] = testCoord[1];
    EXPECT_FALSE( cylinder2.isCoordInObject( testCoord ) );
  }
  for( integer j = 0; j < numSamples; ++j )
  {
    testCoord[1] = startY + j * (endY-startY)/numSamples;
    testCoord[0] = testCoord[1]+0.3;
    EXPECT_TRUE( cylinder2.isCoordInObject( testCoord ) );
  }
}


} /* namespace geos */
