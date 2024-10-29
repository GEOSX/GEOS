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

#include "fieldSpecification/TractionBoundaryCondition.hpp"

#include <gtest/gtest.h>

using namespace geos;

// The `ENUM_STRING` implementation relies on consistency between the order of the `enum`,
// and the order of the `string` array provided. Since this consistency is not enforced, it can be corrupted anytime.
// This unit test aims at preventing from this implicit relationship to bring a bug.

TEST( TractionBoundaryConditionEnums, TractionType )
{
  using EnumType = TractionBoundaryCondition::TractionType;

  ASSERT_EQ( "vector", toString( EnumType::vector ) );
  ASSERT_EQ( "normal", toString( EnumType::normal ) );
  ASSERT_EQ( "stress", toString( EnumType::stress ) );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
