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

#include "constitutive/fluid/singlefluid/ParticleFluid.hpp"
#include "constitutive/ExponentialRelation.hpp"
#include "constitutive/capillaryPressure/JFunctionCapillaryPressure.hpp"

#include <gtest/gtest.h>

using namespace geos;

// The `ENUM_STRING` implementation relies on consistency between the order of the `enum`,
// and the order of the `string` array provided. Since this consistency is not enforced, it can be corrupted anytime.
// This unit test aims at preventing from this implicit relationship to bring a bug.

TEST( ParticleFluidEnums, ParticleSettlingModel )
{
  using EnumType = constitutive::ParticleSettlingModel;

  ASSERT_EQ( "Stokes", toString( EnumType::Stokes ) );
  ASSERT_EQ( "Intermediate", toString( EnumType::Intermediate ) );
  ASSERT_EQ( "Turbulence", toString( EnumType::Turbulence ) );
}


TEST( ParticleFluidEnums, ExponentApproximationType )
{
  using EnumType = constitutive::ExponentApproximationType;

  ASSERT_EQ( "exponential", toString( EnumType::Full ) );
  ASSERT_EQ( "linear", toString( EnumType::Linear ) );
  ASSERT_EQ( "quadratic", toString( EnumType::Quadratic ) );
}


TEST( ParticleFluidEnums, PermeabilityDirection )
{
  using EnumType = constitutive::JFunctionCapillaryPressure::PermeabilityDirection;

  ASSERT_EQ( "XY", toString( EnumType::XY ) );
  ASSERT_EQ( "X", toString( EnumType::X ) );
  ASSERT_EQ( "Y", toString( EnumType::Y ) );
  ASSERT_EQ( "Z", toString( EnumType::Z ) );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
