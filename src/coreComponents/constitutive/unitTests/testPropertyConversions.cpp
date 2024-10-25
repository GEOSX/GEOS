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

#include "gtest/gtest.h"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PropertyConversions.hpp"

using namespace geos;
using namespace ::geos::constitutive;

TEST( PropertyConversionTests, testElasticConversions )
{
  real64 bulk, poisson, young, shear;
  real64 bulk_out, poisson_out, young_out, shear_out;

  // bulk and poisson
  bulk        = 5.4321;
  poisson     = 0.2345;
  young      = conversions::bulkModAndPoissonRatio::toYoungMod( bulk, poisson );
  shear       = conversions::bulkModAndPoissonRatio::toShearMod( bulk, poisson );

  EXPECT_DOUBLE_EQ( young, 3*bulk*(1-2*poisson) );
  EXPECT_DOUBLE_EQ( shear, 3*bulk*(1-2*poisson)/(2+2*poisson) );

  // bulk and shear
  poisson_out = conversions::bulkModAndShearMod::toPoissonRatio( bulk, shear );
  young_out = conversions::bulkModAndShearMod::toYoungMod( bulk, shear );

  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( young, young_out );

  // bulk and young
  shear_out = conversions::bulkModAndYoungMod::toShearMod( bulk, young );
  poisson_out = conversions::bulkModAndYoungMod::toPoissonRatio( bulk, young );

  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // young and shear
  bulk_out    = conversions::shearModAndYoungMod::toBulkMod( shear, young );
  poisson_out = conversions::shearModAndYoungMod::toPoissonRatio( shear, young );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( poisson, poisson_out );

  // young and poisson
  bulk_out = conversions::youngModAndPoissonRatio::toBulkMod( young, poisson );
  shear_out = conversions::youngModAndPoissonRatio::toShearMod( young, poisson );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // shear and poisson
  bulk_out = conversions::shearModAndPoissonRatio::toBulkMod( shear, poisson );
  young_out =conversions::shearModAndPoissonRatio::toYoungMod( shear, poisson );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( young, young_out );

  // first lame (bulk and shear input only)
  real64 lame = conversions::bulkModAndShearMod::toFirstLame( bulk, shear );

  EXPECT_DOUBLE_EQ( lame, bulk-2*shear/3 );
}
