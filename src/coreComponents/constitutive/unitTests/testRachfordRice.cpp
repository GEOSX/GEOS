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


// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/RachfordRice.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

static constexpr real64 relTol = 1e-4;

TEST( RachfordRiceTest, testRachfordRiceTwoComponents )
{
  constexpr integer numComps = 2;

  array1d< real64 > kValues( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  ////////////////////////////////////////

  kValues[0] = 1.12223;
  kValues[1] = 1.12223;

  feed[0] = 0.5;
  feed[1] = 0.5;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction1 = 1;
  real64 const vaporFraction1 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction1, expectedVaporFraction1, relTol );

  ////////////////////////////////////////

  kValues[0] = 56.1091;
  kValues[1] = 56.1091;

  feed[0] = 0.5;
  feed[1] = 0.5;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction2 = 1;
  real64 const vaporFraction2 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction2, expectedVaporFraction2, relTol );

  ////////////////////////////////////////

  kValues[0] = 54.866;
  kValues[1] = 54.866;

  feed[0] = 0.1;
  feed[1] = 0.9;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction3 = 1;
  real64 const vaporFraction3 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction3, expectedVaporFraction3, relTol );

  ////////////////////////////////////////

  kValues[0] = 1.09733;
  kValues[1] = 1.09733;

  feed[0] = 0.1;
  feed[1] = 0.9;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction4 = 1;
  real64 const vaporFraction4 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction4, expectedVaporFraction4, relTol );

  ////////////////////////////////////////

  kValues[0] = 0.9;
  kValues[1] = 1.09733;

  feed[0] = 1.0e-10;
  feed[1] = 1.0 - feed[0];

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction5 = 1;
  real64 const vaporFraction5 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction5, expectedVaporFraction5, relTol );

  ////////////////////////////////////////

  kValues[0] = 1.09733;
  kValues[1] = 0.9;

  feed[0] = 1.0e-10;
  feed[1] = 1.0 - feed[0];

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;

  real64 const expectedVaporFraction6 = 0.0;
  real64 const vaporFraction6 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction6, expectedVaporFraction6, relTol );
}

TEST( RachfordRiceTest, testRachfordRiceFourComponents )
{
  constexpr integer numComps = 4;

  array1d< real64 > kValues( numComps );
  array1d< real64 > feed( numComps );
  array1d< integer > presentComponentIds( numComps );

  ////////////////////////////////////////

  kValues[0] = 537.526;
  kValues[1] = 0.00297518;
  kValues[2] = 6.19123e-08;
  kValues[3] = 1.44696;

  feed[0] = 0.099;
  feed[1] = 0.3;
  feed[2] = 0.6;
  feed[3] = 0.001;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction1 = 0.0975568;
  real64 const vaporFraction1 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction1, expectedVaporFraction1, relTol );

  ////////////////////////////////////////

  kValues[0] = 17.4329;
  kValues[1] = 0.000770753;
  kValues[2] = 1.18694e-06;
  kValues[3] = 0.00968376;

  feed[0] = 0.1;
  feed[1] = 0.1;
  feed[2] = 0.1;
  feed[3] = 0.7;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction2 = 0.045999;
  real64 const vaporFraction2 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction2, expectedVaporFraction2, relTol );

  ////////////////////////////////////////

  kValues[0] = 11.5506;
  kValues[1] = 0.000210682;
  kValues[2] = 1.57686e-08;
  kValues[3] = 0.0442361;

  feed[0] = 0.0984186;
  feed[1] = 0.297297;
  feed[2] = 0.593142;
  feed[3] = 0.0111427;

  presentComponentIds[0] = 0;
  presentComponentIds[1] = 1;
  presentComponentIds[2] = 2;
  presentComponentIds[3] = 3;

  real64 const expectedVaporFraction3 = 0.0130259;
  real64 const vaporFraction3 =
    RachfordRice::solve( kValues.toSliceConst(),
                         feed.toSliceConst(),
                         presentComponentIds.toSliceConst() );

  checkRelativeError( vaporFraction3, expectedVaporFraction3, relTol );

}
