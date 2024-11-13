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
 * @file testComputationalGeometry.cpp
 */
#include "../utilities/ComputationalGeometry.hpp"
#include <gtest/gtest.h>
#include <numeric>

namespace geos
{

TEST( testComputationalGeometry, checkCentroid3DPolygon )
{
  constexpr localIndex n = 12;
  real64 t = 2.0*M_PI/n;
  real64 s = (1.0 + sin( t ))/cos( t ) - 1.0;
  constexpr localIndex nodeDim = 3;
  constexpr localIndex numNodes = n*3;

  real64 e[3] = { -1.0e-6, 0.0, 1.0e-6 };

  array2d< real64, nodes::REFERENCE_POSITION_PERM > points;

  points.resize( numNodes, nodeDim );

  localIndex idx = 0;
  for( localIndex z = 0; z < 3; ++z )
  {
    for( localIndex i = 0; i < n; ++i )
    {
      real64 r = 25.0*(1.0 + (i%2)*s + e[z]);
      real64 x = r*cos( i*t );
      real64 y = r*sin( i*t );
      points( idx, 0 ) = x;
      points( idx, 1 ) = y;
      points( idx, 2 ) = 10.0*z;
      ++idx;
    }
  }

  constexpr real64 areaTolerance = 0.0001;

  for( localIndex ci = 0; ci < 3; ++ci )
  {
    array1d< localIndex > indices;
    for( localIndex i = 0; i < n; ++i )
    {
      indices.emplace_back( ci*n + i );
    }

    real64 faceCenter[ 3 ], faceNormal[ 3 ];
    computationalGeometry::centroid_3DPolygon( indices.toSliceConst(), points.toViewConst(), faceCenter, faceNormal, areaTolerance );

    real64 norm = LvArray::math::square( faceNormal[0] ) + LvArray::math::square( faceNormal[1] ) + LvArray::math::square( faceNormal[2] );
    EXPECT_EQ( abs( sqrt( norm ) - 1.0 ) < 0.0001, true );
  }
}

} /* namespace geos */
