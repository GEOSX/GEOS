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

#ifndef GEOS_FINITEELEMENT_UNITTESTS_TESTFINITEELEMENTHELPERS_HPP_
#define GEOS_FINITEELEMENT_UNITTESTS_TESTFINITEELEMENTHELPERS_HPP_

#include "common/DataTypes.hpp"
#include <cstdlib>

namespace geos
{
inline real64 rando( real64 const min,
                     real64 const max )
{
  return ( rand() % 101 ) / 100.0 *(max-min) + min;
}


template< int NUM_SUPPORT_POINTS >
inline void randomShape( real64 (& N)[NUM_SUPPORT_POINTS] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    N[a] = rando( 0.0, 1.0 );
  }
}

template< int NUM_SUPPORT_POINTS >
inline void randomShapeGradient( real64 (& gradN)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    gradN[a][0] = rando( -0.4, 0.4 );
    gradN[a][1] = rando( -0.4, 0.4 );
    gradN[a][2] = rando( -0.4, 0.4 );
  }
}

template< int NUM_SUPPORT_POINTS >
inline void randomSupportVar( real64 (& v)[NUM_SUPPORT_POINTS] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    v[a] = rando( -1.0, 1.0 );
  }
}

template< int NUM_SUPPORT_POINTS >
inline void randomSupportVar( real64 (& v)[NUM_SUPPORT_POINTS][3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    v[a][0] = rando( -1, 1 );
    v[a][1] = rando( -1, 1 );
    v[a][2] = rando( -1, 1 );
  }
}

inline void randomVar( real64 (& v)[6] )
{
  for( int i=0; i<6; ++i )
  {
    v[i] = rando( -1.0, 1.0 );
  }
}

inline void randomVar( real64 (& v)[3] )
{
  for( int i=0; i<3; ++i )
  {
    v[i] = rando( -1.0, 1.0 );
  }
}

inline void randomVar( real64 (& v)[3][3] )
{
  for( int i=0; i<3; ++i )
  {
    for( int j=0; j<3; ++j )
    {
      v[i][j] = rando( -1.0, 1.0 );
    }
  }
}

/*
   template< typename FE_FORMULATION >
   void createTestElement( real64 (& X)[FE_FORMULATION::numNodes][3] )
   {
   constexpr int numNodes = FE_FORMULATION::numNodes;
   constexpr int numQuadraturePts = FE_FORMULATION::numQuadraturePoints;
   real64 xCoords[numNodes][3];
   real64 pCoords[numNodes][3];

   srand( 1234 );
   for( int a=0; a<numNodes; ++a )
   {
    pCoords[a][0] = FE_FORMULATION::parentCoords0( a );
    pCoords[a][1] = FE_FORMULATION::parentCoords1( a );
    pCoords[a][2] = FE_FORMULATION::parentCoords2( a );

    xCoords[a][0] = FE_FORMULATION::parentCoords0( a ) + rando( -0.3, 0.3 );
    xCoords[a][1] = FE_FORMULATION::parentCoords1( a ) + rando( -0.3, 0.3 );
    xCoords[a][2] = FE_FORMULATION::parentCoords2( a ) + rando( -0.3, 0.3 );
   }
   }
 */

}

#endif
