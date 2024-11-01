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
 * @file PTScotchInterface.cpp
 */

#include "PTScotchInterface.hpp"

#include <ptscotch.h>

#include <numeric>

static_assert( std::is_same< int64_t, SCOTCH_Num >::value,
               "Non-matching index types. Scotch must be configured with 64-bit indices." );

#define GEOS_SCOTCH_CHECK( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_NE_MSG( ierr, 0, "Error in call to:\n" << #call ); \
  } while( false )

namespace geos
{
namespace ptscotch
{

array1d< int64_t >
partition( ArrayOfArraysView< int64_t const, int64_t > const & graph,
           int64_t const numParts,
           MPI_Comm comm )
{
  SCOTCH_Num const numVerts = graph.size();

  array1d< int64_t > part( numVerts ); // all 0 by default
  if( numParts == 1 )
  {
    return part;
  }

  SCOTCH_Dgraph * const gr = SCOTCH_dgraphAlloc();
  GEOS_SCOTCH_CHECK( SCOTCH_dgraphInit( gr, comm ) );

  SCOTCH_Num const numEdges = graph.getOffsets()[numVerts];

  // Technical UB if Scotch writes into these arrays; in practice we discard them right after
  SCOTCH_Num * const offsets = const_cast< SCOTCH_Num * >( graph.getOffsets() );
  SCOTCH_Num * const edges = const_cast< SCOTCH_Num * >( graph.getValues() );

  GEOS_SCOTCH_CHECK( SCOTCH_dgraphBuild( gr,          // graphptr
                                         0,            // baseval
                                         numVerts,     // vertlocnbr
                                         numVerts,     // vertlocmax
                                         offsets,      // vertloctab
                                         offsets + 1,  // vendloctab
                                         nullptr,      // veloloctab
                                         nullptr,      // vlblloctab
                                         numEdges,     // edgelocnbr
                                         numEdges,     // edgelocsiz
                                         edges,        // edgeloctab
                                         nullptr,      // edgegsttab
                                         nullptr       // edloloctab,
                                         ) );

  // TODO: maybe remove?
  GEOS_SCOTCH_CHECK( SCOTCH_dgraphCheck( gr ) );

  SCOTCH_Strat * const strategy = SCOTCH_stratAlloc();
  GEOS_SCOTCH_CHECK( SCOTCH_stratInit( strategy ) );

  GEOS_SCOTCH_CHECK( SCOTCH_dgraphPart( gr, numParts, strategy, part.data() ) );

  SCOTCH_stratExit( strategy );
  SCOTCH_dgraphExit( gr );

  return part;
}

} // namespace ptscotch
} // namespace geos
