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

#ifndef GEOS_MESH_GENERATORS_PRISMUTILITIES_HPP_
#define GEOS_MESH_GENERATORS_PRISMUTILITIES_HPP_

#include "common/DataTypes.hpp"
#include "common/Span.hpp"

namespace geos
{

/**
 * @brief Get the local indices of the nodes in a face of the prism with N-sided polygon base.
 * @tparam N the number of sides in the polygon base
 * @param[in] faceNum the local index of the target face in the element  (this will be 0-numFacesInElement)
 * @param[in] elemNodes Element to nodes mapping.
 * @param[out] faceNodes space to write node indices to, must be of sufficient size
 * @return number of nodes in the face
 *
 * @note The function can be called only for N > 5. For N = 3 and N = 4 function getFaceNodesWedge
 *       and getFaceNodesHex, respectively, should be used.
 */
template< localIndex N >
localIndex getFaceNodesPrism( localIndex const faceNum,
                              arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const & elemNodes,
                              Span< localIndex > const faceNodes )
{
  static_assert( N > 4,
                 "Function getFaceNodePrism can be called for a prism with N-sided polygon base where N > 5." );
  static constexpr auto nodeCountError = "Not enough nodes for {} element (face index = {}).\n";
  GEOS_UNUSED_VAR( nodeCountError ); // Not used in GPU builds.

  if( faceNum == 0 )
  {
    GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, N, faceNum ) << generalMeshErrorAdvice );
    faceNodes[0] = elemNodes[0];
    faceNodes[1] = elemNodes[1];
    faceNodes[2] = elemNodes[N+1];
    faceNodes[3] = elemNodes[N];
    return 4;
  }
  else if( faceNum == 1 )
  {
    GEOS_ERROR_IF_LT_MSG( faceNodes.size(), N, GEOS_FMT( nodeCountError, N, faceNum ) << generalMeshErrorAdvice );
    faceNodes[0] = elemNodes[0];
    for( localIndex i = 1; i <  N; ++i )
    {
      faceNodes[i] = elemNodes[N-i];
    }
    return N;
  }
  else if( faceNum == 2 )
  {
    GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, N, faceNum ) << generalMeshErrorAdvice );
    faceNodes[0] = elemNodes[0];
    faceNodes[1] = elemNodes[N];
    faceNodes[2] = elemNodes[N*2-1];
    faceNodes[3] = elemNodes[N-1];
    return 4;
  }
  else if( faceNum >= 3 && faceNum <= N )
  {
    GEOS_ERROR_IF_LT_MSG( faceNodes.size(), 4, GEOS_FMT( nodeCountError, N, faceNum ) << generalMeshErrorAdvice );
    faceNodes[0] = elemNodes[faceNum-2];
    faceNodes[1] = elemNodes[faceNum-1];
    faceNodes[2] = elemNodes[N+faceNum-1];
    faceNodes[3] = elemNodes[N+faceNum-2];
    return 4;
  }
  else if( faceNum == N + 1 )
  {
    GEOS_ERROR_IF_LT_MSG( faceNodes.size(), N, GEOS_FMT( nodeCountError, N, faceNum ) << generalMeshErrorAdvice );
    for( localIndex i = 0; i <  N; ++i )
    {
      faceNodes[i] = elemNodes[i+N];
    }
    return N;
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Local face index out of range for Prism{} element: face index = {}.\n{}",
                          N, faceNum, generalMeshErrorAdvice ) );
    return 0;
  }

}

}

#endif // GEOS_MESH_GENERATORS_PRISMUTILITIES_HPP_
