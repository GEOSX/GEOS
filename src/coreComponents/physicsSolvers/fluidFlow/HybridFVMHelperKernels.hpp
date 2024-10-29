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

/**
 * @file CompositionalMultiphaseHybridFVMHelperKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_HYBRIDFVMUPWINDINGHELPERKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_HYBRIDFVMUPWINDINGHELPERKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{

namespace hybridFVMKernels
{


/******************************** CellConnectivity ********************************/

struct CellConnectivity
{

  GEOS_HOST_DEVICE
  static bool
  isNeighborFound( localIndex const (&localIds)[3],
                   localIndex const ifaceLoc,
                   arrayView2d< localIndex const > const & elemRegionList,
                   arrayView2d< localIndex const > const & elemSubRegionList,
                   arrayView2d< localIndex const > const & elemList,
                   SortedArrayView< localIndex const > const & regionFilter,
                   arraySlice1d< localIndex const > const & elemToFaces,
                   localIndex ( & neighborIds )[3] )
  {
    neighborIds[0] = localIds[0];
    neighborIds[1] = localIds[1];
    neighborIds[2] = localIds[2];

    // the face has at most two adjacent elements
    // one of these two elements is the current element indexed by er, esr, ei
    // but here we are interested in the indices of the other element
    // this other element is "the neighbor" for this one-sided face
    for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
    {

      localIndex const erNeighbor  = elemRegionList[elemToFaces[ifaceLoc]][k];
      localIndex const esrNeighbor = elemSubRegionList[elemToFaces[ifaceLoc]][k];
      localIndex const eiNeighbor  = elemList[elemToFaces[ifaceLoc]][k];

      // this element is not the current element
      // we have found the neighbor or we are at the boundary
      if( erNeighbor != localIds[0] || esrNeighbor != localIds[1] || eiNeighbor != localIds[2] )
      {
        bool const onBoundary       = (erNeighbor == -1 || esrNeighbor == -1 || eiNeighbor == -1);
        bool const neighborInTarget = regionFilter.contains( erNeighbor );

        // if not on boundary, save the element indices
        if( !onBoundary && neighborInTarget )
        {
          neighborIds[0] = erNeighbor;
          neighborIds[1] = esrNeighbor;
          neighborIds[2] = eiNeighbor;
        }
        // if the face is on the boundary, use the properties of the local elem
      }
    }
    return !( localIds[0] == neighborIds[0] &&
              localIds[1] == neighborIds[1] &&
              localIds[2] == neighborIds[2] );
  }

};


} // namespace hybridFVMUpwindingKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMUPWINDINGHELPERKERNELS_HPP
