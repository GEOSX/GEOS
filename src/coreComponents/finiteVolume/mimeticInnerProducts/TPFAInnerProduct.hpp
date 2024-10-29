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
 * @file TPFAInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_TPFAINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_TPFAINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class TPFAInnerProduct
 *
 * Provides an implementation of a TPFA inner product in the hybrid FVM solvers
 */
class TPFAInnerProduct : public MimeticInnerProductBase
{
public:

  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using TPFA.
   * @param[in] nodePosition the position of the nodes
   * @param[in] transMultiplier the transmissibility multipliers at the mesh faces
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static void
  compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           arrayView1d< real64 const > const & transMultiplier,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const & elemToFaces,
           arraySlice1d< real64 const > const & elemCenter,
           real64 const & elemVolume,
           real64 const (&elemPerm)[ 3 ],
           real64 const & lengthTolerance,
           arraySlice2d< real64 > const & transMatrix );

};

template< localIndex NF >
GEOS_HOST_DEVICE
void
TPFAInnerProduct::compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                           arrayView1d< real64 const > const & transMultiplier,
                           ArrayOfArraysView< localIndex const > const & faceToNodes,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           arraySlice1d< real64 const > const & elemCenter,
                           real64 const & elemVolume,
                           real64 const (&elemPerm)[ 3 ],
                           real64 const & lengthTolerance,
                           arraySlice2d< real64 > const & transMatrix )
{
  GEOS_UNUSED_VAR( elemVolume );

  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance;

  // 0) assemble full coefficient tensor from principal axis/components
  real64 permTensor[ 3 ][ 3 ] = {{ 0 }};
  MimeticInnerProductHelpers::makeFullTensor( elemPerm, permTensor );

  // we are ready to compute the transmissibility matrix
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 const mult = transMultiplier[elemToFaces[ifaceLoc]];

    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      // for now, TPFA trans
      if( ifaceLoc == jfaceLoc )
      {
        real64 faceCenter[ 3 ], faceNormal[ 3 ], faceConormal[ 3 ], cellToFaceVec[ 3 ];
        // 1) compute the face geometry data: center, normal, vector from cell center to face center
        real64 const faceArea =
          computationalGeometry::centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                     nodePosition,
                                                     faceCenter,
                                                     faceNormal,
                                                     areaTolerance );

        LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
        LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );

        if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
        {
          LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
        }

        real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

        LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, elemPerm, faceNormal );

        // 3) compute the one-sided face transmissibility
        transMatrix[ifaceLoc][jfaceLoc]  = LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );
        transMatrix[ifaceLoc][jfaceLoc] *= mult * faceArea / c2fDistance;
        transMatrix[ifaceLoc][jfaceLoc]  = LvArray::math::max( transMatrix[ifaceLoc][jfaceLoc], weightTolerance );
      }
      else
      {
        transMatrix[ifaceLoc][jfaceLoc] = 0;
      }
    }
  }
}

} // end namespace mimeticInnerProduct

} // end namespace geos


#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_TPFAINNERPRODUCT_HPP_
