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
 * @file BdVLMInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_BDVLMINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_BDVLMINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductHelpers.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class BdVLMInnerProduct
 *
 * Provides an implementation of the inner product proposed by Beirao da Veiga, Lipnikov, Manzini in the hybrid FVM solvers
 */
class BdVLMInnerProduct : public MimeticInnerProductBase
{
public:

  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using the inner product of Beirao da Veiga, Lipnikov,
   * Manzini
   * @param[in] nodePosition the position of the nodes
   * @param[in] transMultiplier the transmissibility multipliers at the mesh faces
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   *
   * @details Reference: Beirao da Veiga, Lipnikov, Manzini, "The mimetic finite-difference method for elliptic problems"
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
BdVLMInnerProduct::compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
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
  real64 const weightToleranceInv = 1e30 / lengthTolerance;

  real64 cellToFaceMat[ NF ][ 3 ] = {{ 0 }};
  real64 normalsMat[ NF ][ 3 ] = {{ 0 }};
  real64 permMat[ 3 ][ 3 ] = {{ 0 }};
  real64 faceAreaMat[ NF ][ NF ] = {{ 0 }};

  real64 work_dimByDim[ 3 ][ 3 ] = {{ 0 }};
  real64 work_numFacesByDim[ NF ][ 3 ] = {{ 0 }};
  real64 work_dimByNumFaces[ 3 ][ NF ] = {{ 0 }};
  real64 work_numFacesByNumFaces[ NF ][ NF ] = {{ 0 }};

  real64 tpTransInv[ NF ] = { 0.0 };

  // 0) assemble full coefficient tensor from principal axis/components
  MimeticInnerProductHelpers::makeFullTensor( elemPerm, permMat );

  // 1) fill the matrices cellToFaceMat and normalsMat row by row
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 faceCenter[ 3 ], faceNormal[ 3 ], cellToFaceVec[ 3 ];
    // compute the face geometry data: center, normal, vector from cell center to face center
    faceAreaMat[ ifaceLoc ][ ifaceLoc ] =
      computationalGeometry::centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                 nodePosition,
                                                 faceCenter,
                                                 faceNormal,
                                                 areaTolerance );

    LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
    LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );

    cellToFaceMat[ ifaceLoc ][0] = faceAreaMat[ ifaceLoc ][ ifaceLoc ] * cellToFaceVec[ 0 ];
    cellToFaceMat[ ifaceLoc ][1] = faceAreaMat[ ifaceLoc ][ ifaceLoc ] * cellToFaceVec[ 1 ];
    cellToFaceMat[ ifaceLoc ][2] = faceAreaMat[ ifaceLoc ][ ifaceLoc ] * cellToFaceVec[ 2 ];

    if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    // the two-point transmissibility is computed to computed here because it is needed
    // in the implementation of the transmissibility multiplier (see below)
    // TODO: see what it would take to bring the (harmonically averaged) two-point trans here
    real64 diagEntry = 0.0;
    MimeticInnerProductHelpers::computeInvTPFATransWithMultiplier< NF >( elemPerm,
                                                                         faceNormal,
                                                                         faceAreaMat[ifaceLoc][ifaceLoc],
                                                                         transMultiplier[elemToFaces[ifaceLoc]],
                                                                         weightToleranceInv,
                                                                         cellToFaceVec,
                                                                         diagEntry );
    tpTransInv[ifaceLoc] = diagEntry;

    normalsMat[ ifaceLoc ][ 0 ] = faceNormal[ 0 ];
    normalsMat[ ifaceLoc ][ 1 ] = faceNormal[ 1 ];
    normalsMat[ ifaceLoc ][ 2 ] = faceNormal[ 2 ];

  }

  // 2) compute N of Beirao da Veiga, Lipnikov, Manzini
  LvArray::tensorOps::Rij_eq_AikBkj< NF, 3, 3 >( work_numFacesByDim,
                                                 normalsMat,
                                                 permMat );

  // 3) compute ( N^T R )^-1 of Beirao da Veiga, Lipnikov, Manzini
  LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, NF >( work_dimByDim,
                                                 work_numFacesByDim,
                                                 cellToFaceMat );
  LvArray::tensorOps::invert< 3 >( work_dimByDim );

  // 4) compute N ( N^T R )^-1 N^T
  LvArray::tensorOps::Rij_eq_AikBjk< 3, NF, 3 >( work_dimByNumFaces,
                                                 work_dimByDim,
                                                 work_numFacesByDim );
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, 3 >( transMatrix,
                                                  work_numFacesByDim,
                                                  work_dimByNumFaces );

  // 5) compute the stabilization coefficient \tilde{ \lambda }
  real64 const stabCoef = 2.0 / NF * LvArray::tensorOps::trace< NF >( transMatrix );

  // at this point we have N ( N^T R )^-1 N^T and \tilde{ \lambda }
  // we have to compute I - R ( R^T R )^-1 R^T and sum

  // 6) compute ( R^T R )^-1
  LvArray::tensorOps::Rij_eq_AkiAkj< 3, NF >( work_dimByDim,
                                              cellToFaceMat );
  LvArray::tensorOps::invert< 3 >( work_dimByDim );

  // 7) compute I - R ( R^T R )^-1 R^T
  LvArray::tensorOps::addIdentity< NF >( work_numFacesByNumFaces, -1 );
  LvArray::tensorOps::Rij_eq_AikBjk< 3, NF, 3 >( work_dimByNumFaces,
                                                 work_dimByDim,
                                                 cellToFaceMat );
  LvArray::tensorOps::Rij_add_AikBkj< NF, NF, 3 >( work_numFacesByNumFaces,
                                                   cellToFaceMat,
                                                   work_dimByNumFaces );

  // 8) compute N ( N^T R )^-1 N^T + \tilde{ \lambda } * (I - R ( R^T R )^-1 R^T)
  LvArray::tensorOps::scaledAdd< NF, NF >( transMatrix, work_numFacesByNumFaces, -stabCoef );


  // 9) this IP was designed for velocities,
  //    so we have to pre- and post-multiply by faceAreaMat to get an IP for fluxes
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, NF >( work_numFacesByNumFaces,
                                                   transMatrix,
                                                   faceAreaMat );
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, NF >( transMatrix,
                                                   faceAreaMat,
                                                   work_numFacesByNumFaces );

  // 10) incorporate the transmissbility multipliers
  // Ref: Nilsen, H. M., J. R. Natvig, and K.-A Lie.,
  // "Accurate modeling of faults by multipoint, mimetic, and mixed methods." SPEJ

  if( !isZero( LvArray::tensorOps::l2NormSquared< NF >( tpTransInv ) ) )
  {
    MimeticInnerProductHelpers::computeTransMatrixWithMultipliers< NF >( tpTransInv,
                                                                         transMatrix );
  }

}

} // end namespace mimeticInnerProduct

} // end namespace geos


#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_BDVLMINNERPRODUCT_HPP_
