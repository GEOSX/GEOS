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
 * @file MimeticInnerProductHelpers.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTHELPERS_HPP
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTHELPERS_HPP

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @struct MimeticInnerProductHelpers
 * @brief Helper struct handling inner product for hybrid finite volume schemes.
 */
struct MimeticInnerProductHelpers
{

  /**
   * @brief Create a full tensor from an array.
   * @param[in] values the input array
   * @param[out] result the full tensor
   */
  GEOS_HOST_DEVICE
  static
  void makeFullTensor( real64 const (&values)[ 3 ],
                       real64 (& result)[ 3 ][ 3 ] )
  {
    LvArray::tensorOps::fill< 3, 3 >( result, 0.0 );
    result[ 0 ][ 0 ] = values[ 0 ];
    result[ 1 ][ 1 ] = values[ 1 ];
    result[ 2 ][ 2 ] = values[ 2 ];
  }

  /**
   * @brief Orthonormalize a set of three vectors
   * @tparam NF number of faces in the element
   * @param[in,out] q0 first vector
   * @param[in,out] q1 second vector
   * @param[in,out] q2 third vector
   * @param[out] cellToFaceMat a copy of in/out vectors stacked into a matrix
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static
  void orthonormalize( real64 (& q0)[ NF ],
                       real64 (& q1)[ NF ],
                       real64 (& q2)[ NF ],
                       real64 (& cellToFaceMat)[ NF ][ 3 ] )
  {
    // modified Gram-Schmidt algorithm

    // q0
    LvArray::tensorOps::scale< NF >( q0, 1.0/LvArray::tensorOps::l2Norm< NF >( q0 ) );

    // q1
    real64 const q0Dotq1 = LvArray::tensorOps::AiBi< NF >( q0, q1 );
    LvArray::tensorOps::scaledAdd< NF >( q1, q0, -q0Dotq1 );
    LvArray::tensorOps::scale< NF >( q1, 1.0/LvArray::tensorOps::l2Norm< NF >( q1 ) );

    // q2
    real64 const q0Dotq2 = LvArray::tensorOps::AiBi< NF >( q0, q2 );
    LvArray::tensorOps::scaledAdd< NF >( q2, q0, -q0Dotq2 );
    real64 const q1Dotq2 = LvArray::tensorOps::AiBi< NF >( q1, q2 );
    LvArray::tensorOps::scaledAdd< NF >( q2, q1, -q1Dotq2 );
    LvArray::tensorOps::scale< NF >( q2, 1.0/LvArray::tensorOps::l2Norm< NF >( q2 ) );

    for( localIndex i = 0; i < NF; ++i )
    {
      cellToFaceMat[ i ][ 0 ] = q0[ i ];
      cellToFaceMat[ i ][ 1 ] = q1[ i ];
      cellToFaceMat[ i ][ 2 ] = q2[ i ];
    }
  }

  /**
   * @brief For a given face, compute the TPFA entry incorporating the multiplier
   * @tparam NF number of faces in the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] faceNormal the unit normal vector to the face
   * @param[in] faceArea the area of the face
   * @param[in] transMult the transmissibility multiplier for this face
   * @param[in] weightToleranceInv tolerance used in the trans computation
   * @param[in,out] cellToFaceVec the vector from the element center to the face center
   * @param[out] tpTransInv the TPFA entry incorporating the multiplier
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static void
  computeInvTPFATransWithMultiplier( real64 const (&elemPerm)[ 3 ],
                                     real64 const (&faceNormal)[ 3 ],
                                     real64 const & faceArea,
                                     real64 const & transMult,
                                     real64 const & weightToleranceInv,
                                     real64 (& cellToFaceVec)[ 3 ],
                                     real64 & tpTransInv )
  {
    real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );
    real64 const mult = transMult;
    tpTransInv = c2fDistance / faceArea;

    real64 faceConormal[ 3 ] = { 0.0 };
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, elemPerm, faceNormal );
    real64 halfWeight = LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );
    if( halfWeight < 0.0 )
    {
      LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, elemPerm, cellToFaceVec );
      halfWeight = LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );
    }
    tpTransInv /= halfWeight;
    tpTransInv  = LvArray::math::min( tpTransInv, weightToleranceInv );
    tpTransInv *=  ( 1.0 - mult ) / mult;
  }


  /**
   * @brief Incorporate the transmissibility multiplier into the transmissibility matrix
   * @tparam NF number of faces in the element
   * @param[in,out] tpTransInv inverse of the (diagonal of the) TPFA transmissibility matrix (already accounting for the multiplier
   * @param[in,out] transMatrix transmissibility matrix
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static void
  computeTransMatrixWithMultipliers( real64 const (&tpTransInv)[ NF ],
                                     arraySlice2d< real64 > const & transMatrix )
  {
    // the inverse of the pertubed inverse is computed using the Sherman-Morrison formula
    for( localIndex k = 0; k < NF; ++k )
    {
      real64 const mult = LvArray::math::sqrt( tpTransInv[k] );
      real64 Tmult[ NF ] = { 0.0 };
      for( localIndex i = 0; i < NF; ++i )
      {
        Tmult[i] = transMatrix[k][i] * mult;
      }

      real64 const invDenom = 1.0 / ( 1.0 + Tmult[k] * mult );
      for( localIndex i = 0; i < NF; ++i )
      {
        for( localIndex j = 0; j < NF; ++j )
        {
          transMatrix[i][j] -= Tmult[i]*Tmult[j]*invDenom;
        }
      }
    }
  }

};

} // namespace mimeticInnerProduct

} // namespace geos

#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTHELPERS_HPP
