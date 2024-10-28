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
 * @file SolidMechanicsALMKernelsBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_ROTATIONMATRIXKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_ROTATIONMATRIXKERNEL_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"

namespace geos
{

namespace rotationMatrixKernel
{

/**
 * @brief A struct to compute rotation matrices
 */
struct ComputeRotationMatricesKernel
{

  /**
   * @brief Launch the kernel function to comute rotation matrices
   * @tparam POLICY the type of policy used in the kernel launch
   * @param[in] size the size of the subregion
   * @param[in] faceNormal the array of array containing the face to nodes map
   * @param[in] elemsToFaces the array of array containing the element to faces map
   * @param[out] rotationMatrix the array containing the rotation matrices
   */
  template< typename POLICY >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const > const & faceNormal,
          ArrayOfArraysView< localIndex const > const & elemsToFaces,
          arrayView3d< real64 > const & rotationMatrix )
  {

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      localIndex const & f0 = elemsToFaces[k][0];
      localIndex const & f1 = elemsToFaces[k][1];

      real64 Nbar[3];
      Nbar[0] = faceNormal[f0][0] - faceNormal[f1][0];
      Nbar[1] = faceNormal[f0][1] - faceNormal[f1][1];
      Nbar[2] = faceNormal[f0][2] - faceNormal[f1][2];

      LvArray::tensorOps::normalize< 3 >( Nbar );
      computationalGeometry::RotationMatrix_3D( Nbar, rotationMatrix[k] );

    } );
  }

};


} // namespace solidMechanicsLagrangeContactKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_ROTATIONMATRIXKERNEL_HPP_ */
