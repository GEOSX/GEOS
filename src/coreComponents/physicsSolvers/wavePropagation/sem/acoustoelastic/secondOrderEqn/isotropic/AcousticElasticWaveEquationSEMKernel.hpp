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
 * @file AcousticElasticWaveEquationSEMKernel.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif

#include <utility>

namespace geos
{

namespace acousticElasticWaveEquationSEMKernels
{

template< typename FE_TYPE >
struct CouplingKernel
{
  static constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;


  /**
   * @brief Launches the computation of the coupling vector
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of faces
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] fluidRegionIndex Acoustic region
   * @param[in] fluidSubRegionIndex Acoustic subregion
   * @param[in] faceToSubRegion Array which gives you the subregion on which the face belongs
   * @param[in] faceToRegion Array which gives you the region on which the face belongs
   * @param[in] faceToElement Array which gives you the element on which the face belongs
   * @param[in] facesToNodes Array which gives the nodes on the mesh knowing the face and the local node on an element
   * @param[in] faceNormal Normals of the faces
   * @param[in] faceCenters Centers of the faces
   * @param[in] elemCenters Centers of the elements
   * @param[out] couplingVectorx x-component of the coupling vector
   * @param[out] couplingVectory y-component of the coupling vector
   * @param[out] couplingVectorz z-component of the coupling vector
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const,
                       nodes::REFERENCE_POSITION_USD > const nodeCoords,
          localIndex const fluidRegionIndex,
          localIndex const fluidSubRegionIndex,
          arrayView2d< localIndex const > const faceToSubRegion,
          arrayView2d< localIndex const > const faceToRegion,
          arrayView2d< localIndex const > const faceToElement,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView2d< real64 const > const faceNormals,
          arrayView2d< real64 const > const faceCenters,
          arrayView2d< real64 const > const elemCenters,
          arrayView1d< real32 > const couplingVectorx,
          arrayView1d< real32 > const couplingVectory,
          arrayView1d< real32 > const couplingVectorz )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const f )
    {
      localIndex const e0 = faceToElement( f, 0 ), e1 = faceToElement( f, 1 );
      localIndex const er0 = faceToRegion( f, 0 ), er1 = faceToRegion( f, 1 );
      localIndex const esr0 = faceToSubRegion( f, 0 ), esr1 = faceToSubRegion( f, 1 );

      if( e0 != -1 && e1 != -1 && er0 != er1 )  // an interface is defined as a transition between regions
      {
        // check that one of the region is the fluid subregion for the fluid -> solid coupling term
        bool const e0IsFluid = er0 == fluidRegionIndex && esr0 == fluidSubRegionIndex;
        bool const e1IsFluid = er1 == fluidRegionIndex && esr1 == fluidSubRegionIndex;

        if( e0IsFluid != e1IsFluid )  // xor: a single element must be fluid
        {
          // only the four corners of the mesh face are needed to compute the Jacobian
          real64 xLocal[ 4 ][ 3 ];
          for( localIndex a = 0; a < 4; ++a )
          {
            localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = nodeCoords( nodeIndex, i );
            }
          }

          real64 const nx = faceNormals( f, 0 ), ny = faceNormals( f, 1 ), nz = faceNormals( f, 2 );

          // determine sign to get an outward pointing normal for the fluid -> solid coupling
          localIndex const e = e0IsFluid ? e0 : (e1IsFluid ? e1 : -1);  // fluid element
          localIndex const sgn = (
            (faceCenters( f, 0 ) - elemCenters( e, 0 )) * nx +
            (faceCenters( f, 1 ) - elemCenters( e, 1 )) * ny +
            (faceCenters( f, 2 ) - elemCenters( e, 2 )) * nz
            ) < 0 ? 1 : -1;

          for( localIndex q = 0; q < numNodesPerFace; ++q )
          {
            real64 const aux = FE_TYPE::computeDampingTerm( q, xLocal );

            real32 const localIncrementx = aux * (sgn * nx);
            real32 const localIncrementy = aux * (sgn * ny);
            real32 const localIncrementz = aux * (sgn * nz);

            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectorx[facesToNodes( f, q )], localIncrementx );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectory[facesToNodes( f, q )], localIncrementy );
            RAJA::atomicAdd< ATOMIC_POLICY >( &couplingVectorz[facesToNodes( f, q )], localIncrementz );
          }
        }
      }
    } );

  }
};

} /* namespace acousticElasticWaveEquationSEMKernels */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICELASTICWAVEEQUATIONSEMKERNEL_HPP_ */
