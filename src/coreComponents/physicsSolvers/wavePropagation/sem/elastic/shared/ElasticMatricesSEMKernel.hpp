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
 * @file ElasticMatricesSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICMATRICESSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICMATRICESSEMKERNEL_HPP_

namespace geos
{

struct ElasticMatricesSEM
{
  template< typename FE_TYPE >
  struct MassMatrix
  {

    MassMatrix( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launches the precomputation of the mass matrices
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] numFacesPerElem number of faces per element
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] velocity cell-wise velocity
     * @param[out] mass diagonal of the mass matrix
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeMassMatrix( localIndex const size,
                       arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                       arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                       arrayView1d< real32 const > const density,
                       arrayView1d< real32 > const mass )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {

        // only the eight corners of the mesh cell are needed to compute the Jacobian
        real64 xLocal[ 8 ][ 3 ];
        for( localIndex a = 0; a < 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i = 0; i < 3; ++i )
          {
            xLocal[a][i] = nodeCoords( nodeIndex, i );
          }
        }

        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          real32 const localIncrement = density[e] * m_finiteElement.computeMassTerm( q, xLocal );
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes( e, q )], localIncrement );
        }
      } ); // end loop over element
    }

    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;

  };

  template< typename FE_TYPE >
  struct DampingMatrix
  {

    DampingMatrix( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launches the precomputation of the damping matrices
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToFaces map from elements to faces
     * @param[in] facesToNodes map from face to nodes
     * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
     * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
     * @param[in] faceNormal array containing the normal fo the faces
     * @param[in] density cell-wise density
     * @param[in] velocityVp cell-wise P-wavespeed
     * @param[in] velocityVs cell-wise S-wavespeed
     * @param[out] dampingx diagonal of the damping matrix (x-direction)
     * @param[out] dampingy diagonal of the damping matrix (y-direction)
     * @param[out] dampingz diagonal of the damping matrix (z-direction)
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDampingMatrix( localIndex const size,
                          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                          arrayView2d< localIndex const > const elemsToFaces,
                          ArrayOfArraysView< localIndex const > const facesToNodes,
                          arrayView1d< integer const > const facesDomainBoundaryIndicator,
                          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                          arrayView2d< real64 const > const faceNormal,
                          arrayView1d< real32 const > const density,
                          arrayView1d< real32 const > const velocityVp,
                          arrayView1d< real32 const > const velocityVs,
                          arrayView1d< real32 > const dampingx,
                          arrayView1d< real32 > const dampingy,
                          arrayView1d< real32 > const dampingz )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
        {
          localIndex const f = elemsToFaces( e, i );
          // face on the domain boundary and not on free surface
          if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
          {
            // only the four corners of the mesh face are needed to compute the Jacobian
            real64 xLocal[ 4 ][ 3 ];
            for( localIndex a = 0; a < 4; ++a )
            {
              localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
              for( localIndex d = 0; d < 3; ++d )
              {
                xLocal[a][d] = nodeCoords( nodeIndex, d );
              }
            }

            real32 const nx = faceNormal( f, 0 ), ny = faceNormal( f, 1 ), nz = faceNormal( f, 2 );
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              real32 const aux = density[e] * m_finiteElement.computeDampingTerm( q, xLocal );
              real32 const localIncrementx = aux * ( velocityVp[e] * LvArray::math::abs( nx ) + velocityVs[e] * LvArray::math::sqrt( pow( ny, 2 ) + pow( nz, 2 ) ) );
              real32 const localIncrementy = aux * ( velocityVp[e] * LvArray::math::abs( ny ) + velocityVs[e] * LvArray::math::sqrt( pow( nx, 2 ) + pow( nz, 2 ) ) );
              real32 const localIncrementz = aux * ( velocityVp[e] * LvArray::math::abs( nz ) + velocityVs[e] * LvArray::math::sqrt( pow( nx, 2 ) + pow( ny, 2 ) ) );

              RAJA::atomicAdd< ATOMIC_POLICY >( &dampingx[facesToNodes( f, q )], localIncrementx );
              RAJA::atomicAdd< ATOMIC_POLICY >( &dampingy[facesToNodes( f, q )], localIncrementy );
              RAJA::atomicAdd< ATOMIC_POLICY >( &dampingz[facesToNodes( f, q )], localIncrementz );
            }
          }
        }
      } );
    }

    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;

  };



};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICMATRICESSEMKERNEL_HPP_
