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
 * @file AcousticMatricesSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_

namespace geos
{

struct AcousticMatricesSEM
{
  //Debug
  template< typename FE_TYPE >
  struct DofArrays
  {

    DofArrays( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}
    /**
     * @brief Launches the precomputation of the mass matrices
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @tparam FE_TYPE the type of discretization
     * @param[in] finiteElement The finite element discretization used
     * @param[in] size the number of cells in the subRegion
     * @param[in] numFacesPerElem number of faces per element
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] epsilon cell-wise velocity
     * @param[in] delta cell-wise density
     * @param[out] dofEpsilon Array of epsilon on dof
     * @param[out] dofDelta Array of delta on dof
     * @param[out] dofOrder Number of elements containing the dof
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDofArrays( localIndex const size,
                      arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const GEOS_UNUSED_PARAM( nodeCoords ),
                      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                      arrayView1d< real32 const > const vti_epsilon,
                      arrayView1d< real32 const > const vti_delta,
                      arrayView1d< real32 > const dofEpsilon,
                      arrayView1d< real32 > const dofDelta,
                      arrayView1d< real32 > const dofOrder )

    {
      // First: how many element contains the dofs ?
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          localIndex a = elemsToNodes( e, q ); //global index
          real32 localIncrement = 1.0;
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofOrder[a], localIncrement ); // add one element
        }
      } ); // end loop over element
      // Second: Compute delta and epsilon on dof (approx)
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        real32 localEpsilon = std::fabs( vti_epsilon[e] );
        real32 localDelta = std::fabs( vti_delta[e] );
        if( localEpsilon < 1e-5 )
          localEpsilon = 0.;
        if( localDelta < 1e-5 )
          localDelta = 0.;
        if( localDelta > localEpsilon )
          localDelta = localEpsilon;

        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          localIndex a = elemsToNodes( e, q ); //global index
          real32 const localOrder = dofOrder[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofEpsilon[a], localEpsilon/localOrder );
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofDelta[a], localDelta/localOrder );
        }
      } ); // end loop over element

    }

//end debug

    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeQ1Params( localIndex const size,
                     arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                     arrayView1d< real32 const > const vti_epsilon,
                     arrayView1d< real32 const > const vti_delta,
                     arrayView1d< real32 > const dofEpsilon,
                     arrayView1d< real32 > const dofDelta,
                     arrayView1d< real32 > const dofOrder )

    {
      //First: compute the local mass at each vertices
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        real64 xLocal[ 8 ][ 3 ];
        for( localIndex a = 0; a < 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i = 0; i < 3; ++i )
          {
            xLocal[a][i] = nodeCoords( nodeIndex, i );
          }
        }
        // Compute Jacobian
        real64 J[3][3] = {{0}};
        m_finiteElement.jacobianTransformation( 0, 0, 0, xLocal, J );
        real64 const detJ = std::abs( LvArray::tensorOps::determinant< 3 >( J ));
        // Add contributions
        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          localIndex a = elemsToNodes( e, q ); //global index
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofOrder[a], detJ ); // Add local determinant to every nodes
        }
      } );
      // Second: compute value of parameters: delta[q] = sum(jac*delta[e])/sum(jac)
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        real32 localEpsilon = std::fabs( vti_epsilon[e] );
        real32 localDelta = std::fabs( vti_delta[e] );
        if( localEpsilon < 1e-5 )
          localEpsilon = 0.;
        if( localDelta < 1e-5 )
          localDelta = 0.;
        if( localDelta > localEpsilon )
          localDelta = localEpsilon;
        //Compute Jacobian
        real64 xLocal[ 8 ][ 3 ];
        for( localIndex a = 0; a < 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i = 0; i < 3; ++i )
          {
            xLocal[a][i] = nodeCoords( nodeIndex, i );
          }
        }
        real64 J[3][3] = {{0}};
        m_finiteElement.jacobianTransformation( 0, 0, 0, xLocal, J );
        real64 const detJ = std::abs( LvArray::tensorOps::determinant< 3 >( J ));
        //Loop on the vertices but we keep global numbering of Q_r
        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          localIndex a = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( q ));
          real32 const localMass = dofOrder[a];
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofEpsilon[a], detJ*localEpsilon / localMass );
          RAJA::atomicAdd< ATOMIC_POLICY >( &dofDelta[a], detJ*localDelta / localMass );
        }
      } ); // end loop over element
      // Compute on each nodes?


    }
    FE_TYPE const & m_finiteElement;

  };
  // End debug

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
     * @tparam FE_TYPE the type of discretization
     * @param[in] finiteElement The finite element discretization used
     * @param[in] size the number of cells in the subRegion
     * @param[in] numFacesPerElem number of faces per element
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[out] mass diagonal of the mass matrix
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeMassMatrix( localIndex const size,
                       arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                       arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                       arrayView1d< real32 const > const velocity,
                       arrayView1d< real32 const > const density,
                       arrayView1d< real32 > const mass )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {

        real32 const invC2 = 1.0 / ( density[e] * pow( velocity[e], 2 ) );
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
          real32 const localIncrement = invC2 * m_finiteElement.computeMassTerm( q, xLocal );
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes( e, q )], localIncrement );
        }
      } ); // end loop over element
    }

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
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[out] damping diagonal of the damping matrix
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDampingMatrix( localIndex const size,
                          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                          arrayView2d< localIndex const > const elemsToFaces,
                          ArrayOfArraysView< localIndex const > const facesToNodes,
                          arrayView1d< integer const > const facesDomainBoundaryIndicator,
                          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                          arrayView1d< real32 const > const velocity,
                          arrayView1d< real32 const > const density,
                          arrayView1d< real32 > const damping )
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
            real32 const alpha = 1.0 / (density[e] * velocity[e]);
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              real32 const localIncrement = alpha * m_finiteElement.computeDampingTerm( q, xLocal );
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping[facesToNodes( f, q )], localIncrement );
            }
          }
        }
      } );
    }

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
     * @param[in] lateralSurfaceFaceIndicator flag equal to 1 if the face is on the lateral surface, and to 0 otherwise
     * @param[in] bottomSurfaceFaceIndicator flag equal to 1 if the face is on the bottom surface, and to 0 otherwise
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[in] vti_epsilon cell-wise epsilon (Thomsen parameter)
     * @param[in] vti_delta density cell-wise delta (Thomsen parameter)
     * @param[in] vti_sigma sigma cell-wise parameter
     * @param[out] damping_pp Damping matrix D^{pp}
     * @param[out] damping_pq Damping matrix D^{pq}
     * @param[out] damping_qp Damping matrix D^{qp}
     * @param[out] damping_qq Damping matrix D^{qq}
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeVTIFletcherDampingMatrices( localIndex const size,
                                       arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                       arrayView2d< localIndex const > const elemsToFaces,
                                       ArrayOfArraysView< localIndex const > const facesToNodes,
                                       arrayView1d< integer const > const facesDomainBoundaryIndicator,
                                       arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                                       arrayView1d< localIndex const > const lateralSurfaceFaceIndicator,
                                       arrayView1d< localIndex const > const bottomSurfaceFaceIndicator,
                                       arrayView1d< real32 const > const velocity,
                                       arrayView1d< real32 const > const density,
                                       arrayView1d< real32 const > const vti_epsilon,
                                       arrayView1d< real32 const > const vti_delta,
                                       arrayView1d< real32 const > const vti_sigma,
                                       arrayView1d< real32 > const damping_pp,
                                       arrayView1d< real32 > const damping_pq,
                                       arrayView1d< real32 > const damping_qp,
                                       arrayView1d< real32 > const damping_qq )
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
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            real32 vti_f = 1 - (vti_epsilon[e] - vti_delta[e]) / vti_sigma[e];
            if( lateralSurfaceFaceIndicator[f] == 1 )
            {
              // ABC coefficients
              real32 alpha = 1.0 / (velocity[e] * density[e] * sqrt( 1+2*vti_epsilon[e] ));
              // VTI coefficients
              real32 vti_p_xy = 0;
              real32 vti_q_xy = 0;
              real32 vti_qp_xy= 0;

              vti_p_xy  = (1+2*vti_epsilon[e]);
              vti_q_xy  = -(vti_f - 1);
              vti_qp_xy = (vti_f+2*vti_delta[e]);
              for( localIndex q = 0; q < numNodesPerFace; ++q )
              {
                real32 const aux = m_finiteElement.computeDampingTerm( q, xLocal );
                real32 const localIncrement_p = alpha * vti_p_xy  * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pp[facesToNodes( f, q )], localIncrement_p );

                real32 const localIncrement_q = alpha*vti_q_xy * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qq[facesToNodes( f, q )], localIncrement_q );

                real32 const localIncrement_qp = alpha*vti_qp_xy * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qp[facesToNodes( f, q )], localIncrement_qp );
              }
            }
            if( bottomSurfaceFaceIndicator[f] == 1 )
            {
              // ABC coefficients updated to fit horizontal velocity
              real32 alpha = 1.0 / (velocity[e] * density[e]);
              // VTI coefficients
              real32 vti_p_z  = 0;
              real32 vti_pq_z = 0;
              real32 vti_q_z  = 0;
              vti_p_z  = -(vti_f - 1);
              vti_pq_z = vti_f;
              vti_q_z  = 1;
              for( localIndex q = 0; q < numNodesPerFace; ++q )
              {
                real32 const aux = m_finiteElement.computeDampingTerm( q, xLocal );
                real32 const localIncrement_p = alpha * vti_p_z * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pp[facesToNodes( f, q )], localIncrement_p );

                real32 const localIncrement_pq = alpha*vti_pq_z * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pq[facesToNodes( f, q )], localIncrement_pq );

                real32 const localIncrement_q = alpha * vti_q_z * aux;
                RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qq[facesToNodes( f, q )], localIncrement_q );
              }
            }
          }
        }
      } );
    }

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
     * @param[in] lateralSurfaceFaceIndicator flag equal to 1 if the face is on the lateral surface, and to 0 otherwise
     * @param[in] bottomSurfaceFaceIndicator flag equal to 1 if the face is on the bottom surface, and to 0 otherwise
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[in] vti_epsilon cell-wise epsilon (Thomsen parameter)
     * @param[in] vti_delta density cell-wise delta (Thomsen parameter)
     * @param[out] damping_pp Damping matrix D^{pp}
     * @param[out] damping_pq Damping matrix D^{pq}
     * @param[out] damping_qp Damping matrix D^{qp}
     * @param[out] damping_qq Damping matrix D^{qq}
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeVTIZhangDampingMatrices( localIndex const size,
                                    arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                    arrayView2d< localIndex const > const elemsToFaces,
                                    ArrayOfArraysView< localIndex const > const facesToNodes,
                                    arrayView1d< integer const > const facesDomainBoundaryIndicator,
                                    arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                                    arrayView1d< localIndex const > const lateralSurfaceFaceIndicator,
                                    arrayView1d< localIndex const > const bottomSurfaceFaceIndicator,
                                    arrayView1d< real32 const > const velocity,
                                    arrayView1d< real32 const > const density,
                                    arrayView1d< real32 const > const dofEpsilon,
                                    arrayView1d< real32 const > const dofDelta,
                                    arrayView1d< real32 > const damping_pp,
                                    arrayView1d< real32 > const damping_pq,
                                    arrayView1d< real32 > const damping_qp,
                                    arrayView1d< real32 > const damping_qq )
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
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            // debug

            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              //            real32 epsi = std::fabs( vti_epsilon[e] );
              real32 epsi = std::fabs( dofEpsilon[facesToNodes( f, q )] );
              // end debug

              if( std::fabs( epsi ) < 1e-5 )
                epsi = 0;
              // debug
              //         real32 delt = std::fabs( vti_delta[e] );
              real32 delt = std::fabs( dofDelta[facesToNodes( f, q )] );
              // end debug
              if( std::fabs( delt ) < 1e-5 )
                delt = 0;
              if( delt > epsi )
                delt = epsi;
              real32 sqrtEpsi = sqrt( 1 + 2 * epsi );
              // debug
              //            real32 sqrtDelta = sqrt( 1 + 2 * vti_delta[e] );
              real32 sqrtDelta = sqrt( 1 + 2 * dofDelta[facesToNodes( f, q )] );
              // end debug
              if( lateralSurfaceFaceIndicator[f] == 1 )
              {
                // ABC coefficients updated to fit horizontal velocity
                real32 alpha = 1.0 / (velocity[e] * density[e] * sqrtEpsi);
                // VTI coefficients
                real32 vti_p_xy  = 1 + 2 * epsi;
                real32 vti_qp_xy = sqrtDelta;

                for( localIndex qq = 0; qq < numNodesPerFace; ++qq )
                {
                  real32 const aux = m_finiteElement.computeDampingTerm( qq, xLocal );
                  real32 const localIncrement_p = alpha* vti_p_xy  * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pp[facesToNodes( f, qq )], localIncrement_p );

                  real32 const localIncrement_qp = alpha * vti_qp_xy * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qp[facesToNodes( f, qq )], localIncrement_qp );
                }
              }
              if( bottomSurfaceFaceIndicator[f] == 1 )
              {
                // ABC coefficients updated to fit horizontal velocity
                real32 alpha = 1.0 / (velocity[e] * density[e]);
                // VTI coefficients
                real32 vti_pq_z = sqrtDelta;
                real32 vti_q_z  = 1;
                for( localIndex qq = 0; qq < numNodesPerFace; ++qq )
                {
                  real32 const aux = m_finiteElement.computeDampingTerm( qq, xLocal );

                  real32 const localIncrement_pq = alpha * vti_pq_z * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pq[facesToNodes( f, qq )], localIncrement_pq );

                  real32 const localIncrement_q = alpha * vti_q_z * aux;
                  RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qq[facesToNodes( f, qq )], localIncrement_q );
                }
              }
            }// Debug
          }
        }
      } );
    }

    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;

  };

  template< typename FE_TYPE >
  struct GradientKappaBuoyancy
  {

    GradientKappaBuoyancy( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launch the computation of the 2 gradients relative to the coeff of the wave equation K=1/rho*c2 and b=1/rho
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] q_dt2 second order derivative in time of backward
     * @param[in] q_n current time step of backward
     * @param[in] p_n current time step of forward
     * @param[out] grad first part of gradient vector with respect to K=1/rho*c2
     * @param[out] grad2 second part of gradient vector with respact to b=1/rho
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeGradient( localIndex const size,
                     arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                     arrayView1d< integer const > const elemGhostRank,
                     arrayView1d< real32 const > const q_dt2,
                     arrayView1d< real32 const > const q_n,
                     arrayView1d< real32 const > const p_n,
                     arrayView1d< real32 > const grad,
                     arrayView1d< real32 > const grad2 )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        if( elemGhostRank[e]<0 )
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
            localIndex nodeIdx = elemsToNodes( e, q );
            grad[e] += q_dt2[nodeIdx] * p_n[nodeIdx] * m_finiteElement.computeMassTerm( q, xLocal );
            m_finiteElement.template computeStiffnessTerm( q, xLocal, [&] ( const int i, const int j, const real64 val )
            {
              grad2[e] += val* q_n[elemsToNodes( e, j )] * p_n[elemsToNodes( e, i )];
            } );
          }
        }
      } ); // end loop over element
    }
    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;
  };

};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_
