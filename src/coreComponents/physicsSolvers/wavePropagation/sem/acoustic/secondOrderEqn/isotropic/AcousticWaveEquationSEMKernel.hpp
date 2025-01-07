/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file AcousticWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"

namespace geos
{

/// Namespace to contain the acoustic wave kernels.
namespace acousticWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{

  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] facesToNodes map from faces to nodes
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceNodeConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverNodeConstants constant part of the receiver term
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          localIndex const numFacesPerElem,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView1d< integer const > const elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const faceNormal,
          arrayView2d< real64 const > const faceCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsAccessible,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real32 > const sourceValue,
          real64 const dt,
          real32 const timeSourceFrequency,
          real32 const timeSourceDelay,
          localIndex const rickerOrder )
  {
    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      real64 const center[3] = { elemCenter[k][0],
                                 elemCenter[k][1],
                                 elemCenter[k][2] };

      // Step 1: locate the sources, and precompute the source term

      /// loop over all the source that haven't been found yet
      for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
      {
        if( sourceIsAccessible[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          bool const sourceFound =
            WaveSolverUtils::locateSourceElement( numFacesPerElem,
                                                  center,
                                                  faceNormal,
                                                  faceCenter,
                                                  elemsToFaces[k],
                                                  coords );
          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};


            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              elemsToNodes[k],
                                                                              nodeCoords,
                                                                              coordsOnRefElem );

            sourceIsAccessible[isrc] = 1;
            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes( k, a );
              sourceConstants[isrc][a] = Ntest[a];
            }

            for( localIndex cycle = 0; cycle < sourceValue.size( 0 ); ++cycle )
            {
              sourceValue[cycle][isrc] = WaveSolverUtils::evaluateRicker( cycle * dt, timeSourceFrequency, timeSourceDelay, rickerOrder );
            }
          }
        }
      } // end loop over all sources


      // Step 2: locate the receivers, and precompute the receiver term

      /// loop over all the receivers that haven't been found yet
      for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
      {
        if( receiverIsLocal[ircv] == 0 )
        {
          real64 const coords[3] = { receiverCoordinates[ircv][0],
                                     receiverCoordinates[ircv][1],
                                     receiverCoordinates[ircv][2] };

          real64 coordsOnRefElem[3]{};
          bool const receiverFound =
            WaveSolverUtils::locateSourceElement( numFacesPerElem,
                                                  center,
                                                  faceNormal,
                                                  faceCenter,
                                                  elemsToFaces[k],
                                                  coords );

          if( receiverFound && elemGhostRank[k] < 0 )
          {
            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              elemsToNodes[k],
                                                                              nodeCoords,
                                                                              coordsOnRefElem );

            receiverIsLocal[ircv] = 1;

            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][a] = elemsToNodes( k, a );
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // end loop over receivers

    } );

  }
};

template< typename FE_TYPE >
struct MassMatrixKernel
{

  MassMatrixKernel( FE_TYPE const & finiteElement )
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
   * @param[in] density cell-wise density
   * @param[out] mass diagonal of the mass matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
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

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

template< typename FE_TYPE >
struct DampingMatrixKernel
{

  DampingMatrixKernel( FE_TYPE const & finiteElement )
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
  launch( localIndex const size,
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

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

struct PMLKernelHelper
{
  /**
   * @brief Compute the damping profile for the Perfectly Matched Layer (PML)
   * @param xLocal a given x-y-z coordinates (3-components array)
   * @param xMin coordinate limits of the inner PML boundaries, left-front-top
   * @param xMax coordinate limits of the inner PML boundaries, right-back-bottom
   * @param dMin PML thickness, left-front-top
   * @param dMax PML thickness, right-back-bottom
   * @param cMin PML wave speed, left-front-top
   * @param cMax PML wave speed, right-back-bottom
   * @param r desired reflectivity of the PML
   * @param sigma 3-components array to hold the damping profile in each direction
   */
  GEOS_HOST_DEVICE
  inline
  static void computeDampingProfilePML( real32 const (&xLocal)[3],
                                        real32 const (&xMin)[3],
                                        real32 const (&xMax)[3],
                                        real32 const (&dMin)[3],
                                        real32 const (&dMax)[3],
                                        real32 const (&cMin)[3],
                                        real32 const (&cMax)[3],
                                        real32 const r,
                                        real32 (& sigma)[3] )
  {

    sigma[0] = 0;
    sigma[1] = 0;
    sigma[2] = 0;

    if( xLocal[0] < xMin[0] )
    {
      real32 const factor =  -3.0/2.0*cMin[0]*log( r )/(dMin[0]*dMin[0]*dMin[0]);
      sigma[0] = factor*(xLocal[0]-xMin[0])*(xLocal[0]-xMin[0]);
    }
    else if( xLocal[0] > xMax[0] )
    {
      real32 const factor =  -3.0/2.0*cMax[0]*log( r )/(dMax[0]*dMax[0]*dMax[0]);
      sigma[0] = factor*(xLocal[0]-xMax[0])*(xLocal[0]-xMax[0]);
    }
    if( xLocal[1] < xMin[1] )
    {
      real32 const factor =  -3.0/2.0*cMin[1]*log( r )/(dMin[1]*dMin[1]*dMin[1]);
      sigma[1] = factor*(xLocal[1]-xMin[1])*(xLocal[1]-xMin[1]);
    }
    else if( xLocal[1] > xMax[1] )
    {
      real32 const factor =  -3.0/2.0*cMax[1]*log( r )/(dMax[1]*dMax[1]*dMax[1]);
      sigma[1] = factor*(xLocal[1]-xMax[1])*(xLocal[1]-xMax[1]);
    }
    if( xLocal[2] < xMin[2] )
    {
      real32 const factor =  -3.0/2.0*cMin[2]*log( r )/(dMin[2]*dMin[2]*dMin[2]);
      sigma[2] = factor*(xLocal[2]-xMin[2])*(xLocal[2]-xMin[2]);
    }
    else if( xLocal[2] > xMax[2] )
    {
      real32 const factor =  -3.0/2.0*cMax[2]*log( r )/(dMax[2]*dMax[2]*dMax[2]);
      sigma[2] = factor*(xLocal[2]-xMax[2])*(xLocal[2]-xMax[2]);
    }
  }
};

template< typename FE_TYPE >
struct PMLKernel
{

  PMLKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of field gradients and divergence for PML region
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] targetSet list of cells in the target set
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemToNodes constant array view of map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[in] p_n pressure field at time n
   * @param[in] v_n PML auxiliary field at time n
   * @param[in] u_n PML auxiliary field at time n
   * @param[in] xMin coordinate limits of the inner PML boundaries, left-front-top
   * @param[in] xMax coordinate limits of the inner PML boundaries, right-back-bottom
   * @param[in] dMin PML thickness, left-front-top
   * @param[in] dMax PML thickness, right-back-bottom
   * @param[in] cMin PML wave speed, left-front-top
   * @param[in] cMax PML wave speed, right-back-bottom
   * @param[in] r desired reflectivity of the PML
   * @param[out] grad_n array holding the gradients at time n
   * @param[out] divV_n array holding the divergence at time n
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( SortedArrayView< localIndex const > const targetSet,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodes,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 const > const p_n,
          arrayView2d< real32 const > const v_n,
          arrayView1d< real32 const > const u_n,
          real32 const (&xMin)[3],
          real32 const (&xMax)[3],
          real32 const (&dMin)[3],
          real32 const (&dMax)[3],
          real32 const (&cMin)[3],
          real32 const (&cMax)[3],
          real32 const r,
          arrayView2d< real32 > const grad_n,
          arrayView1d< real32 > const divV_n )
  {
    /// Loop over elements in the subregion, 'l' is the element index within the target set
    forAll< EXEC_POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const l )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

      /// global element index
      localIndex const k = targetSet[l];

      /// wave speed at the element
      real32 const c = velocity[k];

      /// coordinates of the element nodes
      real64 xLocal[ numNodesPerElem ][ 3 ];
      real32 xLocal32[ numNodesPerElem ][ 3 ];

      /// local arrays to store the pressure at all nodes and its gradient at a given node
      real64 pressure[ numNodesPerElem ];
      real64 pressureGrad[ 3 ];

      /// local arrays to store the PML vectorial auxiliary variable at all nodes and its gradient at a given node
      real64 auxV[3][ numNodesPerElem ];
      real64 auxVGrad[3][3];

      /// local arrays to store the PML scalar auxiliary variable at all nodes and its gradient at a given node
      real64 auxU[ numNodesPerElem ];
      real64 auxUGrad[3];

      /// local array to store the PML damping profile
      real32 sigma[ 3 ];

      /// copy from global to local arrays
      for( localIndex i=0; i<numNodesPerElem; ++i )
      {
        pressure[i] = p_n[elemToNodes( k, i )];
        auxU[i] = u_n[elemToNodes( k, i )];
        for( int j=0; j<3; ++j )
        {
          xLocal[i][j]   = nodeCoords[elemToNodes( k, i )][j];
          xLocal32[i][j] = nodeCoords[elemToNodes( k, i )][j];
          auxV[j][i] = v_n[elemToNodes( k, i )][j];
        }
      }

      /// local arrays to store shape functions gradients
      real64 gradN[ numNodesPerElem ][ 3 ];
      using GRADIENT_TYPE = TYPEOFREF( gradN );

      /// loop over the nodes i in the element k
      /// the nodes are implicitly assumed the same as quadrature points
      for( localIndex i=0; i<numNodesPerElem; ++i )
      {

        /// compute the shape functions gradients
        real32 const detJ = m_finiteElement.template getGradN< FE_TYPE >( k, i, xLocal, gradN );
        GEOS_UNUSED_VAR ( detJ );

        /// compute the gradient of the pressure and the PML auxiliary variables at the node
        m_finiteElement.template gradient< numNodesPerElem, GRADIENT_TYPE >( gradN, pressure, pressureGrad );
        m_finiteElement.template gradient< numNodesPerElem, GRADIENT_TYPE >( gradN, auxU, auxUGrad );
        for( int j=0; j<3; ++j )
        {
          m_finiteElement.template gradient< numNodesPerElem, GRADIENT_TYPE >( gradN, auxV[j], auxVGrad[j] );
        }

        /// compute the PML damping profile
        PMLKernelHelper::computeDampingProfilePML(
          xLocal32[i],
          xMin,
          xMax,
          dMin,
          dMax,
          cMin,
          cMax,
          r,
          sigma );

        /// compute B.pressureGrad - C.auxUGrad where B and C are functions of the damping profile
        /// WARNING: the division by 'numNodesPerElem' below is needed because the average of
        /// gradient and divergence at the nodes are sought. It is the number of cells contributing
        /// to each node that is needed. In this case, it is equal to 'numNodesPerElem'. For high-order
        /// SEM, this approach won't work and the average needs to be computed differently (maybe using counters).
        real32 localIncrementArray[3];
        localIncrementArray[0] = (sigma[0]-sigma[1]-sigma[2])*pressureGrad[0] - (sigma[1]*sigma[2])*auxUGrad[0];
        localIncrementArray[1] = (sigma[1]-sigma[0]-sigma[2])*pressureGrad[1] - (sigma[0]*sigma[2])*auxUGrad[1];
        localIncrementArray[2] = (sigma[2]-sigma[0]-sigma[1])*pressureGrad[2] - (sigma[0]*sigma[1])*auxUGrad[2];
        for( int j=0; j<3; ++j )
        {
          RAJA::atomicAdd< ATOMIC_POLICY >( &grad_n[elemToNodes( k, i )][j], localIncrementArray[j]/numNodesPerElem );
        }
        /// compute beta.pressure + gamma.u - c^2 * divV where beta and gamma are functions of the damping profile
        real32 const beta = sigma[0]*sigma[1]+sigma[0]*sigma[2]+sigma[1]*sigma[2];
        real32 const gamma = sigma[0]*sigma[1]*sigma[2];
        real32 const localIncrement = beta*p_n[elemToNodes( k, i )]
                                      + gamma*u_n[elemToNodes( k, i )]
                                      - c*c*( auxVGrad[0][0] + auxVGrad[1][1] + auxVGrad[2][2] );

        RAJA::atomicAdd< ATOMIC_POLICY >( &divV_n[elemToNodes( k, i )], localIncrement/numNodesPerElem );
      }
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};


template< typename FE_TYPE >
struct waveSpeedPMLKernel
{

  waveSpeedPMLKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of average wave speeds in the PML region
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] targetSet list of cells in the target set
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemToNodes constant array view of map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[in] xMin coordinate limits of the inner PML boundaries, left-front-top
   * @param[in] xMax coordinate limits of the inner PML boundaries, right-back-bottom
   * @param[out] cMin PML wave speed, left-front-top
   * @param[out] cMax PML wave speed, right-back-bottom
   * @param[out] counterMin PML wave speed counter, left-front-top
   * @param[out] counterMax PML wave speed counter, left-front-top
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( SortedArrayView< localIndex const > const targetSet,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodes,
          arrayView1d< real32 const > const velocity,
          real32 const (&xMin)[3],
          real32 const (&xMax)[3],
          real32 (& cMin)[3],
          real32 (& cMax)[3],
          int (& counterMin)[3],
          int (& counterMax)[3] )
  {

    RAJA::ReduceSum< parallelDeviceReduce, real32 > subRegionAvgWaveSpeedLeft( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real32 > subRegionAvgWaveSpeedRight( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real32 > subRegionAvgWaveSpeedFront( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real32 > subRegionAvgWaveSpeedBack( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real32 > subRegionAvgWaveSpeedTop( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real32 > subRegionAvgWaveSpeedBottom( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterLeft( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterRight( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterFront( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterBack( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterTop( 0 );
    RAJA::ReduceSum< parallelDeviceReduce, int > subRegionAvgWaveSpeedCounterBottom( 0 );

    /// Loop over elements in the subregion, 'l' is the element index within the target set
    forAll< EXEC_POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const l )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

      /// global element index
      localIndex const k = targetSet[l];

      /// wave speed at the element
      real32 const c = velocity[k];

      /// coordinates of the element center
      real64 xLocal[ 3 ] = {0.0, 0.0, 0.0};

      /// compute the coordinates of the element center
      for( int j=0; j<3; ++j )
      {
        for( localIndex i=0; i<numNodesPerElem; ++i )
        {
          xLocal[j] += nodeCoords[elemToNodes( k, i )][j];
        }
        xLocal[j] /= numNodesPerElem;
      }

      /// check the location of the cell and increment wave speed
      /// and counters accordingly
      if( xLocal[0] < xMin[0]
          && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1]
          && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedLeft += c;
        subRegionAvgWaveSpeedCounterLeft += 1;
      }
      else if( xLocal[0] > xMax[0]
               && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1]
               && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedRight += c;
        subRegionAvgWaveSpeedCounterRight += 1;
      }
      if( xLocal[1] < xMin[1]
          && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
          && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedFront += c;
        subRegionAvgWaveSpeedCounterFront += 1;
      }
      else if( xLocal[1] > xMax[1]
               && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
               && xLocal[2] >= xMin[2] && xLocal[2] <= xMax[2] )
      {
        subRegionAvgWaveSpeedBack += c;
        subRegionAvgWaveSpeedCounterBack += 1;
      }
      if( xLocal[2] < xMin[2]
          && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
          && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1] )
      {
        subRegionAvgWaveSpeedTop += c;
        subRegionAvgWaveSpeedCounterTop += 1;
      }
      else if( xLocal[2] > xMax[2]
               && xLocal[0] >= xMin[0] && xLocal[0] <= xMax[0]
               && xLocal[1] >= xMin[1] && xLocal[1] <= xMax[1] )
      {
        subRegionAvgWaveSpeedBottom += c;
        subRegionAvgWaveSpeedCounterBottom += 1;
      }
    } );

    /// transfer local results to global variables
    cMin[0]+=subRegionAvgWaveSpeedLeft.get();
    cMin[1]+=subRegionAvgWaveSpeedFront.get();
    cMin[2]+=subRegionAvgWaveSpeedTop.get();
    cMax[0]+=subRegionAvgWaveSpeedRight.get();
    cMax[1]+=subRegionAvgWaveSpeedBack.get();
    cMax[2]+=subRegionAvgWaveSpeedBottom.get();
    counterMin[0]+=subRegionAvgWaveSpeedCounterLeft.get();
    counterMin[1]+=subRegionAvgWaveSpeedCounterFront.get();
    counterMin[2]+=subRegionAvgWaveSpeedCounterTop.get();
    counterMax[0]+=subRegionAvgWaveSpeedCounterRight.get();
    counterMax[1]+=subRegionAvgWaveSpeedCounterBack.get();
    counterMax[2]+=subRegionAvgWaveSpeedCounterBottom.get();
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};



/**
 * @brief Implements kernels for solving the acoustic wave equations
 *   explicit central FD method and SEM
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitAcousticSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              FE_TYPE,
                                                              1,
                                                              1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   *   elements to be processed during this kernel launch.
   */
  ExplicitAcousticSEM( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_nodeCoords( nodeManager.getField< fields::referencePosition32 >() ),
    m_p_n( nodeManager.getField< fields::acousticfields::Pressure_n >() ),
    m_stiffnessVector( nodeManager.getField< fields::acousticfields::StiffnessVector >() ),
    m_density( elementSubRegion.template getField< fields::acousticfields::AcousticDensity >() ),
    m_dt( dt )
  {
    GEOS_UNUSED_VAR( edgeManager );
    GEOS_UNUSED_VAR( faceManager );
    GEOS_UNUSED_VAR( targetRegionIndex );
  }

  //*****************************************************************************
  /**
   * @copydoc geos::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitAcousticSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal(),
      stiffnessVectorLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ 8 ][ 3 ];
    real32 stiffnessVectorLocal[ numNodesPerElem ]{};
    real32 invDensity;
  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    stack.invDensity = 1./m_density[k];
    for( localIndex a=0; a< 8; a++ )
    {
      localIndex const nodeIndex =  m_elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    for( int i=0; i<numNodesPerElem; i++ )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector[m_elemsToNodes( k, i )], stack.stiffnessVectorLocal[i] );
    }
    return 0;
  }


  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticSEM Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    m_finiteElementSpace.template computeStiffnessTerm( q, stack.xLocal, [&] ( const int i, const int j, const real64 val )
    {
      real32 const localIncrement = stack.invDensity*val*m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal[ i ] += localIncrement;
    } );
  }

protected:
  /// The array containing the nodal position array.
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_nodeCoords;

  /// The array containing the nodal pressure array.
  arrayView1d< real32 const > const m_p_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure.
  arrayView1d< real32 > const m_stiffnessVector;

  /// The array containing the cell-wise density
  arrayView1d< real32 const > const m_density;

  /// The time increment for this time integration step.
  real64 const m_dt;


};



/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitAcousticSEMFactory = finiteElement::KernelFactory< ExplicitAcousticSEM,
                                                                 real64 >;


} // namespace acousticWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEMKERNEL_HPP_
