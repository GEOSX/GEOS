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
 * @file PrecomputeSourcesAndReceiversKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_PRECOMPUTESOURCESANDRECEIVERSKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_PRECOMPUTESOURCESANDRECEIVERSKERNEL_HPP_

namespace geos
{

struct PreComputeSourcesAndReceivers
{

  using EXEC_POLICY = parallelDevicePolicy< >;

  /**
   * @brief Launches the precomputation of the source and receiver terms for 1D solution (2nd order acoustic)
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] baseFacesToNodes face to node map
   * @param[in] baseNodeCoords coordinates of the nodes
   * @param[in] baseNodeLocalToGlobal local to global index map for nodes
   * @param[in] elementLocalToGlobal local to global index map for elements
   * @param[in] baseNodesToElements node to element map for the base mesh
   * @param[in] baseElemsToNodes element to node map for the base mesh
   * @param[in] elemGhostRank rank of the ghost element
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverConstants constant part of the receiver term
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  Compute1DSourceAndReceiverConstants( localIndex const size,
                                       ArrayOfArraysView< localIndex const > const baseFacesToNodes,
                                       arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const baseNodeCoords,
                                       arrayView1d< globalIndex const > const baseNodeLocalToGlobal,
                                       arrayView1d< globalIndex const > const elementLocalToGlobal,
                                       ArrayOfArraysView< localIndex const > const baseNodesToElements,
                                       arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes,
                                       arrayView1d< integer const > const elemGhostRank,
                                       arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
                                       arrayView2d< localIndex const > const elemsToFaces,
                                       arrayView2d< real64 const > const & elemCenter,
                                       arrayView2d< real64 const > const sourceCoordinates,
                                       arrayView1d< localIndex > const sourceIsAccessible,
                                       arrayView2d< localIndex > const sourceNodeIds,
                                       arrayView2d< real64 > const sourceConstants,
                                       arrayView2d< real64 const > const receiverCoordinates,
                                       arrayView1d< localIndex > const receiverIsLocal,
                                       arrayView2d< localIndex > const receiverNodeIds,
                                       arrayView2d< real64 > const receiverConstants )
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
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );
          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};


            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
                                                                              coordsOnRefElem );

            sourceIsAccessible[isrc] = 1;
            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes( k, a );
              sourceConstants[isrc][a] = Ntest[a];
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
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );

          if( receiverFound && elemGhostRank[k] < 0 )
          {
            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
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


  /**
   * @brief Launches the precomputation of the source and receiver terms with storage of elements and region
   *        in which the receivers and sources are located
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] baseFacesToNodes face to node map of the base mesh
   * @param[in] baseNodeCoords coordinates of the nodes of the base mesh
   * @param[in] baseNodeLocalToGlobal local to global index map for nodes of the base mesh
   * @param[in] elementLocalToGlobal local to global index map for elements (for the base or high order mesh)
   * @param[in] baseNodesToElements local node to element map for the base mesh
   * @param[in] baseElemsToNodes element to node map for the base mesh
   * @param[in] elemGhostRank rank of the ghost element
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceElem element where a source is located
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverElem element where a receiver is located
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverConstants constant part of the receiver term
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  Compute1DSourceAndReceiverConstantsWithElementsAndRegionStorage( localIndex const size,
                                                                   localIndex const regionIndex,
                                                                   ArrayOfArraysView< localIndex const > const baseFacesToNodes,
                                                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const baseNodeCoords,
                                                                   arrayView1d< globalIndex const > const baseNodeLocalToGlobal,
                                                                   arrayView1d< globalIndex const > const elementLocalToGlobal,
                                                                   ArrayOfArraysView< localIndex const > const baseNodesToElements,
                                                                   arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes,
                                                                   arrayView1d< integer const > const elemGhostRank,
                                                                   arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
                                                                   arrayView2d< localIndex const > const elemsToFaces,
                                                                   arrayView2d< real64 const > const & elemCenter,
                                                                   arrayView2d< real64 const > const sourceCoordinates,
                                                                   arrayView1d< localIndex > const sourceIsAccessible,
                                                                   arrayView1d< localIndex > const sourceElem,
                                                                   arrayView2d< localIndex > const sourceNodeIds,
                                                                   arrayView2d< real64 > const sourceConstants,
                                                                   arrayView1d< localIndex > const sourceRegion,
                                                                   arrayView2d< real64 const > const receiverCoordinates,
                                                                   arrayView1d< localIndex > const receiverIsLocal,
                                                                   arrayView1d< localIndex > const receiverElem,
                                                                   arrayView2d< localIndex > const receiverNodeIds,
                                                                   arrayView2d< real64 > const receiverConstants,
                                                                   arrayView1d< localIndex > const receiverRegion )
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
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );
          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};

            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
                                                                              coordsOnRefElem );

            sourceIsAccessible[isrc] = 1;
            sourceElem[isrc] = k;
            sourceRegion[isrc] = regionIndex;
            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
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
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );

          if( receiverFound && elemGhostRank[k] < 0 )
          {
            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
                                                                              coordsOnRefElem );
            receiverIsLocal[ircv] = 1;
            receiverElem[ircv] = k;
            receiverRegion[ircv] = regionIndex;

            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][a] = elemsToNodes[k][a];
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // end loop over receivers

    } );

  }



  /**
   * @brief Launches the precomputation of the source and receiver terms for 3D arrays solution and DAS receiver constants
   *        computation
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] baseFacesToNodes face to node map
   * @param[in] baseNodeCoords coordinates of the nodes
   * @param[in] baseNodeLocalToGlobal local to global index map for nodes
   * @param[in] elementLocalToGlobal local to global index map for elements
   * @param[in] baseNodesToElements node to element map for the base mesh
   * @param[in] baseElemsToNodes element to node map for the base mesh
   * @param[in] elemGhostRank array containing the ghost rank
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceConstantsx constant part of the source terms in x-direction
   * @param[out] sourceConstantsy constant part of the source terms in y-direction
   * @param[out] sourceConstantsz constant part of the source terms in z-direction
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverConstants constant part of the receiver term
   * @param[in] useDAS parameter that determines which kind of receiver needs to be modeled (DAS or not, and which type)
   * @param[in] linearDASSamples parameter that gives the number of integration points to be used when computing the DAS signal via strain
   * integration
   * @param[in] linearDASGeometry geometry of the linear DAS receivers, if needed
   * @param[in] sourceForce force vector of the source
   * @param[in] sourceMoment moment (symmetric rank-2 tensor) of the source
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  Compute3DSourceAndReceiverConstantsWithDAS( localIndex const size,
                                              ArrayOfArraysView< localIndex const > const baseFacesToNodes,
                                              arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const baseNodeCoords,
                                              arrayView1d< globalIndex const > const baseNodeLocalToGlobal,
                                              arrayView1d< globalIndex const > const elementLocalToGlobal,
                                              ArrayOfArraysView< localIndex const > const baseNodesToElements,
                                              arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes,
                                              arrayView1d< integer const > const elemGhostRank,
                                              arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
                                              arrayView2d< localIndex const > const elemsToFaces,
                                              arrayView2d< real64 const > const & elemCenter,
                                              arrayView2d< real64 const > const sourceCoordinates,
                                              arrayView1d< localIndex > const sourceIsAccessible,
                                              arrayView2d< localIndex > const sourceNodeIds,
                                              arrayView2d< real64 > const sourceConstantsx,
                                              arrayView2d< real64 > const sourceConstantsy,
                                              arrayView2d< real64 > const sourceConstantsz,
                                              arrayView2d< real64 const > const receiverCoordinates,
                                              arrayView1d< localIndex > const receiverIsLocal,
                                              arrayView2d< localIndex > const receiverNodeIds,
                                              arrayView2d< real64 > const receiverConstants,
                                              WaveSolverUtils::DASType useDAS,
                                              integer linearDASSamples,
                                              arrayView2d< real64 const > const linearDASGeometry,
                                              R1Tensor const sourceForce,
                                              R2SymTensor const sourceMoment )
  {
    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
    integer nSamples = useDAS == WaveSolverUtils::DASType::none ? 1 : linearDASSamples;
    array1d< real64 > const samplePointLocationsA( nSamples );
    arrayView1d< real64 > const samplePointLocations = samplePointLocationsA.toView();
    array1d< real64 > const sampleIntegrationConstantsA( nSamples );
    arrayView1d< real64 > const sampleIntegrationConstants = sampleIntegrationConstantsA.toView();

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

          real64 xLocal[8][3];

          for( localIndex a = 0; a < 8; ++a )
          {
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = baseNodeCoords( baseElemsToNodes( k, a ), i );
            }
          }


          bool const sourceFound =
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );

          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};


            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
                                                                              coordsOnRefElem );
            sourceIsAccessible[isrc] = 1;

            real64 N[numNodesPerElem];
            real64 gradN[numNodesPerElem][3];
            FE_TYPE::calcN( coordsOnRefElem, N );
            FE_TYPE::calcGradNWithCorners( coordsOnRefElem, xLocal, gradN );
            R2SymTensor moment = sourceMoment;
            for( localIndex q=0; q< numNodesPerElem; ++q )
            {
              real64 inc[3] = { 0, 0, 0 };
              sourceNodeIds[isrc][q] = elemsToNodes( k, q );
              inc[0] += sourceForce[0] * N[q];
              inc[1] += sourceForce[1] * N[q];
              inc[2] += sourceForce[2] * N[q];

              LvArray::tensorOps::Ri_add_symAijBj< 3 >( inc, moment.data, gradN[q] );
              sourceConstantsx[isrc][q] += inc[0];
              sourceConstantsy[isrc][q] += inc[1];
              sourceConstantsz[isrc][q] += inc[2];
            }

          }
        }
      } // end loop over all sources

      // Step 2: locate the receivers, and precompute the receiver term

      // for geophones, we need only a point per receiver.
      // for DAS, we need multiple points

      /// compute locations of samples along receiver
      if( nSamples == 1 )
      {
        samplePointLocations[ 0 ] = 0;
      }
      else
      {
        for( integer i = 0; i < nSamples; ++i )
        {
          samplePointLocations[ i ] = -0.5 + (real64) i / ( linearDASSamples - 1 );
        }
      }

      /// compute integration constants of samples
      /// for displacement difference (dipole) DAS, take the discrete derivative of the pair of geophones
      if( useDAS == WaveSolverUtils::DASType::dipole )
      {
        sampleIntegrationConstants[ 0 ] = -1.0;
        sampleIntegrationConstants[ 1 ] = 1.0;
      }
      /// for strain integration DAS, take the average of strains to average strain data
      else if( nSamples == 1 )
      {
        sampleIntegrationConstants[ 0 ] = 1.0;
      }
      else
      {
        for( integer i = 0; i < linearDASSamples; i++ )
        {
          sampleIntegrationConstants[ i ] = 1.0 / nSamples;
        }
      }

      /// loop over all the receivers
      for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
      {
        R1Tensor receiverCenter = { receiverCoordinates[ ircv ][ 0 ], receiverCoordinates[ ircv ][ 1 ], receiverCoordinates[ ircv ][ 2 ] };
        R1Tensor receiverVector;
        if( useDAS == WaveSolverUtils::DASType::none )
        {
          receiverVector = { 0, 0, 0 };
        }
        else
        {
          receiverVector =  WaveSolverUtils::computeDASVector( linearDASGeometry[ ircv ][ 0 ], linearDASGeometry[ ircv ][ 1 ] );
        }
        real64 receiverLength = useDAS == WaveSolverUtils::DASType::none ? 0 : linearDASGeometry[ ircv ][ 2 ];
        /// loop over samples
        for( integer iSample = 0; iSample < nSamples; ++iSample )
        {
          /// compute sample coordinates and locate the element containing it
          real64 const coords[3] = { receiverCenter[ 0 ] + receiverVector[ 0 ] * receiverLength * samplePointLocations[ iSample ],
                                     receiverCenter[ 1 ] + receiverVector[ 1 ] * receiverLength * samplePointLocations[ iSample ],
                                     receiverCenter[ 2 ] + receiverVector[ 2 ] * receiverLength * samplePointLocations[ iSample ] };
          bool const sampleFound =
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );
          if( sampleFound && elemGhostRank[k] < 0 )
          {
            real64 coordsOnRefElem[3]{};
            real64 xLocal[8][3];

            for( localIndex a = 0; a < 8; ++a )
            {
              for( localIndex i=0; i < 3; ++i )
              {
                xLocal[a][i] = baseNodeCoords( baseElemsToNodes( k, a ), i );
              }
            }

            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
                                                                              coordsOnRefElem );
            real64 N[numNodesPerElem];
            real64 gradN[numNodesPerElem][3];
            FE_TYPE::calcN( coordsOnRefElem, N );
            FE_TYPE::calcGradNWithCorners( coordsOnRefElem, xLocal, gradN );
            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][iSample * numNodesPerElem + a] = elemsToNodes( k,
                                                                                   a );
              if( useDAS == WaveSolverUtils::DASType::strainIntegration )
              {
                receiverConstants[ircv][iSample * numNodesPerElem + a] += ( gradN[a][0] * receiverVector[0] + gradN[a][1] * receiverVector[1] + gradN[a][2] * receiverVector[2] ) *
                                                                          sampleIntegrationConstants[ iSample ];
              }
              else
              {
                receiverConstants[ircv][iSample * numNodesPerElem + a] += N[a] * sampleIntegrationConstants[ iSample ];
              }
            }
            receiverIsLocal[ ircv ] = 2;
          }
        } // end loop over samples
        // determine if the current rank is the owner of this receiver
        real64 const coords[3] = { receiverCenter[ 0 ], receiverCenter[ 1 ], receiverCenter[ 2 ] };
        bool const receiverFound =
          computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                      baseNodeCoords,
                                                                      elemsToFaces,
                                                                      baseFacesToNodes,
                                                                      baseNodesToElements,
                                                                      baseNodeLocalToGlobal,
                                                                      elementLocalToGlobal,
                                                                      center,
                                                                      coords );
        if( receiverFound && elemGhostRank[k] < 0 )
        {
          receiverIsLocal[ ircv ] = 1;
        }
      } // end loop over receivers
    } );

  }

};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_PRECOMPUTESOURCESANDRECEIVERSKERNEL_HPP_
