/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file AcousticFirstOrderWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONDGKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONDGKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"

namespace geos
{

/// Namespace to contain the first order acoustic wave kernels.
namespace AcousticWaveEquationDGKernels
{

struct PrecomputeSourceAndReceiverKernel
{
  using EXEC_POLICY = parallelDevicePolicy< >;

  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] numNodesPerElem number of nodes per element
   * @param[in] numFacesPerElem number of faces per element
   * @param[in] X coordinates of the nodes
   * @param[in] elemGhostRank rank of the ghost element
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] faceNormal normal of each faces
   * @param[in] faceCenter coordinates of the center of a face
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
   * @param[out] sourceValue value of the temporal source (eg. Ricker)
   * @param[in] dt time-step
   * @param[in] timeSourceFrequency the central frequency of the source
   * @param[in] rickerOrder order of the Ricker wavelet
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
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
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView1d< localIndex > const receiverElem,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants )
  {

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
            sourceIsAccessible[isrc] = 1;
            sourceElem[isrc] = k;
            real64 Ntest[FE_TYPE::numNodes];
            constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

            real64 xLocal[4][3];
            for( localIndex a=0; a<4; ++a )
            {
              for( localIndex i=0; i<3; ++i )
              {
                xLocal[a][i] = baseNodeCoords( baseElemsToNodes( k, a ), i );
              }
            }
            FE_TYPE::calcN( xLocal, coords, Ntest );

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
            receiverIsLocal[ircv] = 1;
            receiverElem[ircv] = k;

            real64 Ntest[FE_TYPE::numNodes];
            constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

            real64 xLocal[4][3];
            for( localIndex a=0; a< 4; ++a )
            {
              for( localIndex i=0; i<3; ++i )
              {
                xLocal[a][i] = baseNodeCoords( baseElemsToNodes( k, a ), i );
              }
            }
            FE_TYPE::calcN( xLocal, coords, Ntest );

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
};

struct PrecomputeNeighborhoodKernel
{

  using EXEC_POLICY = parallelDevicePolicy< >;

  /**
   * @brief Launches the precomputation of the element neighborhood information needed by DG
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in]
   * @param[out]
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const & elemsToFaces,
          arrayView2d< localIndex const > const & facesToElems,
          ArrayOfArraysView< localIndex const > const & facesToNodes,
          arrayView1d< localIndex const > const & freeSurfaceFaceIndicator,
          arrayView2d< localIndex > const & elemsToOpposite,
          arrayView2d< integer > const & elemsToOppositePermutation )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k1 )
    {
      printf("sizeeTN=%d\n",elemsToNodes.size());
      localIndex vertices[ 4 ] = { elemsToNodes( k1, 0 ), elemsToNodes( k1, 1 ), elemsToNodes( k1, 2 ), elemsToNodes( k1, 3 ) };
      printf("vertices0=%d\n",elemsToNodes( k1, 0 ));
      printf("vertices1=%d\n",elemsToNodes( k1, 1 ));
      printf("vertices2=%d\n",elemsToNodes( k1, 2 ));
      printf("vertices3=%d\n",elemsToNodes( k1, 3 ));
      for( int i = 0; i < 4; i++ )
      {
        localIndex  k1OrderedVertices[ 3 ];
        localIndex f = elemsToFaces( k1, i );
        localIndex faceVertices[ 3 ] = { facesToNodes( f, 0 ), facesToNodes( f, 1 ), facesToNodes( f, 2 ) };
         printf("index=%d\n",i);
        printf("face=%d\n",f);
        printf("sizefacetonodes=%d\n",facesToNodes.size());
        printf("faceVertices0=%d\n",facesToNodes( f, 0 ));
        printf("faceVertices1=%d\n",facesToNodes( f, 1 ));
        printf("faceVertices2=%d\n",facesToNodes( f, 2 ));
        // find neighboring element, if any
        localIndex k2 = facesToElems( f, 0 );
        if( k2 == k1 )
        {
          k2 = facesToElems( f, 1 );
        }
        // find opposite vertex in first element
        int o1 = -1;
        int indexo1= -1;
        int vertex = -1;
        int count = 0;
        for ( localIndex k=0; k< 4; ++k) {
          vertex = vertices[k];
          bool found = false;
          for ( int j = 0; j < 3; j++ )
          {
            if ( vertex == faceVertices[ j ] )
            {
              found = true;
              break;
            }
          }
          if( !found )
          {
            o1 = vertex;
            indexo1=k;
            printf("o1=%d\n",o1);
          }
          else
          {
            k1OrderedVertices[ count++ ] = vertex;
          }
        }
       
        GEOS_ERROR_IF( o1 < 0, "Topological error in mesh: a face and its adjacent element share all vertices.");
        if( k2 < 0 )
        {
          // boundary element, either free surface, or absorbing boundary
          elemsToOpposite( k1, indexo1 ) = freeSurfaceFaceIndicator( f ) == 1 ? -2 : -1;
          elemsToOppositePermutation( k1, indexo1 ) = 0;
        }
        else
        {
          printf("là\n");
          printf("k1=%d\n",k1);
          printf("o1=%d\n",o1);
          printf("k2=%d\n",k2);
          elemsToOpposite( k1, indexo1 ) = k2;
          printf("apres");
          localIndex oppositeElemVertices[ 4 ] = { elemsToNodes( k2, 0 ), elemsToNodes( k2, 1 ), elemsToNodes( k2, 2 ), elemsToNodes( k2, 3 ) };
          // find opposite vertex in second element
          int o2 = -1;
          int indexo2 = -1;
          count = 0;
          for ( localIndex k=0; k<4; ++k) {
            vertex = vertices[k];
            bool found = false;
            for ( int j = 0; j < 3; j++ )
            {
              if ( vertex == faceVertices[ j ] )
              {
                found = true;
                break;
              }
            }
            if( !found )
            {
              o2 = vertex;
              indexo2 = k;
            }

          }
          GEOS_ERROR_IF( o2 < 0, "Topological error in mesh: a face and its adjacent element share all vertices.");
          // compute permutation
          integer permutation = 0;
          int c = 1;
          for (localIndex k2Vertex : oppositeElemVertices )
          {
            int position = -1;
            for( int j = 0; j < 3; j++ )
            {
              if( k1OrderedVertices[ j ] == k2Vertex )
              {
                position = j;
                break;
              }
            }
            permutation = permutation + c * ( position + 1 );
            c = c * 4;
          }
          elemsToOppositePermutation( k1, indexo1 ) = permutation;
        }
      }
    } );
  }
};

template< typename FE_TYPE >
struct PressureComputation
{
 PressureComputation( FE_TYPE const & finiteElement )
   : m_finiteElement( finiteElement )
 {}

 /**
  * @brief Launches the computation of the pressure for one iteration
  * @tparam EXEC_POLICY the execution policy
  * @tparam ATOMIC_POLICY the atomic policy
  * @param[in] size the number of cells in the subRegion
  * @param[in] X coordinates of the nodes
  * @param[in] p_nm1 pressure  array at time n-1 (only used here)
  * @param[in] p_n pressure array at time n (only used here)
  * @param[in] sourceConstants constant part of the source terms
  * @param[in] sourceValue value of the temporal source (eg. Ricker)
  * @param[in] sourceIsAccessible flag indicating whether the source is accessible or not
  * @param[in] sourceElem element where a source is located
  * @param[in] cycleNumber the number of cycle
  * @param[in] dt time-step
  * @param[out] p_np1 pressure array at time n+1 (updated here)
  */
 //List is not complete, it will need several GEOS maps to add
 template< typename EXEC_POLICY, typename ATOMIC_POLICY >
 void
 launch( localIndex const size,
         localIndex const regionIndex,
         arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
         arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
         arrayView2d< real32 const > const p_n,
         arrayView2d< real32 const > const p_nm1,
         arrayView2d< localIndex > const & elemsToOpposite,
         arrayView2d< integer > const & elemsToOppositePermutation,
         arrayView2d< real64 const > const sourceConstants,
         arrayView1d< localIndex const > const sourceIsAccessible,
         arrayView1d< localIndex const > const sourceElem,
         arrayView1d< localIndex const > const sourceRegion,
         real64 const dt,
         real64 const time_n,
         real32 const timeSourceFrequency,
         real32 const timeSourceDelay,
         localIndex const rickerOrder,
         bool const useSourceWaveletTables,
         arrayView1d< TableFunction::KernelWrapper const > const sourceWaveletTableWrappers,
         arrayView2d< real32 > const p_np1 )

 {

    real64 const rickerValue = useSourceWaveletTables ? 0 : WaveSolverUtils::evaluateRicker( time_n, timeSourceFrequency, timeSourceDelay, rickerOrder );

    //For now lots of comments with ideas  + needed array to add to the method prototype
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      real64 const dt2 = pow(dt,2);

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real32 flowx[numNodesPerElem] = {0.0};

      real64 xLocal[4][3];
      for( localIndex a=0; a< 4; ++a )
      {
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      //Mass matrix declaration and computation
      real32 mass[numNodesPerElem][numNodesPerElem];

      //Multiply by p_{n } by 2*Mass
      m_finiteElement.template computeMassTerm(xLocal, [&] (const int i, const int j, const real64 val)
      {
         flowx[j] = 2.0*val*p_n[k][i];
         flowx[j] -= val*p_nm1[k][i];
      } );




      //First stiffness part (volume)
      m_finiteElement.template computeStiffnessTerm(xLocal, [&] (const int i, const int j, real64 val)
      {
         flowx[j] -= dt2*val*p_n[k][i];
      } );


      // m_finiteElement.template computeSurfaceTerms(xLocal, [&] (const int c1, const int c2, const int f1, const int , const int , const int ,const int i2, const int j2, const int k2, real64 val)
      // {
      //   //We take the neighbour element
      //   const localIndex elemNeigh = elemsToOpposite(k,f1);

      //   // Now we seek the degree of freedom on the neighbour element to use for the computation of the flux (or the penalty)
      //   // First, we compute the four possible values of the permutation of the degrees of freedom depending on the the fixed
      //   // permutation value contained inside elemsToOppositePermutation permutation

      //   const int perm = elemsToOppositePermutation(elemNeigh,f1);

      //   const int p1 = perm%4-1;
      //   const int p2 = (perm/4)%4-1;
      //   const int p3 = (perm/16)%4-1;
      //   const int p4 = (perm/64)-1;

      //   // Then we transform the 3 indices returned by the callback (i2,j2,k2) using the permutations. One of this permutation, will be 0 (depending on which
      //   // degree of freedom is the one at the opposite of the face shared with the neighbour element) and will correspond to the one where p* will be negative

      //   const int Indices[3] = {i2,j2,k2};

      //   const int ii2 = p1 < 0 ? 0 : Indices[p1];
      //   const int jj2 = p2 < 0 ? 0 : Indices[p2];
      //   const int kk2 = p3 < 0 ? 0 : Indices[p3];

      //   // Finally, using the dofIndex function, we compute the number of the global degree of freedom on the element

      //   const int neighDof = m_finiteElement.dofIndex(ii2,jj2,kk2);

      //   //Flux computation

      //   flowx[c1] -= 0.5*dt2*p_n[k][c2];
      //   flowx[c1] += 0.5*dt2*p_n[elemNeigh][neighDof];


      // },
      // [&] (const int c1, const int c2, const int f1, const int i1, const int j1, const int k1, const int i2, const int j2, const int k2, real64 val)
      // {
      //   //We take the neighbour element
      //    const int elemNeigh = elemsToOpposite(k,f1);

      //   // Now we seek the degree of freedom on the neighbour element to use for the computation of the flux (or the penalty)
      //   // First, we compute the four possible values of the permutation of the degrees of freedom depending on the the fixed
      //   // permutation value contained inside elemsToOppositePermutation permutation

      //   const int perm = elemsToOppositePermutation(elemNeigh,f1);

      //   const int p1 = perm%4-1;
      //   const int p2 = (perm/4)%4-1;
      //   const int p3 = (perm/16)%4-1;

      //   // Then we transform the 3 indices returned by the callback (i2,j2,k2) using the permutations. One of this permutation, will be 0 (depending on which
      //   // degree of freedom is the one at the opposite of the face shared with the neighbour element) and will correspond to the one where p* will be negative

      //   const int Indices[3] = {i2,j2,k2};

      //   const int ii2 = p1 < 0 ? 0 : Indices[p1];
      //   const int jj2 = p2 < 0 ? 0 : Indices[p2];
      //   const int kk2 = p3 < 0 ? 0 : Indices[p3];

      //   // Finally, using the dofIndex function, we compute the number of the global degree of freedom on the element


      //   const int neighDof = m_finiteElement.dofIndex(ii2,jj2,kk2);

      //   //Flux computation

      //   flowx[c1] += 0.5*dt2*p_n[elemNeigh][neighDof];
      //   flowx[c1] -= 0.5*dt2*p_n[k][c2];

      //   //Then we need a second time where we take the transpose of the previous values:


      //   const int IndicesTranspose[3] = {i1,j1,k1};

      //   const int ii1 = p1 < 0 ? 0 : IndicesTranspose[p1];
      //   const int jj1 = p2 < 0 ? 0 : IndicesTranspose[p2];
      //   const int kk1 = p3 < 0 ? 0 : IndicesTranspose[p3];

      //   // Finally, using the dofIndex function, we compute the number of the global degree of freedom on the element

      //   const int neighDof2 = m_finiteElement.dofIndex(ii1,jj1,kk1);

      //   //Flux computation

      //   flowx[c2] -= 0.5*dt2*p_n[elemNeigh][neighDof2];
      //   flowx[c2] += 0.5*dt2*p_n[k][c1];


      // } );

      //Source Injection
      for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
      {
        if( sourceIsAccessible[isrc] == 1 )
        {
          if( sourceElem[isrc]==k && sourceRegion[isrc] == regionIndex )
          {
            real64 const srcValue = useSourceWaveletTables ? sourceWaveletTableWrappers[ isrc ].compute( &time_n ) : rickerValue;
            for( localIndex i = 0; i < numNodesPerElem; ++i )
            {
              flowx[i]+= dt2*(sourceConstants[isrc][i]*srcValue);
            }
          }
        }
      }

    } );

  //    //Mass and Damping matrix declaration
  //    real64 mass[numNodesPerelem][numNodesPerElem];
  //    //real64 damping[numNodesPerelem][numNodesPerElem];


  //    //Computation of the mass and damping matrix
  //    computeMassMatrix(mass);

  //    //First steps for time scheme
  //    for (localindex i = 0; i < numNodesPerElem; ++i)
  //    {
  //     p_np1[]
  //    }





  //    real32 flow[numNodesPerElem]  = {0.0};


  //    // Volume  + fluxes computation integration
  //    for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
  //    {


  //      //Stiffness terms
  //      m_finiteElement.template computeStiffnessTerm( q, xLocal, [&] ( int i, int j, real64 val )
  //      {
  //        //Maybe reverse j and i
  //        flow[i] += val * p_n[k][j];
  //      } );

  //      //Fluxes
  //      for( localIndex f = 0; f < numFacesPerElem; ++f )
  //      {
  //        //Possible way:
  //        //Get the global number of face using elemeToFaces :
  //        localIndex face_glob = elemToFaces[k][f];
  //        //Use faceToElemIndex map to know which element shared this global face: faceToElemIndex is a 2d array which knowing a face and a
  //        // index between 0 and 1 can give you the two
  //        // element which share the face and if you get -1 it means that the element will be in the boundary.
  //        // Initialize the storage value for contributions:
  //        real32 fp = 0.0;
  //        for( localIndex m = 0; m < 2; ++m )
  //        {
  //          localIndex elem = faceToElemIndex[face_glob][m];
  //          //We start by the test on the boundaries to skip it directly:
  //          if( elem == -1 )
  //          {
  //            //Nothing we continue the loop
  //          }
  //          else if( elem == k )
  //          {
  //            //Here we compute the fluxes part corresponding to the element itself (the (K,K) part seen in the latex document). We can both
  //            // compute the "classical" flux part + the penalization one:
  //            //PS: Not sure about how to include the normals so I'll just put "normals" (surely missing something with the gradient inside
  //            // the flux matrix)
  //            // Inside the matrix computation: we need the volumic Jacobian (its inverse) and the surface determinant. Due to the
  //            // fact that we take the inverse of the jacobian
  //            // we will have the ratio surface/volume.

  //            m_finiteElement.template computeKKFluxMatrix(q,xLocal,f [&] (int i, int j, real64 val)
  //            {
  //              fp += 0.5* val *  p_n[k][i] + gamma[k]* val * p_n[k][i];
  //            //Needs normals at some point
  //            //flow[j] += fp*faceNormal
  //            } );
  //          }
  //          else
  //          {
  //            //It remains the case where we look at the neighbour element and we need to add the (K,L) contribution.
  //            //It will be transparent here, but inside the mathematical computation we need to be careful on which degrees of freedom we
  //            // send back for the pressure as we get the contribution
  //            //of the neighbour so we need to get the correct dof (can be taken in account inside the math stuff)
  //            m_finiteElement.template computeKLFluxMatrix(q, xLocal,f [&] (int i, int j, real32 val)
  //            {
  //              fp += 0.5* val * p_n[k][i] - gamma[k]* val * p_n[elem][j];
  //              //Again, needs normals at some point
  //              //flow[j] += fp*normals
  //            } );
  //          }
  //        }

  //      }



  //    }

  //  } );


 }

 /// The finite element space/discretization object for the element type in the subRegion
 FE_TYPE const & m_finiteElement;


};

} // namespace AcousticWaveEquationDGKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONDGKERNEL_HPP_
