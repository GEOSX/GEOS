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

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverUtils.hpp"



namespace geos
{

/// Namespace to contain the first order acoustic wave kernels.
namespace AcousticWaveEquationDGKernels
{

struct PrecomputeSourceAndReceiverKernel
{

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
   * @param[out] rcvElem element where a receiver is located
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
          localIndex const numNodesPerElem,
          localIndex const numFacesPerElem,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView1d< integer const > const elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const faceNormal,
          arrayView2d< real64 const > const faceCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsAccessible,
          arrayView1d< localIndex > const sourceElem,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView1d< localIndex > const rcvElem,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real32 > const sourceValue,
          real64 const dt,
          real32 const timeSourceFrequency,
          localIndex const rickerOrder )
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
                                                                              X,
                                                                              coordsOnRefElem );

            sourceIsAccessible[isrc] = 1;
            sourceElem[isrc] = k;
            real64 Ntest[FE_TYPE::numNodes];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }

            for( localIndex cycle = 0; cycle < sourceValue.size( 0 ); ++cycle )
            {
              real64 const time = cycle*dt;
              sourceValue[cycle][isrc] = WaveSolverUtils::evaluateRicker( time, timeSourceFrequency, rickerOrder );
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
                                                                              X,
                                                                              coordsOnRefElem );
            receiverIsLocal[ircv] = 1;
            rcvElem[ircv] = k;

            real64 Ntest[FE_TYPE::numNodes];
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
          localindex const numFacesPerElem
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
          arrayView2d< real32 const > const p_n,
          arrayView2d< real32 const > const p_nm1,
          arrayView2d< real64 const > const sourceConstants,
          arrayView2d< real32 const > const sourceValue,
          arrayView1d< localIndex const > const sourceIsAccessible,
          arrayView1d< localIndex const > const sourceElem,
          real64 const dt,
          integer const cycleNumber,
          arrayView2d< real32 > const p_np1 )

  {

    //For now lots of comments with ideas  + needed array to add to the method prototype
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

      real64 xLocal[numNodesPerElem][3];
      for( localIndex a=0; a< numNodesPerElem; ++a )
      {
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = X( elemsToNodes( k, a ), i );
        }
      }

      real32 flow[numNodesPerElem]  = {0.0};


      // Volume  + fluxes computation integration
      for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
      {


        //Stiffness terms
        m_finiteElement.template computeStiffnessTerm( q, xLocal, [&] ( int i, int j, real64 val )
        {
          //Maybe reverse j and i
          flow[i] += val * p_n[k][j];
        } );

        //Fluxes
        for( localIndex f = 0; f < numFacesPerElem; ++f )
        {
          //Possible way:
          //Get the global number of face using elemeToFaces :
          localIndex face_glob = elemToFaces[k][f];
          //Use faceToElemIndex map to know which element shared this global face: faceToElemIndex is a 2d array which knowing a face and a
          // index between 0 and 1 can give you the two
          // element which share the face and if you get -1 it means that the element will be in the boundary.
          // Initialize the storage value for contributions: 
          real32 fp = 0.0;
          for( localIndex m = 0; m < 2; ++m )
          {
            localIndex elem = faceToElemIndex[face_glob][m];
            //We start by the test on the boundaries to skip it directly:
            if( elem == -1 )
            {
              //Nothing we continue the loop
            }
            else if( elem == k )
            {
              //Here we compute the fluxes part corresponding to the element itself (the (K,K) part seen in the latex document). We can both
              // compute the "classical" flux part + the penalization one:
              //PS: Not sure about how to include the normals so I'll just put "normals" (surely missing something with the gradient inside
              // the flux matrix)
              // Inside the matrix computation: we need the volumic Jacobian (its inverse) and the surface determinant. Due to the
              // fact that we take the inverse of the jacobian
              // we will have the ratio surface/volume.

              m_finiteElement.template computeKKFluxMatrix(q,xLocal,f [&] (int i, int j, real64 val)
              {
                fp += 0.5* val *  p_n[k][i] + gamma[k]* val * p_n[k][i];
              //Needs normals at some point 
              //flow[j] += fp*faceNormal
              } );
            }
            else
            {
              //It remains the case where we look at the neighbour element and we need to add the (K,L) contribution.
              //It will be transparent here, but inside the mathematical computation we need to be careful on which degrees of freedom we
              // send back for the pressure as we get the contribution
              //of the neighbour so we need to get the correct dof (can be taken in account inside the math stuff)
              m_finiteElement.template computeKLFluxMatrix(q, xLocal,f [&] (int i, int j, real32 val)
              {
                fp += 0.5* val * p_n[k][i] - gamma[k]* val * p_n[elem][j];
                //Again, needs normals at some point
                //flow[j] += fp*normals
              } );
            }
          }

        }



      }

    } );


  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

} // namespace AcousticWaveEquationDGKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONDGKERNEL_HPP_
