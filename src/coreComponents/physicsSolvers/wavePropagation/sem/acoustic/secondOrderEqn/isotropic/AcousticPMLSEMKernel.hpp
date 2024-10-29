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
 * @file AcousticPMLSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICPMLSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICPMLSEMKERNEL_HPP_

namespace geos
{

struct AcousticPMLSEM
{
  struct ComputeDamping
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
    static void computeDampingProfile( real32 const (&xLocal)[3],
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
          ComputeDamping::computeDampingProfile(
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


};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICPMLSEMKERNEL_HPP_
