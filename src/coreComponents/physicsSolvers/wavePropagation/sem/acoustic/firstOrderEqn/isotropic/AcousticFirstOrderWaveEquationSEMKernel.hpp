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
 * @file AcousticFirstOrderWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"

namespace geos
{

/// Namespace to contain the first order acoustic wave kernels.
namespace acousticFirstOrderWaveEquationSEMKernels
{

template< typename FE_TYPE >
struct VelocityComputation
{

  VelocityComputation( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of the velocity for one iteration
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] p_np1 pressure array (only used here)
   * @param[in] density cell-wise density
   * @param[in] dt time-step
   * @param[out] velocity_x velocity array in the x direction (updated here)
   * @param[out] velocity_y velocity array in the y direction (updated here)
   * @param[out] velocity_z velocity array in the z direction (updated here)
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const p_np1,
          arrayView1d< real32 const > const density,
          real64 const dt,
          arrayView2d< real32 > const velocity_x,
          arrayView2d< real32 > const velocity_y,
          arrayView2d< real32 > const velocity_z )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // only the eight corners of the mesh cell are needed to compute the Jacobian
      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      real32 uelemx[numNodesPerElem] = {0.0};
      real32 uelemy[numNodesPerElem] = {0.0};
      real32 uelemz[numNodesPerElem] = {0.0};
      real32 flowx[numNodesPerElem] = {0.0};
      real32 flowy[numNodesPerElem] = {0.0};
      real32 flowz[numNodesPerElem] = {0.0};

      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        uelemx[i] = massLoc*velocity_x[k][i];
        uelemy[i] = massLoc*velocity_y[k][i];
        uelemz[i] = massLoc*velocity_z[k][i];
      }

      for( localIndex q = 0; q < numNodesPerElem; q++ )
      {

        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          flowx[j] += dfx1*p_np1[elemsToNodes[k][i]];
          flowy[j] += dfx2*p_np1[elemsToNodes[k][i]];
          flowz[j] += dfx3*p_np1[elemsToNodes[k][i]];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          flowx[j] += dfy1*p_np1[elemsToNodes[k][i]];
          flowy[j] += dfy2*p_np1[elemsToNodes[k][i]];
          flowz[j] += dfy3*p_np1[elemsToNodes[k][i]];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          flowx[j] += dfz1*p_np1[elemsToNodes[k][i]];
          flowy[j] += dfz2*p_np1[elemsToNodes[k][i]];
          flowz[j] += dfz3*p_np1[elemsToNodes[k][i]];

        } );

      }
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        uelemx[i]+=dt*flowx[i]/density[k];
        uelemy[i]+=dt*flowy[i]/density[k];
        uelemz[i]+=dt*flowz[i]/density[k];

        velocity_x[k][i] = uelemx[i]/massLoc;
        velocity_y[k][i] = uelemy[i]/massLoc;
        velocity_z[k][i] = uelemz[i]/massLoc;
      }
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
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
   * @param[in] regionIndex Index of the subregion
   * @param[in] size_node the number of nodes in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[out] velocity_x velocity array in the x direction (only used here)
   * @param[out] velocity_y velocity array in the y direction (only used here)
   * @param[out] velocity_z velocity array in the z direction (only used here)
   * @param[in] mass the mass matrix
   * @param[in] damping the damping matrix
   * @param[in] sourceConstants constant part of the source terms
   * @param[in] sourceValue value of the temporal source (eg. Ricker)
   * @param[in] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[in] sourceElem element where a source is located
   * @param[in] cycleNumber the number of cycle
   * @param[in] dt time-step
   * @param[out] p_np1 pressure array (updated here)
   */

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          localIndex const regionIndex,
          localIndex const size_node,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView2d< real32 const > const velocity_x,
          arrayView2d< real32 const > const velocity_y,
          arrayView2d< real32 const > const velocity_z,
          arrayView1d< real32 const > const mass,
          arrayView1d< real32 const > const damping,
          arrayView2d< real64 const > const sourceConstants,
          arrayView2d< real32 const > const sourceValue,
          arrayView1d< localIndex const > const sourceIsAccessible,
          arrayView1d< localIndex const > const sourceElem,
          arrayView1d< localIndex const > const sourceRegion,
          real64 const dt,
          integer const cycleNumber,
          arrayView1d< real32 > const p_np1 )

  {

    //Pre-mult by the first factor for damping
    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_np1[a] *= 1.0-((dt/2)*(damping[a]/mass[a]));
    } );

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // only the eight corners of the mesh cell are needed to compute the Jacobian
      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      real32 auxx[numNodesPerElem]  = {0.0};
      real32 auyy[numNodesPerElem]  = {0.0};
      real32 auzz[numNodesPerElem]  = {0.0};
      real32 uelemx[numNodesPerElem] = {0.0};


      // Volume integration
      for( localIndex q=0; q < numNodesPerElem; q++ )
      {


        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          auxx[i] -= dfx1*velocity_x[k][j];
          auyy[i] -= dfx2*velocity_y[k][j];
          auzz[i] -= dfx3*velocity_z[k][j];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          auxx[i] -= dfy1*velocity_x[k][j];
          auyy[i] -= dfy2*velocity_y[k][j];
          auzz[i] -= dfy3*velocity_z[k][j];
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          auxx[i] -= dfz1*velocity_x[k][j];
          auyy[i] -= dfz2*velocity_y[k][j];
          auzz[i] -= dfz3*velocity_z[k][j];
        } );



      }

      //Time update + multiplication by inverse of the mass matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 diag=(auxx[i]+auyy[i]+auzz[i]);
        uelemx[i]+=dt*diag;

        real32 const localIncrement = uelemx[i]/mass[elemsToNodes[k][i]];
        RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[elemsToNodes[k][i]], localIncrement );
      }

      //Source Injection
      for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
      {
        if( sourceIsAccessible[isrc] == 1 )
        {
          if( sourceElem[isrc]==k && sourceRegion[isrc] == regionIndex )
          {
            for( localIndex i = 0; i < numNodesPerElem; ++i )
            {
              real32 const localIncrement2 = dt*(sourceConstants[isrc][i]*sourceValue[cycleNumber][isrc])/(mass[elemsToNodes[k][i]]);
              RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[elemsToNodes[k][i]], localIncrement2 );
            }
          }
        }
      }

    } );

    //Pre-mult by the first factor for damping
    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_np1[a] /= 1.0+((dt/2)*(damping[a]/mass[a]));
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

} // namespace AcousticFirstOrderWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIRSTTORDERWAVEEQUATIONSEMKERNEL_HPP_
