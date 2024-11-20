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
 * @file ElasticFirstOrderWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"


namespace geos
{

/// Namespace to contain the elastic wave kernels.
namespace elasticFirstOrderWaveEquationSEMKernels
{

template< typename FE_TYPE >
struct StressComputation
{

  StressComputation( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of the strain tensor for one iteration
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] regionIndex Index of the subregion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] ux_np1 displacement array in the x direction (only used here)
   * @param[in] uy_np1 displacement array in the y direction (only used here)
   * @param[in] uz_np1 displacement array in the z direction (only used here)
   * @param[in] velocityVp P-wavespeed array
   * @param[in] velocityVs S-wavespeed array
   * @param[in] lambda lambda (Lamé parameter) array
   * @param[in] mu mu (Lamé parameter) array
   * @param[in] sourceConstants constant part of the source terms
   * @param[in] sourceIsLocal flag indicating whether the source is accessible or not
   * @param[in] sourceElem element where a source is located
   * @param[in] sourceRegion region where the source is located
   * @param[in] dt time-step
   * @param[in] time_n current time
   * @param[in] timeSourceFrequency the central frequency of the source
   * @param[in] timeSourceDelay the time delay of the source
   * @param[in] rickerOrder order of the Ricker wavelet
   * @param[out] stressxx xx-component of the strain tensor array (updated here)
   * @param[out] stressyy yy-component of the strain tensor array (updated here)
   * @param[out] stresszz zz-component of the strain tensor array (updated here)
   * @param[out] stressxy xy-component of the strain tensor array (updated here)
   * @param[out] stressxz xz-component of the strain tensor array (updated here)
   * @param[out] stressyz yz-component of the strain tensor array (updated here)
   */

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          localIndex const regionIndex,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const ux_np1,
          arrayView1d< real32 const > const uy_np1,
          arrayView1d< real32 const > const uz_np1,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 const > const velocityVp,
          arrayView1d< real32 const > const velocityVs,
          arrayView1d< real32 > const lambda,
          arrayView1d< real32 > const mu,
          arrayView2d< real64 const > const sourceConstants,
          arrayView1d< localIndex const > const sourceIsLocal,
          arrayView1d< localIndex const > const sourceElem,
          arrayView1d< localIndex const > const sourceRegion,
          real64 const dt,
          real64 const time_n,
          real32 const timeSourceFrequency,
          real32 const timeSourceDelay,
          localIndex const rickerOrder,
          bool const useSourceWaveletTables,
          arrayView1d< TableFunction::KernelWrapper const > const sourceWaveletTableWrappers,
          arrayView2d< real32 > const stressxx,
          arrayView2d< real32 > const stressyy,
          arrayView2d< real32 > const stresszz,
          arrayView2d< real32 > const stressxy,
          arrayView2d< real32 > const stressxz,
          arrayView2d< real32 > const stressyz )

  {
    real64 const rickerValue = useSourceWaveletTables ? 0 : WaveSolverUtils::evaluateRicker( time_n, timeSourceFrequency, timeSourceDelay, rickerOrder );
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

      mu[k] = density[k] * pow( velocityVs[k], 2 );
      lambda[k] = density[k] * pow( velocityVp[k], 2 ) - 2.0*mu[k];

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      real32 uelemxx[numNodesPerElem] = {0.0};
      real32 uelemyy[numNodesPerElem] = {0.0};
      real32 uelemzz[numNodesPerElem] = {0.0};
      real32 uelemxy[numNodesPerElem] = {0.0};
      real32 uelemxz[numNodesPerElem] = {0.0};
      real32 uelemyz[numNodesPerElem]= {0.0};
      real32 auxx[numNodesPerElem] = {0.0};
      real32 auyy[numNodesPerElem] = {0.0};
      real32 auzz[numNodesPerElem] = {0.0};
      real32 auxy[numNodesPerElem] = {0.0};
      real32 auxz[numNodesPerElem] = {0.0};
      real32 auyz[numNodesPerElem] = {0.0};


      //Pre-multiplication by mass matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        uelemxx[i] = massLoc*stressxx[k][i];
        uelemyy[i] = massLoc*stressyy[k][i];
        uelemzz[i] = massLoc*stresszz[k][i];
        uelemxy[i] = massLoc*stressxy[k][i];
        uelemxz[i] = massLoc*stressxz[k][i];
        uelemyz[i] = massLoc*stressyz[k][i];
      }

      for( localIndex q = 0; q < numNodesPerElem; q++ )
      {

        //Volume integral
        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int const i, int const j, real32 const dfx1, real32 const dfx2, real32 const dfx3 )
        {
          localIndex const nodeIndex = elemsToNodes[k][i];
          auxx[j] = auxx[j] + dfx1*ux_np1[nodeIndex];
          auyy[j] = auyy[j] + dfx2*uy_np1[nodeIndex];
          auzz[j] = auzz[j] + dfx3*uz_np1[nodeIndex];
          auxy[j] = auxy[j] + dfx1*uy_np1[nodeIndex]+dfx2*ux_np1[nodeIndex];
          auxz[j] = auxz[j] + dfx1*uz_np1[nodeIndex]+dfx3*ux_np1[nodeIndex];
          auyz[j] = auyz[j] + dfx2*uz_np1[nodeIndex]+dfx3*uy_np1[nodeIndex];

        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int const i, int const j, real32 const dfy1, real32 const dfy2, real32 const dfy3 )
        {
          localIndex const nodeIndex = elemsToNodes[k][i];
          auxx[j] = auxx[j] + dfy1*ux_np1[nodeIndex];
          auyy[j] = auyy[j] + dfy2*uy_np1[nodeIndex];
          auzz[j] = auzz[j] + dfy3*uz_np1[nodeIndex];
          auxy[j] = auxy[j] + dfy1*uy_np1[nodeIndex]+dfy2*ux_np1[nodeIndex];
          auxz[j] = auxz[j] + dfy1*uz_np1[nodeIndex]+dfy3*ux_np1[nodeIndex];
          auyz[j] = auyz[j] + dfy2*uz_np1[nodeIndex]+dfy3*uy_np1[nodeIndex];

        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int const i, int const j, real32 const dfz1, real32 const dfz2, real32 const dfz3 )
        {
          localIndex const nodeIndex = elemsToNodes[k][i];
          auxx[j] = auxx[j] + dfz1*ux_np1[nodeIndex];
          auyy[j] = auyy[j] + dfz2*uy_np1[nodeIndex];
          auzz[j] = auzz[j] + dfz3*uz_np1[nodeIndex];
          auxy[j] = auxy[j] + dfz1*uy_np1[nodeIndex]+dfz2*ux_np1[nodeIndex];
          auxz[j] = auxz[j] + dfz1*uz_np1[nodeIndex]+dfz3*ux_np1[nodeIndex];
          auyz[j] = auyz[j] + dfz2*uz_np1[nodeIndex]+dfz3*uy_np1[nodeIndex];

        } );

      }
      //Time integration
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 diag = lambda[k]*(auxx[i]+auyy[i]+auzz[i]);
        uelemxx[i]+= dt*(diag+2*mu[k]*auxx[i]);
        uelemyy[i]+= dt*(diag+2*mu[k]*auyy[i]);
        uelemzz[i]+= dt*(diag+2*mu[k]*auzz[i]);
        uelemxy[i]+= dt*mu[k]*auxy[i];
        uelemxz[i]+= dt*mu[k]*auxz[i];
        uelemyz[i]+= dt*mu[k]*auyz[i];
      }

      // Multiplication by inverse mass matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
        stressxx[k][i] = uelemxx[i]/massLoc;
        stressyy[k][i] = uelemyy[i]/massLoc;
        stresszz[k][i] = uelemzz[i]/massLoc;
        stressxy[k][i] = uelemxy[i]/massLoc;
        stressxz[k][i] = uelemxz[i]/massLoc;
        stressyz[k][i] = uelemyz[i]/massLoc;
      }


      //Source injection
      for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
      {
        if( sourceIsLocal[isrc] == 1 )
        {
          if( sourceElem[isrc]==k && sourceRegion[isrc] == regionIndex )
          {
            real64 const srcValue = useSourceWaveletTables ? sourceWaveletTableWrappers[ isrc ].compute( &time_n ) : rickerValue;
            for( localIndex i = 0; i < numNodesPerElem; ++i )
            {
              real32 massLoc = m_finiteElement.computeMassTerm( i, xLocal );
              real32 const localIncrement = dt*(sourceConstants[isrc][i]*srcValue)/massLoc;
              RAJA::atomicAdd< ATOMIC_POLICY >( &stressxx[k][i], localIncrement );
              RAJA::atomicAdd< ATOMIC_POLICY >( &stressyy[k][i], localIncrement );
              RAJA::atomicAdd< ATOMIC_POLICY >( &stresszz[k][i], localIncrement );
            }

          }

        }
      }

    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion

  FE_TYPE const & m_finiteElement;
};

template< typename FE_TYPE >
struct VelocityComputation
{

  VelocityComputation( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the computation of the displacement for one iteration
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] size_node the number of nodes in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[out] stressxx xx-component of the strain tensor array (only used here)
   * @param[out] stressyy yy-component of the strain tensor array (only used here)
   * @param[out] stresszz zz-component of the strain tensor array (only used here)
   * @param[out] stressxy xy-component of the strain tensor array (only used here)
   * @param[out] stressxz xz-component of the strain tensor array (only used here)
   * @param[out] stressyz yz-component of the strain tensor array (only used here)
   * @param[in] mass the mass matrix
   * @param[in] dampingx the damping (x-component) matrix
   * @param[in] dampingy the damping for (y-component) matrix
   * @param[in] dampingz the damping for (z-component) matrix
   * @param[in] dt time-step
   * @param[in] ux_np1 displacement array in the x direction (updated here)
   * @param[in] uy_np1 displacement array in the y direction (updated here)
   * @param[in] uz_np1 displacement array in the z direction (updated here)
   */

  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          localIndex const size_node,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView2d< real32 const > const stressxx,
          arrayView2d< real32 const > const stressyy,
          arrayView2d< real32 const > const stresszz,
          arrayView2d< real32 const > const stressxy,
          arrayView2d< real32 const > const stressxz,
          arrayView2d< real32 const > const stressyz,
          arrayView1d< const real32 > const mass,
          arrayView1d< real32 const > const dampingx,
          arrayView1d< real32 const > const dampingy,
          arrayView1d< real32 const > const dampingz,
          real64 const dt,
          arrayView1d< real32 > const ux_np1,
          arrayView1d< real32 > const uy_np1,
          arrayView1d< real32 > const uz_np1 )
  {

    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      ux_np1[a] *= 1.0-((dt/2)*(dampingx[a]/mass[a]));
      uy_np1[a] *= 1.0-((dt/2)*(dampingy[a]/mass[a]));
      uz_np1[a] *= 1.0-((dt/2)*(dampingz[a]/mass[a]));
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
      real32 uelemx[numNodesPerElem] = {0.0};
      real32 uelemy[numNodesPerElem] = {0.0};
      real32 uelemz[numNodesPerElem] = {0.0};
      real32 flowx[numNodesPerElem] = {0.0};
      real32 flowy[numNodesPerElem] = {0.0};
      real32 flowz[numNodesPerElem] = {0.0};

      for( localIndex q = 0; q < numNodesPerElem; q++ )
      {


        // Stiffness part
        m_finiteElement.template computeFirstOrderStiffnessTermX( q, xLocal, [&] ( int i, int j, real32 dfx1, real32 dfx2, real32 dfx3 )
        {
          flowx[i] -= stressxx[k][j]*dfx1 + stressxy[k][j]*dfx2 + stressxz[k][j]*dfx3;
          flowy[i] -= stressxy[k][j]*dfx1 + stressyy[k][j]*dfx2 + stressyz[k][j]*dfx3;
          flowz[i] -= stressxz[k][j]*dfx1 + stressyz[k][j]*dfx2 + stresszz[k][j]*dfx3;

        } );

        m_finiteElement.template computeFirstOrderStiffnessTermY( q, xLocal, [&] ( int i, int j, real32 dfy1, real32 dfy2, real32 dfy3 )
        {
          flowx[i] -= stressxx[k][j]*dfy1 + stressxy[k][j]*dfy2 + stressxz[k][j]*dfy3;
          flowy[i] -= stressxy[k][j]*dfy1 + stressyy[k][j]*dfy2 + stressyz[k][j]*dfy3;
          flowz[i] -= stressxz[k][j]*dfy1 + stressyz[k][j]*dfy2 + stresszz[k][j]*dfy3;
        } );

        m_finiteElement.template computeFirstOrderStiffnessTermZ( q, xLocal, [&] ( int i, int j, real32 dfz1, real32 dfz2, real32 dfz3 )
        {
          flowx[i] -= stressxx[k][j]*dfz1 + stressxy[k][j]*dfz2 + stressxz[k][j]*dfz3;
          flowy[i] -= stressxy[k][j]*dfz1 + stressyy[k][j]*dfz2 + stressyz[k][j]*dfz3;
          flowz[i] -= stressxz[k][j]*dfz1 + stressyz[k][j]*dfz2 + stresszz[k][j]*dfz3;
        } );

      }
      // Time update
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        uelemx[i]+=dt*flowx[i];
        uelemy[i]+=dt*flowy[i];
        uelemz[i]+=dt*flowz[i];
      }

      // Mult by inverse mass matrix + damping matrix
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        real32 localIncrement1 = uelemx[i]/mass[elemsToNodes[k][i]];
        real32 localIncrement2 = uelemy[i]/mass[elemsToNodes[k][i]];
        real32 localIncrement3 = uelemz[i]/mass[elemsToNodes[k][i]];
        RAJA::atomicAdd< ATOMIC_POLICY >( &ux_np1[elemsToNodes[k][i]], localIncrement1 );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uy_np1[elemsToNodes[k][i]], localIncrement2 );
        RAJA::atomicAdd< ATOMIC_POLICY >( &uz_np1[elemsToNodes[k][i]], localIncrement3 );
      }

    } );
    forAll< EXEC_POLICY >( size_node, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      ux_np1[a] /= 1.0+((dt/2)*(dampingx[a]/mass[a]));
      uy_np1[a] /= 1.0+((dt/2)*(dampingy[a]/mass[a]));
      uz_np1[a] /= 1.0+((dt/2)*(dampingz[a]/mass[a]));
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};


} // namespace ElasticFirstOrderWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIRSTORDERWAVEEQUATIONSEMKERNEL_HPP_
