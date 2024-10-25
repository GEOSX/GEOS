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
 * @file FlowSolverBaseKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{

namespace flowSolverBaseKernels
{

/// Threshold for the min pore volume (below, a warning is issued)
static constexpr real64 poreVolumeThreshold = 1e-4;

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

/**
 * @struct MinPoreVolumeMaxPorosityKernel
 * @brief Kernel to compute the min pore volume and the max porosity in a subRegion
 */
struct MinPoreVolumeMaxPorosityKernel
{

  /*
   * @brief Kernel computing the min pore volume and the max porosity
   * @param[in] size the number of elements in the subRegion
   * @param[in] ghostRank the ghost ranks
   * @param[in] porosity the current element porosity
   * @param[in] volume the current element volume
   * @param[out] minPoreVolumeInSubRegion the min pore volume
   * @param[out] maxPorosityInSubRegion the max porosity
   * @param[out] numElemsBelowPoreVolumeThresholdInSubRegion the number of elements below the pore volume threshold
   * @param[out] numElemsAbovePorosityThresholdInSubRegion the number of elements with a porosity above 1
   */
  inline static void
  computeMinPoreVolumeMaxPorosity( localIndex const size,
                                   arrayView1d< integer const > const & ghostRank,
                                   arrayView2d< real64 const > const & porosity,
                                   arrayView1d< real64 const > const & volume,
                                   real64 & minPoreVolumeInSubRegion,
                                   real64 & maxPorosityInSubRegion,
                                   localIndex & numElemsBelowPoreVolumeThresholdInSubRegion,
                                   localIndex & numElemsAbovePorosityThresholdInSubRegion )
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > minPoreVolume( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > maxPorosity( -LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, localIndex > numElemsBelowPoreVolumeThreshold( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, localIndex > numElemsAbovePorosityThreshold( 0.0 );

    real64 const pvThreshold = poreVolumeThreshold;

    forAll< parallelDevicePolicy<> >( size, [ghostRank,
                                             porosity,
                                             volume,
                                             pvThreshold,
                                             minPoreVolume,
                                             maxPorosity,
                                             numElemsBelowPoreVolumeThreshold,
                                             numElemsAbovePorosityThreshold] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] >= 0 )
      {
        return;
      }

      real64 const poreVolume = porosity[ei][0] * volume[ei];
      if( poreVolume < pvThreshold )
      {
        numElemsBelowPoreVolumeThreshold += 1;
      }
      if( porosity[ei][0] > 1 )
      {
        numElemsAbovePorosityThreshold += 1;
      }

      minPoreVolume.min( poreVolume );
      maxPorosity.max( porosity[ei][0] );
    } );

    minPoreVolumeInSubRegion = minPoreVolume.get();
    maxPorosityInSubRegion = maxPorosity.get();
    numElemsBelowPoreVolumeThresholdInSubRegion = numElemsBelowPoreVolumeThreshold.get();
    numElemsAbovePorosityThresholdInSubRegion = numElemsAbovePorosityThreshold.get();
  }
};

/**
 * @brief
 *
 * @tparam STENCILWRAPPER
 */
template< typename STENCILWRAPPER >
struct stencilWeightsUpdateKernel
{
  /**
   * @brief
   *
   * @param stencilWrappper
   * @param hydraulicAperture
   */
  inline static void prepareStencilWeights( STENCILWRAPPER & stencilWrapper,
                                            ElementViewConst< arrayView1d< real64 const > > const hydraulicAperture )
  {
    forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      stencilWrapper.removeHydraulicApertureContribution( iconn, hydraulicAperture );
    } );
  }

  /**
   * @brief
   *
   * @param stencilWrappper
   * @param hydraulicAperture
   */
  inline static void updateStencilWeights( STENCILWRAPPER & stencilWrapper,
                                           ElementViewConst< arrayView1d< real64 const > > const hydraulicAperture )
  {
    forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      stencilWrapper.addHydraulicApertureContribution( iconn, hydraulicAperture );
    } );
  }
};

} // namespace flowSolverBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP
