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
 * @file StencilWeightsUpdateKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_STENCILWEIGHTSUPDATEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_STENCILWEIGHTSUPDATEKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{

namespace flowSolverBaseKernels
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

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

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_STENCILWEIGHTSUPDATEKERNEL_HPP
