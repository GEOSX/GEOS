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
 * @file FluidUpdateKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUIDUPDATEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUIDUPDATEKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void launch( FLUID_WRAPPER const & fluidWrapper,
                      arrayView1d< real64 const > const & pres,
                      arrayView1d< real64 const > const & temp )
  {
    forAll< parallelDevicePolicy<> >( fluidWrapper.numElems(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k] );
      }
    } );
  }
};

} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUIDUPDATEKERNEL_HPP
