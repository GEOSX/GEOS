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
 * @file SolutionScalingKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_SOLUTIONSCALINGKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_SOLUTIONSCALINGKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** SolutionScalingKernel ********************************/

struct SolutionScalingKernel
{
  template< typename POLICY >
  static std::pair< real64, real64 > launch( arrayView1d< real64 const > const & localSolution,
                                             globalIndex const rankOffset,
                                             arrayView1d< globalIndex const > const & dofNumber,
                                             arrayView1d< integer const > const & ghostRank,
                                             real64 const maxAbsolutePresChange )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > scalingFactor( 1.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaPres( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei ) mutable
    {
      if( ghostRank[ei] < 0 && dofNumber[ei] >= 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;

        // compute the change in pressure
        real64 const absPresChange = LvArray::math::abs( localSolution[lid] );
        maxDeltaPres.max( absPresChange );

        // maxAbsolutePresChange <= 0.0 means that scaling is disabled, and we are only collecting maxDeltaPres in this kernel
        if( maxAbsolutePresChange > 0.0 && absPresChange > maxAbsolutePresChange )
        {
          real64 const presScalingFactor = maxAbsolutePresChange / absPresChange;
          scalingFactor.min( presScalingFactor );
        }
      }

    } );

    return { scalingFactor.get(), maxDeltaPres.get() };
  }

};

} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_SOLUTIONSCALINGKERNEL_HPP
