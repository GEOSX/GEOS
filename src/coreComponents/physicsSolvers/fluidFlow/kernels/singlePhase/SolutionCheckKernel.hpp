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
 * @file SolutionCheckKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_SOLUTIONCHECKKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_SOLUTIONCHECKKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY >
  static std::pair< integer, real64 > launch( arrayView1d< real64 const > const & localSolution,
                                              globalIndex const rankOffset,
                                              arrayView1d< globalIndex const > const & dofNumber,
                                              arrayView1d< integer const > const & ghostRank,
                                              arrayView1d< real64 const > const & pres,
                                              real64 const scalingFactor )
  {
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegativePressures( 0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPres( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 && dofNumber[ei] >= 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;
        real64 const newPres = pres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          numNegativePressures += 1;
          minPres.min( newPres );
        }
      }

    } );

    return { numNegativePressures.get(), minPres.get() };
  }

};

} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_SOLUTIONCHECKKERNEL_HPP
