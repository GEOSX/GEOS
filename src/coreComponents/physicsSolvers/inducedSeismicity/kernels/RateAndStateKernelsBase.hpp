/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_RATEANDSTATEKERNELSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_RATEANDSTATEKERNELSBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "physicsSolvers/inducedSeismicity/rateAndStateFields.hpp"

namespace geos
{

namespace rateAndStateKernels
{

// TBD: Pass the kernel and add getters for relevant fields to make this function general purpose and avoid
// wrappers?
GEOS_HOST_DEVICE
static void projectSlipRateBase( localIndex const k,
                                 real64 const frictionCoefficient,
                                 real64 const shearImpedance,
                                 arrayView1d< real64 const > const normalTraction,
                                 arrayView2d< real64 const > const shearTraction,
                                 arrayView1d< real64 const > const slipRate,
                                 arrayView2d< real64 > const slipVelocity )
{
  // Project slip rate onto shear traction to get slip velocity components
  real64 const frictionForce = normalTraction[k] * frictionCoefficient;
  real64 const projectionScaling = 1.0 / ( shearImpedance +  frictionForce / slipRate[k] );
  slipVelocity[k][0] = projectionScaling * shearTraction[k][0];
  slipVelocity[k][1] = projectionScaling * shearTraction[k][1];
}

template< typename POLICY, typename KERNEL_TYPE >
static bool newtonSolve( SurfaceElementSubRegion & subRegion,
                         KERNEL_TYPE & kernel,
                         real64 const dt,
                         integer const maxNewtonIter,
                         real64 const newtonTol )
{
  bool allConverged = false;
  for( integer iter = 0; iter < maxNewtonIter; iter++ )
  {
    RAJA::ReduceMin< parallelDeviceReduce, int > converged( 1 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > residualNorm( 0.0 );
    forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;
      kernel.setup( k, dt, stack );
      kernel.solve( k, stack );
      auto const [elementConverged, elementResidualNorm] = kernel.checkConvergence( stack, newtonTol );
      converged.min( elementConverged );
      residualNorm.max( elementResidualNorm );
    } );

    real64 const maxResidualNorm = MpiWrapper::max( residualNorm.get() );
    GEOS_LOG_RANK_0( GEOS_FMT( "   Newton iter {} : residual = {:.10e} ", iter, maxResidualNorm ) );

    if( converged.get() )
    {
      allConverged = true;
      break;
    }
  }
  return allConverged;
}

// template< typename POLICY, typename KERNEL_TYPE >
// static real64 solveRateAndStateEquation( SurfaceElementSubRegion & subRegion,
//                                          KERNEL_TYPE & kernel,
//                                          real64 dt,
//                                          integer const maxNewtonIter,
//                                          real64 const newtonTol )
// {
//   bool converged = false;
//   for( integer attempt = 0; attempt < 5; attempt++ )
//   {
//     if( attempt > 0 )
//     {
//       forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
//       {
//         kernel.resetState( k );
//       } );
//     }
//     GEOS_LOG_RANK_0( GEOS_FMT( "  Attempt {} ", attempt ) );
//     converged = newtonSolve< POLICY >( subRegion, kernel, dt, maxNewtonIter, newtonTol );
//     if( converged )
//     {
//       forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
//       {
//         kernel.udpateVariables( k );
//       } );
//       return dt;
//     }
//     else
//     {
//       GEOS_LOG_RANK_0( GEOS_FMT( "  Attempt {} failed. Halving dt and retrying.", attempt ) );
//       dt *= 0.5;
//     }
//   }
//   if( !converged )
//   {
//     GEOS_ERROR( "Maximum number of attempts reached without convergence." );
//   }
//   return dt;
// }

/**
 * @brief Performs the kernel launch
 * @tparam POLICY the policy used in the RAJA kernels
 */
template< typename KERNEL_TYPE, typename POLICY >
static void
createAndLaunch( SurfaceElementSubRegion & subRegion,
                 string const & frictionLawNameKey,
                 real64 const shearImpedance,
                 integer const maxNewtonIter,
                 real64 const newtonTol,
                 real64 const time_n,
                 real64 const totalDt )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n );

  string const & frictionaLawName = subRegion.getReference< string >( frictionLawNameKey );
  constitutive::RateAndStateFriction const & frictionLaw = subRegion.getConstitutiveModel< constitutive::RateAndStateFriction >( frictionaLawName );
  KERNEL_TYPE kernel( subRegion, frictionLaw, shearImpedance );

  real64 dtRemaining = totalDt;
  real64 dt = totalDt;
  for( integer subStep = 0; subStep < 5 && dtRemaining > 0.0; ++subStep )
  {
    real64 dtAccepted = KERNEL_TYPE::template solveRateAndStateEquation< POLICY >( subRegion, kernel, dt, maxNewtonIter, newtonTol );
    dtRemaining -= dtAccepted;

    if( dtRemaining > 0.0 )
    {
      dt = dtAccepted;
    }
    GEOS_LOG_RANK_0( GEOS_FMT( "  sub-step = {} completed, dt = {}, remaining dt = {}", subStep, dt, dtRemaining ) );
  }
}

} /* namespace rateAndStateKernels */

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_RATEANDSTATEKERNELSBASE_HPP_ */
