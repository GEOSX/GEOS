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
 * @file CO2EOSSolver.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2EOSSolver.hpp"

#include "common/Units.hpp"
#include "common/logger/Logger.hpp"


namespace geos
{

namespace constitutive
{

namespace PVTProps
{

real64
CO2EOSSolver::solve( string const & name,
                     integer const maxNumNewtonIter,
                     integer const maxNumBacktrackIter,
                     real64 const & tolerance,
                     real64 const & minAbsDeriv,
                     real64 const & maxAbsUpdate,
                     real64 const & allowedMinValue,
                     real64 const & initialGuess,
                     real64 const & temp,
                     real64 const & pres,
                     real64 const & presMultiplierForReporting,
                     real64 (* f)( real64 const & t, real64 const & p, real64 const & var ))
{
  // evaluate the residual using the initial guess
  real64 var = initialGuess;
  real64 res = (*f)( temp, pres, var );
  real64 update = LvArray::NumericLimits< real64 >::infinity;

  // follow the standard practice
  // (and multiply by 10, otherwise the FD derivative is often too close to zero.
  //  To remove this fudge factor, we should probably implement analytical derivatives here).
  real64 const dVar = 10 * sqrt( LvArray::NumericLimits< real64 >::epsilon ) * var;

  // start the Newton iterations
  bool newtonHasConverged = false;
  for( integer newtonIter = 0; newtonIter < maxNumNewtonIter; ++newtonIter )
  {

    // chop back if the variable goes below the minimum value
    if( var < allowedMinValue )
    {
      var = allowedMinValue;
      res = (*f)( temp, pres, var );
    }

    // compute finite-difference derivative
    real64 const perturbedVar = var + dVar;
    real64 const perturbedRes = (*f)( temp, pres, perturbedVar );
    real64 deriv = ( perturbedRes - res ) / dVar;

    // we need to make sure that we won't divide by zero
    // this is critical for the Helmholtz energy equation (Span-Wagner density)
    if( deriv > 0 && deriv < minAbsDeriv )
    {
      deriv = minAbsDeriv;
    }
    else if( deriv <= 0 && deriv > -minAbsDeriv )
    {
      deriv = -minAbsDeriv;
    }

    // compute Newton update
    update = -res / deriv;

    // we need to impose a max update value, otherwise Newton struggles when deriv is close to zero
    // this is critical for the Helmholtz energy equation (Span-Wagner density)
    if( update > maxAbsUpdate )
    {
      update = maxAbsUpdate;
    }
    else if( update < -maxAbsUpdate )
    {
      update = -maxAbsUpdate;
    }

    // recompute the residual
    real64 const newVar = var + update;
    real64 const newRes = (*f)( temp, pres, newVar );

    // check convergence based on the magnitude of the residual
    if( LvArray::math::abs( newRes ) < tolerance )
    {
      var = newVar;
      newtonHasConverged = true;
      break;
    }

    // if the new residual is larger than the previous one, do some backtracking steps
    bool backtrackingIsSuccessful = false;
    real64 alpha = 1.0;
    real64 backtrackVar = newVar;
    real64 backtrackRes = LvArray::NumericLimits< real64 >::infinity;
    if( LvArray::math::abs( newRes ) > LvArray::math::abs( res ) )
    {

      // launch the backtracking steps, and break if successful
      for( integer iBacktrackIter = 0; iBacktrackIter < maxNumBacktrackIter; ++iBacktrackIter )
      {
        alpha *= 0.5;
        backtrackVar = var + alpha*update;
        backtrackRes = (*f)( temp, pres, backtrackVar );

        // if we managed to reduce the residual, we can break the loop
        if( LvArray::math::abs( backtrackRes ) < LvArray::math::abs( res ) )
        {
          backtrackingIsSuccessful = true;
          var = backtrackVar;
          res = backtrackRes;
          break;
        }
      }
    }

    // if backtracking failed, keep the original new residual
    if( !backtrackingIsSuccessful )
    {
      var = newVar;
      res = newRes;
    }
  }

  GEOS_THROW_IF( !newtonHasConverged,
                 name << ": Newton's method failed to converge for pair "
                      << "( pressure = " << pres*presMultiplierForReporting << " Pa, temperature = " << units::convertCToK( temp ) << " K) :"
                      << " final residual = " << res << ", final update = " << update << ", tolerance = " << tolerance,
                 InputError );
  return var;
}



} // namespace PVTProps

} // namespace constitutive

} // namespace geos
