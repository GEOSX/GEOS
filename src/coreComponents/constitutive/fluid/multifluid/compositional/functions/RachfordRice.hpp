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
 * @file RachfordRice.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

struct RachfordRice
{
public:
  /// Tolerance of the SSI loop
  static constexpr real64 SSITolerance = MultiFluidConstants::SSITolerance;
  /// Tolerance of the Newton loop
  static constexpr real64 newtonTolerance = MultiFluidConstants::newtonTolerance;
  /// Max number of SSI iterations
  static constexpr integer maxSSIIterations = MultiFluidConstants::maxSSIIterations;
  /// Max number of Newton iterations
  static constexpr integer maxNewtonIterations = MultiFluidConstants::maxNewtonIterations;
  /// Epsilon used in the calculations
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;

  /**
   * @brief Function solving the Rachford-Rice equation
   * @input[in] kValues the array of K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @return the gas mole fraction
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  solve( arraySlice1d< real64 const, USD2 > const & kValues,
         arraySlice1d< real64 const, USD1 > const & feed,
         arraySlice1d< integer const > const & presentComponentIds )
  {
    real64 gasPhaseMoleFraction = 0;

    // min and max Kvalues for non-zero composition
    real64 minK, maxK;
    findKValueRange( kValues, presentComponentIds, minK, maxK );

    // check for trivial solutions.
    // this corresponds to bad Kvalues
    if( maxK < 1.0 )
    {
      return 0.0;
    }
    if( minK > 1.0 )
    {
      return 1.0;
    }

    // start the solver loop

    // step 1: find solution window
    real64 xMin = 1.0 / ( 1 - maxK );
    real64 xMax = 1.0 / ( 1 - minK );
    real64 const sqrtEpsilon = sqrt( epsilon );
    xMin += sqrtEpsilon * ( LvArray::math::abs( xMin ) + sqrtEpsilon );
    xMax -= sqrtEpsilon * ( LvArray::math::abs( xMax ) + sqrtEpsilon );

    real64 currentError = 1 / epsilon;

    // step 2: start the SSI loop
    // Evaluate at the bounds
    real64 funcXMin = evaluate( kValues, feed, presentComponentIds, xMin );
    real64 funcXMax = evaluate( kValues, feed, presentComponentIds, xMax );
    real64 funcXMid = 0.0;

    // If the bound values are the same sign then we have a trivial solution
    if( 0.0 < funcXMin * funcXMax )
    {
      gasPhaseMoleFraction = (0.0 < funcXMin) ? 1.0 : 0.0;
      return gasPhaseMoleFraction;
    }

    integer SSIIteration = 0;

    while( ( currentError > SSITolerance ) && ( SSIIteration < maxSSIIterations ) )
    {
      real64 const xMid = 0.5 * ( xMin + xMax );
      funcXMid = evaluate( kValues, feed, presentComponentIds, xMid );

      if( 0.0 < funcXMax * funcXMid )
      {
        xMax = xMid;
        funcXMax = funcXMid;
      }
      else if( 0.0 < funcXMin * funcXMid )
      {
        xMin = xMid;
        funcXMin = funcXMid;
      }

      currentError = LvArray::math::min( LvArray::math::abs( funcXMid ),
                                         LvArray::math::abs( xMax - xMin ) );
      SSIIteration++;

      // TODO: add warning if max number of SSI iterations is reached
    }

    gasPhaseMoleFraction = 0.5 * ( xMax + xMin );

    // step 3: start the Newton loop
    integer newtonIteration = 0;
    real64 newtonValue = gasPhaseMoleFraction;
    real64 funcNewton = evaluate( kValues, feed, presentComponentIds, newtonValue );

    while( ( currentError > newtonTolerance ) && ( newtonIteration < maxNewtonIterations ) )
    {
      real64 deltaNewton = -funcNewton / evaluateDerivative( kValues, feed, presentComponentIds, newtonValue );

      // test if we are stepping out of the [xMin;xMax] interval
      if( newtonValue + deltaNewton < xMin )
      {
        deltaNewton = 0.5 * ( xMin - newtonValue );
      }
      else if( newtonValue + deltaNewton > xMax )
      {
        deltaNewton = 0.5 * ( xMax - newtonValue );
      }

      newtonValue = newtonValue + deltaNewton;

      funcNewton = evaluate( kValues, feed, presentComponentIds, newtonValue );

      currentError = LvArray::math::min( LvArray::math::abs( funcNewton ),
                                         LvArray::math::abs( deltaNewton ) );
      newtonIteration++;

      // TODO: add warning if max number of Newton iterations is reached
    }
    return gasPhaseMoleFraction = newtonValue;
  }

private:

  /**
   * @brief Function evaluating the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the Rachford-Rice function is evaluated
   * @return the value of the Rachford-Rice function at x
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  evaluate( arraySlice1d< real64 const, USD2 > const & kValues,
            arraySlice1d< real64 const, USD1 > const & feed,
            arraySlice1d< integer const > const & presentComponentIds,
            real64 const & x )
  {
    real64 value = 0.0;
    for( integer i = 0; i < presentComponentIds.size(); ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const k = ( kValues[ic] - 1.0 );
      value += feed[ic] * k / ( 1.0 + x * k );
    }
    return value;
  }

  /**
   * @brief Function evaluating the derivative of the Rachford-Rice function
   * @input[in] kValues the array fo K-values
   * @input[in] feed the component fractions
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[in] x the value at which the derivative of the Rachford-Rice function is evaluated
   * @return the value of the derivative of the Rachford-Rice function at x
   **/
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  real64
  static
  evaluateDerivative( arraySlice1d< real64 const, USD2 > const & kValues,
                      arraySlice1d< real64 const, USD1 > const & feed,
                      arraySlice1d< integer const > const & presentComponentIds,
                      real64 const & x )
  {
    real64 value = 0.0;
    for( integer i = 0; i < presentComponentIds.size(); ++i )
    {
      integer const ic = presentComponentIds[i];
      real64 const k = ( kValues[ic] - 1.0 );
      real64 const r = k / ( 1.0 + x * k );
      value -= feed[ic] * r * r;
    }
    return value;
  }

  /**
   * @brief Calculate the minimum and maximum k-value
   * @input[in] kValues the array fo K-values
   * @input[in] presentComponentIds the indices of components with a non-zero fractions
   * @input[out] minK the minimum k-value for non-zero components
   * @input[out] maxK the maximum k-value for non-zero components
   **/
  template< integer USD >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void
  static
  findKValueRange( arraySlice1d< real64 const, USD > const & kValues,
                   arraySlice1d< integer const > const & presentComponentIds,
                   real64 & minK,
                   real64 & maxK )
  {
    minK = 1.0 / epsilon;
    maxK = 0.0;
    for( integer const ic : presentComponentIds )
    {
      if( kValues[ic] > maxK )
      {
        maxK = kValues[ic];
      }
      if( kValues[ic] < minK )
      {
        minK = kValues[ic];
      }
    }
  }

};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_RACHFORDRICE_HPP_
