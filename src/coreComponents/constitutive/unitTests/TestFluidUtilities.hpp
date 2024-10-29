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
 * @file TestFluidUtilities.hpp
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_

#include "codingUtilities/UnitTestUtilities.hpp"

namespace geos
{

namespace testing
{

namespace internal
{

static constexpr real64 relTol = 1.0e-4;
static constexpr real64 absTol = 1.0e-7;

/**
 * @brief Tests a function against a derivative
 * @details Will calculate the left-sided and the right-sided numerical derivatives of a function
 *          and compare this against a analytically calculated value provided.
 * @param x The value at which the function should be evaluated
 * @param dx The value to use to perturb @c x in the calculation of the numerical derivatives
 * @param derivative The value of the analytically calculated derivative to use for comparison
 * @param function The function which is being tested. This should be a function that returns
          a single real value.
 * @param absTolerance The absolute tolerance to use for the comparison
 * @param relTolerance The relative tolerance to use for the comparison
 */
template< typename FUNCTION >
void testNumericalDerivative( real64 const x, real64 const dx, real64 const derivative,
                              FUNCTION && function,
                              real64 const absTolerance = absTol, real64 const relTolerance = relTol )
{
  real64 const leftValue = function( x-dx );
  real64 const centreValue = function( x );
  real64 const rightValue = function( x+dx );
  real64 const leftDerivative = (centreValue - leftValue) / dx;
  real64 const rightDerivative = (rightValue - centreValue) / dx;
  checkRelativeError( derivative, leftDerivative, relTolerance, absTolerance, "Left derivative" );
  checkRelativeError( derivative, rightDerivative, relTolerance, absTolerance, "Right derivative" );
}

/**
 * @brief Tests a multi-valued function against a derivative
 * @details Will calculate the left-sided and the right-sided numerical derivatives of a function
 *          and compare this against a analytically calculated values provided.
 * @tparam numValues the number of values that the function returns
 * @tparam FUNCTION the type of function (typically a lambda)
 * @param x The value at which the function should be evaluated
 * @param dx The value to use to perturb @c x in the calculation of the numerical derivatives
 * @param derivatives The values of the analytically calculated derivatives to use for comparison
 * @param function The function which is being tested. This should be a function that takes 2 parameters.
 *       The first is the value at which the function is being evaluated (x) and the second is an array
 *       of size @c numValues which is the result of the execution of the function.
 * @param absTolerance The absolute tolerance to use for the comparison
 * @param relTolerance The relative tolerance to use for the comparison
 */
template< integer numValues, typename FUNCTION >
void testNumericalDerivative( real64 const x,
                              real64 const dx,
                              arraySlice1d< real64 const > const & derivatives,
                              FUNCTION && function,
                              real64 const absTolerance = absTol,
                              real64 const relTolerance = relTol )
{
  stackArray1d< real64, numValues > leftValues( numValues );
  stackArray1d< real64, numValues > centreValues( numValues );
  stackArray1d< real64, numValues > rightValues( numValues );
  function( x-dx, leftValues );
  function( x, centreValues );
  function( x+dx, rightValues );

  // Use the same space to calculate the left-sided and right sided derivatives
  for( integer i = 0; i < numValues; ++i )
  {
    // Choose from the left, central and right derivatives, the one that's nearest the analytical value
    real64 minError = LvArray::NumericLimits< real64 >::max;
    real64 selectedDerivative = 0.0;
    for( real64 const distance : {centreValues[i] - leftValues[i],
                                  rightValues[i] - centreValues[i],
                                  0.5*(rightValues[i] - leftValues[i])} )
    {
      real64 const deriv = distance / dx;
      real64 const error = LvArray::math::abs( deriv - derivatives[i] );
      if( error < minError )
      {
        minError = error;
        selectedDerivative = deriv;
      }
    }
    checkRelativeError( derivatives[i], selectedDerivative, relTolerance, absTolerance,
                        GEOS_FMT( "Numerical derivative for component {}", i ) );
  }
}

}// namespace internal

}// namespace testing

}// namespace geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_
