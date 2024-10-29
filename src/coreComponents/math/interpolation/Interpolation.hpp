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

/*
 * File: Interpolation.hpp
 */

#ifndef GEOS_MATH_INTERPOLATION_HPP_
#define GEOS_MATH_INTERPOLATION_HPP_


namespace geos
{

namespace interpolation
{

static inline
real64 parabolicInterpolationThreePoints( real64 const lambdac,
                                          real64 const lambdam,
                                          real64 const ff0,
                                          real64 const ffT,
                                          real64 const ffm )
{
  // Apply three-point safeguarded parabolic model for a line search.
  //
  // C. T. Kelley, April 1, 2003
  //
  // This code comes with no guarantee or warranty of any kind.
  //
  // function lambdap = parab3p(lambdac, lambdam, ff0, ffT, ffm)
  //
  // input:
  //       lambdac = current steplength
  //       lambdam = previous steplength
  //       ff0 = value of \| F(x_c) \|^2
  //       ffc = value of \| F(x_c + \lambdac d) \|^2
  //       ffm = value of \| F(x_c + \lambdam d) \|^2
  //
  // output:
  //       lambdap = new value of lambda given parabolic model
  //
  // internal parameters:
  //       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch

  // Set internal parameters.
  real64 const sigma0 = 0.1;
  real64 const sigma1 = 0.5;

  // Compute coefficients of interpolation polynomial.
  // p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
  // d1 = (lambdac - lambdam)*lambdac*lambdam < 0
  //      so, if c2 > 0 we have negative curvature and default to
  //      lambdap = sigma1 * lambda.
  real64 const c2 = lambdam*(ffT-ff0)-lambdac*(ffm-ff0);
  if( c2 >= 0.0 )
  {
    return ( sigma1*lambdac );
  }
  real64 const c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffT-ff0);
  real64 lambdap = -c1*0.5;
  if( lambdap > sigma0*lambdac*c2 )
  {
    lambdap = sigma0*lambdac*c2;
  }
  if( lambdap < sigma1*lambdac*c2 )
  {
    lambdap = sigma1*lambdac*c2;
  }
  lambdap /= c2;
  return lambdap;
}

/**
 * @brief Perform linear interpolation of function f(x) at x3 on interval [x1,x2]
 * @tparam T the value type
 * @param[in] dx1 length of interval [x1,x3]: x3-x1
 * @param[in] dx2 length of interval [x3,x2]: x2-x3
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @return the interpolated value of function f(x) at x3
 */
template< typename T >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
T linearInterpolation( T const dx1,
                       T const dx2,
                       T const f1,
                       T const f2 )
{
  T const alpha = dx1 / ( dx1 + dx2 );
  return alpha * f2 + ( 1.0 - alpha ) * f1;
}

/**
 * @brief Perform linear interpolation of function f(x) at x3 on interval [x1,x2]
 * @tparam T the value type
 * @param[in] dx1 length of interval [x1,x3]: x3-x1
 * @param[in] dx2 length of interval [x3,x2]: x2-x3
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @param[out] f the interpolated value of function f(x) at x3
 * @param[out] df_dx the derivative (wrt x) of the interpolated value of function f(x) at x3
 */
template< typename T >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void linearInterpolation( T const dx1,
                          T const dx2,
                          T const f1,
                          T const f2,
                          T & f,
                          T & df_dx )
{
  T const alpha = dx1 / ( dx1 + dx2 );
  T const dAlpha_dx = 1.0 / ( dx1 + dx2 );
  f = alpha * f2 + ( 1.0 - alpha ) * f1;
  df_dx = dAlpha_dx * ( f2 - f1 );
}

/**
 * @brief Perform linear interpolation of function f(x) for all the values of a given array
 * @tparam T the value type
 * @param[in] xIn the values at which we know the value of f(x)
 * @param[in] fIn the known values of function f(x) at the values in xIn
 * @param[in] xOut the values where we want to interpolate f(x)
 * @param[out] fOut the interpolated values of f(x) at xOut
 */
template< typename T >
void linearInterpolation( arrayView1d< T const > const & xIn,
                          arrayView1d< T const > const & fIn,
                          arrayView1d< T const > const & xOut,
                          arrayView1d< T > const & fOut )
{
  GEOS_ASSERT_EQ_MSG( xIn.size(), fIn.size(), "linearInterpolation: size mismatch in the xIn and fIn input vectors" );
  GEOS_ASSERT_EQ_MSG( xOut.size(), fOut.size(), "linearInterpolation: size mismatch in the xOut and fOut input vectors" );

  for( localIndex i = 0; i < xOut.size(); ++i )
  {
    integer const idx = LvArray::sortedArrayManipulation::find( xIn.begin(),
                                                                xIn.size(),
                                                                xOut[i] );
    integer const iUp  = LvArray::math::min( LvArray::math::max( idx, 1 ),
                                             LvArray::integerConversion< integer >( xIn.size()-1 ) );
    integer const iLow = iUp-1;
    fOut[i] = linearInterpolation( xOut[i] - xIn[iLow], xIn[iUp] - xOut[i], fIn[iLow], fIn[iUp] );
  }
}


}

}

#endif // GEOS_MATH_INTERPOLATION_HPP_
