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

#ifndef GEOS_MATH_EXTRAPOLATION_HPP_
#define GEOS_MATH_EXTRAPOLATION_HPP_


namespace geos
{

namespace extrapolation
{

/**
 * @brief Perform linear extrapolation of function f(x) at x3 using the known values of f on interval between x1 and x2
 * @tparam T the value type
 * @param[in] x1 first bound of the interval
 * @param[in] x2 second bound of the interval
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @param[in] x3 value at which we want to extrapolate linearly function f(x)
 * @return the value of function f(x) extrapolated at x3
 */
template< typename T >
T linearExtrapolation( T const x1,
                       T const x2,
                       T const f1,
                       T const f2,
                       T const x3 )
{
  return ( f2 - f1 ) / ( x2 - x1 ) * ( x3 - x2 ) + f2;
}

/**
 * @brief Perform log extrapolation of function f(x) at x3 using the known values of f on interval between x1 and x2
 * @tparam T the value type
 * @param[in] x1 first bound of the interval
 * @param[in] x2 second bound of the interval
 * @param[in] f1 value of function f(x) evaluated at x1
 * @param[in] f2 value of function f(x) evaluated at x2
 * @param[in] x3 value at which we want to extrapolate logarithmically function f(x)
 * @return the value of function f(x) extrapolated at x3
 */
template< typename T >
T logExtrapolation( T const x1,
                    T const x2,
                    T const f1,
                    T const f2,
                    T const x3 )
{
  T const lnf3 = ( log( f2 ) - log( f1 ) ) / ( log( x2 ) - log( x1 ) ) * ( log( x3 ) - log( x2 ) ) + log( f2 );
  return exp( lnf3 );
}

}

}

#endif // GEOS_MATH_EXTRAPOLATION_HPP_
