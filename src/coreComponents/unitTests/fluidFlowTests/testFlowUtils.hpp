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

#ifndef GEOS_UNITTESTS_FLUIDFLOWTESTS_TESTFLOWUTILS_HPP
#define GEOS_UNITTESTS_FLUIDFLOWTESTS_TESTFLOWUTILS_HPP

#include "codingUtilities/UnitTestUtilities.hpp"

namespace geos
{

namespace testing
{

void checkDerivative( real64 const valueEps,
                      real64 const value,
                      real64 const deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var )
{
  real64 const numDeriv = (valueEps - value) / eps;
  checkRelativeError( deriv, numDeriv, relTol, absTol, "d(" + name + ")/d(" + var + ")" );
}

void checkDerivative( real64 const valueEps,
                      real64 const value,
                      real64 const deriv,
                      real64 const eps,
                      real64 const relTol,
                      string const & name,
                      string const & var )
{ return checkDerivative( valueEps, value, deriv, eps, relTol, DEFAULT_ABS_TOL, name, var ); }

template< int USD1, int USD2, int USD3 >
void checkDerivative( arraySlice1d< real64 const, USD1 > const & valueEps,
                      arraySlice1d< real64 const, USD2 > const & value,
                      arraySlice1d< real64 const, USD3 > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol, absTol,
                     name + "[" + labels[i] + "]", var );
  }
}

template< int DIM, int USD1, int USD2, int USD3, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM, USD1 > const & valueEps,
                      ArraySlice< real64 const, DIM, USD2 > const & value,
                      ArraySlice< real64 const, DIM, USD3 > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels,
                      Args ... label_lists )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol, absTol,
                     name + "[" + labels[i] + "]", var, label_lists ... );
  }
}

template< int DIM, int USD1, int USD2, int USD3, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM, USD1 > const & valueEps,
                      ArraySlice< real64 const, DIM, USD2 > const & value,
                      ArraySlice< real64 const, DIM, USD3 > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels,
                      Args ... label_lists )
{ return checkDerivative( valueEps, value, deriv, eps, relTol, DEFAULT_ABS_TOL, name, var, labels, label_lists ... ); }

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
template< int USD >
array1d< real64 > invertLayout( arraySlice1d< real64 const, USD > const & input,
                                localIndex N )
{
  array1d< real64 > output( N );
  for( int i = 0; i < N; ++i )
  {
    output[i] = input[i];
  }

  return output;
}

template< int USD >
array2d< real64 > invertLayout( arraySlice2d< real64 const, USD > const & input,
                                localIndex N1,
                                localIndex N2 )
{
  array2d< real64 > output( N2, N1 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      output( j, i ) = input( i, j );
    }
  }

  return output;
}

template< int USD >
array3d< real64 > invertLayout( arraySlice3d< real64 const, USD > const & input,
                                localIndex N1,
                                localIndex N2,
                                localIndex N3 )
{
  array3d< real64 > output( N3, N1, N2 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      for( localIndex k = 0; k < N3; ++k )
      {
        output( k, i, j ) = input( i, j, k );
      }
    }
  }

  return output;
}

} // namespace testing

} // namespace geos

#endif //GEOS_UNITTESTS_FLUIDFLOWTESTS_TESTFLOWUTILS_HPP
