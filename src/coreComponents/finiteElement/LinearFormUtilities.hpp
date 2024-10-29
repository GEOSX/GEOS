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
 * @file BilinearFormUtilities.hpp
 */

#ifndef GEOS_FINITEELEMENT_LINEARFORMUTILITIES_HPP_
#define GEOS_FINITEELEMENT_LINEARFORMUTILITIES_HPP_

#include "finiteElement/PDEUtilities.hpp"

namespace geos
{

namespace LinearFormUtilities
{

template< PDEUtilities::FunctionSpace V,
          PDEUtilities::DifferentialOperator OP_V >
struct Helper
{};

template<>
struct Helper< PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const & Nv,
                       real64 const A,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Nv );
    for( int a = 0; a < numTestDOF; ++a )
    {
      vec[a] = vec[a] + A * weight;
    }
  }

  template< int numTestDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const & Nv,
                       real64 const A[numTestDOF],
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Nv );
    for( int a = 0; a < numTestDOF; ++a )
    {
      vec[a] = vec[a] + A[a] * weight;
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF >
  GEOS_HOST_DEVICE
  void static compute( double (& vec)[numTestDOF],
                       double const (&Nv)[numTestDOF/3],
                       double const (&A)[3],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + Nv[a] * A[0] * weight;
      vec[a*3+1] = vec[a*3+1] + Nv[a] * A[1] * weight;
      vec[a*3+2] = vec[a*3+2] + Nv[a] * A[2] * weight;
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::SymmetricGradient >
{
  // symmetric second-order tensor A
  template< int numTestDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[6],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + ( dNdX[a][0] * A[0] + dNdX[a][1] * A[5] + dNdX[a][2] * A[4] ) * weight;
      vec[a*3+1] = vec[a*3+1] + ( dNdX[a][0] * A[5] + dNdX[a][1] * A[1] + dNdX[a][2] * A[3] ) * weight;
      vec[a*3+2] = vec[a*3+2] + ( dNdX[a][0] * A[4] + dNdX[a][1] * A[3] + dNdX[a][2] * A[2] ) * weight;
    }
  }

  // diagonal second-order tensor
  template< int numTestDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[3],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + dNdX[a][0] * A[0] * weight;
      vec[a*3+1] = vec[a*3+1] + dNdX[a][1] * A[1] * weight;
      vec[a*3+2] = vec[a*3+2] + dNdX[a][2] * A[2] * weight;
    }
  }

  // scalar*identity second-order tensor
  template< int numTestDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const A,
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + dNdX[a][0] * A * weight;
      vec[a*3+1] = vec[a*3+1] + dNdX[a][1] * A * weight;
      vec[a*3+2] = vec[a*3+2] + dNdX[a][2] * A * weight;
    }
  }
};

/**
 * @brief Generic linear form template to assemble elemental vectors.
 *
 * This function computes the elemental vector for a grid cell \f$ \Omega^e \f$
 * with entries defined by the linear form
 *
 * \f[
 *   \int_{\Omega^e}
 *     \texttt{OP}_V ( v^h ) \; \cdot \;
 *     \mathsf{T}  \;
 *     \mathrm{d}\Omega
 * \f]
 *
 * where \f$ v^h \f$ is a finite element function belonging to
 * the test space (@p V), \f$ \texttt{OP}_V \f$ denotes a differential
 * operator (zero- or first-order) applied to the test function, and
 * \f$ \mathsf{T} \f$ indicates a tensor.
 *
 * @tparam V Test function space.
 * @tparam OP_V Differential operator applied to functions in @p V.
 * @tparam VECTOR Derived vector type.
 * @tparam FE_VALUES_V Derived test shape function values or derivatives type.
 * @tparam TENSOR Derived tensor type.
 * @param vector The elemental vector.
 * @param feValuesV Test shape function values or derivatives.
 * @param tensor The tensor defining the linear form.
 * @param weight Quadrature weight.
 *
 */
template< PDEUtilities::FunctionSpace V,
          PDEUtilities::DifferentialOperator OP_V,
          typename VECTOR,
          typename FE_VALUES_V,
          typename TENSOR >
GEOS_HOST_DEVICE
static void compute( VECTOR && vector,
                     FE_VALUES_V const & feValuesV,
                     TENSOR const & tensor,
                     real64 const weight )
{
  Helper< V, OP_V >::compute( vector, feValuesV, tensor, weight );
}

} // namespace LinearFormUtilities

} // namespace geos

#endif //GEOS_FINITEELEMENT_LINEARFORMUTILITIES_HPP_
