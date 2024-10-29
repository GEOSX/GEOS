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

#ifndef GEOS_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_
#define GEOS_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_

#include "finiteElement/PDEUtilities.hpp"

namespace geos
{

namespace BilinearFormUtilities
{

template< PDEUtilities::FunctionSpace V,
          PDEUtilities::FunctionSpace U,
          PDEUtilities::DifferentialOperator OP_V,
          PDEUtilities::DifferentialOperator OP_U >
struct Helper
{};

template<>
struct Helper< PDEUtilities::FunctionSpace::P0,
               PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Divergence >
{
  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const & Nv,
                       real64 const A,
                       real64 const (&dNudX)[numTrialDOF/3][3],
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Nv );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF/3; ++b )
      {
        mat[a][b*3+0] = mat[a][b*3+0] + A * dNudX[b][0] * weight;
        mat[a][b*3+1] = mat[a][b*3+1] + A * dNudX[b][1] * weight;
        mat[a][b*3+2] = mat[a][b*3+2] + A * dNudX[b][2] * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const & Nv,
                       real64 const (&A)[numTestDOF],
                       real64 const (&dNudX)[numTrialDOF/3][3],
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Nv );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF/3; ++b )
      {
        mat[a][b*3+0] = mat[a][b*3+0] + A[a] * dNudX[b][0] * weight;
        mat[a][b*3+1] = mat[a][b*3+1] + A[a] * dNudX[b][1] * weight;
        mat[a][b*3+2] = mat[a][b*3+2] + A[a] * dNudX[b][2] * weight;
      }
    }
  }

};

template<>
struct Helper< PDEUtilities::FunctionSpace::P0,
               PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const & Nv,
                       double const A,
                       double const & Nu,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Nv );
    GEOS_UNUSED_VAR( Nu );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + A * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const & Nv,
                       double const A[numTestDOF],
                       double const & Nu,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Nv );
    GEOS_UNUSED_VAR( Nu );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + A[a] * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const & Nv,
                       double const A[numTestDOF][numTrialDOF],
                       double const & Nu,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Nv );
    GEOS_UNUSED_VAR( Nu );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + A[a][b] * weight;
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&N)[numTestDOF/3],
                       real64 const (&A)[3],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + N[a] * A[0] * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + N[a] * A[1] * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + N[a] * A[2] * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&N)[numTestDOF/3],
                       real64 const (&A)[3][numTrialDOF],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + N[a] * A[0][b] * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + N[a] * A[1][b] * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + N[a] * A[2][b] * weight;
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::SymmetricGradient,
               PDEUtilities::DifferentialOperator::Identity >
{
  // symmetric second-order tensor A
  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[6],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + ( dNdX[a][0] * A[0] + dNdX[a][1] * A[5] + dNdX[a][2] * A[4] ) * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + ( dNdX[a][0] * A[5] + dNdX[a][1] * A[1] + dNdX[a][2] * A[3] ) * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + ( dNdX[a][0] * A[4] + dNdX[a][1] * A[3] + dNdX[a][2] * A[2] ) * weight;
      }
    }
  }

  // diagonal second-order tensor
  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[3],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A[0] * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A[1] * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A[2] * weight;
      }
    }
  }

  // scalar*identity second-order tensor
  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const A,
                       real64 const & Np,
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A * weight;
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::SymmetricGradient,
               PDEUtilities::DifferentialOperator::SymmetricGradient >
{
  // Fourth-order tensor A
  template< int numTestDOF, int numTrialDOF, typename DISCRETIZATION_OPS >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNvdX)[numTestDOF/3][3],
                       DISCRETIZATION_OPS const & A,
                       real64 const (&dNudX)[numTrialDOF/3][3],
                       real64 const weight )
  {
    GEOS_UNUSED_VAR( dNvdX );
    const_cast< DISCRETIZATION_OPS & >(A).template BTDB< numTestDOF/3 >( dNudX, weight, mat );
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Divergence >
{
  template< int numTestDOF, int numTrialDOF >
  GEOS_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&Nv)[numTestDOF/3],
                       real64 const (&A)[3],
                       real64 const (&dNudX)[numTestDOF/3][3],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF/3; ++b )
      {
        mat[a*3+0][b*3+0] = mat[a*3+0][b*3+0] + Nv[a] * A[0] * dNudX[b][0] * weight;
        mat[a*3+0][b*3+1] = mat[a*3+0][b*3+1] + Nv[a] * A[0] * dNudX[b][1] * weight;
        mat[a*3+0][b*3+2] = mat[a*3+0][b*3+2] + Nv[a] * A[0] * dNudX[b][2] * weight;
        mat[a*3+1][b*3+0] = mat[a*3+1][b*3+0] + Nv[a] * A[1] * dNudX[b][0] * weight;
        mat[a*3+1][b*3+1] = mat[a*3+1][b*3+1] + Nv[a] * A[1] * dNudX[b][1] * weight;
        mat[a*3+1][b*3+2] = mat[a*3+1][b*3+2] + Nv[a] * A[1] * dNudX[b][2] * weight;
        mat[a*3+2][b*3+0] = mat[a*3+2][b*3+0] + Nv[a] * A[2] * dNudX[b][0] * weight;
        mat[a*3+2][b*3+1] = mat[a*3+2][b*3+1] + Nv[a] * A[2] * dNudX[b][1] * weight;
        mat[a*3+2][b*3+2] = mat[a*3+2][b*3+2] + Nv[a] * A[2] * dNudX[b][2] * weight;
      }
    }
  }

};

/**
 * @brief Generic bilinear form template to assemble elemental matrices.
 *
 * This function compute the elemental matrix for a grid cell \f$ \Omega^e \f$
 * with entries defined by the bilinear form
 *
 * \f[
 *   \int_{\Omega^e}
 *     \texttt{OP}_V ( v^h ) \; \cdot \;
 *     \mathsf{T}  \; \cdot \;
 *     \texttt{OP}_U ( u^h ) \;
 *     \mathrm{d}\Omega
 * \f]
 *
 * where \f$ v^h \f$ and \f$ u^h \f$ are finite element functions belonging to
 * the test (@p V) and  trial (@p U) space, respectively, \f$ \texttt{OP}_V \f$
 * and \f$ \texttt{OP}_U \f$ denote differential operators (zero- or
 * first-order) applied to the test and trial function, respectively, and
 * \f$ \mathsf{T} \f$ indicates a tensor.
 *
 * @tparam V Test function space.
 * @tparam U Trial function space.
 * @tparam OP_V Differential operator applied to functions in @p V.
 * @tparam OP_U Differential operator applied to functions in @p U.
 * @tparam MATRIX Derived matrix type.
 * @tparam FE_VALUES_V Derived test shape function values or derivatives type.
 * @tparam TENSOR Derived tensor type.
 * @tparam FE_VALUES_U Derived trial shape function values or derivatives type.
 * @param matrix The elemental matrix.
 * @param feValuesV Test shape function values or derivatives.
 * @param tensor The tensor defining the bilinear form.
 * @param feValuesU Trial shape function values or derivatives.
 * @param weight Quadrature weight.
 *
 */
template< PDEUtilities::FunctionSpace V,
          PDEUtilities::FunctionSpace U,
          PDEUtilities::DifferentialOperator OP_V,
          PDEUtilities::DifferentialOperator OP_U,
          typename MATRIX,
          typename FE_VALUES_V,
          typename TENSOR,
          typename FE_VALUES_U >
GEOS_HOST_DEVICE
static void compute( MATRIX && matrix,
                     FE_VALUES_V const & feValuesV,
                     TENSOR const & tensor,
                     FE_VALUES_U const & feValuesU,
                     real64 const weight )
{
  Helper< V, U, OP_V, OP_U >::compute( matrix, feValuesV, tensor, feValuesU, weight );
}

} // namespace BilinearFormUtilities

} // namespace geos

#endif //GEOS_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_
