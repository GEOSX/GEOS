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
 * @file Arnoldi.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_ARNOLDI_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_ARNOLDI_HPP_

#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include "common/LinearOperator.hpp"

namespace geos
{

/**
 * @brief Function implementing the Arnoldi scheme to compute the largest eigenvalue
 * @tparam VECTOR vector type of the linear operator
 * @param op the operator whose largest eigenvalue is required
 * @param m the number of iterations (size of the Krylov subspace)
 * @return the largest eigenvalue
 */
template< typename VECTOR >
real64 ArnoldiLargestEigenvalue( LinearOperator< VECTOR > const & op, localIndex const m = 4 )
{
  localIndex const numGlobalRows = LvArray::integerConversion< localIndex >( op.numGlobalRows() );
  localIndex const numLocalRows = op.numLocalRows();
  localIndex const mInternal = std::min( numGlobalRows, m );

  // Initialize data structure (Hessenberg matrix and Krylov subspace)
  array2d< real64, MatrixLayout::ROW_MAJOR_PERM > H( mInternal+1, mInternal );
  array1d< VECTOR > V( mInternal + 1 );

  // Initial unitary vector
  V[0].create( numLocalRows, op.comm() );
  V[0].set( 1.0 / sqrt( static_cast< real64 >( numGlobalRows ) ) );

  for( localIndex j = 0; j < mInternal; ++j )
  {
    // Apply operator
    V[j+1].create( numLocalRows, op.comm() );
    op.apply( V[j], V[j+1] );
    // Arnoldi process
    for( localIndex i = 0; i <= j; ++i )
    {
      H( i, j ) = V[i].dot( V[j+1] );
      V[j+1].axpy( -H( i, j ), V[i] );
    }
    H( j+1, j ) = V[j+1].norm2();
    V[j+1].scale( 1.0 / H( j+1, j ) );
  }

  // Disregard the last entry and make the matrix square
  // Note: this is ok since we are using ROW_MAJOR_PERM
  H.resize( mInternal, mInternal );

  // Compute the eigenvalues
  array1d< std::complex< real64 > > lambda( mInternal );
  BlasLapackLA::matrixEigenvalues( H, lambda );

  // Find the largest eigenvalues
  real64 lambdaMax = 0.0;
  for( localIndex i = 0; i < mInternal; ++i )
  {
    lambdaMax = std::max( std::abs( lambda[i] ), lambdaMax );
  }

  return lambdaMax;
}

} // namespace geos

#endif //GEOS_LINEARALGEBRA_UTILITIES_ARNOLDI_HPP_
