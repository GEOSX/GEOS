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

#ifndef GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERBLOCKJACOBI_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERBLOCKJACOBI_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

namespace geos
{

/**
 * @brief Common interface for identity preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerBlockJacobi : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param blockSize the size of block diagonal matrices.
   */
  PreconditionerBlockJacobi( localIndex const & blockSize = 0 )
    : m_blockDiag{}
  {
    m_blockSize = blockSize;
  }

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void setup( Matrix const & mat ) override
  {
    GEOS_LAI_ASSERT( mat.ready() );
    GEOS_LAI_ASSERT_GT( m_blockSize, 0 );
    GEOS_LAI_ASSERT_EQ( mat.numLocalRows() % m_blockSize, 0 );
    GEOS_LAI_ASSERT_EQ( mat.numLocalCols() % m_blockSize, 0 );

    PreconditionerBase< LAI >::setup( mat );

    m_blockDiag.createWithLocalSize( mat.numLocalRows(), mat.numLocalCols(), m_blockSize, mat.comm() );
    m_blockDiag.open();

    array1d< globalIndex > idxBlk( m_blockSize );
    array2d< real64 > values( m_blockSize, m_blockSize );
    array2d< real64 > valuesInv( m_blockSize, m_blockSize );
    array1d< globalIndex > cols;
    array1d< real64 > vals;
    for( globalIndex i = mat.ilower(); i < mat.iupper(); i += m_blockSize )
    {
      values.zero();
      for( localIndex j = 0; j < m_blockSize; ++j )
      {
        globalIndex const iRow = i + LvArray::integerConversion< globalIndex >( j );
        idxBlk[j] = iRow;
        localIndex const rowLength = mat.rowLength( iRow );
        cols.resize( rowLength );
        vals.resize( rowLength );
        mat.getRowCopy( iRow, cols, vals );
        for( localIndex k = 0; k < rowLength; ++k )
        {
          localIndex const jCol = LvArray::integerConversion< localIndex >( cols[k]-i );
          if( cols[k] >= i && cols[k] < i+LvArray::integerConversion< globalIndex >( m_blockSize ) )
          {
            values( j, jCol ) = vals[k];
          }
        }
      }
      BlasLapackLA::matrixInverse( values, valuesInv );
      m_blockDiag.insert( idxBlk, idxBlk, valuesInv );
    }
    m_blockDiag.close();
  }

  /**
   * @brief Clean up the preconditioner setup.
   *
   * Releases memory used and allows the matrix to be deleted cleanly.
   * This method should be called before the matrix used to compute the preconditioner
   * goes out of scope or is re-created. Some implementations require the matrix
   * to outlive the preconditioner (for example, Trilinos/ML may crash the program if
   * deleted after the matrix).
   *
   * @note Should be properly overridden in derived classes, which may call this method.
   */
  virtual void clear() override
  {
    m_blockDiag.reset();
  }

  /**
   * @brief Apply operator to a vector.
   *
   * @param src Input vector (src).
   * @param dst Output vector (dst).
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override
  {
    GEOS_LAI_ASSERT( m_blockDiag.ready() );
    GEOS_LAI_ASSERT_EQ( this->numGlobalRows(), dst.globalSize() );
    GEOS_LAI_ASSERT_EQ( this->numGlobalCols(), src.globalSize() );

    m_blockDiag.apply( src, dst );
  }

  /**
   * @brief Whether the preconditioner is available in matrix form
   * @return true: explicit form is available
   */
  virtual bool hasPreconditionerMatrix() const override
  {
    GEOS_LAI_ASSERT( m_blockDiag.ready() );
    return true;
  }

  /**
   * @brief Access the preconditioner in matrix form
   * @return reference to the preconditioner matrix
   */
  virtual Matrix const & preconditionerMatrix() const override
  {
    GEOS_LAI_ASSERT( m_blockDiag.ready() );
    return m_blockDiag;
  }

private:

  /// The preconditioner matrix
  Matrix m_blockDiag;

  /// Block size
  localIndex m_blockSize = 0;
};

}

#endif //GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERBLOCKJACOBI_HPP_
