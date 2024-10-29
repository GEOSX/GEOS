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

#ifndef GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERJACOBI_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERJACOBI_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"

namespace geos
{

/**
 * @brief Common interface for identity preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerJacobi : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void setup( Matrix const & mat ) override
  {
    GEOS_LAI_ASSERT( mat.ready() );
    m_diagInv.createWithLocalSize( mat.numLocalRows(), mat.comm() );
    mat.extractDiagonal( m_diagInv );
    m_diagInv.reciprocal();
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
    m_diagInv.reset();
  }

  /**
   * @brief Get the number of global rows.
   * @return Number of global rows in the operator.
   */
  virtual globalIndex numGlobalRows() const override final
  {
    GEOS_LAI_ASSERT( m_diagInv.ready() );
    return m_diagInv.globalSize();
  }

  /**
   * @brief Get the number of global columns.
   * @return Number of global columns in the operator.
   */
  virtual globalIndex numGlobalCols() const override final
  {
    GEOS_LAI_ASSERT( m_diagInv.ready() );
    return m_diagInv.globalSize();
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
    GEOS_LAI_ASSERT( m_diagInv.ready() );
    GEOS_LAI_ASSERT_EQ( this->numGlobalRows(), dst.globalSize() );
    GEOS_LAI_ASSERT_EQ( this->numGlobalCols(), src.globalSize() );

    m_diagInv.pointwiseProduct( src, dst );
  }

private:

  /// The diagonal of the matrix
  Vector m_diagInv;
};

}

#endif //GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERJACOBI_HPP_
