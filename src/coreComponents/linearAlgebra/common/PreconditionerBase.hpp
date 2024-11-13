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

#ifndef GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_

#include "linearAlgebra/common/common.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"

namespace geos
{

class DofManager;

/**
 * @brief Common interface for preconditioning operators
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerBase : public LinearOperator< typename LAI::ParallelVector >
{
public:

  PreconditionerBase() = default;

  virtual ~PreconditionerBase() = default;

  /// Alias for base type
  using Base = LinearOperator< typename LAI::ParallelVector >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename LAI::ParallelMatrix;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void setup( Matrix const & mat )
  {
    GEOS_LAI_ASSERT( mat.ready() );
    GEOS_LAI_ASSERT_MSG( mat.numLocalRows() == mat.numLocalCols(), "Matrix must be square" );
    m_mat = &mat;
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
  virtual void clear()
  {
    m_mat = nullptr;
  }

  /**
   * @brief Get the number of global rows.
   * @return Number of global rows in the operator.
   */
  virtual globalIndex numGlobalRows() const override
  {
    return m_mat->numGlobalRows();
  }

  /**
   * @brief Get the number of global columns.
   * @return Number of global columns in the operator.
   */
  virtual globalIndex numGlobalCols() const override
  {
    return m_mat->numGlobalCols();
  }

  /**
   * @brief Get the number of local rows.
   * @return Number of local rows in the operator.
   */
  virtual localIndex numLocalRows() const override
  {
    return m_mat->numLocalRows();
  }

  /**
   * @brief Get the number of local columns.
   * @return Number of local columns in the operator.
   */
  virtual localIndex numLocalCols() const override
  {
    return m_mat->numLocalCols();
  }

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator passed in @p create...()
   *
   * @note when build without MPI, may return anything
   *       (MPI_Comm will be a mock type defined in MpiWrapper)
   */
  virtual MPI_Comm comm() const override
  {
    return m_mat->comm();
  }

  /**
   * @brief Chech if preconditioner is ready to use
   * @return @p true if compute() has been called but not clear().
   */
  bool ready() const
  {
    return m_mat != nullptr;
  }

  /**
   * @brief Access the matrix the preconditioner was computed from
   * @return reference to the matrix (user's repsonsibility to ensure it's still valid)
   */
  Matrix const & matrix() const
  {
    GEOS_LAI_ASSERT( ready() );
    return *m_mat;
  }

  /**
   * @brief Check whether the preconditioner is available in matrix (explicit) form.
   * @return if the preconditioner is available in explicit form
   */
  virtual bool hasPreconditionerMatrix() const
  {
    return false;
  }

  /**
   * @brief Access the preconditioner in matrix form (whenever available). It must be
   *        overridden by the specific preconditioner
   * @return reference to the preconditioner matrix
   */
  virtual Matrix const & preconditionerMatrix() const
  {
    GEOS_ERROR( "PreconditionerBase::preconditionerMatrix called. This is not supposed to happen."
                "Check the value of hasPreconditionerMatrix() before accessing this function." );
    // This is here just to be able to compile ...
    return *m_mat;
  }

private:

  /// Pointer to the matrix
  Matrix const * m_mat{};
};

}

#endif //GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_
