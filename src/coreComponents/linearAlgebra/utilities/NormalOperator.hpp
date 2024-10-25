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
 * @file NormalOperator.hpp
 */
#ifndef GEOS_LINEARALGEBRA_NORMALOPERATOR_HPP_
#define GEOS_LINEARALGEBRA_NORMALOPERATOR_HPP_

#include "common/LinearOperator.hpp"

namespace geos
{

/**
 * @brief Wraps a matrix A and represents A^T * A as a linear operator.
 * @tparam LAI the linear algebra interface
 */
template< typename LAI >
class NormalOperator : public LinearOperator< typename LAI::ParallelVector >
{
public:

  /// Alias for base type
  using Base = LinearOperator< typename LAI::ParallelVector >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename LAI::ParallelMatrix;

  /**
   * @brief Constructor
   * @param mat the underlying matrix (must outlive this operator)
   */
  explicit NormalOperator( Matrix const & mat )
    : m_matrix( mat )
  {}

  /**
   * @brief Destructor.
   */
  virtual ~NormalOperator() override = default;

  /**
   * @brief Apply operator to a vector.
   * @param src input vector
   * @param dst output vector
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  void apply( Vector const & src, Vector & dst ) const override
  {
    m_matrix.gemv( 1.0, src, 0.0, dst, false );
    m_matrix.gemv( 1.0, dst, 0.0, dst, true );
  }

  /**
   * @brief @return the global number of rows
   */
  globalIndex numGlobalRows() const override
  {
    return m_matrix.numGlobalCols();
  }

  /**
   * @brief @return the global number of columns
   */
  globalIndex numGlobalCols() const override
  {
    return m_matrix.numGlobalCols();
  }

  /**
   * @brief @return the local number of rows
   */
  localIndex numLocalRows() const override
  {
    return m_matrix.numLocalCols();
  }

  /**
   * @brief @return the local number of columns
   */
  localIndex numLocalCols() const override
  {
    return m_matrix.numLocalCols();
  }

  /**
   * @brief @return the communicator
   */
  MPI_Comm comm() const override
  {
    return m_matrix.comm();
  }

private:

  /// the matrix object
  Matrix const & m_matrix;
};

} // namespace geos

#endif //GEOS_LINEARALGEBRA_NORMALOPERATOR_HPP_
