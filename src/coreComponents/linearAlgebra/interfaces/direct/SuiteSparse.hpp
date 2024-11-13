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
 * @file SuiteSparse.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_SUITESPARSE_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_SUITESPARSE_HPP_

#include "common/DataTypes.hpp"
#include "common/LinearSolverBase.hpp"
#include "common/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"

#include <memory>

namespace geos
{

/// Forward declaration for SuiteSparse data struct
struct SuiteSparseData;

/**
 * @brief Wrapper for UMFPACK direct solver from SuiteSparse package.
 * @tparam LAI type of linear algebra interface
 */
template< typename LAI >
class SuiteSparse final : public LinearSolverBase< LAI >
{
public:

  /// Alias for base type
  using Base = LinearSolverBase< LAI >;

  /// Alias for vector type
  using Vector = typename PreconditionerBase< LAI >::Vector;

  /// Alias for matrix type
  using Matrix = typename PreconditionerBase< LAI >::Matrix;

  /**
   * @brief Constructor with parameters
   * @param[in] params the linear solver parameters
   */
  SuiteSparse( LinearSolverParameters params );

  /**
   * @brief Destructor
   */
  virtual ~SuiteSparse() override;

  using Base::ready;
  using Base::matrix;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void setup( Matrix const & mat ) override;

  /**
   * @brief Apply operator to a vector, <tt>dst = this(src)</tt>.
   * @param src input vector
   * @param dst output vector
   *
   * @warning @p src and @p dst cannot alias the same vector.
   */
  virtual void apply( Vector const & src, Vector & dst ) const override;

  /**
   * @brief Apply transpose operator to a vector, <tt>dst = this^T(src)</tt>.
   * @param src input vector
   * @param dst output vector
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  void applyTranspose( Vector const & src, Vector & dst ) const;

  /**
   * @brief Clean up the preconditioner setup.
   */
  virtual void clear() override;

  /**
   * @brief Solve the system with a particular rhs
   * @param [in] rhs system right hand side
   * @param [out] sol system solution
   */
  virtual void solve( Vector const & rhs, Vector & sol ) const override;

private:

  using Base::m_params;
  using Base::m_result;

  /**
   * @brief Perform the actual solution using LU factors.
   * @param b the rhs vector
   * @param x the solution vector
   * @param transpose whether to do a regular or a transpose solve
   */
  void doSolve( Vector const & b, Vector & x, bool transpose ) const;

  /**
   * @brief Estimates the condition number of the matrix using LU factors.
   * @return the estimated condition number
   */
  real64 estimateConditionNumberBasic() const;

  /**
   * @brief Estimates the condition number of the matrix with Arnoldi iteration.
   * @return the estimated condition number
   */
  real64 estimateConditionNumberAdvanced() const;

  /// Internal SuiteSparse data
  std::unique_ptr< SuiteSparseData > m_data;

  /// MPI rank carrying out the solution
  int m_workingRank;

  /// The export object used to gather/scatter matrices and vectors
  std::unique_ptr< typename Matrix::Export > m_export;

  /// condition number estimation
  mutable real64 m_condEst;
};

}

#endif /*GEOS_LINEARALGEBRA_INTERFACES_SUITESPARSE_HPP_*/
