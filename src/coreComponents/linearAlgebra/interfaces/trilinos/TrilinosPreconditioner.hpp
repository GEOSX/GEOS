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
 * @file TrilinosPreconditioner.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_TRILINOSPRECONDITIONER_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_TRILINOSPRECONDITIONER_HPP_

#include "common/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <memory>

class Epetra_Operator;

namespace geos
{

/**
 * @brief Wrapper around Trilinos-based preconditioners.
 */
class TrilinosPreconditioner final : public PreconditionerBase< TrilinosInterface >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< TrilinosInterface >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   */
  explicit TrilinosPreconditioner( LinearSolverParameters params );

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   * @param nearNullKernel the user-provided near null kernel
   */
  TrilinosPreconditioner( LinearSolverParameters params,
                          array1d< Vector > const & nearNullKernel );

  /**
   * @brief Destructor.
   */
  virtual ~TrilinosPreconditioner() override;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void setup( Matrix const & mat ) override;

  /**
   * @brief Apply operator to a vector
   * @param src Input vector (x).
   * @param dst Output vector (b).
   *
   * @warning @p src and @p dst cannot alias the same vector.
   */
  virtual void apply( Vector const & src, Vector & dst ) const override;

  virtual void clear() override;

  /**
   * @brief Access the underlying Epetra preconditioning operator.
   * @return reference to the Epetra operator
   */
  Epetra_Operator const & unwrapped() const;

  /**
   * @copydoc unwrapped() const
   */
  Epetra_Operator & unwrapped();

private:

  /**
   * @brief Setup additional preconditioning matrix if necessary.
   * @param mat the source matrix
   * @return reference to the matrix that should be used to setup the main preconditioner
   */
  EpetraMatrix const & setupPreconditioningMatrix( EpetraMatrix const & mat );

  /// Parameters for all preconditioners
  LinearSolverParameters m_params;

  /// Preconditioning matrix (if different from input matrix)
  /// Note: must be declared before (destroyed after) m_precond
  EpetraMatrix m_precondMatrix;

  /// Pointer to the Trilinos implementation
  std::unique_ptr< Epetra_Operator > m_precond;

  /// Null space vectors
  array2d< real64 > m_nullSpacePointer;
};

}

#endif //GEOS_LINEARALGEBRA_INTERFACES_TRILINOSPRECONDITIONER_HPP_
