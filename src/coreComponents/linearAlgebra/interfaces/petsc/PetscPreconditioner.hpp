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
 * @file PetscPreconditioner.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_PETSCPRECONDITIONER_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_PETSCPRECONDITIONER_HPP_

#include "common/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/petsc/PetscInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

/**
 * @name PETSc forward declarations.
 *
 * Forward declare PETSc's solver structs and pointer aliases in order
 * to avoid including PETSc headers and leaking into the rest of GEOSX.
 */
///@{

/// PETSc preconditioner struct forward declaration
extern "C" struct _p_PC;
extern "C" struct _p_MatNullSpace;

///@}

namespace geos
{

/**
 * @brief Wrapper around PETSc-based preconditioners.
 */
class PetscPreconditioner final : public PreconditionerBase< PetscInterface >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< PetscInterface >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   */
  explicit PetscPreconditioner( LinearSolverParameters params );

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   * @param nearNullKernel the user-provided near null kernel
   */
  PetscPreconditioner( LinearSolverParameters params, array1d< Vector > const & nearNullKernel );

  /**
   * @brief Destructor.
   */
  virtual ~PetscPreconditioner() override;

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

private:

  /// Preconditioner pointer alias
  using PC = _p_PC *;

public:

  /**
   * @brief Access the underlying implementation.
   * @return the wrapped PETSc preconditioner
   */
  PC const & unwrapped() const;

private:

  /// Near null space pointer alias
  using MatNullSpace = _p_MatNullSpace *;

  /**
   * @brief Setup additional preconditioning matrix if necessary.
   * @param mat the source matrix
   * @return reference to the matrix that should be used to setup the main preconditioner
   */
  PetscMatrix const & setupPreconditioningMatrix( PetscMatrix const & mat );

  /// Parameters for all preconditioners
  LinearSolverParameters m_params;

  /// Pointer to the PETSc implementation
  PC m_precond;

  /// Preconditioning matrix (if different from input matrix)
  PetscMatrix m_precondMatrix;

  /// Pointer to the near null space
  MatNullSpace m_nullsp;
};

}

#endif //GEOS_LINEARALGEBRA_INTERFACES_PETSCPRECONDITIONER_HPP_
