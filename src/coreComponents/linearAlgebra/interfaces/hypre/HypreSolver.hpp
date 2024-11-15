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
 * @file HypreSolver.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPRESOLVER_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPRESOLVER_HPP_

#include "linearAlgebra/interfaces/hypre/HypreInterface.hpp"
#include "linearAlgebra/interfaces/hypre/HyprePreconditioner.hpp"
#include "common/LinearSolverBase.hpp"

namespace geos
{

class DofManager;

/// Forward-declared struct that hosts pointers to preconditioner functions
struct HypreSolverWrapper;

/**
 * @brief This class creates and provides basic support for Hypre solvers.
 */
class HypreSolver final : public LinearSolverBase< HypreInterface >
{
public:

  /// Alias for base type
  using Base = LinearSolverBase< HypreInterface >;

  /**
   * @brief Solver constructor, with parameter list reference
   * @param[in] parameters structure containing linear solver parameters
   */
  explicit HypreSolver( LinearSolverParameters parameters );

  /**
   * @brief Destructor.
   */
  virtual ~HypreSolver() override;

  /**
   * @copydoc PreconditionerBase<PetscInterface>::setup
   */
  virtual void setup( HypreMatrix const & mat ) override;

  /**
   * @copydoc PreconditionerBase<PetscInterface>::apply
   */
  virtual void apply( HypreVector const & src,
                      HypreVector & dst ) const override;

  /**
   * @copydoc LinearSolverBase<PetscInterface>::solve
   */
  virtual void solve( HypreVector const & rhs,
                      HypreVector & sol ) const override;

  /**
   * @copydoc PreconditionerBase<PetscInterface>::clear
   */
  virtual void clear() override;

private:

  /**
   * @brief Perform the solve.
   * @param rhs right-hand side vector
   * @param sol solution vector
   * @return the error code from the hypre call
   */
  int doSolve( HypreVector const & rhs, HypreVector & sol ) const;

  using Base::m_params;
  using Base::m_result;

  /// Preconditioner
  HyprePreconditioner m_precond;

  /// Pointers to hypre functions for the krylov solver
  std::unique_ptr< HypreSolverWrapper > m_solver;
};

} // end geos namespace

#endif /* HYPRESOLVER_HPP_ */
