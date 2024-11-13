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
 * @file PetscSolver.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_PETSCSOLVER_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_PETSCSOLVER_HPP_

#include "linearAlgebra/interfaces/petsc/PetscInterface.hpp"
#include "linearAlgebra/interfaces/petsc/PetscPreconditioner.hpp"
#include "common/LinearSolverBase.hpp"

/// Forward declare PETSC's solver struct
extern "C" struct _p_KSP;

namespace geos
{

/**
 * @brief This class creates and provides basic support for PETSc solvers.
 */
class PetscSolver final : public LinearSolverBase< PetscInterface >
{
public:

  /// Alias for base type
  using Base = LinearSolverBase< PetscInterface >;

  /**
   * @brief Solver constructor, with parameter list reference
   * @param[in] parameters structure containing linear solver parameters
   */
  explicit PetscSolver( LinearSolverParameters parameters );

  /**
   * @brief Destructor.
   */
  virtual ~PetscSolver();

  /**
   * @copydoc PreconditionerBase<PetscInterface>::setup
   */
  virtual void setup( PetscMatrix const & mat ) override;

  /**
   * @copydoc PreconditionerBase<PetscInterface>::apply
   */
  virtual void apply( PetscVector const & src,
                      PetscVector & dst ) const override;

  /**
   * @copydoc LinearSolverBase<PetscInterface>::solve
   */
  virtual void solve( PetscVector const & rhs,
                      PetscVector & sol ) const override;

  /**
   * @copydoc PreconditionerBase<PetscInterface>::clear
   */
  virtual void clear() override;

private:

  using KSP = struct _p_KSP *;

  using Base::m_params;
  using Base::m_result;

  /// Preconditioner
  PetscPreconditioner m_precond;

  /// Krylov solver instance
  KSP m_solver{};
};

} // end geos namespace

#endif /* PETSCSOLVER_HPP_ */
