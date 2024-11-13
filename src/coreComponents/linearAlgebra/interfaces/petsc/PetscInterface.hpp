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
 * @file PetscInterface.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_

#include "linearAlgebra/interfaces/petsc/PetscVector.hpp"
#include "linearAlgebra/interfaces/petsc/PetscMatrix.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"
#include "linearAlgebra/common/LinearSolverBase.hpp"

#include <memory>

namespace geos
{

/**
 * @class PetscInterface
 * @brief This class holds aliases based on the Petsc library.
 */
struct PetscInterface
{
  /**
   * @brief Initializes the Petsc library.
   */
  static void initialize();

  /**
   * @brief Finalizes the Petsc library.
   */
  static void finalize();

  /**
   * @brief Create a petsc-based solver object.
   * @param params the preconditioner parameters
   * @return owning pointer to the newly created solver
   */
  static std::unique_ptr< LinearSolverBase< PetscInterface > >
  createSolver( LinearSolverParameters params );

  /**
   * @brief Create a PETSc-based preconditioner object.
   * @param params the parameters for preconditioner
   * @return owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< PetscInterface > >
  createPreconditioner( LinearSolverParameters params );

  /**
   * @brief Create a PETSc-based preconditioner object.
   * @param params the parameters for preconditioner
   * @param nearNullKernel the user-provided near null kernel
   * @return owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< PetscInterface > >
  createPreconditioner( LinearSolverParameters params,
                        array1d< PetscVector > const & nearNullKernel );

  /// Alias for PetscMatrix
  using ParallelMatrix = PetscMatrix;
  /// Alias for PetscVector
  using ParallelVector = PetscVector;
};

} /* namespace geos */

#endif /*GEOS_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_*/
