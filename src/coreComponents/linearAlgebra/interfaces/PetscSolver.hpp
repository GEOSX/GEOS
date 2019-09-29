/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PetscSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_PETSCSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_PETSCSOLVER_HPP_

namespace geosx
{

class PetscVector;
class PetscSparseMatrix;
class LinearSolverParameters;

/**
 * \class PetscSolver
 * \brief This class creates and provides basic support for PETSc solvers.
 */

class PetscSolver
{
public:

  /**
   * @brief Solver constructor, with parameter list reference
   *
   */
  PetscSolver( LinearSolverParameters const & parameters );

  /**
   * @brief Virtual destructor.
   *
   */
  virtual ~PetscSolver() = default;

  /**
   * @brief Solve system with an iterative solver.
   *
   * Solve Ax=b with A an PetscSparseMatrix, x and b PetscVector.
   */

  void solve( PetscSparseMatrix & mat,
              PetscVector & sol,
              PetscVector & rhs );

private:

  LinearSolverParameters const & m_parameters;

  void solve_direct( PetscSparseMatrix & mat,
                     PetscVector & sol,
                     PetscVector & rhs );

  void solve_krylov( PetscSparseMatrix & mat,
                     PetscVector & sol,
                     PetscVector & rhs );

};

} // end geosx namespace

#endif /* PETSCSOLVER_HPP_ */
