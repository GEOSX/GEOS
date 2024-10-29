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
 * @file TrilinosSolver.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_TRILINOSSOLVER_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_TRILINOSSOLVER_HPP_

#include "linearAlgebra/interfaces/trilinos/TrilinosInterface.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosPreconditioner.hpp"
#include "common/LinearSolverBase.hpp"

#include <memory>

class AztecOO;

namespace geos
{

/**
 * @brief This class creates and provides basic support for AztecOO, Amesos and ML libraries.
 */
class TrilinosSolver final : public LinearSolverBase< TrilinosInterface >
{
public:

  /// Alias for base type
  using Base = LinearSolverBase< TrilinosInterface >;

  /**
   * @brief Solver constructor, with parameter list reference
   * @param[in] parameters structure containing linear solver parameters
   */
  explicit TrilinosSolver( LinearSolverParameters parameters );

  /**
   * @brief Destructor.
   */
  virtual ~TrilinosSolver() override;

  /**
   * @copydoc PreconditionerBase<PetscInterface>::setup
   */
  virtual void setup( EpetraMatrix const & mat ) override;

  /**
   * @copydoc PreconditionerBase<PetscInterface>::apply
   */
  virtual void apply( EpetraVector const & src,
                      EpetraVector & dst ) const override;

  /**
   * @copydoc LinearSolverBase<PetscInterface>::solve
   */
  virtual void solve( EpetraVector const & rhs,
                      EpetraVector & sol ) const override;

private:

  /**
   * @brief Perform the solve.
   * @param rhs right-hand side vector
   * @param sol solution vector
   * @return the error code from the AztecOO call
   */
  int doSolve( EpetraVector const & rhs, EpetraVector & sol ) const;

  using Base::m_params;
  using Base::m_result;

  /// Preconditioner
  TrilinosPreconditioner m_precond;

  /// AztecOO solver instance
  std::unique_ptr< AztecOO > m_solver;
};

} // end geos namespace

#endif /* TRILINOSSOLVER_HPP_ */
