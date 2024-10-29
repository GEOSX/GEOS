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
 * @file LinearSolverBase.hpp
 */
#ifndef GEOS_LINEARSOLVERBASE_HPP
#define GEOS_LINEARSOLVERBASE_HPP

#include "PreconditionerBase.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"
#include "common/Stopwatch.hpp"

namespace geos
{

/**
 * @brief Simple interface for linear solvers that allows to extract solution results.
 * @tparam LAI linear algebra interface type
 */
template< typename LAI >
class LinearSolverBase : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param params solver parameters
   */
  explicit LinearSolverBase( LinearSolverParameters params )
    : m_params( std::move( params ) )
  {}

  /**
   * @brief Solve preconditioned system
   * @param [in] rhs system right hand side.
   * @param [in,out] sol system solution (input = initial guess, output = solution).
   */
  virtual void solve( Vector const & rhs, Vector & sol ) const = 0;

  /**
   * @brief @return parameters of the solver.
   */
  LinearSolverParameters const & parameters() const
  {
    return m_params;
  }

  /**
   * @brief @return result of the most recent solve.
   */
  LinearSolverResult const & result() const
  {
    return m_result;
  }

protected:

  /// Parameters for the solver
  LinearSolverParameters m_params;

  /// Result of most recent solve (status, timings)
  mutable LinearSolverResult m_result;
};

}

#endif //GEOS_LINEARSOLVERBASE_HPP
