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
 * @file LinearSolverResult.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_LINEARSOLVERRESULT_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_LINEARSOLVERRESULT_HPP_

#include "common/DataTypes.hpp"

#include <vector>

namespace geos
{

/**
 * @brief Results/stats of a linear solve
 *
 * Lightweight struct to hold stats of a linear solve, such as
 * convergence flag, number of iterations and residual norms.
 */
struct LinearSolverResult
{
  /**
   * @brief Status of the linear solve.
   */
  enum class Status
  {
    InProgress,
    Success,
    NotConverged,
    Breakdown
  };

  /// convergence flag
  Status status = Status::InProgress;

  /// Number of solver iterations performed
  integer numIterations = 0;

  /// Final relative residual norm
  real64 residualReduction = 0.0;

  /// Setup time (in seconds) for preconditioners and/or direct factorizations
  real64 setupTime = 0.0;

  /// Solve time (in seconds) exclusive of setup costs
  real64 solveTime = 0.0;

  /**
   * @brief Check whether the last solve was successful.
   * @return @p true if last solve was successful, @p false otherwise
   */
  bool success() const
  {
    return status == Status::Success;
  }

  /**
   * @brief Check whether the last solve brokedown.
   * @return @p true if last solve brokedown, @p false otherwise
   */
  bool breakdown() const
  {
    return status == Status::Breakdown;
  }
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( LinearSolverResult::Status,
              "InProgress",
              "Success",
              "NotConverged",
              "Breakdown" );

}

#endif //GEOS_LINEARALGEBRA_UTILITIES_LINEARSOLVERRESULT_HPP_
