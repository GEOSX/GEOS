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
 * @file HyprePreconditioner.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREPRECONDITIONER_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREPRECONDITIONER_HPP_

#include "common/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/hypre/HypreInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <memory>

namespace geos
{

/// Forward-declared struct that hosts pointers to preconditioner functions
struct HyprePrecWrapper;

/// Forward-declared struct that hosts preconditioner MGR data
struct HypreMGRData;

/// Forward-declared struct that hosts null space data
struct HypreNullSpace;

/**
 * @brief Wrapper around hypre-based preconditioners.
 */
class HyprePreconditioner final : public PreconditionerBase< HypreInterface >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< HypreInterface >;

  /// Allow for partial overload of Base::setup()
  using Base::setup;

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   */
  explicit HyprePreconditioner( LinearSolverParameters params );

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   * @param nearNullKernel the user-provided near null kernel
   */
  HyprePreconditioner( LinearSolverParameters params,
                       arrayView1d< HypreVector > const & nearNullKernel );

  /**
   * @brief Destructor.
   */
  virtual ~HyprePreconditioner() override;

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
   * @brief Access the underlying implementation.
   * @return reference to container of preconditioner functions.
   *
   * Intended for use by HypreSolver.
   */
  HyprePrecWrapper const & unwrapped() const;

  /**
   * @brief @return time spent setting up separate component matrix.
   */
  real64 componentFilterTime() const
  {
    return m_componentFilterTime;
  }

  /// @return time to construct restrictor matrix.
  real64 makeRestrictorTime() const
  {
    return m_makeRestrictorTime;
  }

  /// @return time to apply restrictor matrix.
  real64 computeAuuTime() const
  {
    return m_computeAuuTime;
  }

private:

  /**
   * @brief Create the preconditioner object.
   * @param dofManager (optional) pointer to DofManager for the linear system
   */
  void create( DofManager const * const dofManager );

  /**
   * @brief Setup additional preconditioning matrix if necessary.
   * @param mat the source matrix
   * @return reference to the matrix that should be used to setup the main preconditioner
   */
  HypreMatrix const & setupPreconditioningMatrix( HypreMatrix const & mat );

  /// Parameters for all preconditioners
  LinearSolverParameters m_params;

  /// Preconditioning matrix (if different from input matrix)
  HypreMatrix m_precondMatrix;

  /// Pointers to hypre preconditioner and corresponding functions
  std::unique_ptr< HyprePrecWrapper > m_precond;

  /// MGR-specific data
  std::unique_ptr< HypreMGRData > m_mgrData;

  /// Null space vectors
  std::unique_ptr< HypreNullSpace > m_nullSpace;

  /// Timing of separate component matrix construction
  real64 m_componentFilterTime = 0.0;

  /// Timing of the restrictor matrix construction
  real64 m_makeRestrictorTime = 0.0;

  /// Timing of the cost of applying the restrictor matrix to the system
  real64 m_computeAuuTime = 0.0;
};

}

#endif //GEOS_LINEARALGEBRA_INTERFACES_HYPREPRECONDITIONER_HPP_
