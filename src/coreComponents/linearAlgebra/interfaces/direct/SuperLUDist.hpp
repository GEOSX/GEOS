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
 * @file SuperLUDist.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_SUPERLU_DIST_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_SUPERLU_DIST_HPP_

#include "common/DataTypes.hpp"
#include "common/LinearSolverBase.hpp"
#include "common/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"

#include <memory>

namespace geos
{

/// Forward declaration for SuperLU_Dist data struct
struct SuperLUDistData;

/**
 * @brief Wrapper for SuperLU_Dist parallel direct solver.
 * @tparam LAI type of linear algebra interface
 */
template< typename LAI >
class SuperLUDist final : public LinearSolverBase< LAI >
{
public:

  /// Alias for base type
  using Base = LinearSolverBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor with parameters.
   * @param[in] params the linear solver parameters
   */
  explicit SuperLUDist( LinearSolverParameters params );

  /**
   * @brief Destructor
   */
  virtual ~SuperLUDist() override;

  using Base::ready;
  using Base::matrix;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void setup( Matrix const & mat ) override;

  /**
   * @brief Apply operator to a vector, <tt>dst = this(src)</tt>.
   * @param src input vector (x)
   * @param dst output vector (b)
   *
   * @warning @p src and @p dst cannot alias the same vector.
   */
  virtual void apply( Vector const & src, Vector & dst ) const override;

  /**
   * @brief Clean up the preconditioner setup.
   */
  virtual void clear() override;

  /**
   * @brief Solve the system with a particular rhs
   * @param [in] rhs system right hand side
   * @param [out] sol system solution
   */
  virtual void solve( Vector const & rhs, Vector & sol ) const override;

private:

  using Base::m_params;
  using Base::m_result;

  /**
   * @brief Convert input parameters into SuperLU_Dist options.
   */
  void setOptions();

  /**
   * @brief Perform symbolic/numeric factorization of the matrix.
   */
  void factorize();

  /**
   * @brief Estimates the condition number of the matrix using LU factors.
   * @return the estimated condition number
   */
  real64 estimateConditionNumberBasic() const;

  /**
   * @brief Estimates the condition number of the matrix with Arnoldi iteration.
   * @return the estimated condition number
   */
  real64 estimateConditionNumberAdvanced() const;

  /// Internal SuperLU_Dist data
  std::unique_ptr< SuperLUDistData > m_data;

  /// Condition number estimation (cached)
  mutable real64 m_condEst;
};

}

#endif /*GEOS_LINEARALGEBRA_INTERFACES_SUPERLU_DIST_HPP_*/
