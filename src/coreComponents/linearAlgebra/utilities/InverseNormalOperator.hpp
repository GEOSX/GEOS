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
 * @file InverseNormalOperator.hpp
 */
#ifndef GEOS_LINEARALGEBRA_UTILITIES_INVERSENORMALOPERATOR_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_INVERSENORMALOPERATOR_HPP_

#include "common/LinearOperator.hpp"
#include "common/PreconditionerBase.hpp"
#include "common/traits.hpp"

namespace geos
{

namespace internal
{

template< typename LAI, template< typename > class SOLVER >
bool constexpr HasApplyTranspose = traits::VectorBasedTraits< typename LAI::ParallelVector >::template HasMemberFunction_applyTranspose< SOLVER< LAI > >;

template< typename LAI, template< typename > class SOLVER >
class ExplicitTransposeInverse
{
public:

  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;
  using Solver = SOLVER< LAI >;

  ExplicitTransposeInverse( Matrix const & matrix, Solver const & solver ):
    m_transposeSolver( solver.parameters() )
  {
    matrix.transpose( m_transposeMatrix );
    m_transposeSolver.setup( m_transposeMatrix );
  }

  void apply( Solver const &, Vector const & src, Vector & dst ) const
  {
    m_transposeSolver.apply( src, dst );
  }

public:

  Matrix m_transposeMatrix;
  Solver m_transposeSolver;
};

template< typename LAI, template< typename > class SOLVER >
class ImplicitTransposeInverse
{
protected:

  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;
  using Solver = SOLVER< LAI >;

  ImplicitTransposeInverse( Matrix const &, Solver const & )
  {}

  void apply( Solver const & solver, Vector const & src, Vector & dst ) const
  {
    solver.applyTranspose( src, dst );
  }
};

template< typename LAI, template< typename > class SOLVER >
using TransposeInverse = std::conditional_t< HasApplyTranspose< LAI, SOLVER >,
                                             ImplicitTransposeInverse< LAI, SOLVER >,
                                             ExplicitTransposeInverse< LAI, SOLVER > >;

}

/**
 * @brief Wraps a matrix A and represents A^{-1} * A^{-T} as a linear operator.
 * @tparam LAI the linear algebra interface
 * @tparam SOLVER type of solver used to apply inverse of the matrix
 */
template< typename LAI, template< typename > class SOLVER >
class InverseNormalOperator : public LinearOperator< typename LAI::ParallelVector >
{
public:

  /// Alias for base type
  using Base = LinearOperator< typename LAI::ParallelVector >;

  /// Alias for matrix type
  using Matrix = typename LAI::ParallelMatrix;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for solver type
  using Solver = SOLVER< LAI >;

  /**
   * @brief Constructor
   * @param mat the underlying matrix (must outlive this operator)
   * @param solver the solver used to apply inverse of the matrix (must be already setup)
   */
  explicit InverseNormalOperator( Matrix const & mat,
                                  Solver const & solver ):
    m_matrix( mat ),
    m_solver( solver ),
    m_transposeInverse( mat, solver )
  {}

  /**
   * @brief Destructor.
   */
  virtual ~InverseNormalOperator() override = default;

  /**
   * @brief Apply operator to a vector.
   * @param src input vector (x)
   * @param dst output vector (b)
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  void apply( Vector const & src, Vector & dst ) const override
  {
    m_solver.apply( src, dst );
    m_transposeInverse.apply( m_solver, dst, dst );
  }

  /**
   * @brief @return the global number of rows
   */
  globalIndex numGlobalRows() const override
  {
    return m_matrix.numGlobalCols();
  }

  /**
   * @brief @return the global number of columns
   */
  globalIndex numGlobalCols() const override
  {
    return m_matrix.numGlobalCols();
  }

  /**
   * @brief @return the local number of rows
   */
  localIndex numLocalRows() const override
  {
    return m_matrix.numLocalCols();
  }

  /**
   * @brief @return the local number of columns
   */
  localIndex numLocalCols() const override
  {
    return m_matrix.numLocalCols();
  }

  /**
   * @brief @return the communicator
   */
  MPI_Comm comm() const override
  {
    return m_matrix.comm();
  }

private:

  /// Type of transposer helper (either explicit or implicit, depending on solver capabilities)
  using TransposeInverse = internal::TransposeInverse< LAI, SOLVER >;

  /// the matrix object
  Matrix const & m_matrix;

  /// the solver object
  Solver const & m_solver;

  /// the transposer object
  TransposeInverse m_transposeInverse;
};

} // namespace geos

#endif //GEOS_LINEARALGEBRA_UTILITIES_INVERSENORMALOPERATOR_HPP_
