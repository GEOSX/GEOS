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
 * @file SeparateComponentPreconditioner.hpp
 */

#ifndef GEOS_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_

#include "linearAlgebra/common/PreconditionerBase.hpp"

#include <memory>

namespace geos
{

/**
 * @brief Separate component filter implemented as a compound preconditioner.
 * @tparam LAI linear algebra interface to use
 */
template< typename LAI >
class SeparateComponentPreconditioner : public PreconditionerBase< LAI >
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
   * @param numComp number of components in the field
   * @param precond the actual preconditioner to apply to filtered matrix (ownership transferred)
   */
  SeparateComponentPreconditioner( localIndex const numComp,
                                   std::unique_ptr< PreconditionerBase< LAI > > precond );

  /**
   * @brief Destructor.
   */
  virtual ~SeparateComponentPreconditioner() override;

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

  virtual bool hasPreconditionerMatrix() const override
  {
    return m_precond->hasPreconditionerMatrix();
  }

  virtual Matrix const & preconditionerMatrix() const override
  {
    return m_precond->preconditionerMatrix();
  }

  /**
   * @brief @return reference to the filtered matrix
   */
  Matrix const & separateComponentMatrix() const
  {
    GEOS_LAI_ASSERT( Base::ready() );
    return m_matSC;
  }

  /**
   * @brief Access to the nested preconditioner.
   * @return reference to the preconditioner passed at construction
   */
  PreconditionerBase< LAI > const & innerPrecond() const
  {
    return *m_precond;
  }

private:

  /// Number of components in the matrix
  localIndex m_numComp;

  /// Separate component approximation of the original matrix
  Matrix m_matSC;

  /// Actual preconditioner
  std::unique_ptr< PreconditionerBase< LAI > > m_precond;
};

}

#endif //GEOS_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_
