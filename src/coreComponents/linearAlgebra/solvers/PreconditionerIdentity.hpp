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

#ifndef GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/common/PreconditionerBase.hpp"

namespace geos
{

/**
 * @brief Common interface for identity preconditioning operator
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerIdentity : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  virtual ~PreconditionerIdentity() = default;

  /**
   * @brief Apply operator to a vector.
   *
   * @param src Input vector (src).
   * @param dst Output vector (dst).
   */
  virtual void apply( Vector const & src,
                      Vector & dst ) const override
  {
    GEOS_LAI_ASSERT_EQ( this->numGlobalRows(), dst.globalSize() );
    GEOS_LAI_ASSERT_EQ( this->numGlobalCols(), src.globalSize() );
    dst.copy( src );
  }
};

}

#endif //GEOS_LINEARALGEBRA_SOLVERS_PRECONDITIONERIDENTITY_HPP_
