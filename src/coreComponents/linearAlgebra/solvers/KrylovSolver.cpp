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
 * @file KrylovSolver.cpp
 */

#include "KrylovSolver.hpp"
#include "linearAlgebra/solvers/BicgstabSolver.hpp"
#include "linearAlgebra/solvers/CgSolver.hpp"
#include "linearAlgebra/solvers/GmresSolver.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geos
{

template< typename VECTOR >
KrylovSolver< VECTOR >::KrylovSolver( LinearSolverParameters params,
                                      LinearOperator< Vector > const & matrix,
                                      LinearOperator< Vector > const & precond )
  : Base(),
  m_params( std::move( params ) ),
  m_operator( matrix ),
  m_precond( precond )
{
  GEOS_ERROR_IF_LE_MSG( m_params.krylov.maxIterations, 0, "Krylov solver: max number of iteration must be positive." );
  GEOS_LAI_ASSERT_EQ( m_operator.numLocalRows(), m_operator.numLocalCols() );
  GEOS_LAI_ASSERT_EQ( m_operator.numLocalRows(), m_precond.numLocalRows() );
  GEOS_LAI_ASSERT_EQ( m_operator.numLocalCols(), m_precond.numLocalCols() );
  m_residualNorms.reserve( m_params.krylov.maxIterations + 1 );
}

template< typename VECTOR >
std::unique_ptr< KrylovSolver< VECTOR > >
KrylovSolver< VECTOR >::create( LinearSolverParameters const & parameters,
                                LinearOperator< VECTOR > const & matrix,
                                LinearOperator< VECTOR > const & precond )
{
  switch( parameters.solverType )
  {
    case LinearSolverParameters::SolverType::cg:
    {
      return std::make_unique< CgSolver< Vector > >( parameters,
                                                     matrix,
                                                     precond );
    }
    case LinearSolverParameters::SolverType::bicgstab:
    {
      return std::make_unique< BicgstabSolver< Vector > >( parameters,
                                                           matrix,
                                                           precond );
    }
    case LinearSolverParameters::SolverType::gmres:
    {
      return std::make_unique< GmresSolver< Vector > >( parameters,
                                                        matrix,
                                                        precond );
    }
    default:
    {
      GEOS_ERROR( "Unsupported linear solver type: " << parameters.solverType );
    }
  }
  return {};
}

template< typename VECTOR >
void KrylovSolver< VECTOR >::logProgress() const
{
  GEOS_ASSERT( !m_residualNorms.empty() );
  if( m_params.logLevel >= 2 )
  {
    real64 const relNorm = m_residualNorms[0] > 0.0 ? m_residualNorms.back() / m_residualNorms[0] : 0.0;
    GEOS_LOG_RANK_0( GEOS_FMT( "[{}] iteration {}: residual = {:e}", methodName(), m_result.numIterations, relNorm ) );
  }
}

template< typename VECTOR >
void KrylovSolver< VECTOR >::logResult() const
{
  if( m_params.logLevel >= 1 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "[{}] {} in {} iterations ({:.3f} s)", methodName(),
                               m_result.success() ? "converged" : "failed to converge",
                               m_result.numIterations, m_result.solveTime ) );
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class KrylovSolver< TrilinosInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_HYPRE
template class KrylovSolver< HypreInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_PETSC
template class KrylovSolver< PetscInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

}
