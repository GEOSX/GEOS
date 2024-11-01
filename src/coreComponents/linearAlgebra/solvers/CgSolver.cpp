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
 * @file CGsolver.cpp
 *
 */


#include "CgSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "common/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"

namespace geos
{

// BEGIN_RST_NARRATIVE CGsolver.rst
// ==============================
// CG Solver
// ==============================
// Implementation of the Preconditioned Conjugate Gradient (CG) algorithm.
// The notation is consistent with "Iterative Methods for
// Linear and Non-Linear Equations" from C.T. Kelley (1995)
// and "Iterative Methods for Sparse Linear Systems" from Y. Saad (2003).


// ----------------------------
// Constructor
// ----------------------------
template< typename VECTOR >
CgSolver< VECTOR >::CgSolver( LinearSolverParameters params,
                              LinearOperator< Vector > const & A,
                              LinearOperator< Vector > const & M )
  : KrylovSolver< VECTOR >( std::move( params ), A, M )
{
  GEOS_ERROR_IF( !m_params.isSymmetric, "Cannot use CG solver with a non-symmetric system" );
}

// ----------------------------
// Solve method
// ----------------------------
template< typename VECTOR >
void CgSolver< VECTOR >::solve( Vector const & b, Vector & x ) const
{
  Stopwatch watch;

  // Define residual vector
  VectorTemp r = createTempVector( b );

  // Compute initial rk =  b - Ax
  m_operator.residual( x, b, r );

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;

  // Preconditioning
  VectorTemp z = createTempVector( x );

  // Search direction
  VectorTemp p = createTempVector( z );
  VectorTemp Ap = createTempVector( z );
  p.zero();

  // Keep old value of preconditioned residual norm
  real64 tau_old = 0.0;

  // Initialize iteration state
  m_result.status = LinearSolverResult::Status::NotConverged;
  m_residualNorms.clear();

  integer & k = m_result.numIterations;
  for( k = 0; k <= m_params.krylov.maxIterations; ++k )
  {
    real64 const rnorm = r.norm2();
    m_residualNorms.emplace_back( rnorm );
    logProgress();

    // Convergence check on ||rk||/||b||
    if( rnorm <= absTol )
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }

    // Update z = Mr
    m_precond.apply( r, z );

    // Compute beta
    real64 const tau = z.dot( r );
    real64 const beta = k > 0 ? tau / tau_old : 0.0;

    // Update p = z + beta*p
    p.axpby( 1.0, z, beta );

    // Compute Ap
    m_operator.apply( p, Ap );

    // compute alpha
    real64 const pAp = p.dot( Ap );
    GEOS_KRYLOV_BREAKDOWN_IF_ZERO( pAp )
    real64 const alpha = tau / pAp;

    // Update x = x + alpha*p
    x.axpby( alpha, p, 1.0 );

    // Update rk = rk - alpha*Ap
    r.axpby( -alpha, Ap, 1.0 );

    // Keep the old value of tau
    tau_old = tau;
  }

  m_result.residualReduction = rnorm0 > 0.0 ? m_residualNorms.back() / rnorm0 : 0.0;
  m_result.solveTime = watch.elapsedTime();
  logResult();
}

// END_RST_NARRATIVE

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class CgSolver< TrilinosInterface::ParallelVector >;
template class CgSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_HYPRE
template class CgSolver< HypreInterface::ParallelVector >;
template class CgSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_PETSC
template class CgSolver< PetscInterface::ParallelVector >;
template class CgSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

} //namespace geos
