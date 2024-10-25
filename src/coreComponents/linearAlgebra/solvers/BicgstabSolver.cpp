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
 * @file BiCGSTABsolver.cpp
 *
 */

#include "BicgstabSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "common/LinearOperator.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"

namespace geos
{


template< typename VECTOR >
BicgstabSolver< VECTOR >::BicgstabSolver( LinearSolverParameters params,
                                          LinearOperator< Vector > const & A,
                                          LinearOperator< Vector > const & M )
  : KrylovSolver< VECTOR >( std::move( params ), A, M )
{}

template< typename VECTOR >
void BicgstabSolver< VECTOR >::solve( Vector const & b,
                                      Vector & x ) const
{
  Stopwatch watch;

  // Define vectors
  VectorTemp r( x );

  // Compute initial rk
  m_operator.residual( x, b, r );

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;

  // Define vectors
  VectorTemp r0( r );

  // Define scalars and reinitialize some
  real64 rho_old = rnorm0 * rnorm0; // same as r.dot( r0 );
  real64 alpha = 1.0;
  real64 omega = 1.0;

  // Define temporary vectors
  VectorTemp v = createTempVector( r );
  VectorTemp p = createTempVector( r );
  VectorTemp y = createTempVector( r );
  VectorTemp z = createTempVector( r );
  VectorTemp t = createTempVector( r );
  VectorTemp s = createTempVector( r );

  v.zero();
  p.zero();

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

    // Compute r0.rk
    real64 const rho = r.dot( r0 );

    GEOS_KRYLOV_BREAKDOWN_IF_ZERO( rho_old )
    GEOS_KRYLOV_BREAKDOWN_IF_ZERO( omega )

    // Compute beta
    real64 const beta = rho / rho_old * alpha / omega;

    // Update p = r + beta*(p - omega*v)
    p.axpy( -omega, v );
    p.axpby( 1., r, beta );

    // Update vk = MApk
    m_precond.apply( p, y );
    m_operator.apply( y, v );

    // Compute alpha
    real64 const vr0 = v.dot( r0 );
    GEOS_KRYLOV_BREAKDOWN_IF_ZERO( vr0 )
    alpha = rho / vr0;

    // compute x = x + alpha*y
    x.axpy( alpha, y );

    // Compute s = rk - alpha*vk
    s.copy( r );
    s.axpy( -alpha, v );

    // Compute z = Ms
    m_precond.apply( s, z );

    // Compute t = Az
    m_operator.apply( z, t );

    // Update omega
    real64 const t2 = t.dot( t );
    GEOS_KRYLOV_BREAKDOWN_IF_ZERO( t2 )
    omega = t.dot( s ) / t2;

    // Update x = x + omega*z
    x.axpy( omega, z );

    // Update rk = s - omega*t
    s.axpy( -omega, t );
    r.copy( s );

    // Keep the old value of rho
    rho_old = rho;
  }

  m_result.residualReduction = rnorm0 > 0.0 ? m_residualNorms.back() / rnorm0 : 0.0;
  m_result.solveTime = watch.elapsedTime();
  logResult();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class BicgstabSolver< TrilinosInterface::ParallelVector >;
template class BicgstabSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_HYPRE
template class BicgstabSolver< HypreInterface::ParallelVector >;
template class BicgstabSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_PETSC
template class BicgstabSolver< PetscInterface::ParallelVector >;
template class BicgstabSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

} //namespace geos
