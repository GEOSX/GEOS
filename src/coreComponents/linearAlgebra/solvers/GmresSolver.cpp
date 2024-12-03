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
 * @file GMRESsolver.cpp
 */

#include "GmresSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"
#include "denseLinearAlgebra/common/layouts.hpp"

namespace geos
{

template< typename VECTOR >
GmresSolver< VECTOR >::GmresSolver( LinearSolverParameters params,
                                    LinearOperator< Vector > const & A,
                                    LinearOperator< Vector > const & M )
  : KrylovSolver< VECTOR >( std::move( params ), A, M ),
  m_kspace( m_params.krylov.maxRestart + 1 ),
  m_kspaceInitialized( false )
{
  GEOS_ERROR_IF_LE_MSG( m_params.krylov.maxRestart, 0, "GMRES: max number of iterations until restart must be positive." );
}

namespace
{

void ComputeGivensRotation( real64 const x, real64 const y, real64 & c, real64 & s )
{
  if( isZero( y ) )
  {
    c = 1.0;
    s = 0.0;
  }
  else if( std::fabs( y ) > std::fabs( x ) )
  {
    real64 const nu = x / y;
    s = 1.0 / std::sqrt( 1.0 + nu * nu );
    c = nu * s;
  }
  else
  {
    real64 const nu = y / x;
    c = 1.0 / std::sqrt( 1.0 + nu * nu );
    s = nu * c;
  }
}

void ApplyGivensRotation( real64 const c, real64 const s, real64 & dx, real64 & dy )
{
  real64 const temp = c * dx + s * dy;
  dy = -s * dx + c * dy;
  dx = temp;
}

void Backsolve( integer const k,
                arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & H,
                arraySlice1d< real64 > const & g )
{
  for( integer j = k - 1; j >= 0; --j )
  {
    g[j] /= H( j, j );
    for( integer i = j - 1; i >= 0; --i )
    {
      g[i] -= H( i, j ) * g[j];
    }
  }

}

} // namespace

template< typename VECTOR >
void GmresSolver< VECTOR >::solve( Vector const & b,
                                   Vector & x ) const
{
  // We create Krylov subspace vectors once using the size and partitioning of b.
  // On repeated calls to solve() input vectors must have the same size and partitioning.
  if( !m_kspaceInitialized )
  {
    for( VectorTemp & kv : m_kspace )
    {
      kv = createTempVector( b );
    }
    m_kspaceInitialized = true;
  }

  Stopwatch watch;

  // Define vectors
  VectorTemp r = createTempVector( b );
  VectorTemp w = createTempVector( b );
  VectorTemp z = createTempVector( b );

  // Compute initial rk
  m_operator.residual( x, b, r );

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;

  // Create upper Hessenberg matrix
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > H( m_params.krylov.maxRestart + 1, m_params.krylov.maxRestart );

  // Create plane rotation storage
  array1d< real64 > c( m_params.krylov.maxRestart + 1 );
  array1d< real64 > s( m_params.krylov.maxRestart + 1 );
  array1d< real64 > g( m_params.krylov.maxRestart + 1 );

  // Initialize iteration state
  m_result.status = LinearSolverResult::Status::NotConverged;
  m_residualNorms.clear();

  integer & k = m_result.numIterations;
  while( k <= m_params.krylov.maxIterations && m_result.status == LinearSolverResult::Status::NotConverged )
  {
    // Re-initialize Krylov subspace
    g.zero();
    g[0] = k > 0 ? r.norm2() : rnorm0;
    m_kspace[0].copy( r );
    if( g[0] > 0 )
    {
      m_kspace[0].scale( 1.0 / g[0] );
    }

    integer j = 0;
    for(; j < m_params.krylov.maxRestart && k <= m_params.krylov.maxIterations; ++j, ++k )
    {
      // Record iteration progress
      real64 const rnorm = std::fabs( g[j] );
      m_residualNorms.emplace_back( rnorm );
      logProgress();

      // Convergence check
      if( rnorm <= absTol )
      {
        m_result.status = LinearSolverResult::Status::Success;
        break;
      }

      // Compute the new vector
      m_precond.apply( m_kspace[j], z );
      m_operator.apply( z, w );

      // Orthogonalization
      for( integer i = 0; i <= j; ++i )
      {
        H( i, j ) = w.dot( m_kspace[i] );
        w.axpby( -H( i, j ), m_kspace[i], 1.0 );
      }

      H( j+1, j ) = w.norm2();
      GEOS_KRYLOV_BREAKDOWN_IF_ZERO( H( j+1, j ) )
      m_kspace[j+1].axpby( 1.0 / H( j+1, j ), w, 0.0 );

      // Apply all previous rotations to the new column
      for( integer i = 0; i < j; ++i )
      {
        ApplyGivensRotation( c[i], s[i], H( i, j ), H( i+1, j ) );
      }

      // Compute and apply the new rotation to eliminate subdiagonal element
      ComputeGivensRotation( H( j, j ), H( j+1, j ), c[j], s[j] );
      ApplyGivensRotation( c[j], s[j], H( j, j ), H( j+1, j ) );
      ApplyGivensRotation( c[j], s[j], g[j], g[j+1] );
    }

    // Regardless of how we quit out of inner loop, j is the actual size of H
    Backsolve( j, H, g );
    w.zero();
    for( integer i = 0; i < j; ++i )
    {
      w.axpy( g[i], m_kspace[i] );
    }
    m_precond.apply( w, z );

    // Update the solution vector and recompute residual
    x.axpy( 1.0, z );
    m_operator.residual( x, b, r );
  }

  m_result.residualReduction = rnorm0 > 0.0 ? m_residualNorms.back() / rnorm0 : 0.0;
  m_result.solveTime = watch.elapsedTime();
  logResult();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class GmresSolver< TrilinosInterface::ParallelVector >;
template class GmresSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_HYPRE
template class GmresSolver< HypreInterface::ParallelVector >;
template class GmresSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOS_USE_PETSC
template class GmresSolver< PetscInterface::ParallelVector >;
template class GmresSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

} // namespace geos
