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
 * @file SuiteSparse.cpp
 */

#include "SuiteSparse.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/common/common.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/Arnoldi.hpp"
#include "linearAlgebra/utilities/NormalOperator.hpp"
#include "linearAlgebra/utilities/InverseNormalOperator.hpp"
#include "common/MpiWrapper.hpp"

// Pre-define some suitesparse variables since they are not properly defined
// in the header for alternate index types.
//#if GEOS_GLOBALINDEX_TYPE_FLAG==0
#if 0
/// Set alias for SuiteSparse_long
#define SuiteSparse_long int

/// Set value for SuiteSparse_long_max
#define SuiteSparse_long_max 2147483647

/// Set string key for SuiteSparse_long
#define SuiteSparse_long_idd "d"

/// Set printf key for SuiteSparse_long
#define SuiteSparse_long_id "%d"
#endif

#if !defined(DLONG)
#define DLONG
#endif

#include <umfpack.h>

namespace geos
{

/**
 * SuiteSparse integer definition
 */
using SSlong = SuiteSparse_long;

static_assert( std::is_same< double, real64 >::value,
               "SuiteSparse real and geos::real64 must be the same type" );

struct SuiteSparseData
{
  array1d< SSlong > rowPtr{};          /// row pointers
  array1d< SSlong > colIndices{};      /// column indices
  array1d< double > values{};          /// values
  array1d< double > rhs{};             /// right-hand side vector
  array1d< double > sol{};             /// solution vector
  real64 info[UMFPACK_INFO]{};       /// data structure to gather various info
  real64 control[UMFPACK_CONTROL]{}; /// SuiteSparse options
  void * symbolic{};                 /// pointer to the symbolic factorization
  void * numeric{};                  /// pointer to the numeric factorization

  SuiteSparseData( SSlong const numRows,
                   SSlong const numNonzeros )
    : rowPtr( numRows + 1 ),
    colIndices( numNonzeros ),
    values( numNonzeros ),
    rhs( numRows ),
    sol( numRows )
  {}

  ~SuiteSparseData()
  {
    if( symbolic )
    {
      umfpack_dl_free_symbolic( &symbolic );
    }
    if( numeric )
    {
      umfpack_dl_free_numeric( &numeric );
    }
  }
};

template< typename LAI >
SuiteSparse< LAI >::SuiteSparse( LinearSolverParameters params )
  : Base( std::move( params ) ),
  m_workingRank( -1 ),
  m_condEst( -1.0 )
{}

template< typename LAI >
SuiteSparse< LAI >::~SuiteSparse() = default;

namespace
{

void factorize( SuiteSparseData & data, LinearSolverParameters const & params )
{
  // To be able to use UMFPACK direct solver we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  SSlong status;
  SSlong const numRows = data.rowPtr.size() - 1;

  data.rowPtr.move( hostMemorySpace, false );
  data.colIndices.move( hostMemorySpace, false );
  data.values.move( hostMemorySpace, false );

  // symbolic factorization
  status = umfpack_dl_symbolic( numRows,
                                numRows,
                                data.rowPtr.data(),
                                data.colIndices.data(),
                                data.values.data(),
                                &data.symbolic,
                                data.control,
                                data.info );
  if( status < 0 )
  {
    umfpack_dl_report_info( data.control, data.info );
    umfpack_dl_report_status( data.control, status );
    GEOS_ERROR( "SuiteSparse: umfpack_dl_symbolic failed." );
  }

  // print the symbolic factorization
  if( params.logLevel > 1 )
  {
    umfpack_dl_report_symbolic( data.symbolic, data.control );
  }

  // numeric factorization
  status = umfpack_dl_numeric( data.rowPtr.data(),
                               data.colIndices.data(),
                               data.values.data(),
                               data.symbolic,
                               &data.numeric,
                               data.control,
                               data.info );

  if( status < 0 )
  {
    umfpack_dl_report_info( data.control, data.info );
    umfpack_dl_report_status( data.control, status );
    GEOS_ERROR( "SuiteSparse: umfpack_dl_numeric failed." );
  }

  // print the numeric factorization
  if( params.logLevel > 1 )
  {
    umfpack_dl_report_numeric( data.symbolic, data.control );
  }
}

void setOptions( SuiteSparseData & data, LinearSolverParameters const & params )
{
  // Get the default control parameters
  umfpack_dl_defaults( data.control );
  data.control[UMFPACK_PRL] = params.logLevel > 1 ? 6 : 1;
  data.control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
}

} // namespace

template< typename LAI >
void SuiteSparse< LAI >::setup( Matrix const & mat )
{
  clear();
  PreconditionerBase< LAI >::setup( mat );

  // Choose working rank that will carry out the solve
  int const rank = MpiWrapper::commRank( mat.comm() );
  m_workingRank = MpiWrapper::min( mat.numLocalRows() > 0 ? rank : std::numeric_limits< int >::max(), mat.comm() );

  SSlong const numGR = LvArray::integerConversion< SSlong >( mat.numGlobalRows() );
  SSlong const numNZ = LvArray::integerConversion< SSlong >( mat.numGlobalNonzeros() );

  // Allocate memory and control structures on working rank only
  m_data = std::make_unique< SuiteSparseData >( rank == m_workingRank ? numGR : 0,
                                                rank == m_workingRank ? numNZ : 0 );
  setOptions( *m_data, m_params );

  // Export needs to be carried collectively on all ranks
  m_export = std::make_unique< typename Matrix::Export >( mat, m_workingRank );

  m_export->exportCRS( mat,
                       m_data->rowPtr,
                       m_data->colIndices,
                       m_data->values );

  // Perform matrix factorization on working rank
  if( rank == m_workingRank )
  {
    Stopwatch timer( m_result.setupTime );
    factorize( *m_data, m_params );
  }

  // Sync timer to all ranks
  MpiWrapper::bcast( &m_result.setupTime, 1, m_workingRank, mat.comm() );
}

template< typename LAI >
void SuiteSparse< LAI >::apply( Vector const & src,
                                Vector & dst ) const
{
  doSolve( src, dst, false );
}

template< typename LAI >
void SuiteSparse< LAI >::applyTranspose( Vector const & src,
                                         Vector & dst ) const
{
  doSolve( src, dst, true );
}

template< typename LAI >
void SuiteSparse< LAI >::clear()
{
  PreconditionerBase< LAI >::clear();
  m_data.reset();
  m_condEst = -1.0;
}

template< typename LAI >
void SuiteSparse< LAI >::solve( Vector const & rhs,
                                Vector & sol ) const
{
  real64 const bnorm = rhs.norm2();
  if( isZero( bnorm, 0.0 ) )
  {
    sol.zero();
    m_result.numIterations = 0;
    m_result.residualReduction = 0.0;
    m_result.solveTime = 0.0;
    m_result.status = LinearSolverResult::Status::Success;
    return;
  }

  {
    Stopwatch timer( m_result.solveTime );
    apply( rhs, sol );
  }
  MpiWrapper::bcast( &m_result.solveTime, 1, m_workingRank, rhs.comm() );

  Vector r( rhs );
  matrix().residual( sol, r, r );
  m_result.residualReduction = r.norm2() / bnorm;

  m_result.status = LinearSolverResult::Status::Success;
  m_result.numIterations = 1;

  if( m_params.direct.checkResidual )
  {
    real64 constexpr precTol = 100.0 * std::numeric_limits< real64 >::epsilon();
    real64 condEst = estimateConditionNumberBasic();
    if( m_result.residualReduction > condEst * precTol )
    {
      condEst = estimateConditionNumberAdvanced();
      if( m_result.residualReduction > condEst * precTol )
      {
        if( m_params.logLevel > 0 )
        {
          GEOS_WARNING( "SuiteSparse: failed to reduce residual below tolerance.\n"
                        "Condition number estimate: " << condEst );
        }
        m_result.status = LinearSolverResult::Status::Breakdown;
      }
    }
  }

  if( m_params.logLevel >= 1 )
  {
    GEOS_LOG_RANK_0( "        Linear Solver | " << m_result.status <<
                     " | Iterations: " << m_result.numIterations <<
                     " | Final Rel Res: " << m_result.residualReduction <<
                     " | Setup Time: " << m_result.setupTime << " s" <<
                     " | Solve Time: " << m_result.solveTime << " s" );
  }
}

template< typename LAI >
void SuiteSparse< LAI >::doSolve( Vector const & b, Vector & x, bool transpose ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( b.ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT_EQ( b.localSize(), x.localSize() );
  GEOS_LAI_ASSERT_EQ( b.localSize(), matrix().numLocalRows() );

  m_export->exportVector( b, m_data->rhs );

  if( MpiWrapper::commRank( b.comm() ) == m_workingRank )
  {
    m_data->rhs.move( hostMemorySpace, false );
    m_data->sol.move( hostMemorySpace, true );

    // To be able to use UMFPACK direct solver we need to disable floating point exceptions
    LvArray::system::FloatingPointExceptionGuard guard;

    // Note: UMFPACK expects column-sparse matrix, but we have row-sparse, so we flip the transpose flag
    SSlong const status = umfpack_dl_solve( transpose ? UMFPACK_A : UMFPACK_At,
                                            m_data->rowPtr.data(),
                                            m_data->colIndices.data(),
                                            m_data->values.data(),
                                            m_data->sol.data(),
                                            m_data->rhs.data(),
                                            m_data->numeric,
                                            m_data->control,
                                            m_data->info );

    if( status < 0 )
    {
      umfpack_dl_report_info( m_data->control, m_data->info );
      umfpack_dl_report_status( m_data->control, status );
      GEOS_ERROR( "SuiteSparse interface: umfpack_dl_solve failed." );
    }
  }

  m_export->importVector( m_data->sol, x );
}

template< typename LAI >
real64 SuiteSparse< LAI >::estimateConditionNumberBasic() const
{
  GEOS_LAI_ASSERT( ready() );
  if( m_condEst >= 0.0 )
  {
    return m_condEst; // used cached result, possibly more accurate
  }

  int const rank = MpiWrapper::commRank( matrix().comm() );
  if( rank == m_workingRank )
  {
    m_condEst = 1.0 / m_data->info[UMFPACK_RCOND];
  }
  MpiWrapper::bcast( &m_condEst, 1, m_workingRank, matrix().comm() );

  return m_condEst;
}

template< typename LAI >
real64 SuiteSparse< LAI >::estimateConditionNumberAdvanced() const
{
  GEOS_LAI_ASSERT( ready() );
  localIndex constexpr numIterations = 4;

  NormalOperator< LAI > const normalOperator( matrix() );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue( normalOperator, numIterations );

  InverseNormalOperator< LAI, SuiteSparse > const inverseNormalOperator( matrix(), *this );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue( inverseNormalOperator, numIterations );

  m_condEst = sqrt( lambdaDirect * lambdaInverse );
  return m_condEst;
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class SuiteSparse< TrilinosInterface >;
#endif

#ifdef GEOS_USE_HYPRE
template class SuiteSparse< HypreInterface >;
#endif

#ifdef GEOS_USE_PETSC
template class SuiteSparse< PetscInterface >;
#endif

}
