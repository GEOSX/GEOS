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
 * @file SuperLUDist.cpp
 */

#include "SuperLUDist.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/MpiWrapper.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/common/common.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/Arnoldi.hpp"
#include "linearAlgebra/utilities/NormalOperator.hpp"
#include "linearAlgebra/utilities/InverseNormalOperator.hpp"

#include <superlu_ddefs.h>

namespace geos
{

// Check matching requirements on index/value types between GEOSX and SuperLU_Dist

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CPU
static_assert( sizeof( int_t ) == sizeof( globalIndex ),
               "SuperLU_Dist int_t and geos::globalIndex must have the same size" );

static_assert( std::is_signed< int_t >::value == std::is_signed< globalIndex >::value,
               "SuperLU_Dist int_t and geos::globalIndex must both be signed or unsigned" );
#endif

static_assert( std::is_same< double, real64 >::value,
               "SuperLU_Dist real and geos::real64 must be the same type" );

namespace
{

/**
 * @brief Converts from GEOSX to SuperLU_Dist columns permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
colperm_t const & getColPermType( LinearSolverParameters::Direct::ColPerm const & value )
{
  static std::map< LinearSolverParameters::Direct::ColPerm, colperm_t > const optionMap =
  {
    { LinearSolverParameters::Direct::ColPerm::none, NATURAL },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtplusA, MMD_AT_PLUS_A },
    { LinearSolverParameters::Direct::ColPerm::MMD_AtA, MMD_ATA },
    { LinearSolverParameters::Direct::ColPerm::colAMD, COLAMD },
    { LinearSolverParameters::Direct::ColPerm::metis, METIS_AT_PLUS_A },
    { LinearSolverParameters::Direct::ColPerm::parmetis, PARMETIS },
  };

  GEOS_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist columns permutation option: " << value );
  return optionMap.at( value );
}

/**
 * @brief Converts from GEOSX to SuperLU_Dist rows permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
rowperm_t const & getRowPermType( LinearSolverParameters::Direct::RowPerm const & value )
{
  static std::map< LinearSolverParameters::Direct::RowPerm, rowperm_t > const optionMap =
  {
    { LinearSolverParameters::Direct::RowPerm::none, NOROWPERM },
    { LinearSolverParameters::Direct::RowPerm::mc64, LargeDiag_MC64 },
  };

  GEOS_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported SuperLU_Dist rows permutation option: " << value );
  return optionMap.at( value );
}

/**
 * @brief Attempts to arrange N entities (e.g. procs) into a grid as square as possible.
 * @tparam T index type (should be an integral type)
 * @param N total number of entities
 * @return a pair (x,y) of grid dimensions
 */
template< typename T >
std::pair< T, T > makeNearSquareGrid( T const N )
{
  static_assert( std::is_integral< T >::value, "T should be an integral type" );
  T x = static_cast< T >( std::sqrt( N ) );
  while( N % x )
  {
    --x;
  }
  T const y = N / x;
  return std::minmax( x, y );
}

} // namespace

struct SuperLUDistData
{
  array1d< int_t > rowPtr{};          ///< row pointers
  array1d< int_t > colIndices{};      ///< column indices
  array1d< double > values{};         ///< values
  array1d< double > rhs{};            ///< rhs/solution vector values
  SuperMatrix mat{};                  ///< SuperLU_Dist matrix format
  dScalePermstruct_t scalePerm{};     ///< data structure to scale and permute the matrix
  dLUstruct_t lu{};                   ///< data structure to store the LU factorization
  SuperLUStat_t stat{};               ///< data structure to gather some statistics
  gridinfo_t grid{};                  ///< SuperLU_Dist MPI subdivision of load
  dSOLVEstruct_t solve{};             ///< data structure to solve the matrix
  superlu_dist_options_t options{};   ///< SuperLU_Dist options

  SuperLUDistData( int_t const numGlobalRows,
                   int_t const numLocalRows,
                   int_t const numLocalNonzeros,
                   MPI_Comm const & comm )
  {
    rowPtr.resize( numLocalRows + 1 );
    colIndices.resize( numLocalNonzeros );
    values.resize( numLocalNonzeros );
    rhs.resize( numLocalRows );
    dScalePermstructInit( numGlobalRows, numGlobalRows, &scalePerm );
    dLUstructInit( numGlobalRows, &lu );
    PStatInit( &stat );

    // Create process grid: the goal is to have the process grid as square as possible
    std::pair< int, int > const gridsize = makeNearSquareGrid( MpiWrapper::commSize( comm ) );
    superlu_gridinit( comm, gridsize.first, gridsize.second, &grid );
  }

  SuperLUDistData( SuperLUDistData const & ) = delete;
  SuperLUDistData( SuperLUDistData && ) = delete;

  ~SuperLUDistData()
  {
    // Only cleanup memory taken by the NRformat_loc struct
    // The actual CSR memory is managed by array1d<> objects
    SUPERLU_FREE( (NRformat_loc *)mat.Store );

    dScalePermstructFree( &scalePerm );
    dDestroy_LU( mat.nrow, &grid, &lu );
    dLUstructFree( &lu );
    PStatFree( &stat );
    superlu_gridexit( &grid );
    if( options.SolveInitialized )
    {
      dSolveFinalize( &options, &solve );
    }
  }
};

template< typename LAI >
SuperLUDist< LAI >::SuperLUDist( LinearSolverParameters params )
  : Base( std::move( params ) ),
  m_condEst( -1.0 )
{}

template< typename LAI >
SuperLUDist< LAI >::~SuperLUDist() = default;

template< typename LAI >
void SuperLUDist< LAI >::setup( Matrix const & mat )
{
  clear();
  Base::setup( mat );

  int_t const numGR = LvArray::integerConversion< int_t >( mat.numGlobalRows() );
  int_t const numLR = LvArray::integerConversion< int_t >( mat.numLocalRows() );
  int_t const numNZ = LvArray::integerConversion< int_t >( mat.numLocalNonzeros() );

  m_data = std::make_unique< SuperLUDistData >( numGR, numLR, numNZ, mat.comm() );
  setOptions();

  typename Matrix::Export matExport;
  matExport.exportCRS( mat, m_data->rowPtr, m_data->colIndices, m_data->values );
  m_data->rowPtr.move( hostMemorySpace, false );
  m_data->colIndices.move( hostMemorySpace, false );
  m_data->values.move( hostMemorySpace, false );

  dCreate_CompRowLoc_Matrix_dist( &m_data->mat,
                                  numGR,
                                  numGR,
                                  numNZ,
                                  numLR,
                                  LvArray::integerConversion< int_t >( mat.ilower() ),
                                  m_data->values.data(),
                                  m_data->colIndices.data(),
                                  m_data->rowPtr.data(),
                                  SLU_NR_loc,
                                  SLU_D,
                                  SLU_GE );

  {
    Stopwatch timer( m_result.setupTime );
    factorize();
  }
}

template< typename LAI >
void SuperLUDist< LAI >::apply( Vector const & src,
                                Vector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( src.localSize(), dst.localSize() );
  GEOS_LAI_ASSERT_EQ( src.localSize(), matrix().numLocalRows() );

  // To be able to use SuperLU_Dist solver we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  // Export the rhs to a host-based array (this is required when vector is on GPU)
  typename Matrix::Export vecExport;
  vecExport.exportVector( src, m_data->rhs );
  m_data->rhs.move( hostMemorySpace, true );

  // Call the linear equation solver to solve the matrix.
  real64 berr = 0.0;
  int info = 0;

  m_data->options.Fact = FACTORED;
  pdgssvx( &m_data->options,
           &m_data->mat,
           &m_data->scalePerm,
           m_data->rhs.data(),
           m_data->rhs.size(),
           1,
           &m_data->grid,
           &m_data->lu,
           &m_data->solve,
           &berr,
           &m_data->stat,
           &info );

  GEOS_LAI_ASSERT_EQ( info, 0 );
  GEOS_LAI_ASSERT( !std::isnan( berr ) );
  GEOS_LAI_ASSERT( !std::isinf( berr ) );

  // Import the solution back into the vector
  vecExport.importVector( m_data->rhs, dst );
}

template< typename LAI >
void SuperLUDist< LAI >::solve( Vector const & rhs,
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
          GEOS_WARNING( "SuperLUDist: failed to reduce residual below tolerance.\n"
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
void SuperLUDist< LAI >::clear()
{
  Base::clear();
  m_data.reset();
  m_condEst = -1.0;
}

template< typename LAI >
void SuperLUDist< LAI >::setOptions()
{
  // Initialize options.
  set_default_options_dist( &m_data->options );
  m_data->options.PrintStat = m_params.logLevel > 1 ? YES : NO;
  m_data->options.Equil = m_params.direct.equilibrate ? YES : NO;
  m_data->options.ColPerm = getColPermType( m_params.direct.colPerm );
  m_data->options.RowPerm = getRowPermType( m_params.direct.rowPerm );
  m_data->options.ParSymbFact = m_params.direct.colPerm == LinearSolverParameters::Direct::ColPerm::parmetis ? YES : NO;
  m_data->options.ReplaceTinyPivot = m_params.direct.replaceTinyPivot ? YES : NO;
  m_data->options.IterRefine = m_params.direct.iterativeRefine ? SLU_DOUBLE : NOREFINE;

  if( m_params.logLevel > 0 )
  {
    print_sp_ienv_dist( &m_data->options );
    print_options_dist( &m_data->options );
  }
}

template< typename LAI >
void SuperLUDist< LAI >::factorize()
{
  // To be able to use SuperLU_Dist solver we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  // Call the linear equation solver to factorize the matrix.
  int info = 0;
  m_data->options.Fact = DOFACT;
  pdgssvx( &m_data->options,
           &m_data->mat,
           &m_data->scalePerm,
           nullptr,
           matrix().numLocalRows(),
           0,
           &m_data->grid,
           &m_data->lu,
           &m_data->solve,
           nullptr,
           &m_data->stat,
           &info );
  m_data->options.Fact = FACTORED;

  if( m_data->options.PrintStat == YES )
  {
    // Print the statistics.
    PStatPrint( &m_data->options, &m_data->stat, &m_data->grid );
  }

  GEOS_LAI_ASSERT_EQ( info, 0 );
}

template< typename LAI >
real64 SuperLUDist< LAI >::estimateConditionNumberBasic() const
{
  GEOS_LAI_ASSERT( ready() );
  if( m_condEst >= 0 )
  {
    return m_condEst; // used cached result, possibly more accurate
  }

  array1d< real64 > diagU( m_data->mat.nrow );
  pdGetDiagU( m_data->mat.nrow,
              &m_data->lu,
              &m_data->grid,
              diagU.data() );

  if( m_data->options.Equil == YES )
  {
    real64 const * const R = m_data->scalePerm.R;
    real64 const * const C = m_data->scalePerm.C;
    int_t const * const perm_c = m_data->scalePerm.perm_c;

    real64 minU = std::numeric_limits< real64 >::lowest();
    real64 maxU = std::numeric_limits< real64 >::max();
    real64 minL = std::numeric_limits< real64 >::lowest();
    real64 maxL = std::numeric_limits< real64 >::max();

    for( int_t i = 0; i < m_data->mat.nrow; ++i )
    {
      real64 const u = std::abs( diagU[perm_c[i]] / C[i] );
      minU = std::min( u, minU );
      maxU = std::max( u, maxU );
      real64 const l = std::abs( 1.0 / R[i] );
      minL = std::min( l, minL );
      maxL = std::max( l, maxL );
    }
    m_condEst = ( maxU / minU ) * ( maxL / minL );
  }
  else
  {
    real64 minU = std::numeric_limits< real64 >::lowest();
    real64 maxU = std::numeric_limits< real64 >::max();
    for( int_t i = 0; i < m_data->mat.nrow; ++i )
    {
      minU = std::min( std::abs( diagU[i] ), minU );
      maxU = std::max( std::abs( diagU[i] ), maxU );
    }
    m_condEst = maxU / minU;
  }
  return m_condEst;
}

template< typename LAI >
real64 SuperLUDist< LAI >::estimateConditionNumberAdvanced() const
{
  GEOS_LAI_ASSERT( ready() );
  localIndex constexpr numIterations = 4;

  NormalOperator< LAI > const normalOperator( matrix() );
  real64 const lambdaDirect = ArnoldiLargestEigenvalue( normalOperator, numIterations );

  InverseNormalOperator< LAI, SuperLUDist > const inverseNormalOperator( matrix(), *this );
  real64 const lambdaInverse = ArnoldiLargestEigenvalue( inverseNormalOperator, numIterations );

  m_condEst = sqrt( lambdaDirect * lambdaInverse );
  return m_condEst;
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class SuperLUDist< TrilinosInterface >;
#endif

#ifdef GEOS_USE_HYPRE
template class SuperLUDist< HypreInterface >;
#endif

#ifdef GEOS_USE_PETSC
template class SuperLUDist< PetscInterface >;
#endif

} // namespace geos
