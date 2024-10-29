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
 * @file HypreExport.cpp
 */

#include "HypreExport.hpp"

#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "LvArray/src/sortedArrayManipulation.hpp"

#include <_hypre_parcsr_mv.h>
#include <_hypre_IJ_mv.h>

namespace geos
{

HypreExport::HypreExport() = default;

HypreExport::HypreExport( HypreMatrix const & mat,
                          integer const targetRank )
  : m_targetRank( targetRank )
{
  // make a sub-communicator for scatter and ensure target rank is mapped to 0 in new comm
  int const rank = MpiWrapper::commRank( mat.comm() );
  int const color = ( mat.numLocalRows() > 0 ) ? 0 : MPI_UNDEFINED;
  int const key = ( rank == m_targetRank ) ? 0 : ( rank < m_targetRank ) ? rank + 1 : rank;
  m_subComm = MpiWrapper::commSplit( mat.comm(), color, key );
}

HypreExport::~HypreExport()
{
  if( m_subComm != MPI_COMM_NULL )
  {
    MpiWrapper::commFree( m_subComm );
  }
}

namespace
{

template< typename T >
void exportArray( HYPRE_MemoryLocation const location,
                  T const * const src,
                  arrayView1d< T > const & dst )
{
  // hypre does not maintain const-correctness of its APIs, hence the const_cast<>
  dst.move( hypre::getLvArrayMemorySpace( location ), true );
  hypre_TMemcpy( dst.data(), const_cast< T * >( src ), T, dst.size(), location, location );
}

template< typename T >
void exportArray( HYPRE_MemoryLocation const location,
                  arrayView1d< T const > const & src,
                  T * const dst )
{
  // hypre does not maintain const-correctness of its APIs, hence the const_cast<>
  src.move( hypre::getLvArrayMemorySpace( location ), false );
  hypre_TMemcpy( dst, const_cast< T * >( src.data() ), T, src.size(), location, location );
}

template< typename T, typename U >
void exportArray( HYPRE_MemoryLocation const location,
                  T const * const src,
                  arrayView1d< U > const & dst )
{
  if( location == HYPRE_MEMORY_HOST )
  {
    dst.move( hostMemorySpace, true );
    std::transform( src, src + dst.size(), dst.begin(),
                    []( T const v ) { return static_cast< U >( v ); } );
  }
  else // src is on device
  {
    forAll< hypre::execPolicy >( dst.size(), [dst, src] GEOS_HYPRE_DEVICE ( localIndex const i )
    {
      dst[i] = static_cast< U >( src[i] );
    } );
  }
}

template< typename T, typename U >
void exportArray( HYPRE_MemoryLocation const location,
                  arrayView1d< T const > const & src,
                  U * const dst )
{
  if( location == HYPRE_MEMORY_HOST )
  {
    src.move( hostMemorySpace, false );
    std::transform( src.begin(), src.end(), dst,
                    []( T const v ) { return static_cast< U >( v ); } );
  }
  else // dst is on device
  {
    forAll< hypre::execPolicy >( src.size(), [dst, src] GEOS_HYPRE_DEVICE ( localIndex const i )
    {
      dst[i] = static_cast< U >( src[i] );
    } );
  }
}

} // namespace

template< typename OFFSET_TYPE, typename COLUMN_TYPE >
void HypreExport::exportCRS( HypreMatrix const & mat,
                             arrayView1d< OFFSET_TYPE > const & rowOffsets,
                             arrayView1d< COLUMN_TYPE > const & colIndices,
                             arrayView1d< real64 > const & values ) const
{
  int const rank = MpiWrapper::commRank( mat.comm() );

  // import on target rank if needed, or extract diag+offd part in each rank
  hypre_CSRMatrix * const localMatrix = m_targetRank < 0
                                      ? hypre_MergeDiagAndOffd( mat.unwrapped() )
                                      : hypre_ParCSRMatrixToCSRMatrixAll( mat.unwrapped() );
  GEOS_ERROR_IF( rank == m_targetRank && !localMatrix, "HypreExport: matrix is empty on target rank" );

  if( m_targetRank < 0 || m_targetRank == rank )
  {
    HYPRE_Int const numRow = hypre_CSRMatrixNumRows( localMatrix );
    HYPRE_Int const numNz  = hypre_CSRMatrixNumNonzeros( localMatrix );

    GEOS_LAI_ASSERT_EQ( rowOffsets.size(), numRow + 1 );
    GEOS_LAI_ASSERT_EQ( colIndices.size(), numNz );
    GEOS_LAI_ASSERT_EQ( values.size(), numNz );

    HYPRE_MemoryLocation const location = hypre_CSRMatrixMemoryLocation( localMatrix );

    exportArray( location, hypre_CSRMatrixI( localMatrix ), rowOffsets );
    exportArray( location, hypre_CSRMatrixData( localMatrix ), values );

    // We have to handle two cases differently because hypre uses two different struct members
    // (j/big_j) to store the column indices depending on how we obtained the local matrix.
    if( m_targetRank < 0 )
    {
      exportArray( location, hypre_CSRMatrixBigJ( localMatrix ), colIndices );
    }
    else
    {
      exportArray( location, hypre_CSRMatrixJ( localMatrix ), colIndices );
    }

    // Sort the values by column index after copying (some solvers expect this)
    forAll< hypre::execPolicy >( numRow, [rowOffsets, colIndices, values] GEOS_HYPRE_DEVICE ( HYPRE_Int const i )
    {
      using LvArray::sortedArrayManipulation::dualSort;
      dualSort( colIndices.data() + rowOffsets[i],
                colIndices.data() + rowOffsets[i + 1],
                values.data() + rowOffsets[i] );
    } );
  }

  GEOS_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( localMatrix ) );
}

void HypreExport::exportVector( HypreVector const & vec,
                                arrayView1d< real64 > const & values ) const
{
  int const rank = MpiWrapper::commRank( vec.comm() );

  // Gather vector on target rank, or just get the local part
  hypre_Vector * const localVector = m_targetRank < 0
                                   ? hypre_ParVectorLocalVector( vec.unwrapped() )
                                   : (hypre_Vector *)hypre::parVectorToVectorAll( vec.unwrapped() );
  GEOS_ERROR_IF( rank == m_targetRank && !localVector, "HypreExport: vector is empty on target rank" );

  if( m_targetRank < 0 || m_targetRank == rank )
  {
    GEOS_LAI_ASSERT_EQ( values.size(), hypre_VectorSize( localVector ) );
    exportArray( hypre_VectorMemoryLocation( localVector ),
                 hypre_VectorData( localVector ),
                 values );
  }

  if( m_targetRank >= 0 )
  {
    GEOS_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( localVector ) );
  }
}

void HypreExport::importVector( arrayView1d< real64 const > const & values,
                                HypreVector & vec ) const
{
  hypre_Vector * const localVector = hypre_ParVectorLocalVector( vec.unwrapped() );
  if( m_subComm != MPI_COMM_NULL )
  {
    hypre_Vector * wrapperVector{};
    if( MpiWrapper::commRank( vec.comm() ) == m_targetRank )
    {
      GEOS_LAI_ASSERT_EQ( values.size(), vec.globalSize() );
      values.move( hostMemorySpace, false );

      // HACK: create a hypre vector that points to local data; we have to use const_cast,
      //       but this is ok because we don't modify the values, only scatter the vector.
      wrapperVector = hypre_SeqVectorCreate( LvArray::integerConversion< HYPRE_Int >( values.size() ) );
      hypre_VectorOwnsData( wrapperVector ) = false;
      hypre_VectorData( wrapperVector ) = const_cast< real64 * >( values.data() );
      hypre_SeqVectorInitialize_v2( wrapperVector, HYPRE_MEMORY_HOST );
    }

    // scatter the data
    hypre_ParVector * const parVector = hypre_VectorToParVector( m_subComm,
                                                                 wrapperVector,
                                                                 hypre_ParVectorPartitioning( vec.unwrapped() ) );
    // copy local part of the data over to the output vector
    hypre_TMemcpy( hypre_VectorData( localVector ),
                   hypre_VectorData( hypre_ParVectorLocalVector( parVector ) ),
                   HYPRE_Real,
                   vec.localSize(),
                   hypre_VectorMemoryLocation( localVector ),
                   hypre_ParVectorMemoryLocation( parVector ) );

    GEOS_LAI_CHECK_ERROR( hypre_ParVectorDestroy( parVector ) );
    GEOS_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( wrapperVector ) );
  }
  else
  {
    exportArray( hypre_VectorMemoryLocation( localVector ),
                 values,
                 hypre_VectorData( localVector ) );
  }
  vec.touch();
}

/**
 * Explicit template instantiation macro for HypreExport::exportCRS.
 * We need to explicitly instantiate this function template because:
 * - we want CRS consumers to specify their own destination buffer types;
 * - we're "hiding" Hypre headers from being included by consumer code.
 */
#define INST_HYPRE_EXPORT_CRS( OFFSET_TYPE, COLUMN_TYPE ) \
  template void \
  HypreExport::exportCRS< OFFSET_TYPE, COLUMN_TYPE >( HypreMatrix const &, \
                                                      arrayView1d< OFFSET_TYPE > const &, \
                                                      arrayView1d< COLUMN_TYPE > const &, \
                                                      arrayView1d< real64 > const & ) const

// Add other instantiations as needed (only use built-in types)
INST_HYPRE_EXPORT_CRS( int, int );
INST_HYPRE_EXPORT_CRS( long, long );
INST_HYPRE_EXPORT_CRS( long long, long long );

#undef INST_HYPRE_EXPORT_CRS

} // namespace geos
