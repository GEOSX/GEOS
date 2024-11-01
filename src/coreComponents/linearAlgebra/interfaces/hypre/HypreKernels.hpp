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
 * @file HypreKernels.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREKERNELS_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREKERNELS_HPP_

#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/common/common.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#include <_hypre_parcsr_mv.h>

namespace geos
{
namespace hypre
{

/// @cond DO_NOT_DOCUMENT

namespace ops
{

template< typename T >
GEOS_HYPRE_HOST_DEVICE
constexpr T identity( T const v )
{
  return v;
}

template< typename T >
GEOS_HYPRE_HOST_DEVICE
constexpr T plus( T const lhs, T const rhs )
{
  return lhs + rhs;
}

}

template< bool CONST >
struct CSRData
{
  HYPRE_Int const * rowptr;
  HYPRE_Int const * colind;
  add_const_if_t< HYPRE_Real, CONST > * values;
  HYPRE_Int nrow;
  HYPRE_Int ncol;
  HYPRE_Int nnz;

  explicit CSRData( add_const_if_t< hypre_CSRMatrix, CONST > * const mat )
    : rowptr( hypre_CSRMatrixI( mat ) ),
    colind( hypre_CSRMatrixJ( mat ) ),
    values( hypre_CSRMatrixData( mat ) ),
    nrow( hypre_CSRMatrixNumRows( mat ) ),
    ncol( hypre_CSRMatrixNumCols( mat ) ),
    nnz( hypre_CSRMatrixNumNonzeros( mat ) )
  {}
};

HYPRE_BigInt const * getOffdColumnMap( hypre_ParCSRMatrix const * const mat );

void scaleMatrixValues( hypre_CSRMatrix * const mat,
                        real64 const factor );

void scaleMatrixRows( hypre_CSRMatrix * const mat,
                      hypre_Vector const * const vec );

void clampMatrixEntries( hypre_CSRMatrix * const mat,
                         real64 const lo,
                         real64 const hi,
                         bool const skip_diag );

real64 computeMaxNorm( hypre_CSRMatrix const * const mat );

real64 computeMaxNorm( hypre_CSRMatrix const * const mat,
                       arrayView1d< globalIndex const > const & rowIndices,
                       globalIndex const firstLocalRow );

namespace internal
{

/// This type is needed because CUDA does not allow a local type to be captured in a device lambda.
template< typename F, typename R >
struct RowReducer
{
  F transform;
  R reduce;

  auto GEOS_HYPRE_HOST_DEVICE
  operator()( double acc, double v ) const
  {
    return reduce( acc, transform( v ) );
  }
};

} // namespace internal

template< typename F, typename R >
void rescaleMatrixRows( hypre_ParCSRMatrix * const mat,
                        arrayView1d< globalIndex const > const & rowIndices,
                        F transform,
                        R reduce )
{
  CSRData< false > diag{ hypre_ParCSRMatrixDiag( mat ) };
  CSRData< false > offd{ hypre_ParCSRMatrixOffd( mat ) };
  HYPRE_BigInt const firstLocalRow = hypre_ParCSRMatrixFirstRowIndex( mat );
  internal::RowReducer< F, R > reducer{ std::move( transform ), std::move( reduce ) };

  forAll< execPolicy >( rowIndices.size(), [diag, offd, reducer, rowIndices, firstLocalRow] GEOS_HYPRE_HOST_DEVICE ( localIndex const i )
      {
        HYPRE_Int const localRow = LvArray::integerConversion< HYPRE_Int >( rowIndices[i] - firstLocalRow );
        GEOS_ASSERT( 0 <= localRow && localRow < diag.nrow );

        HYPRE_Real scale = 0.0;
        for( HYPRE_Int k = diag.rowptr[localRow]; k < diag.rowptr[localRow + 1]; ++k )
        {
          scale = reducer( scale, diag.values[k] );
        }
        if( offd.ncol > 0 )
        {
          for( HYPRE_Int k = offd.rowptr[localRow]; k < offd.rowptr[localRow + 1]; ++k )
          {
            scale = reducer( scale, offd.values[k] );
          }
        }

        GEOS_ASSERT_MSG( !isZero( scale ), "Zero row sum in row " << rowIndices[i] );
        scale = 1.0 / scale;
        for( HYPRE_Int k = diag.rowptr[localRow]; k < diag.rowptr[localRow + 1]; ++k )
        {
          diag.values[k] *= scale;
        }
        if( offd.ncol > 0 )
        {
          for( HYPRE_Int k = offd.rowptr[localRow]; k < offd.rowptr[localRow + 1]; ++k )
          {
            offd.values[k] *= scale;
          }
        }
      } );
}

template< typename F, typename R >
void computeRowsSums( hypre_ParCSRMatrix const * const mat,
                      hypre_ParVector * const vec,
                      F transform,
                      R reduce )
{
  CSRData< true > const diag{ hypre_ParCSRMatrixDiag( mat ) };
  CSRData< true > const offd{ hypre_ParCSRMatrixOffd( mat ) };
  HYPRE_Real * const values = hypre_VectorData( hypre_ParVectorLocalVector( vec ) );
  internal::RowReducer< F, R > reducer{ std::move( transform ), std::move( reduce ) };

  forAll< execPolicy >( diag.nrow, [diag, offd, reducer, values] GEOS_HYPRE_HOST_DEVICE ( HYPRE_Int const localRow )
      {
        HYPRE_Real sum = 0.0;
        for( HYPRE_Int k = diag.rowptr[localRow]; k < diag.rowptr[localRow + 1]; ++k )
        {
          sum = reducer( sum, diag.values[k] );
        }
        if( offd.ncol )
        {
          for( HYPRE_Int k = offd.rowptr[localRow]; k < offd.rowptr[localRow + 1]; ++k )
          {
            sum = reducer( sum, offd.values[k] );
          }
        }
        values[localRow] = sum;
      } );
}

namespace internal
{

template< typename MAP >
void GEOS_HYPRE_HOST_DEVICE
makeSortedPermutation( HYPRE_Int const * const indices,
                       HYPRE_Int const size,
                       HYPRE_Int * const perm,
                       MAP map )
{
  for( HYPRE_Int i = 0; i < size; ++i )
  {
    perm[i] = i; // std::iota
  }
  auto const comp = [indices, map] GEOS_HYPRE_HOST_DEVICE ( HYPRE_Int i, HYPRE_Int j )
  {
    return map( indices[i] ) < map( indices[j] );
  };
  LvArray::sortedArrayManipulation::makeSorted( perm, perm + size, comp );
}

} // namespace internal

template< typename SRC_COLMAP, typename DST_COLMAP >
void addEntriesRestricted( hypre_CSRMatrix const * const src_mat,
                           SRC_COLMAP const src_colmap,
                           hypre_CSRMatrix * const dst_mat,
                           DST_COLMAP const dst_colmap,
                           real64 const scale )
{
  GEOS_LAI_ASSERT( src_mat != nullptr );
  GEOS_LAI_ASSERT( dst_mat != nullptr );

  CSRData< true > src{ src_mat };
  CSRData< false > dst{ dst_mat };
  GEOS_LAI_ASSERT_EQ( src.nrow, dst.nrow );

  if( src.ncol == 0 || isZero( scale ) )
  {
    return;
  }

  // Allocate contiguous memory to store sorted column permutations of each row
  array1d< HYPRE_Int > const src_permutation_arr( hypre_CSRMatrixNumNonzeros( src_mat ) );
  array1d< HYPRE_Int > const dst_permutation_arr( hypre_CSRMatrixNumNonzeros( dst_mat ) );

  arrayView1d< HYPRE_Int > const src_permutation = src_permutation_arr.toView();
  arrayView1d< HYPRE_Int > const dst_permutation = dst_permutation_arr.toView();
  // Each thread adds one row of src into dst
  forAll< hypre::execPolicy >( dst.nrow,
                               [ src,
                                 src_colmap,
                                 dst,
                                 dst_colmap,
                                 scale,
                                 src_permutation,
                                 dst_permutation] GEOS_HYPRE_DEVICE ( HYPRE_Int const localRow )
  {
    HYPRE_Int const src_offset = src.rowptr[localRow];
    HYPRE_Int const src_length = src.rowptr[localRow + 1] - src_offset;
    HYPRE_Int const * const src_indices = src.colind + src_offset;
    HYPRE_Real const * const src_values = src.values + src_offset;
    HYPRE_Int * const src_perm = src_permutation.data() + src_offset;

    HYPRE_Int const dst_offset = dst.rowptr[localRow];
    HYPRE_Int const dst_length = dst.rowptr[localRow + 1] - dst_offset;
    HYPRE_Int const * const dst_indices = dst.colind + dst_offset;
    HYPRE_Real * const dst_values = dst.values + dst_offset;
    HYPRE_Int * const dst_perm = dst_permutation.data() + dst_offset;

    // Since hypre does not store columns in sorted order, create a sorted "view" of src and dst rows
    // TODO: it would be nice to cache the permutation arrays somewhere to avoid recomputing
    internal::makeSortedPermutation( src_indices, src_length, src_perm, src_colmap );
    internal::makeSortedPermutation( dst_indices, dst_length, dst_perm, dst_colmap );

    // Add entries looping through them in sorted column order, skipping src entries not in dst
    for( HYPRE_Int i = 0, j = 0; i < dst_length && j < src_length; ++i )
    {
      while( j < src_length && src_colmap( src_indices[src_perm[j]] ) < dst_colmap( dst_indices[dst_perm[i]] ) )
      {
        ++j;
      }
      if( j < src_length && src_colmap( src_indices[src_perm[j]] ) == dst_colmap( dst_indices[dst_perm[i]] ) )
      {
        dst_values[dst_perm[i]] += scale * src_values[src_perm[j++]];
      }
    }
  } );
}

/// @endcond

} // namespace hypre
} // namespace geos

#endif //GEOS_LINEARALGEBRA_INTERFACES_HYPREKERNELS_HPP_
