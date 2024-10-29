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
 * @file HypreKernels.cpp
 */

#include "HypreKernels.hpp"

namespace geos
{
namespace hypre
{

HYPRE_BigInt const * getOffdColumnMap( hypre_ParCSRMatrix const * const mat )
{
  return ( hypre_CSRMatrixMemoryLocation( hypre_ParCSRMatrixDiag( mat ) ) == HYPRE_MEMORY_HOST )
         ? hypre_ParCSRMatrixColMapOffd( mat )
         : hypre_ParCSRMatrixDeviceColMapOffd( mat );
}

void scaleMatrixValues( hypre_CSRMatrix * const mat,
                        real64 const factor )
{
  GEOS_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return;
  }
  if( isEqual( factor, 1.0 ) )
  {
    return;
  }
  HYPRE_Real * const va = hypre_CSRMatrixData( mat );
  forAll< execPolicy >( hypre_CSRMatrixNumNonzeros( mat ), [=] GEOS_HYPRE_DEVICE ( HYPRE_Int const i )
  {
    va[i] *= factor;
  } );
}

void scaleMatrixRows( hypre_CSRMatrix * const mat,
                      hypre_Vector const * const vec )
{
  GEOS_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return;
  }

  CSRData< false > const csr{ mat };
  HYPRE_Real const * const scalingFactors = hypre_VectorData( vec );

  forAll< execPolicy >( csr.nrow, [=] GEOS_HYPRE_DEVICE ( HYPRE_Int const localRow )
  {
    real64 const factor = scalingFactors[localRow];
    if( !isEqual( factor, 1.0 ) )
    {
      for( HYPRE_Int j = csr.rowptr[localRow]; j < csr.rowptr[localRow + 1]; ++j )
      {
        csr.values[j] *= factor;
      }
    }
  } );
}

void clampMatrixEntries( hypre_CSRMatrix * const mat,
                         real64 const lo,
                         real64 const hi,
                         bool const skip_diag )
{
  GEOS_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return;
  }

  CSRData< false > const csr{ mat };
  forAll< execPolicy >( csr.nrow, [=] GEOS_HYPRE_DEVICE ( HYPRE_Int const localRow )
  {
    // Hypre stores diagonal element at the beginning of each row, we assume it's always present
    for( HYPRE_Int k = csr.rowptr[localRow] + skip_diag; k < csr.rowptr[localRow+1]; ++k )
    {
      csr.values[k] = LvArray::math::min( hi, LvArray::math::max( lo, csr.values[k] ) );
    }
  } );
}

real64 computeMaxNorm( hypre_CSRMatrix const * const mat )
{
  GEOS_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return 0;
  }

  HYPRE_Real const * const va = hypre_CSRMatrixData( mat );
  RAJA::ReduceMax< ReducePolicy< execPolicy >, real64 > maxAbsElement( 0.0 );
  forAll< execPolicy >( hypre_CSRMatrixNumNonzeros( mat ), [=] GEOS_HYPRE_DEVICE ( HYPRE_Int const k )
  {
    maxAbsElement.max( LvArray::math::abs( va[k] ) );
  } );
  return maxAbsElement.get();
}

real64 computeMaxNorm( hypre_CSRMatrix const * const mat,
                       arrayView1d< globalIndex const > const & rowIndices,
                       globalIndex const firstLocalRow )
{
  GEOS_ASSERT( mat != nullptr );
  if( hypre_CSRMatrixNumCols( mat ) == 0 )
  {
    return 0;
  }

  CSRData< true > const csr{ mat };
  HYPRE_Int const numRows = hypre_CSRMatrixNumRows( mat );
  GEOS_DEBUG_VAR( numRows );

  RAJA::ReduceMax< ReducePolicy< execPolicy >, real64 > maxAbsElement( 0.0 );
  forAll< execPolicy >( rowIndices.size(), [=] GEOS_HYPRE_DEVICE ( localIndex const i )
  {
    localIndex const localRow = rowIndices[i] - firstLocalRow;
    GEOS_ASSERT( 0 <= localRow && localRow < numRows );
    for( HYPRE_Int j = csr.rowptr[localRow]; j < csr.rowptr[localRow + 1]; ++j )
    {
      maxAbsElement.max( LvArray::math::abs( csr.values[j] ) );
    }
  } );
  return maxAbsElement.get();
}

} // namespace hypre
} // namespace geos
