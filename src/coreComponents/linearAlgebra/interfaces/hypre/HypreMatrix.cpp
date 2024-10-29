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
 * @file HypreMatrix.cpp
 */

#include "HypreMatrix.hpp"

#include "common/TimingMacros.hpp"
#include "common/GeosxConfig.hpp"
#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/hypre/HypreKernels.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "LvArray/src/output.hpp"

#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"

#include <iomanip>
#include <numeric>

namespace geos
{

// Helper function that performs the following sequence of IJMatrix
// call: Create, SetObjectType, Initialize.
static void initialize( MPI_Comm const & comm,
                        HYPRE_BigInt const & ilower,
                        HYPRE_BigInt const & iupper,
                        HYPRE_BigInt const & jlower,
                        HYPRE_BigInt const & jupper,
                        arrayView1d< HYPRE_Int const > const & ncols,
                        HYPRE_IJMatrix & ij_matrix )
{
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixCreate( comm,
                                              ilower,
                                              iupper,
                                              jlower,
                                              jupper,
                                              &ij_matrix ) );

  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetObjectType( ij_matrix, HYPRE_PARCSR ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetRowSizes( ij_matrix, ncols.data() ) );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( ij_matrix ) );
}

HypreMatrix::HypreMatrix()
  : LinearOperator(),
  MatrixBase()
{}

HypreMatrix::HypreMatrix( HypreMatrix const & src )
  : HypreMatrix()
{
  *this = src;
}

HypreMatrix::HypreMatrix( HypreMatrix && src ) noexcept
  : HypreMatrix()
{
  *this = std::move( src );
}

HypreMatrix & HypreMatrix::operator=( HypreMatrix const & src )
{
  if( &src != this )
  {
    reset();
    if( src.ready() )
    {
      // Copy parcsr matrix
      HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParCSRMatrixClone( src.m_parcsr_mat, 1 );
      // Create IJ layer (with matrix closed)
      parCSRtoIJ( dst_parcsr );
    }
    m_dofManager = src.dofManager();
  }
  return *this;
}

HypreMatrix & HypreMatrix::operator=( HypreMatrix && src ) noexcept
{
  if( &src != this )
  {
    std::swap( m_ij_mat, src.m_ij_mat );
    std::swap( m_parcsr_mat, src.m_parcsr_mat );
    MatrixBase::operator=( std::move( src ) );
  }
  return *this;
}

HypreMatrix::~HypreMatrix()
{
  reset();
}

void HypreMatrix::createWithGlobalSize( globalIndex const globalRows,
                                        globalIndex const globalCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  GEOS_LAI_ASSERT( closed() );
  GEOS_LAI_ASSERT_GE( globalRows, 0 );
  GEOS_LAI_ASSERT_GE( globalCols, 0 );
  GEOS_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  HYPRE_Int const rank  = LvArray::integerConversion< HYPRE_Int >( MpiWrapper::commRank( comm ) );
  HYPRE_Int const nproc = LvArray::integerConversion< HYPRE_Int >( MpiWrapper::commSize( comm ) );

  HYPRE_Int const localRowSize = LvArray::integerConversion< HYPRE_Int >( globalRows / nproc );
  HYPRE_Int const rowResidual = LvArray::integerConversion< HYPRE_Int >( globalRows % nproc );

  HYPRE_Int const localColSize = LvArray::integerConversion< HYPRE_Int >( globalCols / nproc );
  HYPRE_Int const colResidual = LvArray::integerConversion< HYPRE_Int >( globalCols % nproc );

  HYPRE_BigInt const ilower = rank * localRowSize + ( rank == 0 ? 0 : rowResidual );
  HYPRE_BigInt const iupper = ilower + localRowSize + ( rank == 0 ? rowResidual : 0 ) - 1;
  HYPRE_BigInt const jlower = rank * localColSize + ( rank == 0 ? 0 : colResidual );
  HYPRE_BigInt const jupper = jlower + localColSize + ( rank == 0 ? colResidual : 0 ) - 1;

  array1d< HYPRE_Int > row_sizes;
  row_sizes.resizeDefault( LvArray::integerConversion< localIndex >( iupper - ilower + 1 ),
                           LvArray::integerConversion< HYPRE_Int >( maxEntriesPerRow ) );

  initialize( comm,
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );
}

void HypreMatrix::create( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
                          localIndex const numLocalColumns,
                          MPI_Comm const & comm )
{
  GEOS_MARK_FUNCTION;

  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, localIndex > maxRowEntries( 0 );
  forAll< hypre::execPolicy >( localMatrix.numRows(),
                               [localMatrix, maxRowEntries] GEOS_HYPRE_DEVICE ( localIndex const row )
  {
    maxRowEntries.max( localMatrix.numNonZeros( row ) );
  } );

  createWithLocalSize( localMatrix.numRows(), numLocalColumns, maxRowEntries.get(), comm );
  globalIndex const rankOffset = ilower();

  array1d< HYPRE_BigInt > rows;
  rows.resizeWithoutInitializationOrDestruction( hypre::memorySpace, localMatrix.numRows() );

  array1d< HYPRE_Int > sizes;
  sizes.resizeWithoutInitializationOrDestruction( hypre::memorySpace, localMatrix.numRows() );

  array1d< HYPRE_Int > offsets;
  offsets.resizeWithoutInitializationOrDestruction( hypre::memorySpace, localMatrix.numRows() );

  forAll< hypre::execPolicy >( localMatrix.numRows(),
                               [localMatrix, rankOffset,
                                rowsView = rows.toView(),
                                sizesView = sizes.toView(),
                                offsetsView = offsets.toView()] GEOS_HYPRE_DEVICE ( localIndex const row )
  {
    rowsView[row] = LvArray::integerConversion< HYPRE_BigInt >( row + rankOffset );
    sizesView[row] = LvArray::integerConversion< HYPRE_Int >( localMatrix.numNonZeros( row ) );
    offsetsView[row] = LvArray::integerConversion< HYPRE_Int >( localMatrix.getOffsets()[row] );
  } );

  // This is necessary so that localMatrix.getColumns() and localMatrix.getEntries() return device pointers
  localMatrix.move( hypre::memorySpace, false );

  open();
  GEOS_HYPRE_CHECK_DEVICE_ERRORS( "before HYPRE_IJMatrixAddToValues2" );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues2( m_ij_mat,
                                                    localMatrix.numRows(),
                                                    sizes.data(),
                                                    rows.data(),
                                                    offsets.data(),
                                                    localMatrix.getColumns(),
                                                    localMatrix.getEntries() ) );
  close();
}

void HypreMatrix::createWithLocalSize( localIndex const localRows,
                                       localIndex const localCols,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  GEOS_LAI_ASSERT_GE( localRows, 0 );
  GEOS_LAI_ASSERT_GE( localCols, 0 );
  GEOS_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  HYPRE_BigInt const ilower = MpiWrapper::prefixSum< HYPRE_BigInt >( localRows );
  HYPRE_BigInt const iupper = ilower + localRows - 1;

  HYPRE_BigInt const jlower = MpiWrapper::prefixSum< HYPRE_BigInt >( localCols );
  HYPRE_BigInt const jupper = jlower + localCols - 1;

  array1d< HYPRE_Int > row_sizes;
  row_sizes.resizeDefault( localRows, LvArray::integerConversion< HYPRE_Int >( maxEntriesPerRow ) );

  initialize( comm,
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );
}

void HypreMatrix::set( real64 const value )
{
  GEOS_LAI_ASSERT( ready() );
  open();
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetConstantValues( m_ij_mat, value ) );
  close();
}

void HypreMatrix::reset()
{
  MatrixBase::reset();
  if( m_ij_mat )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixDestroy( m_ij_mat ) );
    m_ij_mat = nullptr;
    m_parcsr_mat = nullptr;
  }
}

void HypreMatrix::zero()
{
  set( 0.0 );
}

void HypreMatrix::open()
{
  GEOS_LAI_ASSERT( created() && closed() );
  if( m_assembled )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( m_ij_mat ) );
  }
  m_closed = false;
}

void HypreMatrix::close()
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( !closed() );

  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAssemble( m_ij_mat ) );

  // Get a reference to the constructed matrix object. Done only on the first
  // assembly call when the sparsity pattern of the matrix is defined.
  if( !m_assembled )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixGetObject( m_ij_mat, (void * *) &m_parcsr_mat ) );
    if( !hypre_ParCSRMatrixCommPkg( m_parcsr_mat ) )
    {
      GEOS_LAI_CHECK_ERROR( hypre_MatvecCommPkgCreate( m_parcsr_mat ) );
    }
  }

  m_closed = true;
  m_assembled = true;
}

bool HypreMatrix::created() const
{
  return m_ij_mat != nullptr;
}

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOS_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = 1;
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                   1,
                                                   &ncols,
                                                   hypre::toHypreBigInt( &rowIndex ),
                                                   hypre::toHypreBigInt( &colIndex ),
                                                   &value ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOS_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = 1;
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                 1,
                                                 &ncols,
                                                 hypre::toHypreBigInt( &rowIndex ),
                                                 hypre::toHypreBigInt( &colIndex ),
                                                 &value ) );

}

void HypreMatrix::insert( globalIndex const rowIndex0,
                          globalIndex const colIndex0,
                          real64 const value0 )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( insertable() );

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  array1d< HYPRE_BigInt > rowIndexDevice( 1 );
  array1d< HYPRE_BigInt > colIndexDevice( 1 );
  array1d< HYPRE_Int > ncolsDevice( 1 );
  array1d< real64 > valueDevice( 1 );

  rowIndexDevice[0] = rowIndex0;
  colIndexDevice[0] = colIndex0;
  ncolsDevice[0] = 1;
  valueDevice[0] = value0;

  rowIndexDevice.move( parallelDeviceMemorySpace, false );
  colIndexDevice.move( parallelDeviceMemorySpace, false );
  ncolsDevice.move( parallelDeviceMemorySpace, false );
  valueDevice.move( parallelDeviceMemorySpace, false );

  HYPRE_Int * const ncols = ncolsDevice.data();
  HYPRE_BigInt const * const rowIndex = rowIndexDevice.data();
  HYPRE_BigInt const * const colIndex = colIndexDevice.data();
  real64 * const value = valueDevice.data();
#else
  HYPRE_Int one = 1;
  HYPRE_Int * const ncols = &one;
  HYPRE_BigInt const rowIndexData = rowIndex0;
  HYPRE_BigInt const colIndexData = colIndex0;
  HYPRE_BigInt const * const rowIndex = &rowIndexData;
  HYPRE_BigInt const * const colIndex = &colIndexData;
  real64 const * const value = &value0;
#endif

  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                   1,
                                                   ncols,
                                                   rowIndex,
                                                   colIndex,
                                                   value ) );


}

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( size );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                   1,
                                                   &ncols,
                                                   hypre::toHypreBigInt( &rowIndex ),
                                                   hypre::toHypreBigInt( colIndices ),
                                                   values ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOS_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( size );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                 1,
                                                 &ncols,
                                                 hypre::toHypreBigInt( &rowIndex ),
                                                 hypre::toHypreBigInt( colIndices ),
                                                 values ) );
}

void HypreMatrix::insert( globalIndex const rowIndex0,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( insertable() );

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  array1d< globalIndex > rowIndexDevice( 1 );
  array1d< HYPRE_Int > ncolsDevice( 1 );

  rowIndexDevice[0] = rowIndex0;
  ncolsDevice[0] = LvArray::integerConversion< HYPRE_Int >( size );

  rowIndexDevice.move( parallelDeviceMemorySpace, false );
  ncolsDevice.move( parallelDeviceMemorySpace, false );

  globalIndex const * const rowIndex = rowIndexDevice.data();
  HYPRE_Int * const ncols = ncolsDevice.data();
#else
  globalIndex const * const rowIndex = &rowIndex0;
  HYPRE_Int hypreSize = size;
  HYPRE_Int * const ncols = &hypreSize;
#endif

  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                   1,
                                                   ncols,
                                                   rowIndex,
                                                   hypre::toHypreBigInt( colIndices ),
                                                   values ) );
}

void HypreMatrix::add( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( colIndices.size() );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                   1,
                                                   &ncols,
                                                   hypre::toHypreBigInt( &rowIndex ),
                                                   hypre::toHypreBigInt( colIndices ),
                                                   values ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOS_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( colIndices.size() );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                 1,
                                                 &ncols,
                                                 hypre::toHypreBigInt( &rowIndex ),
                                                 hypre::toHypreBigInt( colIndices ),
                                                 values ) );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice1d< real64 const > const & values )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT( insertable() );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( colIndices.size() );
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                   1,
                                                   &ncols,
                                                   hypre::toHypreBigInt( &rowIndex ),
                                                   hypre::toHypreBigInt( colIndices ),
                                                   values ) );
}

void HypreMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values )
{
  for( localIndex i = 0; i < rowIndices.size(); ++i )
  {
    add( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values )
{
  for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( rowIndices.size() ); ++i )
  {
    set( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice2d< real64 const > const & values )
{
  for( localIndex i = 0; i < rowIndices.size(); ++i )
  {
    insert( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::add( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    add( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::set( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    set( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::insert( globalIndex const * rowIndices,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex const numRows,
                          localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    insert( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::insert( arrayView1d< globalIndex const > const & rowIndices,
                          arrayView1d< globalIndex const > const & colIndices,
                          arrayView1d< real64 const > const & values )
{
  GEOS_MARK_FUNCTION;

  GEOS_LAI_ASSERT_EQ( rowIndices.size(), colIndices.size() );
  GEOS_LAI_ASSERT_GE( rowIndices.size(), values.size() );

  HYPRE_BigInt const numRows = rowIndices.size();

  array1d< HYPRE_Int > nCols( numRows );
  for( int i=0; i<numRows; ++i )
  {
    nCols[i] = 1;
  }
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  rowIndices.move( parallelDeviceMemorySpace, false );
  colIndices.move( parallelDeviceMemorySpace, false );
  values.move( parallelDeviceMemorySpace, false );
  nCols.move( parallelDeviceMemorySpace, false );
#endif
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                   numRows,
                                                   nCols.data(),
                                                   hypre::toHypreBigInt( rowIndices.data() ),
                                                   hypre::toHypreBigInt( colIndices.data() ),
                                                   values.data()
                                                   ) );
}


void HypreMatrix::apply( HypreVector const & src,
                         HypreVector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( numLocalRows(), dst.localSize() );
  GEOS_LAI_ASSERT_EQ( numLocalCols(), src.localSize() );

  GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvec( 1.0,
                                                  m_parcsr_mat,
                                                  src.unwrapped(),
                                                  0.0,
                                                  dst.unwrapped() ) );
  dst.touch();
}

void HypreMatrix::applyTranspose( HypreVector const & src,
                                  HypreVector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( numLocalCols(), dst.localSize() );
  GEOS_LAI_ASSERT_EQ( numLocalRows(), src.localSize() );

  GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvecT( 1.0,
                                                   m_parcsr_mat,
                                                   src.unwrapped(),
                                                   0.0,
                                                   dst.unwrapped() ) );
  dst.touch();
}

void HypreMatrix::multiply( HypreMatrix const & src,
                            HypreMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT_EQ( numLocalCols(), src.numLocalRows() );

  // Compute product
  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParMatmul( m_parcsr_mat, src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );
}

void HypreMatrix::leftMultiplyTranspose( HypreMatrix const & src,
                                         HypreMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT_EQ( numLocalRows(), src.numLocalRows() );

  // Compute product
  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParTMatmul( m_parcsr_mat, src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );
}

void HypreMatrix::rightMultiplyTranspose( HypreMatrix const & src,
                                          HypreMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT_EQ( numLocalCols(), src.numLocalCols() );

  // Transpose this
  HypreMatrix tmp;
  transpose( tmp );

  // Compute product
  src.multiply( tmp, dst );
}

void HypreMatrix::multiplyRAP( HypreMatrix const & R,
                               HypreMatrix const & P,
                               HypreMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( R.ready() );
  GEOS_LAI_ASSERT( P.ready() );
  GEOS_LAI_ASSERT_EQ( numLocalRows(), R.numLocalCols() );
  GEOS_LAI_ASSERT_EQ( numLocalCols(), P.numLocalRows() );

  HypreMatrix Rt;
  R.transpose( Rt );

  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParCSRMatrixRAP( Rt.unwrapped(),
                                                               m_parcsr_mat,
                                                               P.unwrapped() );
  dst.parCSRtoIJ( dst_parcsr );
}

void HypreMatrix::multiplyPtAP( HypreMatrix const & P,
                                HypreMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( P.ready() );
  GEOS_LAI_ASSERT_EQ( numLocalRows(), P.numLocalRows() );
  GEOS_LAI_ASSERT_EQ( numLocalCols(), P.numLocalRows() );

  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParCSRMatrixRAPKT( P.unwrapped(),
                                                                 m_parcsr_mat,
                                                                 P.unwrapped(),
                                                                 0 );

  dst.parCSRtoIJ( dst_parcsr );
}

void HypreMatrix::parCSRtoIJ( HYPRE_ParCSRMatrix const & parCSRMatrix )
{
  reset();
  m_closed = false;

  hypre_IJMatrix * const ijmatrix = hypre_CTAlloc( hypre_IJMatrix, 1, HYPRE_MEMORY_HOST );

  hypre_IJMatrixComm( ijmatrix ) = hypre_ParCSRMatrixComm( parCSRMatrix );
  hypre_IJMatrixObject( ijmatrix ) = parCSRMatrix;
  hypre_IJMatrixTranslator( ijmatrix ) = nullptr;
  hypre_IJMatrixAssumedPart( ijmatrix ) = hypre_ParCSRMatrixAssumedPartition( parCSRMatrix );
  hypre_ParCSRMatrixOwnsAssumedPartition( parCSRMatrix ) = 0;

  hypre_IJMatrixAssembleFlag( ijmatrix ) = 1;

  hypre_IJMatrixObjectType( ijmatrix ) = HYPRE_PARCSR;
#ifdef HYPRE_USING_OPENMP
  hypre_IJMatrixOMPFlag( ijmatrix ) = 1;
#else
  hypre_IJMatrixOMPFlag( ijmatrix ) = 0;
#endif
  hypre_IJMatrixPrintLevel( ijmatrix ) = 0;

  array1d< HYPRE_BigInt > info( 2 );
  if( MpiWrapper::commRank( hypre_IJMatrixComm( ijmatrix ) ) == 0 )
  {
    info( 0 ) = hypre_ParCSRMatrixFirstRowIndex( parCSRMatrix );
    info( 1 ) = hypre_ParCSRMatrixFirstColDiag( parCSRMatrix );
  }
  MpiWrapper::bcast( info.data(), 2, 0, hypre_IJMatrixComm( ijmatrix ) );
  hypre_IJMatrixGlobalFirstRow( ijmatrix ) = info( 0 );
  hypre_IJMatrixGlobalFirstCol( ijmatrix ) = info( 1 );

  hypre_IJMatrixGlobalNumRows( ijmatrix ) = hypre_ParCSRMatrixGlobalNumRows( parCSRMatrix );
  hypre_IJMatrixGlobalNumCols( ijmatrix ) = hypre_ParCSRMatrixGlobalNumCols( parCSRMatrix );

  // Set row partitioning
  hypre_IJMatrixRowPartitioning( ijmatrix )[0] = hypre_ParCSRMatrixRowStarts( parCSRMatrix )[0];
  hypre_IJMatrixRowPartitioning( ijmatrix )[1] = hypre_ParCSRMatrixRowStarts( parCSRMatrix )[1];

  // Set column partitioning
  if( hypre_IJMatrixGlobalNumRows( ijmatrix ) != hypre_IJMatrixGlobalNumCols( ijmatrix ) )
  {
    // Rectangular matrix
    hypre_IJMatrixColPartitioning( ijmatrix )[0] = hypre_ParCSRMatrixColStarts( parCSRMatrix )[0];
    hypre_IJMatrixColPartitioning( ijmatrix )[1] = hypre_ParCSRMatrixColStarts( parCSRMatrix )[1];
  }
  else
  {
    // Square matrix
    hypre_IJMatrixColPartitioning( ijmatrix )[0] = hypre_IJMatrixRowPartitioning( ijmatrix )[0];
    hypre_IJMatrixColPartitioning( ijmatrix )[1] = hypre_IJMatrixRowPartitioning( ijmatrix )[1];
  }

  m_ij_mat = (HYPRE_IJMatrix) ijmatrix;
  close();
}

void HypreMatrix::gemv( real64 const alpha,
                        HypreVector const & x,
                        real64 const beta,
                        HypreVector & y,
                        bool useTranspose ) const
{
  GEOS_LAI_ASSERT( ready() );

  if( !useTranspose )
  {
    GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvec( alpha,
                                                    m_parcsr_mat,
                                                    x.unwrapped(),
                                                    beta,
                                                    y.unwrapped() ) );
  }
  else
  {
    GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvecT( alpha,
                                                     m_parcsr_mat,
                                                     x.unwrapped(),
                                                     beta,
                                                     y.unwrapped() ) );
  }
  y.touch();
}

void HypreMatrix::scale( real64 const scalingFactor )
{
  GEOS_LAI_ASSERT( ready() );
  hypre::scaleMatrixValues( hypre_ParCSRMatrixDiag( m_parcsr_mat ), scalingFactor );
  hypre::scaleMatrixValues( hypre_ParCSRMatrixOffd( m_parcsr_mat ), scalingFactor );
}

void HypreMatrix::leftScale( HypreVector const & vec )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( vec.ready() );
  GEOS_LAI_ASSERT_EQ( vec.localSize(), numLocalRows() );
  hypre::scaleMatrixRows( hypre_ParCSRMatrixDiag( m_parcsr_mat ), hypre_ParVectorLocalVector( vec.unwrapped() ) );
  hypre::scaleMatrixRows( hypre_ParCSRMatrixOffd( m_parcsr_mat ), hypre_ParVectorLocalVector( vec.unwrapped() ) );
}

void HypreMatrix::rescaleRows( arrayView1d< globalIndex const > const & rowIndices,
                               RowSumType const rowSumType )
{
  GEOS_LAI_ASSERT( ready() );

  switch( rowSumType )
  {
    case RowSumType::SumValues:
    {
      hypre::rescaleMatrixRows( unwrapped(), rowIndices, hypre::ops::identity< HYPRE_Real >, hypre::ops::plus< HYPRE_Real > );
      break;
    }
    case RowSumType::SumAbsValues:
    {
      hypre::rescaleMatrixRows( unwrapped(), rowIndices, LvArray::math::abs< HYPRE_Real >, hypre::ops::plus< HYPRE_Real > );
      break;
    }
    case RowSumType::SumSqrValues:
    {
      hypre::rescaleMatrixRows( unwrapped(), rowIndices, LvArray::math::square< HYPRE_Real >, hypre::ops::plus< HYPRE_Real > );
      break;
    }
    case RowSumType::MaxAbsValues:
    {
      hypre::rescaleMatrixRows( unwrapped(), rowIndices, LvArray::math::abs< HYPRE_Real >, LvArray::math::max< HYPRE_Real > );
      break;
    }
  }
}

void HypreMatrix::separateComponentFilter( HypreMatrix & dst,
                                           integer const dofsPerNode ) const
{
  GEOS_MARK_FUNCTION;

  localIndex const maxRowEntries = maxRowLength();
  integer const temp = maxRowEntries % dofsPerNode;
  GEOS_LAI_ASSERT_EQ( temp, 0 );

  CRSMatrix< real64 > tempMat;

  tempMat.resize( numLocalRows(), numGlobalCols(), maxRowEntries / dofsPerNode );
  CRSMatrixView< real64 > const tempMatView = tempMat.toView();

  globalIndex const firstLocalRow = ilower();
  globalIndex const firstLocalCol = jlower();

  hypre::CSRData< true > const diag{ hypre_ParCSRMatrixDiag( unwrapped() ) };

  hypre::CSRData< true > const offd{ hypre_ParCSRMatrixOffd( unwrapped() ) };

  HYPRE_BigInt const * const colMap = hypre::getOffdColumnMap( unwrapped() );

  auto const getComponent = [dofsPerNode] GEOS_HYPRE_DEVICE ( auto const i )
  {
    return LvArray::integerConversion< integer >( i % dofsPerNode );
  };

  forAll< hypre::execPolicy >( numLocalRows(), [diag, offd, tempMatView, getComponent,
                                                firstLocalRow, firstLocalCol, colMap] GEOS_HYPRE_DEVICE ( localIndex const localRow )
  {
    integer const rowComponent = getComponent( firstLocalRow + localRow );
    for( HYPRE_Int k = diag.rowptr[localRow]; k < diag.rowptr[localRow + 1]; ++k )
    {
      HYPRE_BigInt const globalCol = firstLocalCol + diag.colind[k];
      if( getComponent( globalCol ) == rowComponent )
      {
        tempMatView.insertNonZero( localRow, globalCol, diag.values[k] );
      }
    }
    if( offd.ncol > 0 )
    {
      for( HYPRE_Int k = offd.rowptr[localRow]; k < offd.rowptr[localRow + 1]; ++k )
      {
        HYPRE_BigInt const globalCol = colMap[offd.colind[k]];
        if( getComponent( globalCol ) == rowComponent )
        {
          tempMatView.insertNonZero( localRow, globalCol, offd.values[k] );
        }
      }
    }
  } );

  dst.create( tempMatView.toViewConst(), numLocalCols(), comm() );

  dst.setDofManager( dofManager() );
}

void HypreMatrix::addEntries( HypreMatrix const & src,
                              MatrixPatternOp const op,
                              real64 const scale )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( numLocalRows() == src.numLocalRows() );
  GEOS_LAI_ASSERT( numLocalCols() == src.numLocalCols() );

  switch( op )
  {
    case MatrixPatternOp::Restrict:
    {
      hypre::addEntriesRestricted( hypre_ParCSRMatrixDiag( src.unwrapped() ),
                                   hypre::ops::identity< HYPRE_Int >,
                                   hypre_ParCSRMatrixDiag( unwrapped() ),
                                   hypre::ops::identity< HYPRE_Int >,
                                   scale );
      if( hypre_CSRMatrixNumCols( hypre_ParCSRMatrixOffd( unwrapped() ) ) > 0 )
      {
        HYPRE_BigInt const * const src_colmap = hypre::getOffdColumnMap( src.unwrapped() );
        HYPRE_BigInt const * const dst_colmap = hypre::getOffdColumnMap( unwrapped() );
        hypre::addEntriesRestricted( hypre_ParCSRMatrixOffd( src.unwrapped() ),
                                     [src_colmap] GEOS_HYPRE_DEVICE ( auto i ) { return src_colmap[i]; },
                                     hypre_ParCSRMatrixOffd( unwrapped() ),
                                     [dst_colmap] GEOS_HYPRE_DEVICE ( auto i ) { return dst_colmap[i]; },
                                     scale );
      }
      break;
    }
    case MatrixPatternOp::Same:
    case MatrixPatternOp::Subset:
    case MatrixPatternOp::Extend:
    {
      HYPRE_ParCSRMatrix sumMat;
      GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixAdd( 1.0, unwrapped(), scale, src.unwrapped(), &sumMat ) );
      parCSRtoIJ( sumMat );
      break;
    }
  }
}

void HypreMatrix::addDiagonal( HypreVector const & src,
                               real64 const scale )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( numGlobalRows() == numGlobalCols() );
  GEOS_LAI_ASSERT( numLocalRows() == src.localSize() );

  hypre::CSRData< false > const csr{ hypre_ParCSRMatrixDiag( m_parcsr_mat ) };
  arrayView1d< real64 const > const values = src.values();

  if( isEqual( scale, 1.0 ) )
  {
    forAll< hypre::execPolicy >( numLocalRows(), [=] GEOS_HYPRE_DEVICE ( localIndex const localRow )
    {
      // Hypre stores diagonal element at the beginning of each row, we assume it's always present
      csr.values[csr.rowptr[localRow]] += values[localRow];
    } );
  }
  else
  {
    forAll< hypre::execPolicy >( numLocalRows(), [=] GEOS_HYPRE_DEVICE ( localIndex const localRow )
    {
      // Hypre stores diagonal element at the beginning of each row, we assume it's always present
      csr.values[csr.rowptr[localRow]] += scale * values[localRow];
    } );
  }
}

void HypreMatrix::clampEntries( real64 const lo,
                                real64 const hi,
                                bool const excludeDiag )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_ERROR_IF( excludeDiag && numGlobalRows() != numGlobalCols(), "excludeDiag = true, but matrix is not square" );

  hypre::clampMatrixEntries( hypre_ParCSRMatrixDiag( m_parcsr_mat ), lo, hi, excludeDiag );
  hypre::clampMatrixEntries( hypre_ParCSRMatrixOffd( m_parcsr_mat ), lo, hi, false );
}

localIndex HypreMatrix::maxRowLength() const
{
  GEOS_LAI_ASSERT( assembled() );

  HYPRE_Int const * const ia_diag = hypre_CSRMatrixI( hypre_ParCSRMatrixDiag( m_parcsr_mat ) );
  HYPRE_Int const * const ia_offd = hypre_CSRMatrixI( hypre_ParCSRMatrixOffd( m_parcsr_mat ) );

  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, localIndex > localMaxRowLength( 0 );
  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOS_HYPRE_DEVICE ( localIndex const localRow )
  {
    localMaxRowLength.max( (ia_diag[localRow + 1] - ia_diag[localRow]) + (ia_offd[localRow + 1] - ia_offd[localRow] ) );
  } );

  return MpiWrapper::max( localMaxRowLength.get(), comm() );
}

localIndex HypreMatrix::rowLength( globalIndex const globalRowIndex ) const
{
  GEOS_LAI_ASSERT( assembled() );

  localIndex const localRow = LvArray::integerConversion< localIndex >( globalRowIndex - ilower() );
  GEOS_ASSERT( 0 <= localRow && localRow < numLocalRows() );

  HYPRE_Int const * const ia_diag = hypre_CSRMatrixI( hypre_ParCSRMatrixDiag( unwrapped() ) );
  HYPRE_Int const * const ia_offd = hypre_CSRMatrixI( hypre_ParCSRMatrixOffd( unwrapped() ) );
  HYPRE_Int ia_diag_h[2];
  HYPRE_Int ia_offd_h[2];

#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA
  // Don't know if this is faster or slower than launching a kernel. We should deprecate this function in any case.
  cudaMemcpy( ia_diag_h, ia_diag + localRow, 2 * sizeof( HYPRE_Int ), cudaMemcpyDeviceToHost );
  cudaMemcpy( ia_offd_h, ia_offd + localRow, 2 * sizeof( HYPRE_Int ), cudaMemcpyDeviceToHost );
#elif GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  // Don't know if this is faster or slower than launching a kernel. We should deprecate this function in any case.
  hipError_t err = hipMemcpy( ia_diag_h, ia_diag + localRow, 2 * sizeof( HYPRE_Int ), hipMemcpyDeviceToHost );
  GEOS_ERROR_IF( err != 0, GEOS_FMT( "{}", hipGetErrorString( err ) ) );
  err = hipMemcpy( ia_offd_h, ia_offd + localRow, 2 * sizeof( HYPRE_Int ), hipMemcpyDeviceToHost );
  GEOS_ERROR_IF( err != 0, GEOS_FMT( "{}", hipGetErrorString( err ) ) );
#else
  ia_diag_h[0] = ia_diag[localRow]; ia_diag_h[1] = ia_diag[localRow + 1];
  ia_offd_h[0] = ia_offd[localRow]; ia_offd_h[1] = ia_offd[localRow + 1];
#endif

  return LvArray::integerConversion< localIndex >( ( ia_diag_h[1] - ia_diag_h[0] ) + ( ia_offd_h[1] - ia_offd_h[0] ) );
}

void HypreMatrix::getRowLengths( arrayView1d< localIndex > const & lengths ) const
{
  GEOS_LAI_ASSERT( assembled() );
  HYPRE_Int const * const ia_diag = hypre_CSRMatrixI( hypre_ParCSRMatrixDiag( m_parcsr_mat ) );
  HYPRE_Int const * const ia_offd = hypre_CSRMatrixI( hypre_ParCSRMatrixOffd( m_parcsr_mat ) );
  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOS_HYPRE_DEVICE ( localIndex const localRow )
  {
    lengths[localRow] = (ia_diag[localRow + 1] - ia_diag[localRow]) + (ia_offd[localRow + 1] - ia_offd[localRow]);
  } );
}

void HypreMatrix::getRowCopy( globalIndex const globalRowIndex,
                              arraySlice1d< globalIndex > const & colIndices,
                              arraySlice1d< real64 > const & values ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT_GE( globalRowIndex, ilower() );
  GEOS_LAI_ASSERT_GT( iupper(), globalRowIndex );

  HYPRE_BigInt row = LvArray::integerConversion< HYPRE_BigInt >( globalRowIndex );
  HYPRE_Int numEntries = LvArray::integerConversion< HYPRE_Int >( rowLength( globalRowIndex ) );

  GEOS_LAI_ASSERT_GE( colIndices.size(), numEntries );
  GEOS_LAI_ASSERT_GE( values.size(), numEntries );

  // XXX: this is only correct on host! We should deprecate row-wise functions.
  GEOS_LAI_CHECK_ERROR( hypre_IJMatrixGetValuesParCSR( m_ij_mat,
                                                       -1,
                                                       &numEntries,
                                                       &row,
                                                       hypre::toHypreBigInt( colIndices ),
                                                       values ) );
}

void HypreMatrix::extractDiagonal( HypreVector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  HYPRE_Real * const data = hypre_VectorData( hypre_ParVectorLocalVector( dst.unwrapped() ) );
  hypre_CSRMatrixExtractDiagonal( hypre_ParCSRMatrixDiag( m_parcsr_mat ), data, 0 );
  dst.touch();
}

namespace
{

constexpr HYPRE_Int getHypreRowSumType( RowSumType const rowSumType )
{
  switch( rowSumType )
  {
    case RowSumType::SumValues: return 0;
    case RowSumType::SumAbsValues: return 1;
    case RowSumType::SumSqrValues: return 2;
    default: return -1;
  }
}

}

void HypreMatrix::getRowSums( HypreVector & dst,
                              RowSumType const rowSumType ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  switch( rowSumType )
  {
    case RowSumType::SumValues:
    case RowSumType::SumAbsValues:
    case RowSumType::SumSqrValues:
    {
      HYPRE_Real * const values = hypre_VectorData( hypre_ParVectorLocalVector( dst.unwrapped() ) );
      HYPRE_Int const type = getHypreRowSumType( rowSumType );
      hypre_CSRMatrixComputeRowSum( hypre_ParCSRMatrixDiag( m_parcsr_mat ), nullptr, nullptr, values, type, 1.0, "set" );
      if( hypre_CSRMatrixNumCols( hypre_ParCSRMatrixOffd( m_parcsr_mat ) ) > 0 )
      {
        hypre_CSRMatrixComputeRowSum( hypre_ParCSRMatrixOffd( m_parcsr_mat ), nullptr, nullptr, values, type, 1.0, "add" );
      }
      break;
    }
    case RowSumType::MaxAbsValues:
    {
      hypre::computeRowsSums( m_parcsr_mat, dst.unwrapped(), LvArray::math::abs< HYPRE_Real >, LvArray::math::max< HYPRE_Real > );
      break;
    }
  }
  dst.touch();
}

real64 HypreMatrix::clearRow( globalIndex const globalRow,
                              bool const keepDiag,
                              real64 const diagValue )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_ASSERT_GE( globalRow, ilower() );
  GEOS_LAI_ASSERT_GT( iupper(), globalRow );

  // Get local row index
  HYPRE_Int const localRow = LvArray::integerConversion< HYPRE_Int >( getLocalRowID( globalRow ) );

  // Clear row in diagonal block
  hypre::CSRData< false > const diag{ hypre_ParCSRMatrixDiag( m_parcsr_mat ) };
  bool const square = numGlobalRows() == numGlobalCols();
  real64 oldDiag = 0.0;

  for( HYPRE_Int j = diag.rowptr[localRow]; j < diag.rowptr[localRow + 1]; ++j )
  {
    if( square && diag.colind[j] == localRow )
    {
      oldDiag = diag.values[j];
    }
    diag.values[j] = 0.0;
  }

  // Clear row in off-diagonal block
  hypre::CSRData< false > const offd{ hypre_ParCSRMatrixOffd( m_parcsr_mat ) };
  for( HYPRE_Int j = offd.rowptr[localRow]; j < offd.rowptr[localRow + 1]; ++j )
  {
    offd.values[j] = 0.0;
  }

  // Set diagonal value
  real64 const newDiag = keepDiag ? oldDiag : diagValue;
  if( square && std::fabs( newDiag ) > 0.0 )
  {
    set( globalRow, globalRow, newDiag );
  }
  return oldDiag;
}

HYPRE_ParCSRMatrix const & HypreMatrix::unwrapped() const
{
  return m_parcsr_mat;
}

HYPRE_IJMatrix const & HypreMatrix::unwrappedIJ() const
{
  return m_ij_mat;
}

localIndex HypreMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOS_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOS_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return (index >= ilower && index <= iupper ) ? LvArray::integerConversion< localIndex >( index - ilower ) : -1;
}

globalIndex HypreMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOS_LAI_ASSERT( created() );
  GEOS_LAI_ASSERT_GE( index, 0 );
  GEOS_LAI_ASSERT_GT( numLocalRows(), index );
  return ilower() + index;
}

globalIndex HypreMatrix::numGlobalRows() const
{
  GEOS_LAI_ASSERT( created() );
  return hypre_IJMatrixGlobalNumRows( m_ij_mat );
}

globalIndex HypreMatrix::numGlobalCols() const
{
  GEOS_LAI_ASSERT( created() );
  return hypre_IJMatrixGlobalNumCols( m_ij_mat );
}

localIndex HypreMatrix::numLocalRows() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< localIndex >( iupper() - ilower() );
}

localIndex HypreMatrix::numLocalCols() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< localIndex >( jupper() - jlower() );
}

globalIndex HypreMatrix::ilower() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_IJMatrixRowPartitioning( m_ij_mat )[0] );
}

globalIndex HypreMatrix::iupper() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_IJMatrixRowPartitioning( m_ij_mat )[1] );
}

globalIndex HypreMatrix::jlower() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_IJMatrixColPartitioning( m_ij_mat )[0] );
}

globalIndex HypreMatrix::jupper() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_IJMatrixColPartitioning( m_ij_mat )[1] );
}

localIndex HypreMatrix::numLocalNonzeros() const
{
  GEOS_LAI_ASSERT( assembled() );
  HYPRE_Int const nnz_diag = hypre_CSRMatrixNumNonzeros( hypre_ParCSRMatrixDiag( m_parcsr_mat ) );
  HYPRE_Int const nnz_offd = hypre_CSRMatrixNumNonzeros( hypre_ParCSRMatrixOffd( m_parcsr_mat ) );
  return LvArray::integerConversion< localIndex >( nnz_diag + nnz_offd );
}

globalIndex HypreMatrix::numGlobalNonzeros() const
{
  return MpiWrapper::sum( LvArray::integerConversion< globalIndex >( numLocalNonzeros() ), comm() );
}

void HypreMatrix::print( std::ostream & os ) const
{
  GEOS_LAI_ASSERT( ready() );

  int const myRank = MpiWrapper::commRank( comm() );
  int const numProcs = MpiWrapper::commSize( comm() );
  char str[77];

  constexpr char const lineFormat[] = "{:>11}{:>18}{:>18}{:>28.16e}\n";
  constexpr char const headFormat[] = "{:>11}{:>18}{:>18}{:>28}\n";

  if( myRank == 0 )
  {
    GEOS_FMT_TO( str, sizeof( str ), headFormat, "MPI_Process", "GlobalRowID", "GlobalColID", "Value" );
    os << str;
  }

  for( int rank = 0; rank < numProcs; ++rank )
  {
    MpiWrapper::barrier( comm() );
    if( rank == myRank )
    {
      globalIndex const firstRowID = ilower();
      globalIndex const firstDiagColID = jlower();

      hypre::CSRData< true > diag{ hypre_ParCSRMatrixDiag( m_parcsr_mat ) };
      hypre::CSRData< true > offd{ hypre_ParCSRMatrixOffd( m_parcsr_mat ) };
      HYPRE_BigInt const * const colMapOffd = hypre_ParCSRMatrixColMapOffd( m_parcsr_mat );

      for( HYPRE_Int i = 0; i < diag.nrow; ++i )
      {
        for( HYPRE_Int k = diag.rowptr[i]; k < diag.rowptr[i + 1]; ++k )
        {
          GEOS_FMT_TO( str, sizeof( str ), lineFormat,
                       rank,
                       firstRowID + i,
                       firstDiagColID + diag.colind[k],
                       diag.values[k] );
          os << str;
        }
        if( offd.ncol > 0 )
        {
          for( HYPRE_Int k = offd.rowptr[i]; k < offd.rowptr[i + 1]; ++k )
          {
            GEOS_FMT_TO( str, sizeof( str ), lineFormat,
                         rank,
                         firstRowID + i,
                         colMapOffd[offd.colind[k]],
                         offd.values[k] );
            os << str;
          }
        }
      }
    }
  }
}

void HypreMatrix::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOS_LAI_ASSERT( ready() );

  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixPrintIJ( m_parcsr_mat, 1, 1, filename.c_str() ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      int const rank = MpiWrapper::commRank( comm() );

      // Write MatrixMarket header
      if( rank == 0 )
      {
        std::ofstream os( filename );
        GEOS_ERROR_IF( !os, GEOS_FMT( "Unable to open file for writing: {}", filename ) );
        os << "%%MatrixMarket matrix coordinate real general\n";
        os << GEOS_FMT( "{} {} {}\n", numGlobalRows(), numGlobalCols(), numGlobalNonzeros() );
      }

      // Write matrix values
      if( numGlobalRows() > 0 && numGlobalCols() > 0 )
      {
        // Copy distributed parcsr matrix in a local CSR matrix on every process with at least one row
        // Warning: works for a parcsr matrix that is smaller than 2^31-1
        hypre_CSRMatrix * const fullMatrix = hypre_ParCSRMatrixToCSRMatrixAll( m_parcsr_mat );

        // Identify the smallest process where CSRmatrix exists
        int const printRank = MpiWrapper::min( fullMatrix ? rank : MpiWrapper::commSize( comm() ), comm() );

        // Write to file CSRmatrix on one rank
        if( rank == printRank )
        {
          hypre::CSRData< true > csr{ fullMatrix };
          std::ofstream os( filename, std::ios_base::app );
          GEOS_ERROR_IF( !os, GEOS_FMT( "Unable to open file for writing on rank {}: {}", rank, filename ) );
          char str[64];

          for( HYPRE_Int i = 0; i < csr.nrow; i++ )
          {
            for( HYPRE_Int k = csr.rowptr[i]; k < csr.rowptr[i + 1]; k++ )
            {
              // MatrixMarket row/col indices are 1-based
              GEOS_FMT_TO( str, sizeof( str ), "{} {} {:>28.16e}\n", i + 1, csr.colind[k] + 1, csr.values[k] );
              os << str;
            }
          }
        }

        // Destroy temporary matrix
        GEOS_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( fullMatrix ) );
      }
      break;
    }
    default:
    {
      GEOS_ERROR( "Unsupported matrix output format" );
    }
  }
}

real64 HypreMatrix::norm1() const
{
  GEOS_LAI_ASSERT( ready() );

  HypreMatrix matT;
  transpose( matT );
  return matT.normInf();
}

real64 HypreMatrix::normInf() const
{
  GEOS_LAI_ASSERT( ready() );

  hypre::CSRData< true > diag{ hypre_ParCSRMatrixDiag( m_parcsr_mat ) };
  hypre::CSRData< true > offd{ hypre_ParCSRMatrixOffd( m_parcsr_mat ) };

  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, HYPRE_Real > maxRowAbsSum( 0.0 );
  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOS_HYPRE_DEVICE ( localIndex const localRow )
  {
    HYPRE_Real rowAbsSum = 0.0;
    for( HYPRE_Int j = diag.rowptr[localRow]; j < diag.rowptr[localRow + 1]; ++j )
    {
      rowAbsSum += LvArray::math::abs( diag.values[j] );
    }
    for( HYPRE_Int j = offd.rowptr[localRow]; j < offd.rowptr[localRow + 1]; ++j )
    {
      rowAbsSum += LvArray::math::abs( offd.values[j] );
    }
    maxRowAbsSum.max( rowAbsSum );
  } );

  return MpiWrapper::max( maxRowAbsSum.get(), comm() );

}

real64 HypreMatrix::normFrobenius() const
{
  GEOS_LAI_ASSERT( ready() );
  return hypre_ParCSRMatrixFnorm( m_parcsr_mat );
}

real64 HypreMatrix::normMax() const
{
  GEOS_LAI_ASSERT( ready() );
  real64 const maxNorm = std::max( hypre::computeMaxNorm( hypre_ParCSRMatrixDiag( m_parcsr_mat ) ),
                                   hypre::computeMaxNorm( hypre_ParCSRMatrixOffd( m_parcsr_mat ) ) );
  return MpiWrapper::max( maxNorm, comm() );
}

real64 HypreMatrix::normMax( arrayView1d< globalIndex const > const & rowIndices ) const
{
  GEOS_LAI_ASSERT( ready() );
  real64 const maxNorm = std::max( hypre::computeMaxNorm( hypre_ParCSRMatrixDiag( m_parcsr_mat ), rowIndices, ilower() ),
                                   hypre::computeMaxNorm( hypre_ParCSRMatrixOffd( m_parcsr_mat ), rowIndices, ilower() ) );
  return MpiWrapper::max( maxNorm, comm() );
}

void HypreMatrix::rightScale( HypreVector const & vec )
{
  GEOS_LAI_ASSERT( ready() );
  HypreMatrix t;
  transpose( t );
  t.leftScale( vec );
  t.transpose( *this );
}

void HypreMatrix::leftRightScale( HypreVector const & vecLeft,
                                  HypreVector const & vecRight )
{
  leftScale( vecLeft );
  rightScale( vecRight );
}

void HypreMatrix::transpose( HypreMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );

  // Transpose this->m_parcsr_mat
  HYPRE_ParCSRMatrix dst_parcsr;
  GEOS_LAI_CHECK_ERROR( hypre_ParCSRMatrixTranspose( m_parcsr_mat, &dst_parcsr, 1 ) );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );
}

MPI_Comm HypreMatrix::comm() const
{
  GEOS_LAI_ASSERT( created() );
  return hypre_IJMatrixComm( m_ij_mat );
}

}// end namespace geos
