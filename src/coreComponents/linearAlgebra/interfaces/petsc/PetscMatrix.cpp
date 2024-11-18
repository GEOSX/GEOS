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
 * @file PetscMatrix.cpp
 */

#include "PetscMatrix.hpp"

#if !defined(PETSC_USE_64BIT_INDICES)
#define PETSC_USE_64BIT_INDICES
#endif

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/petsc/PetscUtils.hpp"
#include "common/MpiWrapper.hpp"

#include <petscvec.h>
#include <petscmat.h>

#include <cfenv>
#include <numeric>

namespace geos
{

PetscMatrix::PetscMatrix()
  : LinearOperator(),
  MatrixBase()
{}

PetscMatrix::PetscMatrix( PetscMatrix const & src )
  : PetscMatrix()
{
  *this = src;
}

PetscMatrix::PetscMatrix( PetscMatrix && src ) noexcept
  : PetscMatrix()
{
  *this = std::move( src );
}

PetscMatrix & PetscMatrix::operator=( PetscMatrix const & src )
{
  if( &src != this )
  {
    reset();
    if( src.ready() )
    {
      GEOS_LAI_CHECK_ERROR( MatDuplicate( src.m_mat, MAT_COPY_VALUES, &m_mat ) );
      m_assembled = true;
      m_closed = true;
    }
    m_dofManager = src.dofManager();
  }
  return *this;
}

PetscMatrix & PetscMatrix::operator=( PetscMatrix && src ) noexcept
{
  if( &src != this )
  {
    std::swap( m_mat, src.m_mat );
    MatrixBase::operator=( std::move( src ) );
  }
  return *this;
}

PetscMatrix::~PetscMatrix()
{
  reset();
}

void PetscMatrix::createWithLocalSize( localIndex const localRows,
                                       localIndex const localCols,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  GEOS_LAI_ASSERT( closed() );
  GEOS_LAI_ASSERT_GE( localRows, 0 );
  GEOS_LAI_ASSERT_GE( localCols, 0 );
  GEOS_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  // set up matrix
  GEOS_LAI_CHECK_ERROR( MatCreate( comm, &m_mat ) );
  GEOS_LAI_CHECK_ERROR( MatSetType( m_mat, MATMPIAIJ ) );
  GEOS_LAI_CHECK_ERROR( MatSetSizes( m_mat, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE ) );
  GEOS_LAI_CHECK_ERROR( MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr ) );
  GEOS_LAI_CHECK_ERROR( MatSetUp( m_mat ) );
  GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) );
}

void PetscMatrix::createWithGlobalSize( globalIndex const globalRows,
                                        globalIndex const globalCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  GEOS_LAI_ASSERT( closed() );
  GEOS_LAI_ASSERT_GE( globalRows, 0 );
  GEOS_LAI_ASSERT_GE( globalCols, 0 );
  GEOS_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  // set up matrix
  GEOS_LAI_CHECK_ERROR( MatCreate( comm, &m_mat ) );
  GEOS_LAI_CHECK_ERROR( MatSetType( m_mat, MATMPIAIJ ) );
  GEOS_LAI_CHECK_ERROR( MatSetSizes( m_mat, PETSC_DECIDE, PETSC_DECIDE, globalRows, globalCols ) );
  GEOS_LAI_CHECK_ERROR( MatMPIAIJSetPreallocation( m_mat, maxEntriesPerRow, nullptr, maxEntriesPerRow, nullptr ) );
  GEOS_LAI_CHECK_ERROR( MatSetUp( m_mat ) );
  GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) );
}

bool PetscMatrix::created() const
{
  return m_mat != nullptr;
}

void PetscMatrix::reset()
{
  MatrixBase::reset();
  GEOS_LAI_CHECK_ERROR( MatDestroy( &m_mat ) );
}

namespace
{

/// Abstracts a general looping pattern to modify entries in PETSc matrix
template< typename LAMBDA >
void modifyMatrixRows( Mat mat,
                       LAMBDA && lambda )
{
  PetscInt firstRow, lastRow;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( mat, &firstRow, &lastRow ) );

  MPI_Comm comm;
  GEOS_LAI_CHECK_ERROR( PetscObjectGetComm( reinterpret_cast< PetscObject >( mat ), &comm ) );
  PetscInt const maxNumRows = MpiWrapper::max( lastRow - firstRow, comm );

  // Disable global reductions during local calls to MatSetValues
  PetscBool flag;
  GEOS_LAI_CHECK_ERROR( MatGetOption( mat, MAT_NO_OFF_PROC_ENTRIES, &flag ) );
  GEOS_LAI_CHECK_ERROR( MatSetOption( mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE ) );

  array1d< PetscInt > cols;
  array1d< PetscScalar > vals;

  // loop over local rows
  for( PetscInt globalRow = firstRow; globalRow < lastRow; ++globalRow )
  {
    PetscInt numEntries;
    PetscInt const * colsPtr;
    PetscScalar const * valsPtr;

    // Pointers returned from MatGetRow point to an internal buffer in Mat that persists
    // until the next call to MatGetRow. Theoretically we could use it directly to modify
    // entries. We would have to cast away const-ness though, which is not a good idea.
    // Therefore we have to pay for copying column indices and values into another buffer.

    // get entries in row and copy into local buffer
    GEOS_LAI_CHECK_ERROR( MatGetRow( mat, globalRow, &numEntries, &colsPtr, &valsPtr ) );
    cols.resize( numEntries );
    std::copy( colsPtr, colsPtr + numEntries, cols.data() );
    vals.resize( numEntries );
    std::copy( colsPtr, colsPtr + numEntries, cols.data() );
    GEOS_LAI_CHECK_ERROR( MatRestoreRow( mat, globalRow, &numEntries, &colsPtr, &valsPtr ) );

    // call user functions to modify values
    lambda( globalRow, cols.toSliceConst(), vals.toSlice() );

    // set entries to value
    GEOS_LAI_CHECK_ERROR( MatSetValues( mat, 1, &globalRow, cols.size(), cols.data(), vals.data(), INSERT_VALUES ) );
    GEOS_LAI_CHECK_ERROR( MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY ) );
    GEOS_LAI_CHECK_ERROR( MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY ) );
  }

  // ensure all ranks call MatAssemblyEnd the same number of times
  PetscInt const numExtra = maxNumRows - ( lastRow - firstRow );
  for( PetscInt i = 0; i < numExtra; ++i )
  {
    GEOS_LAI_CHECK_ERROR( MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY ) );
    GEOS_LAI_CHECK_ERROR( MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY ) );
  }

  // restore off-proc option
  GEOS_LAI_CHECK_ERROR( MatSetOption( mat, MAT_NO_OFF_PROC_ENTRIES, flag ) );
}

/// Abstracts a general looping pattern to modify entries in PETSc matrix
template< typename LAMBDA >
void modifyMatrixRows( Mat mat,
                       arrayView1d< globalIndex const > const & rowIndices,
                       LAMBDA && lambda )
{
  PetscInt firstRow, lastRow;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( mat, &firstRow, &lastRow ) );

  MPI_Comm comm;
  GEOS_LAI_CHECK_ERROR( PetscObjectGetComm( reinterpret_cast< PetscObject >( mat ), &comm ) );
  PetscInt const maxNumRows = MpiWrapper::max( lastRow - firstRow, comm );

  // Disable global reductions during local calls to MatSetValues
  PetscBool flag;
  GEOS_LAI_CHECK_ERROR( MatGetOption( mat, MAT_NO_OFF_PROC_ENTRIES, &flag ) );
  GEOS_LAI_CHECK_ERROR( MatSetOption( mat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE ) );

  array1d< PetscInt > cols;
  array1d< PetscScalar > vals;

  // loop selected rows
  forAll< serialPolicy >( rowIndices.size(), [=, &cols, &vals]( localIndex const i )
  {
    PetscInt const globalRow = LvArray::integerConversion< PetscInt >( rowIndices[i] );
    GEOS_ASSERT( firstRow <= globalRow && globalRow < lastRow );

    PetscInt numEntries;
    PetscInt const * colsPtr;
    PetscScalar const * valsPtr;

    // Pointers returned from MatGetRow point to an internal buffer in Mat that persists
    // until the next call to MatGetRow. Theoretically we could use it directly to modify
    // entries. We would have to cast away const-ness though, which is not a good idea.
    // Therefore we have to pay for copying column indices and values into another buffer.

    // get entries in row and copy into local buffer
    GEOS_LAI_CHECK_ERROR( MatGetRow( mat, globalRow, &numEntries, &colsPtr, &valsPtr ) );
    cols.resize( numEntries );
    std::copy( colsPtr, colsPtr + numEntries, cols.data() );
    vals.resize( numEntries );
    std::copy( colsPtr, colsPtr + numEntries, cols.data() );
    GEOS_LAI_CHECK_ERROR( MatRestoreRow( mat, globalRow, &numEntries, &colsPtr, &valsPtr ) );

    // call user functions to modify values
    lambda( globalRow, cols.toSliceConst(), vals.toSlice() );

    // set entries to value
    GEOS_LAI_CHECK_ERROR( MatSetValues( mat, 1, &globalRow, cols.size(), cols.data(), vals.data(), INSERT_VALUES ) );
    GEOS_LAI_CHECK_ERROR( MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY ) );
    GEOS_LAI_CHECK_ERROR( MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY ) );
  } );

  // ensure all ranks call MatAssemblyEnd the same number of times
  PetscInt const numExtra = maxNumRows - ( lastRow - firstRow );
  for( PetscInt i = 0; i < numExtra; ++i )
  {
    GEOS_LAI_CHECK_ERROR( MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY ) );
    GEOS_LAI_CHECK_ERROR( MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY ) );
  }

  // restore off-proc option
  GEOS_LAI_CHECK_ERROR( MatSetOption( mat, MAT_NO_OFF_PROC_ENTRIES, flag ) );
}

} // namespace

void PetscMatrix::set( real64 const value )
{
  GEOS_LAI_ASSERT( ready() );

  modifyMatrixRows( m_mat, [=]( PetscInt const globalRow,
                                arraySlice1d< PetscInt const > const cols,
                                arraySlice1d< PetscScalar > const vals )
  {
    GEOS_UNUSED_VAR( globalRow, cols );
    std::fill( vals.begin(), vals.end(), value );
  } );
}

void PetscMatrix::zero()
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( MatZeroEntries( m_mat ) );
}

void PetscMatrix::open()
{
  GEOS_LAI_ASSERT( created() && closed() );
  m_closed = false;
}

void PetscMatrix::close()
{
  GEOS_LAI_ASSERT( !closed() );

  // Initiate global assembly
  GEOS_LAI_CHECK_ERROR( MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY ) );
  GEOS_LAI_CHECK_ERROR( MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY ) );

  // Clear any rows that have been marked for clearing
  if( assembled() )
  {
    GEOS_LAI_CHECK_ERROR( MatZeroRows( m_mat, m_rowsToClear.size(), m_rowsToClear.data(), 0.0, nullptr, nullptr ) );

    for( localIndex i = 0; i < m_rowsToClear.size(); ++i )
    {
      if( std::fabs( m_diagValues[i] ) > 0.0 )
      {
        set( m_rowsToClear[i], m_rowsToClear[i], m_diagValues[i] );
      }
    }
    GEOS_LAI_CHECK_ERROR( MatAssemblyBegin( m_mat, MAT_FINAL_ASSEMBLY ) );
    GEOS_LAI_CHECK_ERROR( MatAssemblyEnd( m_mat, MAT_FINAL_ASSEMBLY ) );
    m_rowsToClear.clear();
    m_diagValues.clear();
  }

  // If this is initial assembly, set some useful options on the matrix
  if( !assembled() )
  {
    GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE ) );
    GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE ) );
    GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE ) );
    GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE ) );
  }
  m_assembled = true;
  m_closed = true;
}

void PetscMatrix::add( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_CHECK_ERROR( MatSetValue( m_mat, rowIndex, colIndex, value, ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_CHECK_ERROR( MatSetValue( m_mat, rowIndex, colIndex, value, INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const rowIndex,
                          globalIndex const colIndex,
                          real64 const value )
{
  GEOS_LAI_ASSERT( insertable() );
  GEOS_LAI_CHECK_ERROR( MatSetValue( m_mat, rowIndex, colIndex, value, ADD_VALUES ) );
}

void PetscMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat, 1, rows, size, petsc::toPetscInt( colIndices ), values, ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOS_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat, 1, rows, size, petsc::toPetscInt( colIndices ), values, INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const rowIndex,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
  GEOS_LAI_ASSERT( insertable() );
  PetscInt rows[1] = {rowIndex};
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat, 1, rows, size, petsc::toPetscInt( colIndices ), values, ADD_VALUES ) );
}

void PetscMatrix::add( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOS_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      1,
                                      rows,
                                      colIndices.size(),
                                      petsc::toPetscInt( colIndices ),
                                      values,
                                      ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOS_LAI_ASSERT( modifiable() );
  PetscInt rows[1] = {rowIndex};
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      1,
                                      rows,
                                      colIndices.size(),
                                      petsc::toPetscInt( colIndices ),
                                      values,
                                      INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const rowIndex,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice1d< real64 const > const & values )
{
  GEOS_LAI_ASSERT( insertable() );
  PetscInt rows[1] = {rowIndex};
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      1,
                                      rows,
                                      colIndices.size(),
                                      petsc::toPetscInt( colIndices ),
                                      values,
                                      ADD_VALUES ) );
}

void PetscMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      rowIndices.size(),
                                      petsc::toPetscInt( rowIndices ),
                                      colIndices.size(),
                                      petsc::toPetscInt( colIndices ),
                                      values.dataIfContiguous(),
                                      ADD_VALUES ) );
}

void PetscMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      rowIndices.size(),
                                      petsc::toPetscInt( rowIndices ),
                                      colIndices.size(),
                                      petsc::toPetscInt( colIndices ),
                                      values.dataIfContiguous(),
                                      INSERT_VALUES ) );
}

void PetscMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice2d< real64 const > const & values )
{
  GEOS_LAI_ASSERT( insertable() );
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      rowIndices.size(),
                                      petsc::toPetscInt( rowIndices ),
                                      colIndices.size(),
                                      petsc::toPetscInt( colIndices ),
                                      values.dataIfContiguous(),
                                      ADD_VALUES ) );
}

void PetscMatrix::add( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      numRows,
                                      petsc::toPetscInt( rowIndices ),
                                      numCols,
                                      petsc::toPetscInt( colIndices ),
                                      values,
                                      ADD_VALUES ) );
}

void PetscMatrix::set( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      numRows,
                                      petsc::toPetscInt( rowIndices ),
                                      numCols,
                                      petsc::toPetscInt( colIndices ),
                                      values,
                                      INSERT_VALUES ) );
}

void PetscMatrix::insert( globalIndex const * rowIndices,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex const numRows,
                          localIndex const numCols )
{
  GEOS_LAI_ASSERT( insertable() );
  GEOS_LAI_CHECK_ERROR( MatSetValues( m_mat,
                                      numRows,
                                      petsc::toPetscInt( rowIndices ),
                                      numCols,
                                      petsc::toPetscInt( colIndices ),
                                      values,
                                      ADD_VALUES ) );
}

void PetscMatrix::insert( arrayView1d< globalIndex const > const & rowIndices,
                          arrayView1d< globalIndex const > const & colIndices,
                          arrayView1d< real64 const > const & values )
{
  GEOS_LAI_ASSERT( insertable() );
  for( localIndex a = 0; a < rowIndices.size(); ++a )
  {
    insert( rowIndices[a], colIndices[a], values[a] );
  }
}

void PetscMatrix::apply( PetscVector const & src,
                         PetscVector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( numGlobalRows(), dst.globalSize() );
  GEOS_LAI_ASSERT_EQ( numGlobalCols(), src.globalSize() );

  GEOS_LAI_CHECK_ERROR( MatMult( m_mat, src.unwrapped(), dst.unwrapped() ) );
  dst.touch();
}

void PetscMatrix::applyTranspose( Vector const & src,
                                  Vector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( numGlobalCols(), dst.globalSize() );
  GEOS_LAI_ASSERT_EQ( numGlobalRows(), src.globalSize() );

  GEOS_LAI_CHECK_ERROR( MatMultTranspose( m_mat, src.unwrapped(), dst.unwrapped() ) );
  dst.touch();
}

void PetscMatrix::multiply( PetscMatrix const & src,
                            PetscMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalRows() );

  dst.reset();
  GEOS_LAI_CHECK_ERROR( MatMatMult( m_mat, src.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  dst.m_assembled = true;
  dst.m_closed = true;
}

void PetscMatrix::leftMultiplyTranspose( PetscMatrix const & src,
                                         PetscMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT_EQ( numGlobalRows(), src.numGlobalRows() );

  dst.reset();
  GEOS_LAI_CHECK_ERROR( MatTransposeMatMult( m_mat, src.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  dst.m_assembled = true;
  dst.m_closed = true;
}

void PetscMatrix::rightMultiplyTranspose( PetscMatrix const & src,
                                          PetscMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalCols() );

  dst.reset();
  GEOS_LAI_CHECK_ERROR( MatMatTransposeMult( m_mat, src.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  dst.m_assembled = true;
  dst.m_closed = true;
}

void PetscMatrix::gemv( real64 const alpha,
                        PetscVector const & x,
                        real64 const beta,
                        PetscVector & y,
                        bool useTranspose ) const
{
  GEOS_LAI_ASSERT( ready() );

  PetscVector x_( x );
  PetscVector b_( y );

  x_.scale( alpha ); // alpha*x_
  y.scale( beta ); // beta*y

  if( useTranspose )
  {
    GEOS_LAI_CHECK_ERROR( MatMultTranspose( m_mat, x_.unwrapped(), b_.unwrapped() ) );
  }
  else
  {
    GEOS_LAI_CHECK_ERROR( MatMult( m_mat, x_.unwrapped(), b_.unwrapped() ) ); // alpha*A*x_ = b_
  }
  GEOS_LAI_CHECK_ERROR( VecAXPY( y.unwrapped(), 1, b_.unwrapped() ) ); // alpha*A*x_ + beta*y = y
  y.touch();
}

void PetscMatrix::scale( real64 const scalingFactor )
{
  GEOS_LAI_ASSERT( ready() );

  if( isEqual( scalingFactor, 1.0 ) )
  {
    return;
  }

  GEOS_LAI_CHECK_ERROR( MatScale( m_mat, scalingFactor ) );
}

void PetscMatrix::leftScale( PetscVector const & vec )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( MatDiagonalScale( m_mat, vec.unwrapped(), nullptr ) );
}

void PetscMatrix::rightScale( PetscVector const & vec )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( MatDiagonalScale( m_mat, nullptr, vec.unwrapped() ) );
}

void PetscMatrix::leftRightScale( PetscVector const & vecLeft,
                                  PetscVector const & vecRight )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( MatDiagonalScale( m_mat, vecLeft.unwrapped(), vecRight.unwrapped() ) );
}

void PetscMatrix::multiplyRAP( PetscMatrix const & R,
                               PetscMatrix const & P,
                               PetscMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( P.ready() );
  GEOS_LAI_ASSERT( R.ready() );
  GEOS_LAI_ASSERT_EQ( R.numGlobalCols(), numGlobalRows() );
  GEOS_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

  dst.reset();

  GEOS_LAI_CHECK_ERROR( MatMatMatMult( R.m_mat, m_mat, P.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );

  dst.m_assembled = true;
}

void PetscMatrix::multiplyPtAP( PetscMatrix const & P,
                                PetscMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( P.ready() );
  GEOS_LAI_ASSERT_EQ( numGlobalRows(), P.numGlobalRows() );
  GEOS_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

  dst.reset();

  // To be able to use the PtAP product in some cases, we need to disable floating point exceptions
  {
    LvArray::system::FloatingPointExceptionGuard guard( FE_ALL_EXCEPT );

    GEOS_LAI_CHECK_ERROR( MatPtAP( m_mat, P.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dst.m_mat ) );
  }

  dst.m_assembled = true;
}

void PetscMatrix::transpose( PetscMatrix & dst ) const
{
  GEOS_LAI_ASSERT( ready() );

  dst.reset();
  GEOS_LAI_CHECK_ERROR( MatTranspose( m_mat, MAT_INITIAL_MATRIX, &dst.m_mat ) );
  dst.m_assembled = true;
}

void PetscMatrix::separateComponentFilter( PetscMatrix & dst,
                                           integer const dofsPerNode ) const
{
  localIndex const maxRowEntries = maxRowLength();
  GEOS_LAI_ASSERT_EQ( maxRowEntries % dofsPerNode, 0 );

  CRSMatrix< real64 > tempMat;
  tempMat.resize( numLocalRows(), numGlobalCols(), maxRowEntries / dofsPerNode );
  CRSMatrixView< real64 > const tempMatView = tempMat.toView();

  PetscInt firstRow, lastRow;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstRow, &lastRow ) );

  auto const getComponent = [dofsPerNode] ( auto i )
  {
    return LvArray::integerConversion< integer >( i % dofsPerNode );
  };

  for( PetscInt row = firstRow; row < lastRow; ++row )
  {
    integer const rowComponent = getComponent( row );
    PetscInt numEntries;
    PetscInt const * cols;
    PetscReal const * vals;
    GEOS_LAI_CHECK_ERROR( MatGetRow( m_mat, row, &numEntries, &cols, &vals ) );
    for( int k = 0; k < numEntries; ++k )
    {
      if( getComponent( cols[k] ) == rowComponent )
      {
        tempMatView.insertNonZero( row - firstRow, cols[k], vals[k] );
      }
    }
    GEOS_LAI_CHECK_ERROR( MatRestoreRow( m_mat, row, &numEntries, nullptr, &vals ) );
  }

  dst.create( tempMatView.toViewConst(), numLocalCols(), MPI_COMM_GEOS );
  dst.setDofManager( dofManager() );
}

real64 PetscMatrix::clearRow( globalIndex const globalRow,
                              bool const keepDiag,
                              real64 const diagValue )
{
  GEOS_LAI_ASSERT( modifiable() );
  GEOS_LAI_ASSERT_GE( globalRow, ilower() );
  GEOS_LAI_ASSERT_GT( iupper(), globalRow );

  PetscInt numEntries = 0;
  PetscInt const * inds;
  PetscReal const * vals;
  GEOS_LAI_CHECK_ERROR( MatGetRow( m_mat, globalRow, &numEntries, &inds, &vals ) );

  bool const square = numGlobalRows() == numGlobalCols();

  real64 oldDiag = 0.0;
  if( square )
  {
    for( localIndex i = 0; i < numEntries; ++i )
    {
      if( inds[i] == globalRow )
      {
        oldDiag = vals[i];
      }
    }
  }
  GEOS_LAI_CHECK_ERROR( MatRestoreRow( m_mat, globalRow, &numEntries, &inds, &vals ) );

  // There is no legitimate way in PETSc to zero out a row in a local fashion,
  // without any collective calls. Therefore, we store the rows to be cleared
  // and perform the clearing at the next close() call.
  m_rowsToClear.emplace_back( globalRow );
  m_diagValues.emplace_back( keepDiag ? oldDiag : diagValue );

  return oldDiag;
}

void PetscMatrix::addEntries( PetscMatrix const & src,
                              MatrixPatternOp const op,
                              real64 const scale )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( numGlobalRows() == src.numGlobalRows() );
  GEOS_LAI_ASSERT( numGlobalCols() == src.numGlobalCols() );

  switch( op )
  {
    case MatrixPatternOp::Same:
    {
      GEOS_LAI_CHECK_ERROR( MatAXPY( m_mat, scale, src.m_mat, SAME_NONZERO_PATTERN ) );
      break;
    }
    case MatrixPatternOp::Subset:
    {
      GEOS_LAI_CHECK_ERROR( MatAXPY( m_mat, scale, src.m_mat, SUBSET_NONZERO_PATTERN ) );
      break;
    }
    case MatrixPatternOp::Extend:
    {
      GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE ) );
      GEOS_LAI_CHECK_ERROR( MatAXPY( m_mat, scale, src.m_mat, DIFFERENT_NONZERO_PATTERN ) );
      GEOS_LAI_CHECK_ERROR( MatSetOption( m_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE ) );
      break;
    }
    case MatrixPatternOp::Restrict:
    {
      GEOS_ERROR( "Not implemented" );
      break;
    }
  }
}

void PetscMatrix::addDiagonal( PetscVector const & src,
                               real64 const scale )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( src.ready() );
  GEOS_LAI_ASSERT( numGlobalRows() == numGlobalCols() );
  GEOS_LAI_ASSERT( numLocalRows() == src.localSize() );

  if( isEqual( scale, 1.0 ) )
  {
    GEOS_LAI_CHECK_ERROR( MatDiagonalSet( m_mat, src.unwrapped(), ADD_VALUES ) );
  }
  else
  {
    PetscVector tmp( src );
    tmp.scale( scale );
    GEOS_LAI_CHECK_ERROR( MatDiagonalSet( m_mat, tmp.unwrapped(), ADD_VALUES ) );
  }
}

void PetscMatrix::clampEntries( real64 const lo,
                                real64 const hi,
                                bool const excludeDiag )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_ERROR_IF( excludeDiag && numGlobalRows() != numGlobalCols(), "excludeDiag = true, but matrix is not square" );

  modifyMatrixRows( m_mat, [=]( PetscInt const globalRow,
                                arraySlice1d< PetscInt const > const cols,
                                arraySlice1d< PetscScalar > const vals )
  {
    for( localIndex k = 0; k < cols.size(); ++k )
    {
      if( !( excludeDiag && cols[k] == globalRow ) )
      {
        vals[k] = LvArray::math::min( hi, LvArray::math::max( lo, vals[k] ) );
      }
    }
  } );
}

localIndex PetscMatrix::maxRowLength() const
{
  GEOS_LAI_ASSERT( assembled() );
  RAJA::ReduceMax< parallelHostReduce, localIndex > maxLocalLength( 0 );
  // Can't use parallel policy here, because of PETSc's single row checkout policy
  forRange< serialPolicy >( ilower(), iupper(), [=]( globalIndex const globalRow )
  {
    maxLocalLength.max( rowLength( globalRow ) );
  } );
  return MpiWrapper::max( maxLocalLength.get(), comm() );
}

localIndex PetscMatrix::rowLength( globalIndex const globalRowIndex ) const
{
  GEOS_LAI_ASSERT( assembled() );
  PetscInt ncols;
  GEOS_LAI_CHECK_ERROR( MatGetRow( m_mat, LvArray::integerConversion< PetscInt >( globalRowIndex ), &ncols, nullptr, nullptr ) );
  localIndex const nnz = ncols;
  GEOS_LAI_CHECK_ERROR( MatRestoreRow( m_mat, LvArray::integerConversion< PetscInt >( globalRowIndex ), &ncols, nullptr, nullptr ) );
  return nnz;
}

void PetscMatrix::getRowLengths( arrayView1d< localIndex > const & lengths ) const
{
  GEOS_LAI_ASSERT( assembled() );
  globalIndex const rowOffset = ilower();
  // Can't use parallel policy here, because of PETSc's single row checkout policy
  forAll< serialPolicy >( numLocalRows(), [=]( localIndex const localRow )
  {
    lengths[localRow] = rowLength( rowOffset + localRow );
  } );
}

void PetscMatrix::getRowCopy( globalIndex const globalRow,
                              arraySlice1d< globalIndex > const & colIndices,
                              arraySlice1d< real64 > const & values ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT_GE( globalRow, ilower() );
  GEOS_LAI_ASSERT_GT( iupper(), globalRow );

  PetscScalar const * vals;
  PetscInt const * inds;
  PetscInt numEntries;
  GEOS_LAI_CHECK_ERROR( MatGetRow( m_mat, LvArray::integerConversion< PetscInt >( globalRow ), &numEntries, &inds, &vals ) );

  GEOS_LAI_ASSERT_GE( colIndices.size(), numEntries );
  GEOS_LAI_ASSERT_GE( values.size(), numEntries );
  std::copy( inds, inds + numEntries, colIndices.dataIfContiguous() );
  std::copy( vals, vals + numEntries, values.dataIfContiguous() );

  GEOS_LAI_CHECK_ERROR( MatRestoreRow( m_mat, LvArray::integerConversion< PetscInt >( globalRow ), &numEntries, &inds, &vals ) );
}

void PetscMatrix::extractDiagonal( PetscVector & dst ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  GEOS_LAI_CHECK_ERROR( MatGetDiagonal( m_mat, dst.unwrapped() ) );
  dst.touch();
}

namespace
{

template< typename R >
double reduceRow( Mat mat,
                  PetscInt const globalRow,
                  R reducer )
{
  PetscScalar const * vals;
  PetscInt const * inds;
  PetscInt numEntries;
  GEOS_LAI_CHECK_ERROR( MatGetRow( mat, globalRow, &numEntries, &inds, &vals ) );
  PetscScalar const res = std::accumulate( vals, vals + numEntries, 0.0, reducer );
  GEOS_LAI_CHECK_ERROR( MatRestoreRow( mat, LvArray::integerConversion< PetscInt >( globalRow ), &numEntries, &inds, &vals ) );
  return res;
}

template< typename F, typename R >
void getRowSumsImpl( Mat const & mat,
                     Vec & vec,
                     F transform,
                     R reduce )
{
  PetscScalar * values;
  GEOS_LAI_CHECK_ERROR( VecGetArray( vec, &values ) );

  PetscInt numLocalRows, firstLocalRow;
  GEOS_LAI_CHECK_ERROR( MatGetLocalSize( mat, &numLocalRows, nullptr ) );
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( mat, &firstLocalRow, nullptr ) );
  auto const reducer = [=]( double acc, double v ){ return reduce( acc, transform( v ) ); };
  for( PetscInt localRow = 0; localRow < numLocalRows; ++localRow )
  {
    values[localRow] = reduceRow( mat, firstLocalRow + localRow, reducer );
  }

  GEOS_LAI_CHECK_ERROR( VecRestoreArray( vec, &values ) );
}

template< typename F, typename R >
void rescaleRowsImpl( Mat mat,
                      arrayView1d< globalIndex const > const & rowIndices,
                      F transform,
                      R reduce )
{
  auto const reducer = [=]( double acc, double v ){ return reduce( acc, transform( v ) ); };
  modifyMatrixRows( mat, rowIndices, [=]( PetscInt const globalRow,
                                          arraySlice1d< PetscInt const > const cols,
                                          arraySlice1d< PetscScalar > const vals )
  {
    GEOS_UNUSED_VAR( globalRow, cols );
    PetscScalar const scale = std::accumulate( vals.begin(), vals.end(), 0.0, reducer );
    std::transform( vals.begin(), vals.end(), vals.begin(), [scale]( double const v ){ return v / scale; } );
  } );
}

} // namespace

void PetscMatrix::getRowSums( PetscMatrix::Vector & dst,
                              RowSumType const rowSumType ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( dst.ready() );
  GEOS_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  switch( rowSumType )
  {
    case RowSumType::SumValues:
    {
      getRowSumsImpl( unwrapped(), dst.unwrapped(), []( auto v ){ return v; }, std::plus<>{} );
      break;
    }
    case RowSumType::SumAbsValues:
    {
      getRowSumsImpl( unwrapped(), dst.unwrapped(), LvArray::math::abs< double >, std::plus<>{} );
      break;
    }
    case RowSumType::SumSqrValues:
    {
      getRowSumsImpl( unwrapped(), dst.unwrapped(), LvArray::math::square< double >, std::plus<>{} );
      break;
    }
    case RowSumType::MaxAbsValues:
    {
      getRowSumsImpl( unwrapped(), dst.unwrapped(), LvArray::math::abs< double >, LvArray::math::max< double > );
      break;
    }
  }
}

void PetscMatrix::rescaleRows( arrayView1d< globalIndex const > const & rowIndices,
                               RowSumType const rowSumType )
{
  GEOS_LAI_ASSERT( ready() );

  switch( rowSumType )
  {
    case RowSumType::SumValues:
    {
      rescaleRowsImpl( unwrapped(), rowIndices, []( auto v ){ return v; }, std::plus<>{} );
      break;
    }
    case RowSumType::SumAbsValues:
    {
      rescaleRowsImpl( unwrapped(), rowIndices, LvArray::math::abs< double >, std::plus<>{} );
      break;
    }
    case RowSumType::SumSqrValues:
    {
      rescaleRowsImpl( unwrapped(), rowIndices, LvArray::math::square< double >, std::plus<>{} );
      break;
    }
    case RowSumType::MaxAbsValues:
    {
      rescaleRowsImpl( unwrapped(), rowIndices, LvArray::math::abs< double >, LvArray::math::max< double > );
      break;
    }
  }
}

Mat & PetscMatrix::unwrapped()
{
  GEOS_LAI_ASSERT( created() );
  return m_mat;
}

const Mat & PetscMatrix::unwrapped() const
{
  GEOS_LAI_ASSERT( created() );
  return m_mat;
}

globalIndex PetscMatrix::numGlobalRows() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt num_rows;
  PetscInt num_cols;
  GEOS_LAI_CHECK_ERROR( MatGetSize( m_mat, &num_rows, &num_cols ) );
  return LvArray::integerConversion< globalIndex >( num_rows );
}

globalIndex PetscMatrix::numGlobalCols() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt num_rows;
  PetscInt num_cols;
  GEOS_LAI_CHECK_ERROR( MatGetSize( m_mat, &num_rows, &num_cols ) );
  return LvArray::integerConversion< globalIndex >( num_cols );
}

globalIndex PetscMatrix::ilower() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt firstRow;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstRow, nullptr ) );
  return LvArray::integerConversion< globalIndex >( firstRow );
}

globalIndex PetscMatrix::iupper() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt lastRow;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, nullptr, &lastRow ) );
  return LvArray::integerConversion< globalIndex >( lastRow );
}

globalIndex PetscMatrix::jlower() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt firstCol;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRangeColumn( m_mat, &firstCol, nullptr ) );
  return LvArray::integerConversion< globalIndex >( firstCol );
}

globalIndex PetscMatrix::jupper() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt lastCol;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRangeColumn( m_mat, nullptr, &lastCol ) );
  return LvArray::integerConversion< globalIndex >( lastCol );
}

localIndex PetscMatrix::numLocalNonzeros() const
{
  GEOS_LAI_ASSERT( assembled() );
  MatInfo info;
  GEOS_LAI_CHECK_ERROR( MatGetInfo( m_mat, MAT_LOCAL, &info ) );
  return static_cast< localIndex >( info.nz_used );
}

globalIndex PetscMatrix::numGlobalNonzeros() const
{
  GEOS_LAI_ASSERT( assembled() );
  MatInfo info;
  GEOS_LAI_CHECK_ERROR( MatGetInfo( m_mat, MAT_GLOBAL_SUM, &info ) );
  return static_cast< globalIndex >( info.nz_used );
}

real64 PetscMatrix::normInf() const
{
  GEOS_LAI_ASSERT( ready() );
  real64 normInf;
  GEOS_LAI_CHECK_ERROR( MatNorm( m_mat, NORM_INFINITY, &normInf ) );
  return normInf;
}

real64 PetscMatrix::norm1() const
{
  GEOS_LAI_ASSERT( ready() );
  real64 norm1;
  GEOS_LAI_CHECK_ERROR( MatNorm( m_mat, NORM_1, &norm1 ) );
  return norm1;
}

real64 PetscMatrix::normFrobenius() const
{
  GEOS_LAI_ASSERT( ready() );
  real64 normFrob;
  GEOS_LAI_CHECK_ERROR( MatNorm( m_mat, NORM_FROBENIUS, &normFrob ) );
  return normFrob;
}

namespace
{
auto const maxAbsReduce = []( real64 const m, PetscReal const v ){ return std::max( m, std::abs( v ) ); };
}

real64 PetscMatrix::normMax() const
{
  GEOS_LAI_ASSERT( ready() );

  PetscInt firstRow, lastRow;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstRow, &lastRow ) );

  real64 norm = 0.0;
  for( PetscInt row = firstRow; row < lastRow; ++row )
  {
    PetscInt numEntries;
    PetscReal const * vals;
    GEOS_LAI_CHECK_ERROR( MatGetRow( m_mat, row, &numEntries, nullptr, &vals ) );
    norm = std::accumulate( vals, vals + numEntries, norm, maxAbsReduce );
    GEOS_LAI_CHECK_ERROR( MatRestoreRow( m_mat, row, &numEntries, nullptr, &vals ) );
  }

  return MpiWrapper::max( norm, comm() );
}

real64 PetscMatrix::normMax( arrayView1d< globalIndex const > const & rowIndices ) const
{
  GEOS_LAI_ASSERT( ready() );

  PetscInt firstRow, lastRow;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &firstRow, &lastRow ) );

  real64 norm = 0.0;
  for( globalIndex const globalRow : rowIndices )
  {
    PetscInt const row = LvArray::integerConversion< PetscInt >( globalRow );
    GEOS_ASSERT( firstRow <= row && row < lastRow );
    PetscInt numEntries;
    PetscReal const * vals;
    GEOS_LAI_CHECK_ERROR( MatGetRow( m_mat, row, &numEntries, nullptr, &vals ) );
    norm = std::accumulate( vals, vals + numEntries, norm, maxAbsReduce );
    GEOS_LAI_CHECK_ERROR( MatRestoreRow( m_mat, row, &numEntries, nullptr, &vals ) );
  }

  return MpiWrapper::max( norm, comm() );
}

localIndex PetscMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &low, &high ) );
  return (index >= low && index < high) ? LvArray::integerConversion< localIndex >( index - low ) : -1;
}

globalIndex PetscMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOS_LAI_ASSERT( created() );
  GEOS_LAI_ASSERT_GE( index, 0 );
  GEOS_LAI_ASSERT_GT( numLocalRows(), index );
  PetscInt low, high;
  GEOS_LAI_CHECK_ERROR( MatGetOwnershipRange( m_mat, &low, &high ) );
  return LvArray::integerConversion< globalIndex >( index + low );
}

localIndex PetscMatrix::numLocalCols() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt cols;
  GEOS_LAI_CHECK_ERROR( MatGetLocalSize( m_mat, nullptr, &cols ) );
  return LvArray::integerConversion< localIndex >( cols );
}

localIndex PetscMatrix::numLocalRows() const
{
  GEOS_LAI_ASSERT( created() );
  PetscInt rows;
  GEOS_LAI_CHECK_ERROR( MatGetLocalSize( m_mat, &rows, nullptr ) );
  return LvArray::integerConversion< localIndex >( rows );
}

MPI_Comm PetscMatrix::comm() const
{
  GEOS_LAI_ASSERT( created() );
  MPI_Comm comm;
  GEOS_LAI_CHECK_ERROR( PetscObjectGetComm( reinterpret_cast< PetscObject >( m_mat ), &comm ) );
  return comm;
}

void PetscMatrix::print( std::ostream & os ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_ERROR_IF( &os != &std::cout, "Only output to stdout currently supported" );
  GEOS_LAI_CHECK_ERROR( MatView( m_mat, PETSC_VIEWER_STDOUT_( comm() ) ) );
}

void PetscMatrix::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOS_LAI_ASSERT( ready() );

  bool ASCIIfile = true;
  PetscViewer viewer;
  PetscViewerFormat petscFormat = PETSC_VIEWER_DEFAULT;

  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_COMMON;
    }
    break;
    case LAIOutputFormat::NATIVE_BINARY:
    {
      ASCIIfile = false;
      petscFormat = PETSC_VIEWER_NATIVE;
    }
    break;
    case LAIOutputFormat::MATLAB_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_MATLAB;
    }
    break;
    case LAIOutputFormat::MATLAB_BINARY:
    {
      ASCIIfile = false;
      petscFormat = PETSC_VIEWER_BINARY_MATLAB;
    }
    break;
    case LAIOutputFormat::MATRIX_MARKET:
    {
      petscFormat = PETSC_VIEWER_ASCII_MATRIXMARKET;
    }
    break;
    default:
      GEOS_ERROR( "Unsupported matrix output format" );
  }

  if( ASCIIfile )
  {
    GEOS_LAI_CHECK_ERROR( PetscViewerASCIIOpen( comm(), filename.c_str(), &viewer ) );
  }
  else
  {
    GEOS_LAI_CHECK_ERROR( PetscViewerBinaryOpen( comm(), filename.c_str(), FILE_MODE_WRITE, &viewer ) );
  }
  GEOS_LAI_CHECK_ERROR( PetscViewerPushFormat( viewer, petscFormat ) );
  GEOS_LAI_CHECK_ERROR( MatView( m_mat, viewer ) );
  GEOS_LAI_CHECK_ERROR( PetscViewerDestroy( &viewer ) );
}

} // end geos namespace
