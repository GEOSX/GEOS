/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/HypreVector.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

// Just a placeholder to avoid to include two HYPRE header files
//#include "_hypre_IJ_mv.h"
//#include "_hypre_parcsr_mv.h"

// IJMatrix definition
struct hypre_IJMatrix_struct;
typedef struct hypre_IJMatrix_struct *HYPRE_IJMatrix;

// ParCSRMatrix definition
struct hypre_ParCSRMatrix_struct;
typedef struct hypre_ParCSRMatrix_struct *HYPRE_ParCSRMatrix;

namespace geosx
{

/**
 * \class HypreMatrix
 * \brief This class ...
 */
class HypreMatrix final : public LinearOperator<HypreVector>,
                          private MatrixBase<HypreMatrix, HypreVector>
{
public:

  /// @name Constructor/Destructor Methods
  ///@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */

  HypreMatrix();

  /**
   * @brief Copy constructor.
   *
   * Create new matrix from matrix <tt>src</tt>.
   */
  HypreMatrix( HypreMatrix const & src );

  /**
   * @brief Virtual destructor.
   */
  ~HypreMatrix() override;

  ///@}

  using MatrixBase::createWithLocalSize;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::closed;
  using MatrixBase::assembled;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;

  void createWithLocalSize( localIndex const localRows,
                            localIndex const localCols,
                            localIndex const maxEntriesPerRow,
                            MPI_Comm const & comm = MPI_COMM_WORLD ) override;

  void createWithGlobalSize( globalIndex const globalRows,
                             globalIndex const globalCols,
                             localIndex const maxEntriesPerRow,
                             MPI_Comm const & comm = MPI_COMM_WORLD ) override;

  void open() override;

  void close() override;

  bool created() const override;

  void reset() override;

  void set( real64 const value ) override;

  void zero() override;

  void add( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) override;

  void set( globalIndex const rowIndex,
            globalIndex const colIndex,
            real64 const value ) override;

  void insert( globalIndex const rowIndex,
               globalIndex const colIndex,
               real64 const value ) override;

  void add( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) override;

  void set( globalIndex const rowIndex,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const size ) override;

  void insert( globalIndex const rowIndex,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const size ) override;

  void add( globalIndex const rowIndex,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice1d< real64 const > const & values ) override;

  void set( globalIndex const rowIndex,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice1d< real64 const > const & values ) override;

  void insert( globalIndex const rowIndex,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice1d< real64 const > const & values ) override;

  void add( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  void set( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  void insert( arraySlice1d< globalIndex const > const & rowIndices,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  void add( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  void set( arraySlice1d< globalIndex const > const & rowIndices,
            arraySlice1d< globalIndex const > const & colIndices,
            arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  void insert( arraySlice1d< globalIndex const > const & rowIndices,
               arraySlice1d< globalIndex const > const & colIndices,
               arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  void add( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) override;

  void set( globalIndex const * rowIndices,
            globalIndex const * colIndices,
            real64 const * values,
            localIndex const numRows,
            localIndex const numCols ) override;

  void insert( globalIndex const * rowIndices,
               globalIndex const * colIndices,
               real64 const * values,
               localIndex const numRows,
               localIndex const numCols ) override;

  void multiply( HypreVector const & src,
                 HypreVector & dst ) const override;

  void multiply( HypreMatrix const & src,
                 HypreMatrix & dst,
                 bool const closeResult = true ) const override;

  void leftMultiplyTranspose( HypreMatrix const & src,
                              HypreMatrix & dst,
                              bool const closeResult = true ) const override;

  void rightMultiplyTranspose( HypreMatrix const & src,
                               HypreMatrix & dst,
                               bool const closeResult = true ) const override;

  void gemv( real64 const alpha,
             HypreVector const & x,
             real64 const beta,
             HypreVector & y,
             bool useTranspose=false ) const override;

  void scale( real64 const scalingFactor ) override;

  void leftScale( HypreVector const & vec ) override;

  void rightScale( HypreVector const &vec ) override;

  void leftRightScale( HypreVector const &vecLeft,
                       HypreVector const &vecRight ) override;

  void clearRow( globalIndex const row,
                 real64 const diagValue = 0 ) override;

  localIndex maxRowLength() const override;

  localIndex localRowLength( localIndex localRowIndex ) const override;

  localIndex globalRowLength( globalIndex globalRowIndex ) const override;

  void getRowCopy( globalIndex globalRow,
                   array1d< globalIndex > & colIndices,
                   array1d< real64 > & values ) const override;

  real64 getDiagValue( globalIndex globalRow ) const override;

  globalIndex numGlobalRows() const override;

  globalIndex numGlobalCols() const override;

  localIndex numLocalRows() const override;

  localIndex numLocalCols() const override;

  globalIndex ilower() const override;

  globalIndex iupper() const override;

  localIndex numLocalNonzeros() const override;

  globalIndex numGlobalNonzeros() const override;

  real64 normInf() const override;

  real64 norm1() const override;

  real64 normFrobenius() const override;

  localIndex getLocalRowID( globalIndex const index ) const override;

  globalIndex getGlobalRowID( localIndex const index ) const override;

  virtual MPI_Comm getComm() const override;

  void print( std::ostream & os = std::cout ) const override;

  void write( string const & filename,
              MatrixOutputFormat const format ) const override;

  ///@}

  /**
   * @brief Returns a pointer to the underlying HYPRE_IJMatrix object.
   */
  HYPRE_IJMatrix const & unwrapped() const;

  HYPRE_IJMatrix & unwrapped();

  HYPRE_ParCSRMatrix const & unwrappedParCSR() const;

  HYPRE_ParCSRMatrix & unwrappedParCSR();

private:

  /**
   * @brief Returns the index of the first global col owned by that processor.
   */
  globalIndex jlower() const;

  /**
   * @brief Returns the next index after last global col owned by that processor.
   *
   * @note The intention is for [jlower; jupper) to be used as a half-open index range
   */
  globalIndex jupper() const;

  /**
   * @brief Perform a matrix matrix product with Parallel Matrix
   */
  void parCSRtoIJ( HYPRE_ParCSRMatrix const & parCSRMatrix );

  /**
   * Pointer to underlying HYPRE_IJMatrix type.
   */
  HYPRE_IJMatrix m_ij_mat = nullptr;

  /**
   * Pointer to underlying HYPRE_ParCSRMatrix type.
   */
  HYPRE_ParCSRMatrix m_parcsr_mat = nullptr;

};

/**
 * @brief Stream insertion operator for EpetraMatrix
 * @param os the output stream
 * @param matrix the matrix to be printed
 * @return reference to the output stream
 */
std::ostream & operator<<( std::ostream & os,
                           HypreMatrix const & matrix );

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_*/
