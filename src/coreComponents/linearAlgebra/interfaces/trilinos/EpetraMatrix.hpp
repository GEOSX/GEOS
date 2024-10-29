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
 * @file EpetraMatrix.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraVector.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraExport.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

class Epetra_Map;
class Epetra_CrsMatrix;
class Epetra_FECrsMatrix;

namespace geos
{

/**
 * @brief Wrapper class for Epetra's CrsMatrix.
 */
class EpetraMatrix final : public virtual LinearOperator< EpetraVector >,
  private MatrixBase< EpetraMatrix, EpetraVector >
{
public:

  /// Compatible vector type
  using Vector = EpetraVector;

  /// Associated exporter type
  using Export = EpetraExport;

  /**
   * @name Constructor/Destructor methods
   */
  ///@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  EpetraMatrix();

  /**
   * @brief Copy constructor.
   * @param[in] src the matrix to be copied
   */
  EpetraMatrix( EpetraMatrix const & src );

  /**
   * @brief Move constructor.
   * @param[in] src the matrix to be copied
   */
  EpetraMatrix( EpetraMatrix && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src matrix to be copied.
   * @return the new vector.
   */
  EpetraMatrix & operator=( EpetraMatrix const & src );

  /**
   * @brief Move assignment.
   * @param src matrix to be moved from.
   * @return the new matrix.
   */
  EpetraMatrix & operator=( EpetraMatrix && src ) noexcept;

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

  using MatrixBase::createWithLocalSize;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::create;
  using MatrixBase::closed;
  using MatrixBase::assembled;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;
  using MatrixBase::residual;
  using MatrixBase::setDofManager;
  using MatrixBase::dofManager;

  virtual void createWithLocalSize( localIndex const localRows,
                                    localIndex const localCols,
                                    localIndex const maxEntriesPerRow,
                                    MPI_Comm const & comm ) override;

  virtual void createWithGlobalSize( globalIndex const globalRows,
                                     globalIndex const globalCols,
                                     localIndex const maxEntriesPerRow,
                                     MPI_Comm const & comm ) override;

  virtual void open() override;

  virtual void close() override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::created
   */
  virtual bool created() const override;

  virtual void reset() override;

  virtual void set( real64 const value ) override;

  virtual void zero() override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value ) override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const size ) override;

  virtual void add( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void set( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void insert( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values ) override;

  virtual void add( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const > const & values ) override;

  virtual void set( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const > const & values ) override;

  virtual void insert( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values ) override;

  virtual void add( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void set( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void insert( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols ) override;

  virtual void insert( arrayView1d< globalIndex const > const & rowIndices,
                       arrayView1d< globalIndex const > const & colIndices,
                       arrayView1d< real64 const > const & values ) override;

  virtual void apply( EpetraVector const & src,
                      EpetraVector & dst ) const override;

  virtual void applyTranspose( EpetraVector const & src,
                               EpetraVector & dst ) const override;

  virtual void multiply( EpetraMatrix const & src,
                         EpetraMatrix & dst ) const override;

  virtual void leftMultiplyTranspose( EpetraMatrix const & src,
                                      EpetraMatrix & dst ) const override;

  virtual void rightMultiplyTranspose( EpetraMatrix const & src,
                                       EpetraMatrix & dst ) const override;

  virtual void multiplyRAP( EpetraMatrix const & R,
                            EpetraMatrix const & P,
                            EpetraMatrix & dst ) const override;

  virtual void multiplyPtAP( EpetraMatrix const & P,
                             EpetraMatrix & dst ) const override;

  virtual void gemv( real64 const alpha,
                     EpetraVector const & x,
                     real64 const beta,
                     EpetraVector & y,
                     bool useTranspose = false ) const override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void leftScale( EpetraVector const & vec ) override;

  virtual void rightScale( EpetraVector const & vec ) override;

  virtual void leftRightScale( EpetraVector const & vecLeft,
                               EpetraVector const & vecRight ) override;

  virtual void rescaleRows( arrayView1d< globalIndex const > const & rowIndices,
                            RowSumType const rowSumType ) override;

  virtual void transpose( EpetraMatrix & dst ) const override;

  virtual void separateComponentFilter( EpetraMatrix & dst,
                                        integer const dofsPerNode ) const override;

  virtual real64 clearRow( globalIndex const row,
                           bool const keepDiag = false,
                           real64 const diagValue = 0.0 ) override;

  virtual void addEntries( EpetraMatrix const & src,
                           MatrixPatternOp const op,
                           real64 const scale ) override;

  virtual void addDiagonal( EpetraVector const & src,
                            real64 const scale ) override;

  virtual void clampEntries( real64 const lo,
                             real64 const hi,
                             bool const excludeDiag ) override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::maxRowLength
   */
  virtual localIndex maxRowLength() const override;

  virtual localIndex rowLength( globalIndex const globalRowIndex ) const override;

  virtual void getRowLengths( arrayView1d< localIndex > const & lengths ) const override;

  virtual void getRowCopy( globalIndex globalRow,
                           arraySlice1d< globalIndex > const & colIndices,
                           arraySlice1d< real64 > const & values ) const override;

  virtual void extractDiagonal( EpetraVector & dst ) const override;

  virtual void getRowSums( EpetraVector & dst,
                           RowSumType const rowSumType ) const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::numGlobalRows
   */
  virtual globalIndex numGlobalRows() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::numGlobalCols
   */
  virtual globalIndex numGlobalCols() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::numLocalRows
   */
  virtual localIndex numLocalRows() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::numLocalCols
   */
  virtual localIndex numLocalCols() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::ilower
   */
  virtual globalIndex ilower() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::iupper
   */
  virtual globalIndex iupper() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::jlower
   */
  virtual globalIndex jlower() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::jupper
   */
  virtual globalIndex jupper() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::numLocalNonzeros
   */
  virtual localIndex numLocalNonzeros() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::numGlobalNonzeros
   */
  virtual globalIndex numGlobalNonzeros() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::normInf
   */
  virtual real64 normInf() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::norm1
   */
  virtual real64 norm1() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::normFrobenius
   */
  virtual real64 normFrobenius() const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::normMax
   */
  virtual real64 normMax() const override;

  virtual real64 normMax( arrayView1d< globalIndex const > const & rowIndices ) const override;

  virtual localIndex getLocalRowID( globalIndex const index ) const override;

  virtual globalIndex getGlobalRowID( localIndex const index ) const override;

  /**
   * @copydoc MatrixBase<EpetraMatrix,EpetraVector>::comm
   */
  virtual MPI_Comm comm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Returns a const pointer to the underlying matrix.
   * @return const pointer to the underlying matrix
   */
  Epetra_FECrsMatrix const & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying matrix.
   * @return non-const pointer to the underlying matrix
   */
  Epetra_FECrsMatrix & unwrapped();

private:

  /**
   * @brief Perform a matrix matrix product with Parallel Matrix
   */
  void multiply( bool const transA,
                 EpetraMatrix const & B,
                 bool const transB,
                 EpetraMatrix & C ) const;

  /**
   * @brief Create the matrix by copying data from an Epetra_CrsMatrix
   * @param src the source matrix
   */
  void create( Epetra_CrsMatrix const & src );

  /// Pointer to the underlying Epetra_CrsMatrix.
  std::unique_ptr< Epetra_FECrsMatrix > m_matrix;

  /// Map representing the parallel partitioning of a source vector (x in y=Ax)
  std::unique_ptr< Epetra_Map > m_src_map;

  /// Map representing the parallel partitioning of a destination vector (y in y=Ax)
  std::unique_ptr< Epetra_Map > m_dst_map;
};

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_*/
