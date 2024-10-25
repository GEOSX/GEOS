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
 * @file EpetraVector.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

class Epetra_Vector;

namespace geos
{

/**
 * @brief Wrapper around Trilinos' Epetra_Vector object.
 */
class EpetraVector final : private VectorBase< EpetraVector >
{
public:

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Empty vector constructor.
   * Create an empty (distributed) vector.
   */
  EpetraVector();

  /**
   * @brief Copy constructor.
   * @param src EpetraVector to be copied.
   */
  EpetraVector( EpetraVector const & src );

  /**
   * @brief Move constructor
   * @param src EpetraVector to move from
   */
  EpetraVector( EpetraVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src EpetraVector to be copied.
   * @return the new vector
   */
  EpetraVector & operator=( EpetraVector const & src );

  /**
   * @brief Move assignment.
   * @param src EpetraVector to be moved from.
   * @return the new vector
   */
  EpetraVector & operator=( EpetraVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~EpetraVector();

  ///@}

  /**
   * @name VectorBase interface
   */
  ///@{

  using VectorBase::setName;
  using VectorBase::closed;
  using VectorBase::ready;
  using VectorBase::open;
  using VectorBase::zero;
  using VectorBase::values;

  /**
   * @copydoc VectorBase<EpetraVector>::created
   */
  virtual bool created() const override;

  virtual void create( localIndex const localSize,
                       MPI_Comm const & comm ) override;

  virtual void close() override;

  virtual void touch() override;

  virtual void reset() override;

  virtual void set( real64 const value ) override;

  virtual void rand( unsigned const seed ) override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void reciprocal() override;

  virtual real64 dot( EpetraVector const & vec ) const override;

  virtual void copy( EpetraVector const & x ) override;

  virtual void axpy( real64 const alpha,
                     EpetraVector const & x ) override;

  virtual void axpby( real64 const alpha,
                      EpetraVector const & x,
                      real64 const beta ) override;

  virtual void pointwiseProduct( EpetraVector const & x,
                                 EpetraVector & y ) const override;

  /**
   * @copydoc VectorBase<EpetraVector>::norm1
   */
  virtual real64 norm1() const override;

  /**
   * @copydoc VectorBase<EpetraVector>::norm2
   */
  virtual real64 norm2() const override;

  /**
   * @copydoc VectorBase<EpetraVector>::normInf
   */
  virtual real64 normInf() const override;

  /**
   * @copydoc VectorBase<EpetraVector>::globalSize
   */
  virtual globalIndex globalSize() const override;

  /**
   * @copydoc VectorBase<EpetraVector>::localSize
   */
  virtual localIndex localSize() const override;

  /**
   * @copydoc VectorBase<EpetraVector>::ilower
   */
  virtual globalIndex ilower() const override;

  /**
   * @copydoc VectorBase<EpetraVector>::iupper
   */
  virtual globalIndex iupper() const override;

  /**
   * @copydoc VectorBase<EpetraVector>::comm
   */
  virtual MPI_Comm comm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Returns a const pointer to the underlying Epetra object.
   * @return const pointer to the underlying Epetra object
   */
  Epetra_Vector const & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying Epetra object.
   * @return non-const pointer to the underlying Epetra object
   */
  Epetra_Vector & unwrapped();

private:

  /// Unique pointer to underlying Epetra_FEVector object.
  std::unique_ptr< Epetra_Vector > m_vec;
};

} // end geos namespace

#endif /*GEOS_LINEARALGEBRA_INTERFACES_EPETRAVECTOR_HPP_*/
