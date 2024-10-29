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
 * @file PetscVector.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_PETSCVECTOR_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_PETSCVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

/**
 * @name PETSc forward declarations.
 *
 * Forward declare PETSc's vector struct and pointer aliases in order
 * to avoid including PETSc headers and leaking into the rest of GEOSX.
 */
///@{

/// Vec struct forward declaration
extern "C" struct _p_Vec;

///@}

namespace geos
{

/**
 * @brief This class creates and provides basic support for Vec
 *        vector object type used in PETSc.
 */
class PetscVector final : private VectorBase< PetscVector >
{
public:

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Empty vector constructor.
   */
  PetscVector();

  /**
   * @brief Copy constructor.
   * @param src PetscVector to be copied.
   */
  PetscVector( PetscVector const & src );

  /**
   * @brief Move constructor
   * @param src PetscVector to move from
   */
  PetscVector( PetscVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src PetscVector to be copied.
   * @return the new vector.
   */
  PetscVector & operator=( PetscVector const & src );

  /**
   * @brief Move assignment.
   * @param src PetscVector to be moved from.
   * @return the new vector.
   */
  PetscVector & operator=( PetscVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~PetscVector();

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
   * @copydoc VectorBase<PetscVector>::created
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

  virtual real64 dot( PetscVector const & vec ) const override;

  virtual void copy( PetscVector const & x ) override;

  virtual void axpy( real64 const alpha,
                     PetscVector const & x ) override;

  virtual void axpby( real64 const alpha,
                      PetscVector const & x,
                      real64 const beta ) override;

  virtual void pointwiseProduct( PetscVector const & x,
                                 PetscVector & y ) const override;

  /**
   * @copydoc VectorBase<PetscVector>::norm1
   */
  virtual real64 norm1() const override;

  /**
   * @copydoc VectorBase<PetscVector>::norm2
   */
  virtual real64 norm2() const override;

  /**
   * @copydoc VectorBase<PetscVector>::normInf
   */
  virtual real64 normInf() const override;

  /**
   * @copydoc VectorBase<PetscVector>::globalSize
   */
  virtual globalIndex globalSize() const override;

  /**
   * @copydoc VectorBase<PetscVector>::localSize
   */
  virtual localIndex localSize() const override;

  /**
   * @copydoc VectorBase<PetscVector>::ilower
   */
  virtual globalIndex ilower() const override;

  /**
   * @copydoc VectorBase<PetscVector>::iupper
   */
  virtual globalIndex iupper() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  /**
   * @copydoc VectorBase<PetscVector>::comm
   */
  virtual MPI_Comm comm() const override;

  ///@}

private:

  /// Alias for PETSc vector struct pointer
  using Vec = struct _p_Vec *;

public:

  /**
   * @brief Returns a const pointer to the underlying Vec.
   * @return the const pointer to the underlying Vec.
   */
  const Vec & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying Vec.
   * @return the non-const pointer to the underlying Vec.
   */
  Vec & unwrapped();

private:

  /**
   * Pointer to underlying PETSc Vec
   */
  Vec m_vec;
};

} // end geos namespace

#endif /*GEOS_LINEARALGEBRA_INTERFACES_PETSCVECTOR_HPP_*/
