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
 * @file HypreVector.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_

#include "linearAlgebra/interfaces/VectorBase.hpp"

/**
 * @name Hypre forward declarations.
 *
 * Forward declare hypre's vector structs and pointer aliases in order
 * to avoid including hypre headers and leaking into the rest of GEOSX.
 */
///@{

extern "C"
{
/// ParVector struct forward declaration
struct hypre_ParVector_struct;
}

///@}

namespace geos
{

/**
 * @brief Wrapper class for hypre's ParVector.
 *
 * This class creates and provides basic support for the HYPRE_ParVector object
 * type used in Hypre using the linear-algebraic system interface (IJ interface).
 */
class HypreVector final : private VectorBase< HypreVector >
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
  HypreVector();

  /**
   * @brief Copy constructor.
   * @param src vector to be copied
   */
  HypreVector( HypreVector const & src );

  /**
   * @brief Move constructor.
   * @param src vector to be moved
   */
  HypreVector( HypreVector && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src HypreVector to be copied.
   * @return the new vector.
   */
  HypreVector & operator=( HypreVector const & src );

  /**
   * @brief Move assignment.
   * @param src HypreVector to be moved.
   * @return the new vector.
   */
  HypreVector & operator=( HypreVector && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~HypreVector();

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
   * @copydoc VectorBase<HypreVector>::created
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

  virtual real64 dot( HypreVector const & vec ) const override;

  virtual void copy( HypreVector const & x ) override;

  virtual void axpy( real64 const alpha,
                     HypreVector const & x ) override;

  virtual void axpby( real64 const alpha,
                      HypreVector const & x,
                      real64 const beta ) override;

  virtual void pointwiseProduct( HypreVector const & x,
                                 HypreVector & y ) const override;

  /**
   * @copydoc VectorBase<HypreVector>::norm1
   */
  virtual real64 norm1() const override;

  /**
   * @copydoc VectorBase<HypreVector>::norm2
   */
  virtual real64 norm2() const override;

  /**
   * @copydoc VectorBase<HypreVector>::normInf
   */
  virtual real64 normInf() const override;

  /**
   * @copydoc VectorBase<HypreVector>::globalSize
   */
  virtual globalIndex globalSize() const override;

  /**
   * @copydoc VectorBase<HypreVector>::localSize
   */
  virtual localIndex localSize() const override;

  /**
   * @copydoc VectorBase<HypreVector>::ilower
   */
  virtual globalIndex ilower() const override;

  /**
   * @copydoc VectorBase<HypreVector>::iupper
   */
  virtual globalIndex iupper() const override;

  /**
   * @copydoc VectorBase<HypreVector>::comm
   */
  virtual MPI_Comm comm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

private:

  /// ParVector pointer alias
  using HYPRE_ParVector = struct hypre_ParVector_struct *;

public:

  /**
   * @brief Returns a pointer to the implementation.
   * @return the underlying HYPRE_ParVector object.
   */
  HYPRE_ParVector const & unwrapped() const;

private:

  /**
   * Pointer to underlying HYPRE_ParVector type.
   */
  HYPRE_ParVector m_vec;

};

}// end namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREVECTOR_HPP_*/
