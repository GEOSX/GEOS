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
 * @file VectorBase.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_VECTORBASE_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_VECTORBASE_HPP_

#include "linearAlgebra/common/common.hpp"
#include "common/MpiWrapper.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

/**
 * @brief Common base template for all vector wrapper types.
 * @tparam VECTOR derived vector type
 *
 * This class template provides a common interface for all derived vector
 * wrapper types. Most methods are pure abstract in order to get the compiler
 * to enforce a common interface in all derived classes; however there is no
 * runtime polymorphism or virtual dispatch - derived classes are not related
 * between each other, and pointers/references to base should never be used -
 * the destructor is protected to enforce that.
 *
 * As an added benefit, the documentation for vector interface also lives here
 * and does not need to be duplicated across vector wrappers. Derived classes
 * should still document specific functions in case a particular LA package's
 * behavior deviates from expectations or has unexpected performance impacts.
 * In that case, @c \@copydoc tag can be used to copy over the documentation.
 */
template< typename VECTOR >
class VectorBase
{
protected:

  /// Alias for VECTOR
  using Vector = VECTOR;

  /**
   * @name Status query methods
   */
  ///@{

  /**
   * @brief Query vector closed status
   * @return @p true if vector has been opened and has not been closed since; @p false otherwise
   */
  inline bool closed() const { return m_closed; }

  /**
   * @brief Query vector creation status
   * @return @p true if vector has been created
   */
  virtual bool created() const = 0;

  /**
   * @brief Query vector ready status
   * @return @p true if vector has been created and is currently closed
   */
  inline bool ready() const { return created() && closed(); }

  ///@}

  /**
   * @name Create Methods
   */
  ///@{

  /**
   * @brief Create a vector based on local number of elements.
   * @param localSize local number of elements
   * @param comm MPI communicator to use
   *
   * Create a vector based on local number of elements.
   * Global size is the sum across processors.
   */
  virtual void create( localIndex const localSize, MPI_Comm const & comm )
  {
    GEOS_UNUSED_VAR( comm );
    GEOS_LAI_ASSERT( closed() );
    GEOS_LAI_ASSERT_GE( localSize, 0 );
    reset();

    // Ideally, resizing to the same size should be a no-op.
    // But bufferManipulation::resize() always forces a touch in host memory.
    // We want to avoid moving values to device again in that case.
    if( m_values.size() != localSize )
    {
      m_values.resize( localSize );
    }
  }

  /**
   * @brief Set a name for the vector (mainly used during various logging).
   * @param name the name
   */
  void setName( string const & name )
  {
    m_values.setName( name );
  }

  ///@}

  /**
   * @name Open/close methods
   */
  ///@{

  /**
   * @brief Open the vector for modifying entries.
   * @return an array view to assemble local values into
   */
  virtual arrayView1d< real64 > open()
  {
    GEOS_LAI_ASSERT( ready() );
    m_closed = false;
    return m_values.toView();
  }

  /**
   * @brief Close vector for modification.
   *
   * After calling this method, the view obtained via open() should not be used.
   */
  virtual void close() = 0;

  /**
   * @brief Notify the vector about external modification through direct data pointer.
   *
   * This method MUST be called after any changes to vector values performed through direct pointer access,
   * so that the underlying Array object can be made aware of external changes in a specific memory space.
   */
  virtual void touch() = 0;

  /**
   * @brief Reset the vector to default state
   */
  virtual void reset()
  {
    // Clearing the array causes issues on GPU when it is later resized again
    // (the bug is in bufferManipulation::resize(), which calls buf.data() without
    // moving the buffer to host, if a capacity increase does not occur, i.e. the new
    // array size is exactly the same as the one prior to clearing).

    //m_values.clear();
    m_closed = true;
  };

  ///@}

  /**
   * @name Modification methods
   */
  ///@{

  /**
   * @brief Set all elements to a constant value.
   * @param value value to set vector elements to
   */
  virtual void set( real64 const value ) = 0;

  /**
   * @brief Set vector elements to zero.
   */
  virtual void zero()
  {
    set( 0.0 );
  }

  /**
   * @brief Set vector elements to random entries.
   * @param seed the random number seed to use
   */
  virtual void rand( unsigned const seed ) = 0;

  ///@}

  /**
   * @name Algebraic Operations
   */
  ///@{

  /**
   * @brief Multiply all elements by factor.
   * @param factor scaling factor
   */
  virtual void scale( real64 const factor ) = 0;

  /**
   * @brief Replace vector elements by their reciprocals
   * @note No guarding is done against division by zero.
   */
  virtual void reciprocal() = 0;

  /**
   * @brief Dot product with the vector vec.
   * @param vec vector to dot-product with
   * @return dot product
   */
  virtual real64 dot( Vector const & vec ) const = 0;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   * @param x vector to copy
   * @note Unlike copy assignment operator, this method expects both vectors to be created
   *       and have identical parallel distributions, and never reallocates memory.
   */
  virtual void copy( Vector const & x ) = 0;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + y</tt>.
   * @param alpha scaling factor for added vector
   * @param x vector to add
   */
  virtual void axpy( real64 const alpha,
                     Vector const & x ) = 0;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + beta*y</tt>.
   * @param alpha scaling factor for added vector
   * @param x vector to add
   * @param beta scaling factor for self vector
   */
  virtual void axpby( real64 const alpha,
                      Vector const & x,
                      real64 const beta ) = 0;

  /**
   * @brief Compute the component-wise multiplication <tt>y</tt> = <tt>v * x</tt>.
   * @param x first vector (input)
   * @param y second vector (output)
   */
  virtual void pointwiseProduct( Vector const & x,
                                 Vector & y ) const = 0;

  /**
   * @brief 1-norm of the vector.
   * @return the 1-norm value
   */
  virtual real64 norm1() const = 0;

  /**
   * @brief 2-norm of the vector.
   * @return the 2-norm value
   */
  virtual real64 norm2() const = 0;

  /**
   * @brief Infinity-norm of the vector.
   * @return the inf-norm value
   */
  virtual real64 normInf() const = 0;

  ///@}

  /**
   * @name Accessor Methods
   */
  ///@{

  /**
   * @brief Returns the global of the vector.
   * @return the global size
   */
  virtual globalIndex globalSize() const = 0;

  /**
   * @brief Returns the local size of the vector.
   * @return the local size (on this processor)
   */
  virtual localIndex localSize() const = 0;

  /**
   * @brief Get lower bound of local partition
   * @return index of the first global row owned by this processor
   */
  virtual globalIndex ilower() const = 0;

  /**
   * @brief Get upper bound of local partition
   * @return next index after last global row owned by that processor
   * @note [ v.ilower(); v.iupper() ) is a half-open index range
   */
  virtual globalIndex iupper() const = 0;

  /**
   * @brief @return a const access view to local vector values
   */
  arrayView1d< real64 const > values() const
  {
    GEOS_LAI_ASSERT( ready() );
    return m_values.toViewConst();
  }

  /**
   * @brief Get the communicator used by this vector
   * @return the MPI communicator
   */
  virtual MPI_Comm comm() const = 0;

  ///@}

  /**
   * @name I/O Methods
   */
  ///@{

  /**
   * @brief Print the vector in Trilinos format to the terminal.
   * @param os the output stream to print to
   */
  virtual void print( std::ostream & os = std::cout ) const = 0;

  /**
   * @brief Write the vector to a file.
   * @param filename name of the output file
   * @param[in] format output format
   */
  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const = 0;

  ///@}

  /**
   * @brief Stream insertion operator for all vector types
   * @param os the output stream to print to
   * @param vec the vector to print
   * @return reference to @p os
   */
  friend std::ostream & operator<<( std::ostream & os, Vector const & vec )
  {
    vec.print( os );
    return os;
  }

  /// Flag indicating whether the vector is closed
  bool m_closed = true;

  /// Actual storage for the local vector values
  array1d< real64 > m_values;
};

} // namespace geos

#endif //GEOS_LINEARALGEBRA_INTERFACES_VECTORBASE_HPP_
