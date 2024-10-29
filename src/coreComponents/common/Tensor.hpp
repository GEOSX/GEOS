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

#ifndef GEOS_COMMON_TENSOR_HPP_
#define GEOS_COMMON_TENSOR_HPP_

#include "GeosxMacros.hpp"

namespace geos
{

/**
 * @brief Lightweight wrapper around a c-array
 * @tparam T data type
 * @tparam SIZE_TPARAM number of values
 */
template< typename T, int SIZE_TPARAM >
class Tensor
{
public:

  static_assert( SIZE_TPARAM > 0, "Tensor size must be a positive value" );

  /// Alias for type template parameter
  using value_type = T;

  /// Alias for size template parameter
  static constexpr int SIZE = SIZE_TPARAM;

  /**
   * @brief Const element access.
   * @param i element index
   * @return const reference to the i-th element
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  T const & operator[]( std::ptrdiff_t const i ) const
  {
    return data[i];
  }

  /**
   * @brief Non-const element access.
   * @param i element index
   * @return reference to the i-th element
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  T & operator[]( std::ptrdiff_t const i )
  {
    return data[i];
  }

  /**
   * @brief Equality comparison operator.
   * @tparam U dummy parameter to enable SFINAE, do not provide explicitly
   * @param rhs the right-hand side value to compare to
   * @return true iff @code data[i] == rhs.data[i] @endcode
   * @note Version for floating point types that avoids direct equality comparison.
   */
  template< typename U = T >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  std::enable_if_t< std::is_floating_point< U >::value, bool >
  operator==( Tensor< U, SIZE > const & rhs ) const
  {
    for( int i = 0; i < SIZE; ++i )
    {
      if( (data[i] > rhs.data[i]) || (data[i] < rhs.data[i]) )
      {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Equality comparison operator.
   * @tparam U dummy parameter to enable SFINAE, do not provide explicitly
   * @param rhs the right-hand side value to compare to
   * @return true iff @code data[i] == rhs.data[i] @endcode
   * @note Version for all types except floating point.
   */
  template< typename U = T >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  std::enable_if_t< !std::is_floating_point< U >::value, bool >
  operator==( Tensor< U, SIZE > const & rhs ) const
  {
    for( int i = 0; i < SIZE; ++i )
    {
      if( data[i] != rhs.data[i] )
      {
        return false;
      }
    }
    return true;
  }

  /**
   * @brief Returns the size of the tensor
   * @param junk Unused
   * @return The value of the template parameter SIZE.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr int size( int junk ) const
  {
    GEOS_UNUSED_VAR( junk );
    return SIZE;
  }

  /// Declare assignment operators
  GEOS_HOST_DEVICE Tensor< T, SIZE_TPARAM > & operator=( const double & rhs );
  GEOS_HOST_DEVICE Tensor< T, SIZE_TPARAM > & operator+=( const double & rhs );
  GEOS_HOST_DEVICE Tensor< T, SIZE_TPARAM > & operator+=( const Tensor< T, SIZE_TPARAM > & rhs );

  /// Define dot product. TODO: Check compatibility of lhs and rhs.
  friend double GEOS_HOST_DEVICE operator*( const Tensor< T, SIZE_TPARAM > & lhs, const Tensor< T, SIZE_TPARAM > & rhs )
  {
    double result = 0;
    for( int i = 0; i < SIZE_TPARAM; ++i )
    {
      result += lhs.data[i] * rhs.data[i];
    }
    return result;
  };

  /// Underlying array
  T data[SIZE] = {};

private:

  /**
   * @brief Stream insertion operator for Tensor.
   * @param os the output stream
   * @param t the tensor value
   * @return reference to @p os
   */
  friend inline std::ostream & operator<<( std::ostream & os, Tensor< T, SIZE > const & t )
  {
    os << '{' << t.data[0];
    for( int i = 1; i < SIZE; ++i )
    {
      os << ',' << t.data[i];
    }
    os << '}';
    return os;
  }
};

// *****************************************************************************
// ****************************** END DECLARATION ******************************
// *****************************************************************************

/**
 * @brief Assigns all tensor components to a single input value
 * @param rhs the input value
 * @return The updated tensor.
 */
template< typename T, int SIZE_TPARAM >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE Tensor< T, SIZE_TPARAM > &
Tensor< T, SIZE_TPARAM >::operator=( const double & rhs )
{
  for( int i = 0; i < SIZE_TPARAM; ++i )
  {
    data[i] = rhs;
  }
  return *this;
}

/**
 * @brief Adds a single value to all tensor components
 * @param rhs the input value
 * @return The updated tensor.
 */
template< typename T, int SIZE_TPARAM >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE Tensor< T, SIZE_TPARAM > &
Tensor< T, SIZE_TPARAM >::operator+=( const double & rhs )
{
  for( int i = 0; i < SIZE_TPARAM; ++i )
  {
    data[i] += rhs;
  }
  return *this;
}

/**
 * @brief Component-wise addition of two tensors
 * @param rhs the tensor being added to 'this' one
 * @return The updated tensor.
 */
template< typename T, int SIZE_TPARAM >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE Tensor< T, SIZE_TPARAM > &
Tensor< T, SIZE_TPARAM >::operator+=( const Tensor< T, SIZE_TPARAM > & rhs )
{
  for( int i = 0; i < SIZE_TPARAM; ++i )
  {
    data[i] += rhs.data[i];
  }
  return *this;
}

}

#endif //GEOS_COMMON_TENSOR_HPP_
