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
 * @file Span.hpp
 */

#ifndef GEOS_COMMON_SPAN_HPP
#define GEOS_COMMON_SPAN_HPP

#include "codingUtilities/traits.hpp"
#include "common/logger/Logger.hpp"

#include <memory>
#include <iterator>
#include <type_traits>

namespace geos
{

/**
 * @brief Lightweight non-owning wrapper over a contiguous range of elements.
 * @tparam T type of range element
 *
 * This is a simplified version of C++17's std::span<> that doesn't support
 * compile-time extents.
 */
template< typename T >
class Span
{
public:

  /// Type of range element
  using element_type = T;

  /// Type of underlying value
  using value_type = std::remove_cv_t< T >;

  /// Type used for indexing the range
  using size_type = std::size_t;

  /// Type used for indexing the range
  using difference_type = std::ptrdiff_t;

  /**
   * @brief Construct an empty span.
   */
  constexpr Span() noexcept = default;

  /**
   * @brief Construct a span from pointer and size.
   * @param ptr pointer to beginning of contiguous range
   * @param size number of elements
   */
  Span( T * const ptr, size_type const size ) noexcept
    : m_data( ptr ),
    m_size( size )
  {}

  /**
   * @brief Construct a span from pair of iterators.
   * @tparam ITER type of iterators
   * @param begin iterator to start of the range
   * @param end iterator to end of the range
   */
  template< typename ITER >
  constexpr Span( ITER const begin, ITER const end ) noexcept
    : Span( std::addressof( *begin ), std::distance( begin, end ) )
  {}

  /**
   * @brief Construct a span from a c-array.
   * @tparam N compile-time size of array
   * @param arr reference to array
   */
  template< int N >
  constexpr Span( T (& arr)[N] ) noexcept
    : Span( std::addressof( arr[0] ), N )
  {}

  /**
   * @brief Construct a span from a range-like object (anything that has begin() and end()).
   * @tparam R the range type
   * @param range reference to the range
   * @note The user must guarantee that the range is contiguous in memory.
   */
  template< typename R, typename std::enable_if_t< traits::is_range_like< R > > * = nullptr >
  constexpr Span( R const & range )
    : Span( range.begin(), range.end() )
  {}

  /**
   * @brief @return size of the range
   */
  constexpr size_type size() const noexcept
  {
    return m_size;
  }

  /**
   * @brief @return size of the range in bytes
   */
  constexpr size_type size_bytes() const noexcept
  {
    return size() * sizeof( element_type );
  }

  /**
   * @brief @return @p true iff range is empty
   */
  constexpr bool empty() const noexcept
  {
    return m_size == 0;
  }

  /**
   * @brief @return raw pointer to the data
   */
  constexpr T * data() const noexcept
  {
    return m_data;
  }

  /**
   * @brief @return iterator to start of the range
   */
  constexpr T * begin() const noexcept
  {
    return m_data;
  }

  /**
   * @brief @return iterator past-the-end of the range
   */
  constexpr T * end() const noexcept
  {
    return m_data + m_size;
  }

  /**
   * @brief @return reverse iterator to start or the range
   */
  constexpr std::reverse_iterator< T * > rbegin() const noexcept
  {
    return std::reverse_iterator< T * >( m_data + m_size - 1 );
  }

  /**
   * @brief @return reverse iterator to end of the range
   */
  constexpr std::reverse_iterator< T * > rend() const noexcept
  {
    return std::reverse_iterator< T * >( m_data - 1 );
  }

  /**
   * @brief @return reference to the first element
   * @pre !empty()
   */
  T & front() const
  {
    GEOS_ASSERT_GT( m_size, 0 );
    return m_data[0];
  }

  /**
   * @brief @return reference to the last element
   * @pre !empty()
   */
  T & back() const
  {
    GEOS_ASSERT_GT( m_size, 0 );
    return m_data[m_size-1];
  }

  /**
   * @brief @return reference to @p i-th element
   * @param i element index
   * @pre 0 <= i < size()
   */
  T & operator[]( size_type const i ) const
  {
    GEOS_ASSERT_GT( m_size, i );
    return m_data[i];
  }

  /**
   * @brief @return a new span of @p count starting elements
   * @param count
   */
  Span< element_type > first( size_type const count ) const
  {
    GEOS_ASSERT_GE( m_size, count );
    return { m_data, count };
  }

  /**
   * @brief @return a new span of @p count trailing elements
   * @param count number of elements
   */
  Span< element_type > last( size_type const count ) const
  {
    GEOS_ASSERT_GE( m_size, count );
    return { m_data + (m_size - count), count };
  }

  /**
   * @brief @return a new span of @p count elements starting at @p offset
   * @param offset starting index
   * @param count number of elements
   */
  Span< element_type > subspan( size_type const offset, size_type const count ) const
  {
    GEOS_ASSERT_GE( m_size, offset + count );
    return { m_data + offset, count };
  }

private:

  /// Pointer to contiguous range of elements
  T * m_data{};

  /// Number of elements represented by this span
  size_type m_size{};

};

}

#endif //GEOS_COMMON_SPAN_HPP
