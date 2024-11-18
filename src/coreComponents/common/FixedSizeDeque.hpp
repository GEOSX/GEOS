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

#ifndef FIXEDSIZEDEQUE_HPP
#define FIXEDSIZEDEQUE_HPP

#include "LvArray/src/Array.hpp"
#include "LvArray/src/memcpy.hpp"
#include "LvArray/src/ChaiBuffer.hpp"
#include "common/logger/Logger.hpp"

/// Get the positive value of a module b
#define POSITIVE_MODULO( a, b ) ( ( ( a ) % ( b ) ) + b ) % ( b )
namespace geos
{
template< typename T, typename INDEX_TYPE >
/// Implement a double ended queue with fixed number of fixed size buffer to be stored
class FixedSizeDeque
{
  /// The integer type used for indexing.
  using IndexType = INDEX_TYPE;
  /// 1D array slice
  using ArraySlice1DLarge = LvArray::ArraySlice< T const, 1, 0, std::ptrdiff_t >;
  /// 2D array type. See LvArray:Array for details.
  using Array2D = LvArray::Array< T, 2, camp::make_idx_seq_t< 2 >, std::ptrdiff_t, LvArray::ChaiBuffer >;
public:
  /**
   * Create a fixed size double ended queue.
   *
   * @param maxEntries     Maximum number of array to store in the queue.
   * @param valuesPerEntry Number of values in each array of the deque.
   * @param space          Space used to store que queue.
   * @param stream         Camp resource to perform the copies.
   */
  FixedSizeDeque( IndexType maxEntries, IndexType valuesPerEntry, LvArray::MemorySpace space, camp::resources::Resource stream ):
    m_stream( stream )
  {
    GEOS_ERROR_IF( maxEntries < 0, "Fixed sized queue size must be positive" );
    GEOS_ERROR_IF( valuesPerEntry < 0, "Fixed sized queue array size must be positive" );
    m_storage.resizeWithoutInitializationOrDestruction( space, maxEntries, valuesPerEntry );
  }

  /// @returns true if the queue is empty
  bool empty() const
  {
    return m_begin > m_end;
  }

  /// @returns true if the queue is full
  bool full() const
  {
    return size() == (size_t)m_storage.size( 0 );
  }

  /// @returns the number of arrays stores in the queue
  size_t size() const
  {
    return (size_t)( m_end - m_begin + 1 );
  }

  /// @returns the maximum number of array that can be store in the queue
  size_t capacity() const
  {
    return m_storage.size( 0 );
  }

  /// @returns the first array in the queue
  ArraySlice1DLarge front() const
  {
    GEOS_ERROR_IF( empty(), "Can't get front from empty queue" );
    return m_storage[ POSITIVE_MODULO( m_begin, m_storage.size( 0 ) ) ];
  }

  /// @returns the future first array in the queue after inc_front will be called
  ArraySlice1DLarge next_front() const
  {
    GEOS_ERROR_IF( full(), "Can't increase in a full queue" );
    return m_storage[ POSITIVE_MODULO( m_begin-1, m_storage.size( 0 ) ) ];
  }

  /// @returns the last array of the queue
  ArraySlice1DLarge back() const
  {
    GEOS_ERROR_IF( empty(), "Can't get back from empty queue" );
    return m_storage[ POSITIVE_MODULO( m_end, m_storage.size( 0 ) ) ];
  }

  /// @returns the future last array of the queue when inc_back will be called
  ArraySlice1DLarge next_back() const
  {
    GEOS_ERROR_IF( full(), "Can't increase in a full queue" );
    return m_storage[ POSITIVE_MODULO( m_end+1, m_storage.size( 0 ) ) ];
  }

  /// Removes first array of the queue
  void pop_front()
  {
    GEOS_ERROR_IF( empty(), "Can't pop front from empty queue" );
    m_begin++;
  }

  /// Removes last array of the queue
  void pop_back()
  {
    GEOS_ERROR_IF( empty(), "Can't pop back from empty queue" );
    m_end--;
  }

  /// Add one array (uninitialized) at the front of the queue
  void inc_front()
  {
    GEOS_ERROR_IF( full(), "Can't increase in a full queue" );
    m_begin--;
  }

  /// Add one array (uninitialized) at the end of the queue
  void inc_back()
  {
    GEOS_ERROR_IF( full(), "Can't increase in a full queue" );
    m_end++;
  }

  /**
   * Add one array (copy of src) at the front of the queue
   *
   * @param src Array to emplace at the front of the queue
   * @return Event associated to the copy.
   */
  template< typename INDEX_TYPE2 >
  camp::resources::Event emplace_front( const LvArray::ArraySlice< T const, 1, 0, INDEX_TYPE2 > & src )
  {
    GEOS_ERROR_IF( full(), "Can't emplace in a full  queue" );
    camp::resources::Event e = LvArray::memcpy( m_stream, m_storage[ POSITIVE_MODULO( m_begin-1, m_storage.size( 0 ) ) ], src );
    --m_begin;
    return e;
  }

  /**
   * Add one array (copy of src) at the end of the queue
   *
   * @param src Array to emplace at the end of the queue
   * @return Event associated to the copy.
   */
  template< typename INDEX_TYPE2 >
  camp::resources::Event emplace_back( const LvArray::ArraySlice< T const, 1, 0, INDEX_TYPE2 > & src )
  {
    GEOS_ERROR_IF( full(), "Can't emplace in a full queue" );
    camp::resources::Event e = LvArray::memcpy( m_stream, m_storage[ POSITIVE_MODULO( m_end+1, m_storage.size( 0 ) ) ], src );
    ++m_end;
    return e;

  }

  /**
   * Getter for the associated stream.
   *
   * @return the associated stream.
   */
  camp::resources::Resource getStream()
  {
    return m_stream;
  }


private:
  Array2D m_storage;
  IndexType m_begin = 0;
  IndexType m_end = -1;
  camp::resources::Resource m_stream;
};
}
#endif // FIXEDSIZEDEQUE_HPP
