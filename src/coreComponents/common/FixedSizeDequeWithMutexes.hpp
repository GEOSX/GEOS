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

#ifndef FIXEDSIZEDEQUEWITHMUTEXES_HPP
#define FIXEDSIZEDEQUEWITHMUTEXES_HPP

#include <future>
#include <mutex>
#include <condition_variable>
#include <camp/camp.hpp>

#include "common/FixedSizeDeque.hpp"
#include "common/MultiMutexesLock.hpp"

namespace geos
{
/// Associate mutexes with the fixedSizeDeque
template< typename T, typename INDEX_TYPE >
class FixedSizeDequeWithMutexes : public FixedSizeDeque< T, INDEX_TYPE >
{
public:
  /// Mutex to protect access to the front
  std::mutex m_frontMutex;
  /// Mutex to protect access to the back
  std::mutex m_backMutex;
  /// Mutex to prevent two simulteaneous pop (can be an issue for last one)
  std::mutex m_popMutex;
  /// Mutex to prevent two simulteaneous emplace (can be an issue for last one)
  std::mutex m_emplaceMutex;
  /// Condition used to notify when queue is not full
  std::condition_variable_any m_notFullCond;
  /// Condition used to notify when queue is not empty
  std::condition_variable_any m_notEmptyCond;

  /**
   * Create a fixed size double ended queue with associated mutexes and condition variables.
   *
   * @param maxEntries     Maximum number of array to store in the queue.
   * @param valuesPerEntry Number of values in each array of the deque.
   * @param space          Space used to store que queue.
   */
  FixedSizeDequeWithMutexes( int maxEntries, int valuesPerEntry, LvArray::MemorySpace space ): FixedSizeDeque< T, INDEX_TYPE >( maxEntries, valuesPerEntry, space,
#ifdef GEOS_USE_CUDA
                                                                                                                                camp::resources::Resource{ camp::resources::Cuda{} }
#else
                                                                                                                                camp::resources::Resource{ camp::resources::Host{} }
#endif
                                                                                                                                ) {}

  /**
   * Emplace on front from array with locks.
   *
   * @param array The array to emplace
   * @returns Event to sync with the memcpy.
   */
  camp::resources::Event emplaceFront( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    camp::resources::Event e;
    {
      auto lock = make_multilock( m_emplaceMutex, m_frontMutex );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notFullCond.wait( lock, [ this ]  { return !this->full(); } );
      }
      {
        LIFO_MARK_SCOPE( copy );
        e = FixedSizeDeque< T, INDEX_TYPE >::emplace_front( array.toSliceConst() );
      }
    }
    m_notEmptyCond.notify_all();
    return e;
  }

  /**
   * Pop from front to array with locks
   *
   * @param array The array to copy data into
   * @returns Event to sync with the memcpy.
   */
  camp::resources::Event popFront( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    camp::resources::Event e;
    {
      auto lock = make_multilock( m_popMutex, m_frontMutex );
      {
        LIFO_MARK_SCOPE( waitingForBuffer );
        m_notEmptyCond.wait( lock, [ this ]  { return !this->empty(); } );
      }
      // deadlock can occur if frontMutex is taken after an
      // emplaceMutex (inside pushAsync) but this is prevented by the
      // pushWait() in popAsync.
      {
        LIFO_MARK_SCOPE( copy );
        camp::resources::Resource r = this->getStream();
        e = LvArray::memcpy( r, array.toSlice(), this->front() );
        this->pop_front();
      }
    }
    m_notFullCond.notify_all();
    return e;
  }

  /**
   * Emplace front from back of given queue
   *
   * @param q2 The queue to copy data from.
   */
  void emplaceFrontFromBack( FixedSizeDequeWithMutexes< T, INDEX_TYPE > & q2 )
  {
    LIFO_MARK_FUNCTION;
    {
      auto lock = make_multilock( m_emplaceMutex, q2.m_popMutex, m_frontMutex, q2.m_backMutex );
      while( this->full() || q2.empty() )
      {
        {
          LIFO_MARK_SCOPE( WaitForBufferToEmplace );
          m_notFullCond.wait( lock, [ this ]  { return !this->full(); } );
        }
        {
          LIFO_MARK_SCOPE( WaitForBufferToPop );
          q2.m_notEmptyCond.wait( lock, [ &q2 ] { return !q2.empty(); } );
        }
      }
      LIFO_MARK_SCOPE( Transfert );
      this->emplace_front( q2.back() ).wait();
      q2.pop_back();
    }
    q2.m_notFullCond.notify_all();
    m_notEmptyCond.notify_all();
  }

  /**
   * Emplace back from front of given queue
   *
   * @param q2 The queue to copy data from.
   */
  void emplaceBackFromFront( FixedSizeDequeWithMutexes< T, INDEX_TYPE > & q2 )
  {
    LIFO_MARK_FUNCTION;
    {
      auto lock = make_multilock( m_emplaceMutex, q2.m_popMutex, m_backMutex, q2.m_frontMutex );
      while( this->full() || q2.empty() )
      {
        m_notFullCond.wait( lock, [ this ]  { return !this->full(); } );
        q2.m_notEmptyCond.wait( lock, [ &q2 ] { return !q2.empty(); } );
      }
      this->emplace_back( q2.front() ).wait();
      q2.pop_front();
    }
    m_notEmptyCond.notify_all();
    q2.m_notFullCond.notify_all();
  }
};
}
#endif // FIXEDSIZEDEQUEWITHMUTEXES_HPP
