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

#ifndef MULTIMUTEXESLOCK_HPP
#define MULTIMUTEXESLOCK_HPP

#include <mutex>
#include <tuple>

#include "codingUtilities/Utilities.hpp"

namespace geos
{

/**
 * @brief Class to handle locks using 2 mutexes
 */
template< typename ... Mutexes >
class MultiMutexesLock
{
public:
  /**
   * @brief Construct a multi mutexes lock and lock the mutexes.
   * @param mutexes The mutexes associated with the lock.
   */
  MultiMutexesLock( Mutexes &... mutexes )
    : m_islocked( false ), m_mutexes( mutexes ... )
  {
    lock();
  }

  /**
   * @brief Unlock the mutexes and destroy the locks.
   */
  ~MultiMutexesLock()
  {
    unlock();
  }

  /**
   * @brief Lock the two mutexes using std::lock is not already owning lock
   */
  void lock()
  {
    if( m_islocked ) return;
    apply( []( auto && ... mutexes ){ std::lock( mutexes ... ); }, m_mutexes );
    m_islocked = true;
  }

  /**
   * @brief Unlock the two mutexes is owning them.
   */
  void unlock()
  {
    if( !m_islocked ) return;
    forEachArgInTuple( m_mutexes, []( auto & mutex, auto ){ mutex.unlock(); } );
    m_islocked = false;
  }

private:
  /// Indicate if the mutexes are owned by the lock
  bool m_islocked;
  /// Array of references to the mutexes
  std::tuple< Mutexes &... > m_mutexes;
};

/**
 * @brief Helper to construct MultiMutexesLock (usage auto lock = make_multilock( mutex1, mutex2, ... ))
 * @param mutexes List of mutex parameters
 * @return A MultiMutexesLock using the mutexes given in parameter.
 */
template< typename ... Mutexes >
auto make_multilock( Mutexes && ... mutexes )
{
  return MultiMutexesLock< Mutexes... >( std::forward< Mutexes >( mutexes )... );
}
}
#endif // MULTIMUTEXESLOCK_HPP
