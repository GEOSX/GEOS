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

#ifndef LIFOSTORAGEHOST_HPP
#define LIFOSTORAGEHOST_HPP

#include <deque>
#include <future>
#include <mutex>
#include <condition_variable>
#include <camp/camp.hpp>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"
#include "common/LifoStorageCommon.hpp"

namespace geos
{
/**
 * This class is used to store in a LIFO way buffers, first on device, then on host, then on disk.
 */
template< typename T, typename INDEX_TYPE >
class LifoStorageHost : public LifoStorageCommon< T, INDEX_TYPE >
{
  typedef LifoStorageCommon< T, INDEX_TYPE > baseLifo;

public:


  /**
   * A LIFO storage will store numberOfBuffersToStoreDevice buffer on
   * deevice, numberOfBuffersToStoreHost on host and the rest on disk.
   *
   * @param name                           Prefix of the files used to save the occurenncy of the saved buffer on disk.
   * @param elemCnt                        Number of elments in the LvArray we want to store in the LIFO storage.
   * @param numberOfBuffersToStoreOnHost   Maximum number of array to store on host memory ( -1 = use 80% of remaining memory ).
   * @param maxNumberOfBuffers             Number of arrays expected to be stores in the LIFO.
   */
  LifoStorageHost( std::string name, size_t elemCnt, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    LifoStorageCommon< T, INDEX_TYPE >( name, elemCnt, numberOfBuffersToStoreOnHost, maxNumberOfBuffers ),
    m_pushToHostFutures( maxNumberOfBuffers ),
    m_popFromHostFutures( maxNumberOfBuffers )
  {}

  /**
   * Asynchroneously push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void pushAsync( arrayView1d< T > array ) override final
  {
    //To be sure 2 pushes are not mixed
    pushWait();
    int id = baseLifo::m_bufferCount++;
    GEOS_ERROR_IF( baseLifo::m_hostDeque.capacity() == 0,
                   "Cannot save on a Lifo without host storage (please set lifoSize, lifoOnDevice and lifoOnHost in xml file)" );

    std::packaged_task< void() > task( std::bind( [ this ] ( int pushId, arrayView1d< T > pushedArray ) {
      baseLifo::m_hostDeque.emplaceFront( pushedArray );

      if( baseLifo::m_maxNumberOfBuffers - pushId > (int)baseLifo::m_hostDeque.capacity() )
      {
        LIFO_MARK_SCOPE( geos::lifoStorage::pushAddTasks );
        // This buffer will go to host memory, and maybe on disk
        std::packaged_task< void() > t2( std::bind( &LifoStorageHost< T, INDEX_TYPE >::hostToDisk, this, baseLifo::m_bufferToDiskCount++ ) );
        {
          std::unique_lock< std::mutex > l2( baseLifo::m_task_queue_mutex[1] );
          baseLifo::m_task_queue[1].emplace_back( std::move( t2 ) );
        }
        baseLifo::m_task_queue_not_empty_cond[1].notify_all();
      }
    }, id, array ) );
    m_pushToHostFutures[id] = task.get_future();
    {
      std::unique_lock< std::mutex > lock( baseLifo::m_task_queue_mutex[0] );
      baseLifo::m_task_queue[0].emplace_back( std::move( task ) );
    }
    baseLifo::m_task_queue_not_empty_cond[0].notify_all();
  }

  /**
   * Waits for last push to be terminated
   */
  void pushWait() override final
  {
    if( baseLifo::m_bufferCount > 0 )
    {
      m_pushToHostFutures[baseLifo::m_bufferCount-1].wait();
    }
  }

  /**
   * Asynchroneously copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void popAsync( arrayView1d< T > array ) override final
  {
    int id = --baseLifo::m_bufferCount;

    std::packaged_task< void() > task( std::bind ( [ this ] ( arrayView1d< T > poppedArray ) {
      baseLifo::m_hostDeque.popFront( poppedArray );

      if( baseLifo::m_bufferToDiskCount > 0 )
      {
        LIFO_MARK_SCOPE( geos::LifoStorageHost::popAddTasks );
        // Trigger pull one buffer from host, and maybe from disk
        std::packaged_task< void() > task2( std::bind( &LifoStorageHost< T, INDEX_TYPE >::diskToHost, this, --baseLifo::m_bufferToDiskCount ) );
        {
          std::unique_lock< std::mutex > lock2( baseLifo::m_task_queue_mutex[1] );
          baseLifo::m_task_queue[1].emplace_back( std::move( task2 ) );
        }
        baseLifo::m_task_queue_not_empty_cond[1].notify_all();
      }
    }, array ) );
    m_popFromHostFutures[id] = task.get_future();
    {
      std::unique_lock< std::mutex > lock( baseLifo::m_task_queue_mutex[0] );
      baseLifo::m_task_queue[0].emplace_back( std::move( task ) );
    }
    baseLifo::m_task_queue_not_empty_cond[0].notify_all();
  }

  /**
   * Waits for last pop to be terminated
   */
  void popWait() override final
  {
    if( baseLifo::m_bufferCount < baseLifo::m_maxNumberOfBuffers )
    {
      m_popFromHostFutures[baseLifo::m_bufferCount].wait();
    }
  }

  /**
   * Compute the number of arrays that can be stored on device
   * @param percent Percentage of the remaining device memory that can be dedicated to the LIFO storage.
   * @param bufferSize Size of one buffer
   * @param maxNumberOfBuffers Maximum number of buffers to store in the LIFO storage
   * @return The maximum number of buffer to allocate to fit in the percentage of the available memory.
   */
  static int computeNumberOfBufferOnDevice( int percent, size_t bufferSize, int maxNumberOfBuffers )
  {
    GEOS_UNUSED_VAR( percent, bufferSize, maxNumberOfBuffers );
    return 0;
  }

private:

  // Futures associated to push to host in case we have no device buffers
  std::vector< std::future< void > > m_pushToHostFutures;
  // Futures associated to pop from host in case we have no device buffers
  std::vector< std::future< void > > m_popFromHostFutures;
};
}
#endif // LIFOSTORAGEHOST_HPP
