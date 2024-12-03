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

#ifndef LIFOSTORAGECUDA_HPP
#define LIFOSTORAGECUDA_HPP

#include <deque>
#include <future>
#include <mutex>
#include <condition_variable>
#include <camp/camp.hpp>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>

#include "common/FixedSizeDeque.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"
#include "common/FixedSizeDequeWithMutexes.hpp"
#include "common/MultiMutexesLock.hpp"
#include "common/LifoStorageCommon.hpp"

namespace geos
{
/**
 * This class is used to store in a LIFO way buffers, first on device, then on host, then on disk.
 */

template< typename T, typename INDEX_TYPE >
class LifoStorageCuda : public LifoStorageCommon< T, INDEX_TYPE >
{
  typedef LifoStorageCommon< T, INDEX_TYPE > baseLifo;
public:


  /**
   * A LIFO storage will store numberOfBuffersToStoreDevice buffer on
   * deevice, numberOfBuffersToStoreHost on host and the rest on disk.
   *
   * @param name                           Prefix of the files used to save the occurenncy of the saved buffer on disk.
   * @param elemCnt                        Number of elments in the LvArray we want to store in the LIFO storage.
   * @param numberOfBuffersToStoreOnDevice Maximum number of array to store on device memory ( -1 = use 80% of remaining memory ).
   * @param numberOfBuffersToStoreOnHost   Maximum number of array to store on host memory ( -1 = use 80% of remaining memory ).
   * @param maxNumberOfBuffers             Number of arrays expected to be stores in the LIFO.
   */
  LifoStorageCuda( std::string name, size_t elemCnt, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    LifoStorageCommon< T, INDEX_TYPE >( name, elemCnt, numberOfBuffersToStoreOnHost, maxNumberOfBuffers ),
    m_deviceDeque( numberOfBuffersToStoreOnDevice, elemCnt, LvArray::MemorySpace::cuda ),
    m_pushToDeviceEvents( maxNumberOfBuffers ),
    m_popFromDeviceEvents( maxNumberOfBuffers )
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
    GEOS_ERROR_IF( baseLifo::m_hostDeque.capacity() == 0 && m_deviceDeque.capacity() < baseLifo::m_maxNumberOfBuffers,
                   "Cannot save on a Lifo without host storage (please set lifoSize, lifoOnDevice and lifoOnHost in xml file)" );

    m_pushToDeviceEvents[id] = m_deviceDeque.emplaceFront( array );

    if( baseLifo::m_maxNumberOfBuffers - id > (int)m_deviceDeque.capacity() )
    {
      LIFO_MARK_SCOPE( geos::lifoStorage::pushAddTasks );
      // This buffer will go to host memory, and maybe on disk
      std::packaged_task< void() > task( std::bind( &LifoStorageCuda< T, INDEX_TYPE >::deviceToHost, this, baseLifo::m_bufferToHostCount++ ) );
      {
        std::unique_lock< std::mutex > lock( baseLifo::m_task_queue_mutex[0] );
        baseLifo::m_task_queue[0].emplace_back( std::move( task ) );
      }
      baseLifo::m_task_queue_not_empty_cond[0].notify_all();
    }
  }

  /**
   * Waits for last push to be terminated
   */
  void pushWait() override final
  {
    if( baseLifo::m_bufferCount > 0 )
    {
      m_pushToDeviceEvents[baseLifo::m_bufferCount-1].wait();
    }
  }

  /**
   * Asynchroneously copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void popAsync( arrayView1d< T > array ) override final
  {
    LIFO_MARK_FUNCTION;
    int id = --baseLifo::m_bufferCount;

    if( baseLifo::m_bufferToHostCount > 0 )
    {
      LIFO_MARK_SCOPE( geos::LifoStorageCuda::popAddTasks );
      // Trigger pull one buffer from host, and maybe from disk
      std::packaged_task< void() > task( std::bind( &LifoStorageCuda< T, INDEX_TYPE >::hostToDevice, this, --baseLifo::m_bufferToHostCount, id ) );
      {
        std::unique_lock< std::mutex > lock( baseLifo::m_task_queue_mutex[0] );
        baseLifo::m_task_queue[0].emplace_back( std::move( task ) );
      }
      baseLifo::m_task_queue_not_empty_cond[0].notify_all();
    }
    m_popFromDeviceEvents[id] = m_deviceDeque.popFront( array );
  }

  /**
   * Waits for last pop to be terminated
   */
  void popWait() override final
  {
    if( baseLifo::m_bufferCount < baseLifo::m_maxNumberOfBuffers )
    {
      int bufferCount = baseLifo::m_bufferCount;
      m_popFromDeviceEvents[bufferCount].wait();
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
    GEOS_ERROR_IF( percent > 100, "Error, percentage of memory should be smaller than 100, check lifoOnDevice (should be greater than -100)" );
    size_t free, total;
    GEOS_ERROR_IF( cudaSuccess != cudaMemGetInfo( &free, &total ), "Error getting CUDA device available memory" );
    double freeGB = ( ( double ) free ) / ( 1024.0 * 1024.0 * 1024.0 );
    LIFO_LOG_RANK( " LIFO : available memory on device " << freeGB << " GB" );
    return std::min( ( int )( 0.01 * percent * free / bufferSize ), maxNumberOfBuffers );
  }

private:

  /**
   * Copy data from device memory to host memory
   *
   * @param ID of the buffer to copy on host.
   */
  void deviceToHost( int id )
  {
    LIFO_MARK_FUNCTION;
    // The copy to host will only start when the data is copied on device buffer
    baseLifo::m_hostDeque.getStream().wait_for( const_cast< camp::resources::Event * >( &m_pushToDeviceEvents[id] ) );
    baseLifo::m_hostDeque.emplaceFrontFromBack( m_deviceDeque );

    if( baseLifo::m_maxNumberOfBuffers - id > (int)(m_deviceDeque.capacity() + baseLifo::m_hostDeque.capacity()) )
    {
      // This buffer will go to host then maybe to disk
      std::packaged_task< void() > task( std::bind( &LifoStorageCuda< T, INDEX_TYPE >::hostToDisk, this, baseLifo::m_bufferToDiskCount++ ) );
      {
        std::unique_lock< std::mutex > lock( baseLifo::m_task_queue_mutex[1] );
        baseLifo::m_task_queue[1].emplace_back( std::move( task ) );
      }
      baseLifo::m_task_queue_not_empty_cond[1].notify_all();
    }
  }

  /**
   * Copy data from host memory to device memory
   *
   * @param id ID of the buffer to load from host memory.
   * @param id_pop ID of the last popped buffer from device
   */
  void hostToDevice( int id, int id_pop )
  {
    LIFO_MARK_FUNCTION;
    // enqueue diskToHost on worker #2 if needed
    if( baseLifo::m_bufferToDiskCount > 0 )
    {
      // This buffer will go to host then to disk
      std::packaged_task< void() > task( std::bind( &LifoStorageCuda< T, INDEX_TYPE >::diskToHost, this, --baseLifo::m_bufferToDiskCount ) );
      {
        std::unique_lock< std::mutex > lock( baseLifo::m_task_queue_mutex[1] );
        baseLifo::m_task_queue[1].emplace_back( std::move( task ) );
      }
      baseLifo::m_task_queue_not_empty_cond[1].notify_all();
    }

    m_deviceDeque.emplaceBackFromFront( baseLifo::m_hostDeque );
  }

  /// ueue of data stored on device
  FixedSizeDequeWithMutexes< T, INDEX_TYPE > m_deviceDeque;
  // Events associated to ith  copies to device buffer
  std::vector< camp::resources::Event > m_pushToDeviceEvents;
  // Events associated to ith  copies from device buffer
  std::vector< camp::resources::Event > m_popFromDeviceEvents;
};
}
#endif // LIFOSTORAGE_HPP
