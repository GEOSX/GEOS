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

#ifndef LIFOSTORAGECOMMON_HPP
#define LIFOSTORAGECOMMON_HPP

#include <deque>
#include <future>
#include <mutex>
#include <condition_variable>
#include <camp/camp.hpp>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>

#ifdef LIFO_DISABLE_CALIPER
#define LIFO_MARK_FUNCTION
#define LIFO_MARK_SCOPE( a )
#define LIFO_LOG_RANK( a ) std::cerr << a << std::endl;
#else
#define LIFO_MARK_FUNCTION GEOS_MARK_FUNCTION
#define LIFO_MARK_SCOPE( a ) GEOS_MARK_SCOPE( a )
#define LIFO_LOG_RANK( a ) GEOS_LOG_RANK( a )
#endif

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"
#include "common/FixedSizeDequeWithMutexes.hpp"
#include "common/MultiMutexesLock.hpp"


namespace geos
{
/**
 * This class is used to store in a LIFO way buffers, first on device, then on host, then on disk.
 */
template< typename T, typename INDEX_TYPE >
class LifoStorageCommon
{
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
  LifoStorageCommon( std::string name, size_t elemCnt, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    m_maxNumberOfBuffers( maxNumberOfBuffers ),
    m_bufferSize( elemCnt*sizeof( T ) ),
    m_name( name ),
    m_hostDeque( numberOfBuffersToStoreOnHost, elemCnt, LvArray::MemorySpace::host ),
    m_bufferCount( 0 ), m_bufferToHostCount( 0 ), m_bufferToDiskCount( 0 ),
    m_continue( true ),
    m_hasPoppedBefore( false )
  {
    m_worker[0] = std::thread( &LifoStorageCommon< T, INDEX_TYPE >::wait_and_consume_tasks, this, 0 );
    m_worker[1] = std::thread( &LifoStorageCommon< T, INDEX_TYPE >::wait_and_consume_tasks, this, 1 );
  }

  virtual ~LifoStorageCommon()
  {
    m_continue = false;
    m_task_queue_not_empty_cond[0].notify_all();
    m_task_queue_not_empty_cond[1].notify_all();
    m_worker[0].join();
    m_worker[1].join();
  }

  /**
   * Asynchroneously push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  virtual void pushAsync( arrayView1d< T > array ) = 0;

  /**
   * Waits for last push to be terminated
   */
  virtual void pushWait() = 0;

  /**
   * Asynchroneously copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  virtual void popAsync( arrayView1d< T > array ) = 0;

  /**
   * Prelude for pop async
   */
  void popAsyncPrelude( )
  {
    LIFO_MARK_FUNCTION;
    //wait the last push to avoid race condition
    pushWait();
    if( m_hasPoppedBefore )
    {
      // Ensure last pop is finished
      popWait();
    }
    else
    {
      if( m_maxNumberOfBuffers != m_bufferCount )
        LIFO_LOG_RANK( " LIFO : warning number of entered buffered (" << m_bufferCount
                                                                      << ") != max LIFO size (" << m_maxNumberOfBuffers << ") !" );
      // Ensure that all push step are ended
      for( int queueId = 0; queueId < 2; queueId++ )
      {
        std::unique_lock< std::mutex > lock( m_task_queue_mutex[queueId] );
        m_task_queue_not_empty_cond[queueId].wait( lock, [ this, &queueId ] { return m_task_queue[queueId].empty(); } );
      }
    }
    m_hasPoppedBefore = true;
  }

  /**
   * Waits for last pop to be terminated
   */
  virtual void popWait() = 0;

  /**
   * Check if the LIFO is empty
   *
   * @return true if the LIFO does not contain a buffer.
   */
  bool empty()
  {
    return m_bufferCount == 0;
  }

  /**
   * Compute the number of arrays that can be stored on host
   * @param percent Percentage of the remaining device memory that can be dedicated to the LIFO storage.
   * @param bufferSize Size of one buffer
   * @param maxNumberOfBuffers Maximum number of buffers to store in the LIFO storage
   * @param numberOfBuffersToStoreOnDevice The number of buffer that will be stored on device by the LIFO.
   * @return The maximum number of buffer to allocate to fit in the percentage of the available memory.
   */
  static int computeNumberOfBufferOnHost( int percent, size_t bufferSize, int maxNumberOfBuffers, int numberOfBuffersToStoreOnDevice )
  {
    GEOS_ERROR_IF( percent > 100, "Error, percentage of memory should be smallerer than -100, check lifoOnHost (should be greater that -100)" );
#if defined( _SC_AVPHYS_PAGES ) && defined( _SC_PAGESIZE )
    size_t const free = sysconf( _SC_AVPHYS_PAGES ) * sysconf( _SC_PAGESIZE );
#else
    size_t const free = 0;
    GEOS_ERROR( "To use LifoStorage, both _SC_AVPHYS_PAGES and _SC_PAGESIZE must be defined." );
#endif
    int numberOfBuffersToStoreOnHost = std::max( 1, std::min( ( int )( 0.01 * percent * free / bufferSize ), maxNumberOfBuffers - numberOfBuffersToStoreOnDevice ) );
    double freeGB = ( ( double ) free ) / ( 1024.0 * 1024.0 * 1024.0 ) / MpiWrapper::nodeCommSize();
    LIFO_LOG_RANK( " LIFO : available memory on host " << freeGB << " GB" );
    return numberOfBuffersToStoreOnHost;
  }

protected:

  /// number of buffers to be inserted into the LIFO
  int m_maxNumberOfBuffers;
  /// size of one buffer in bytes
  size_t m_bufferSize;
  /// name used to store data on disk
  std::string m_name;
  /// Queue of data stored on host memory
  FixedSizeDequeWithMutexes< T, INDEX_TYPE > m_hostDeque;

  /// counter of buffer stored in LIFO
  int m_bufferCount;
  /// counter of buffer pushed to host
  int m_bufferToHostCount;
  /// counter of buffer pushed to disk
  int m_bufferToDiskCount;


  /// condition used to tell m_worker queue has been filled or processed is stopped.
  std::condition_variable m_task_queue_not_empty_cond[2];
  /// mutex to protect access to m_task_queue.
  std::mutex m_task_queue_mutex[2];
  /// queue of task to be executed by m_worker.
  std::deque< std::packaged_task< void() > > m_task_queue[2];
  /// thread to execute tasks.
  std::thread m_worker[2];
  /// boolean to keep m_worker alive.
  bool m_continue;
  /// marker to detect first pop
  bool m_hasPoppedBefore;

  /**
   * Copy data from host memory to disk
   *
   * @param id ID of the buffer to store on disk.
   */
  void hostToDisk( int id )
  {
    LIFO_MARK_FUNCTION;
    {
      auto lock = make_multilock( m_hostDeque.m_popMutex, m_hostDeque.m_backMutex );
      writeOnDisk( m_hostDeque.back().dataIfContiguous(), id );
      m_hostDeque.pop_back();
    }
    m_hostDeque.m_notFullCond.notify_all();
  }

  /**
   * Copy data from disk to host memory
   *
   * @param id ID of the buffer to read on disk.
   */
  void diskToHost( int id )
  {
    LIFO_MARK_FUNCTION;
    {
      auto lock = make_multilock( m_hostDeque.m_emplaceMutex, m_hostDeque.m_backMutex );
      m_hostDeque.m_notFullCond.wait( lock, [ this ]  { return !( m_hostDeque.full() ); } );
      readOnDisk( const_cast< T * >(m_hostDeque.next_back().dataIfContiguous()), id );
      m_hostDeque.inc_back();
    }
    m_hostDeque.m_notEmptyCond.notify_all();
  }
  /**
   * Checks if a directory exists.
   *
   * @param dirName Directory name to check existence of.
   * @return true is dirName exists and is a directory.
   */
  bool dirExists( const std::string & dirName )
  {
    struct stat buffer;
    return stat( dirName.c_str(), &buffer ) == 0;
  }

  /**
   * Write data on disk
   *
   * @param d Data to store on disk.
   * @param id ID of the buffer to read on disk
   */
  void writeOnDisk( const T * d, int id )
  {
    LIFO_MARK_FUNCTION;
    std::string fileName = GEOS_FMT( "{}_{:08}.dat", m_name, id );
    int lastDirSeparator = fileName.find_last_of( "/\\" );
    std::string dirName = fileName.substr( 0, lastDirSeparator );
    if( string::npos != (size_t)lastDirSeparator && !dirExists( dirName ))
      makeDirsForPath( dirName );

    std::ofstream wf( fileName, std::ios::out | std::ios::binary );
    GEOS_ERROR_IF( !wf || wf.fail() || !wf.is_open(),
                   "Could not open file "<< fileName << " for writting" );
    wf.write( (const char *)d, m_bufferSize );
    GEOS_ERROR_IF( wf.bad() || wf.fail(),
                   "An error occured while writting "<< fileName );
    wf.close();
  }

  /**
   * Read data from disk
   *
   * @param d  Buffer to store data read from disk.
   * @param id ID of the buffer on disk.
   */
  void readOnDisk( T * d, int id )
  {
    LIFO_MARK_FUNCTION;
    std::string fileName = GEOS_FMT( "{}_{:08}.dat", m_name, id );
    std::ifstream wf( fileName, std::ios::in | std::ios::binary );
    GEOS_ERROR_IF( !wf,
                   "Could not open file "<< fileName << " for reading" );
    wf.read( (char *)d, m_bufferSize );
    wf.close();
    remove( fileName.c_str() );
  }


private:
  /**
   * Execute the tasks pushed in m_task_queue in the given order.
   *
   * @param queueId index of the queue (0: queue to handle device/host transfers; 1: host/disk transfers)
   */
  void wait_and_consume_tasks( int queueId )
  {
    LIFO_MARK_FUNCTION;
    while( m_continue )
    {
      std::unique_lock< std::mutex > lock( m_task_queue_mutex[queueId] );
      {
        LIFO_MARK_SCOPE( waitForTask );
        m_task_queue_not_empty_cond[queueId].wait( lock, [ this, &queueId ] { return !( m_task_queue[queueId].empty()  && m_continue ); } );
      }
      if( m_continue == false ) break;
      std::packaged_task< void() > task( std::move( m_task_queue[queueId].front() ) );
      m_task_queue[queueId].pop_front();
      lock.unlock();
      m_task_queue_not_empty_cond[queueId].notify_all();
      {
        LIFO_MARK_SCOPE( runningTask );
        task();
      }
    }
  }
};
}
#endif // LIFOSTORAGECOMMON_HPP
