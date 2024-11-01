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

#ifndef LIFOSTORAGE_HPP
#define LIFOSTORAGE_HPP

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
#include "common/LifoStorageHost.hpp"
#ifdef GEOS_USE_CUDA
#include "common/LifoStorageCuda.hpp"
#endif

namespace geos
{
/**
 * This class is used to store in a LIFO way buffers, first on device, then on host, then on disk.
 */
template< typename T, typename INDEX_TYPE >
class LifoStorage
{

public:


  /**
   * A LIFO storage will store numberOfBuffersToStoreDevice buffer on
   * deevice, numberOfBuffersToStoreHost on host and the rest on disk.
   *
   * @param name                           Prefix of the files used to save the occurenncy of the saved buffer on disk.
   * @param elemCnt                        Number of elments in the LvArray we want to store in the LIFO storage.
   * @param numberOfBuffersToStoreOnDevice Maximum number of array to store on device memory. If negative opposite of the percent of left
   * memory we want to use( -80 = use 80% of remaining memory ).
   * @param numberOfBuffersToStoreOnHost   Maximum number of array to store on host memory . If negative opposite of the percent of left
   * memory we want to use( -80 = use 80% of remaining memory ).
   * @param maxNumberOfBuffers             Number of arrays expected to be stores in the LIFO.
   */
  LifoStorage( std::string name, size_t elemCnt, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    m_maxNumberOfBuffers( maxNumberOfBuffers ),
    m_bufferSize( elemCnt*sizeof( T ) ),
    m_bufferCount( 0 )
  {
    LIFO_LOG_RANK( " LIFO : maximum size "<< m_maxNumberOfBuffers << " buffers " );
    LIFO_LOG_RANK( " LIFO : buffer size " << m_bufferSize / ( 1024.0 * 1024.0 ) << "MB" );
    if( numberOfBuffersToStoreOnDevice < 0 )
    {
#ifdef GEOS_USE_CUDA
      numberOfBuffersToStoreOnDevice = LifoStorageCuda< T, INDEX_TYPE >::computeNumberOfBufferOnDevice( -numberOfBuffersToStoreOnDevice, m_bufferSize, m_maxNumberOfBuffers );
#else
      numberOfBuffersToStoreOnDevice = 0;
#endif
    }
    if( numberOfBuffersToStoreOnHost < 0 )
    {
      numberOfBuffersToStoreOnHost =
        LifoStorageCommon< T, INDEX_TYPE >::computeNumberOfBufferOnHost( -numberOfBuffersToStoreOnHost, m_bufferSize, m_maxNumberOfBuffers, numberOfBuffersToStoreOnDevice );
    }
    LIFO_LOG_RANK( " LIFO : allocating "<< numberOfBuffersToStoreOnHost <<" buffers on host" );
    LIFO_LOG_RANK( " LIFO : allocating "<< numberOfBuffersToStoreOnDevice <<" buffers on device" );
#ifdef GEOS_USE_CUDA
    if( numberOfBuffersToStoreOnDevice > 0 )
    {
      m_lifo = std::make_unique< LifoStorageCuda< T, INDEX_TYPE > >( name, elemCnt, numberOfBuffersToStoreOnDevice, numberOfBuffersToStoreOnHost, maxNumberOfBuffers );
    }
    else
#endif
    {
      m_lifo = std::make_unique< LifoStorageHost< T, INDEX_TYPE > >( name, elemCnt, numberOfBuffersToStoreOnHost, maxNumberOfBuffers );
    }

  }

  /**
   * Build a LIFO storage for a given LvArray array.
   *
   * @param name                           Prefix of the files used to save the occurenncy of the saved buffer on disk.
   * @param array                          The LvArray that will be store in the LIFO.
   * @param numberOfBuffersToStoreOnDevice Maximum number of array to store on device memory.
   * @param numberOfBuffersToStoreOnHost   Maximum number of array to store on host memory.
   * @param maxNumberOfBuffers             Number of arrays expected to be stores in the LIFO.
   */
  LifoStorage( std::string name, arrayView1d< T > array, int numberOfBuffersToStoreOnDevice, int numberOfBuffersToStoreOnHost, int maxNumberOfBuffers ):
    LifoStorage( name, array.size(), numberOfBuffersToStoreOnDevice, numberOfBuffersToStoreOnHost, maxNumberOfBuffers ) {}

  /**
   * Asynchroneously push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void pushAsync( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    m_lifo->pushAsync( array );
  }

  /**
   * Waits for last push to be terminated
   */
  void pushWait()
  {
    LIFO_MARK_FUNCTION;
    m_lifo->pushWait();
  }

  /**
   * Push a copy of the given LvArray into the LIFO
   *
   * @param array The LvArray to store in the LIFO, should match the size of the data used in constructor.
   */
  void push( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    pushAsync( array );
    pushWait();
  }

  /**
   * Asynchroneously copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void popAsync( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    m_lifo->popAsyncPrelude();
    m_lifo->popAsync( array );
  }

  /**
   * Waits for last pop to be terminated
   */
  void popWait()
  {
    LIFO_MARK_FUNCTION;
    m_lifo->popWait();
  }

  /**
   * Copy last data from the LIFO into the LvArray.
   *
   * @param array LvArray to store data from the LIFO into it.
   */
  void pop( arrayView1d< T > array )
  {
    LIFO_MARK_FUNCTION;
    popAsync( array );
    popWait();
  }

  /**
   * Check if the LIFO is empty
   *
   * @return true if the LIFO does not contain a buffer.
   */
  bool empty()
  {
    return m_lifo->empty();
  }

private:
  /// number of buffers to be inserted into the LIFO
  int m_maxNumberOfBuffers;
  /// size of one buffer in bytes
  size_t m_bufferSize;
  /// counter of buffer stored in LIFO
  int m_bufferCount;

  /// pointer either to CUDA aware or host only LifoStorage
  std::unique_ptr< LifoStorageCommon< T, INDEX_TYPE > > m_lifo;

};
}
#endif // LIFOSTORAGE_HPP
