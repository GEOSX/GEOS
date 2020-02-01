/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file NeighborCommunicator.hpp
 */

#ifndef GEOSX_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_
#define GEOSX_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_

#include "MpiWrapper.hpp"

#include "common/DataTypes.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "cxx-utilities/src/IntegerConversion.hpp"

namespace geosx
{
inline int CommTag( int const GEOSX_UNUSED_ARG( senderRank ),
                    int const GEOSX_UNUSED_ARG( receiverRank ),
                    int const comm )
{
//  int m_size;
//  MPI_Comm_size( MPI_COMM_GEOSX, &m_size );
//  return senderRank * m_size + receiverRank + m_size * m_size * comm;
  return comm;
}

class MeshLevel;
class ObjectManagerBase;
class NeighborCommunicator
{
public:

  NeighborCommunicator();

  void MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                         int const sendSize,
                         MPI_Request & sendRequest,
                         buffer_unit_type * const receiveBuffer,
                         int const receiveSize,
                         MPI_Request & receiveRequest,
                         int const commID,
                         MPI_Comm mpiComm );


  void MPI_iSendReceiveBufferSizes( int const commID,
                                    MPI_Comm mpiComm );

  void MPI_iSendReceiveBufferSizes( int const commID,
                                    MPI_Request & mpiSendRequest,
                                    MPI_Request & mpiRecvRequest,
                                    MPI_Comm mpiComm );

  void MPI_iSendReceiveBuffers( int const commID,
                                MPI_Comm mpiComm );

  void MPI_iSendReceiveBuffers( int const commID,
                                MPI_Request & mpiSendRequest,
                                MPI_Request & mpiRecvRequest,
                                MPI_Comm mpiComm );

  void MPI_iSendReceive( int const commID,
                         MPI_Comm mpiComm );

  void MPI_iSendReceive( int const commID,
                         MPI_Request & mpiSendRequest,
                         MPI_Request & mpiRecvRequest,
                         MPI_Comm mpiComm );

  template< typename T >
  void MPI_iSendReceive( T const * const sendBuffer,
                         int const sendSize,
                         MPI_Request & sendReq,
                         T * const recvBuffer,
                         int const recvSize,
                         MPI_Request & recvReq,
                         int const commID,
                         MPI_Comm mpiComm )
  {
    MPI_iSendReceive( reinterpret_cast< buffer_unit_type const * >( sendBuffer ),
                      sendSize * sizeof(T),
                      sendReq,
                      reinterpret_cast< buffer_unit_type * >( recvBuffer ),
                      recvSize * sizeof(T),
                      recvReq,
                      commID,
                      mpiComm );
  }


  template< typename T >
  void MPI_iSendReceive( array1d< T > const & sendBuffer,
                         MPI_Request & sendReq,
                         array1d< T > & recvBuffer,
                         MPI_Request & recvReq,
                         int const commID,
                         MPI_Comm mpiComm );


  template< typename T >
  void MPI_iSendReceive( array1d< T > const & sendBuffer,
                         array1d< T > & recvBuffer,
                         int const commID,
                         MPI_Comm mpiComm )
  {
    MPI_iSendReceive( sendBuffer,
                      m_mpiSendBufferRequest[commID],
                      recvBuffer,
                      m_mpiRecvBufferRequest[commID],
                      commID,
                      mpiComm );
  }

  template< typename T >
  void MPI_iSendReceiveWait( array1d< T > const & sendBuffer,
                             array1d< T > & recvBuffer,
                             int const commID,
                             MPI_Comm mpiComm )
  {
    MPI_iSendReceive( sendBuffer, recvBuffer, commID, mpiComm );
    MPI_WaitAll( commID );
  }


  void MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                         int const sendSize,
                         int const commID,
                         MPI_Comm mpiComm );

  void MPI_WaitAll( int const commID,
                    MPI_Request & mpiSendRequest,
                    MPI_Status & mpiSendStatus,
                    MPI_Request & mpiRecvRequest,
                    MPI_Status & mpiReceiveStatus );

  void MPI_WaitAll( int const commID );

  int PostSizeRecv( int const commID );

  MPI_Request GetSizeRecvRequest( int const commID );

  int PostSizeSend( int const commID );

  int PostRecv( int const commID );

  MPI_Request GetRecvRequest( int const commID );

  int PostSend( int const commID );

  /**
   * Posts non-blocking sends to m_neighborRank for
   *  both the size and regular communication buffers
   *  to exchange ghost information.
   *  Additionally posts a non-blocking recv for size
   *  information from m_neighborRank, this size recv
   *  must be completed before PostRecv is called in order
   *  to correctly resize the receive buffer.
   */
  void PrepareAndSendGhosts( bool const contactActive,
                             int const depth,
                             MeshLevel * const meshLevel,
                             int const commID );

  /**
   * Unpack the receive buffer and process ghosting
   *  information recieved from m_neighborRank.
   *  This must be called after PostRecv is called, and
   *  the request associated with that recv has
   *  completed (retrieve the request using GetRecvRequest)
   */
  void UnpackGhosts( MeshLevel * const meshLevel,
                     int const commID );

  /**
   * Posts non-blocking sends to m_neighborRank for
   *  both the size and regular communication buffers
   *  to exchange synchronization lists.
   *  Additionally posts a non-blocking recv for size
   *  information from m_neighborRank, this size recv
   *  must be completed before PostRecv is called in order
   *  to correctly resize the receive buffer.
   */
  void PrepareAndSendSyncLists( MeshLevel * const meshLevel,
                                int const commID );

  /**
   * Unpack the receive buffer and process synchronization
   *  list information recieved from m_neighborRank.
   *  This must be called after PostRecv is called, and
   *  the request associated with that recv has
   *  completed (retrieve the request using GetRecvRequest)
   */
  void UnpackAndRebuildSyncLists( MeshLevel * const meshLevel,
                                  int const CommID );

  void PackCommBufferForSync( std::map< string, string_array > const & fieldNames,
                              MeshLevel * const meshLevel,
                              int const commID,
                              bool on_device = false );

  int PackCommSizeForSync( std::map< string, string_array > const & fieldNames,
                           MeshLevel * const meshLevel,
                           int const commID,
                           bool on_device = false );

  void SendRecvBuffers( int const commID );

  void UnpackBufferForSync( std::map< string, string_array > const & fieldNames,
                            MeshLevel * const meshLevel,
                            int const commID,
                            bool on_device = false );

  static int Rank();
  static int MPISize();

  void SetNeighborRank( int const rank ) { m_neighborRank = rank; }
  int NeighborRank() const { return m_neighborRank; }

  void Clear();

  static int constexpr maxComm = 100;

  buffer_type const & ReceiveBuffer( int commID ) const
  {
    return m_receiveBuffer[commID];
  }
  buffer_type & ReceiveBuffer( int commID )
  {
    return m_receiveBuffer[commID];
  }

  int const & ReceiveBufferSize( int commID ) const
  {
    return m_receiveBufferSize[commID];
  }
  int & ReceiveBufferSize( int commID )
  {
    return m_receiveBufferSize[commID];
  }


  buffer_type const & SendBuffer( int commID ) const
  {
    return m_sendBuffer[commID];
  }
  buffer_type & SendBuffer( int commID )
  {
    return m_sendBuffer[commID];
  }

  void resizeSendBuffer( int const commID, int const newSize )
  {
    m_sendBufferSize[commID] = newSize;
    m_sendBuffer[commID].resize( newSize );
  }

  void resizeRecvBuffer( int const commID, int const newSize )
  {
    m_receiveBufferSize[commID] = newSize;
    m_receiveBuffer[commID].resize( newSize );
  }

  void AddNeighborGroupToMesh( MeshLevel * const mesh ) const;

private:
  int m_neighborRank;

  int m_sendBufferSize[maxComm];
  int m_receiveBufferSize[maxComm];

  std::vector< buffer_type > m_sendBuffer;
  std::vector< buffer_type > m_receiveBuffer;

  MPI_Request m_mpiSendBufferRequest[maxComm];
  MPI_Request m_mpiRecvBufferRequest[maxComm];

  MPI_Request m_mpiSendSizeRequest[maxComm];
  MPI_Request m_mpiRecvSizeRequest[maxComm];

  MPI_Status m_mpiSendBufferStatus[maxComm];
  MPI_Status m_mpiRecvBufferStatus[maxComm];
};



template< typename T >
void NeighborCommunicator::MPI_iSendReceive( array1d< T > const & sendBuffer,
                                             MPI_Request & sendReq,
                                             array1d< T > & recvBuffer,
                                             MPI_Request & recvReq,
                                             int const commID,
                                             MPI_Comm mpiComm )
{
  m_sendBufferSize[commID] = integer_conversion< int >( sendBuffer.size());

  MPI_iSendReceive( &m_sendBufferSize[commID],
                    1,
                    sendReq,
                    &m_receiveBufferSize[commID],
                    1,
                    recvReq,
                    commID,
                    mpiComm );

  MpiWrapper::Wait( &( recvReq ), &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Wait( &( sendReq ), &( m_mpiSendBufferStatus[commID] ) );

  recvBuffer.resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( sendBuffer.data(),
                    m_sendBufferSize[commID],
                    sendReq,
                    recvBuffer.data(),
                    m_receiveBufferSize[commID],
                    recvReq,
                    commID,
                    mpiComm );
}


} /* namespace geosx */

#endif /* GEOSX_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_ */
