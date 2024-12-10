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

#include "MPI_iCommData.hpp"
#include "CommunicationTools.hpp"

namespace geos
{

MPI_iCommData::MPI_iCommData():
  m_size( 0 ),
  m_commID( CommunicationTools::getInstance().getCommID() ),
  m_mpiSendBufferRequest(),
  m_mpiRecvBufferRequest(),
  m_mpiSendBufferStatus(),
  m_mpiRecvBufferStatus(),
  m_mpiSendBufferSizeRequest(),
  m_mpiRecvBufferSizeRequest(),
  m_mpiSendBufferSizeStatus(),
  m_mpiRecvBufferSizeStatus()
{}


MPI_iCommData::~MPI_iCommData()
{
  for( int neighbor=0; neighbor<m_size; ++neighbor )
  {
    GEOS_ERROR_IF_NE_MSG( m_mpiSendBufferRequest[neighbor], MPI_REQUEST_NULL,
                          "Destroying MPI_iCommData with uncompleted m_mpiSendBufferRequest for neighbor "<<neighbor );
    GEOS_ERROR_IF_NE_MSG( m_mpiRecvBufferRequest[neighbor], MPI_REQUEST_NULL,
                          "Destroying MPI_iCommData with uncompleted m_mpiRecvBufferRequest for neighbor "<<neighbor );
    GEOS_ERROR_IF_NE_MSG( m_mpiSendBufferSizeRequest[neighbor], MPI_REQUEST_NULL,
                          "Destroying MPI_iCommData with uncompleted m_mpiSendBufferSizeRequest for neighbor "<<neighbor );
    GEOS_ERROR_IF_NE_MSG( m_mpiRecvBufferSizeRequest[neighbor], MPI_REQUEST_NULL,
                          "Destroying MPI_iCommData with uncompleted m_mpiRecvBufferSizeRequest for neighbor "<<neighbor );
  }
}


void MPI_iCommData::resize( localIndex numMessages )
{
  for( int neighbor=0; neighbor<m_size; ++neighbor )
  {
    GEOS_ERROR_IF_NE_MSG( m_mpiSendBufferRequest[neighbor], MPI_REQUEST_NULL,
                          "resize(localIndex) called on MPI_iCommData with uncompleted m_mpiSendBufferRequest for neighbor "<<neighbor );
    GEOS_ERROR_IF_NE_MSG( m_mpiRecvBufferRequest[neighbor], MPI_REQUEST_NULL,
                          "resize(localIndex) called on MPI_iCommData with uncompleted m_mpiRecvBufferRequest for neighbor "<<neighbor );
    GEOS_ERROR_IF_NE_MSG( m_mpiSendBufferSizeRequest[neighbor], MPI_REQUEST_NULL,
                          "resize(localIndex) called on MPI_iCommData with uncompleted m_mpiSendBufferSizeRequest for neighbor "<<neighbor );
    GEOS_ERROR_IF_NE_MSG( m_mpiRecvBufferSizeRequest[neighbor], MPI_REQUEST_NULL,
                          "resize(localIndex) called on MPI_iCommData with uncompleted m_mpiRecvBufferSizeRequest for neighbor "<<neighbor );
  }

  m_mpiSendBufferRequest.resize( numMessages );
  m_mpiRecvBufferRequest.resize( numMessages );
  m_mpiSendBufferStatus.resize( numMessages );
  m_mpiRecvBufferStatus.resize( numMessages );
  m_mpiSendBufferSizeRequest.resize( numMessages );
  m_mpiRecvBufferSizeRequest.resize( numMessages );
  m_mpiSendBufferSizeStatus.resize( numMessages );
  m_mpiRecvBufferSizeStatus.resize( numMessages );
  m_size = static_cast< int >(numMessages);

  for( int neighbor=0; neighbor<numMessages; ++neighbor )
  {
    m_mpiSendBufferRequest[neighbor] = MPI_REQUEST_NULL;
    m_mpiRecvBufferRequest[neighbor] = MPI_REQUEST_NULL;
    m_mpiSendBufferSizeRequest[neighbor] = MPI_REQUEST_NULL;
    m_mpiRecvBufferSizeRequest[neighbor] = MPI_REQUEST_NULL;
  }

}

} /* namespace geos */
