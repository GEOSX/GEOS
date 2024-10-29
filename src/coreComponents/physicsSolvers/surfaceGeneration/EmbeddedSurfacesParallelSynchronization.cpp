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

/**
 * @file EmebeddedSurfacesParallelSynchronization.cpp
 */

#include "EmbeddedSurfacesParallelSynchronization.hpp"

#include "common/GeosxMacros.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/MeshFields.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"


namespace geos
{

namespace  embeddedSurfacesParallelSynchronization
{

using namespace dataRepository;

namespace parallelSynchronizationHelpers
{

void packNewNodes( NeighborCommunicator * const neighbor,
                   int commID,
                   MeshLevel & mesh,
                   NewObjectLists & newObjects )
{
  EmbeddedSurfaceNodeManager & nodeManager = mesh.getEmbSurfNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();

  int neighborRank = neighbor->neighborRank();

  localIndex_array newNodesToSend;

  arrayView1d< localIndex const > const & edgeGhostsToSend = edgeManager.getNeighborData( neighbor->neighborRank() ).ghostsToSend();
  arrayView1d< localIndex > const & parentIndex = nodeManager.getField< fields::parentEdgeIndex >();
  arrayView1d< integer const > const & nodeGhostRank = nodeManager.ghostRank();

  for( auto const ni : newObjects.newNodes )
  {
    if( nodeGhostRank[ni] == neighborRank )
    {
      // a node is sent if it is a ghost on neighborRank
      newNodesToSend.emplace_back( ni );
    }
    else
    {
      forAll< serialPolicy >( edgeGhostsToSend.size(), [=, &newNodesToSend] ( localIndex const a )
      {
        if( edgeGhostsToSend[a] == parentIndex[ni] )
        {
          // a node is sent if the edge on which it was created had to be sent.
          newNodesToSend.emplace_back( ni );
        }
      } );
    }
  }

  int bufferSize = 0;

  bufferSize += nodeManager.packNewNodesGlobalMapsSize( newNodesToSend );

  neighbor->resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor->sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = 0;

  packedSize += nodeManager.packNewNodesGlobalMaps( sendBufferPtr, newNodesToSend );

  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );
}

void unpackNewNodes( NeighborCommunicator * const neighbor,
                     int commID,
                     MeshLevel & mesh )
{

  EmbeddedSurfaceNodeManager & nodeManager = mesh.getEmbSurfNodeManager();

  buffer_type const & receiveBuffer = neighbor->receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  localIndex_array newGhostNodes;

  parallelDeviceEvents events;

  nodeManager.unpackNewNodesGlobalMaps( receiveBufferPtr, newGhostNodes );

  waitAllDeviceEvents( events );
}

void packNewObjectsToGhosts( NeighborCommunicator * const neighbor,
                             int commID,
                             MeshLevel & mesh,
                             NewObjectLists & newObjects )
{
  EmbeddedSurfaceNodeManager & nodeManager = mesh.getEmbSurfNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  int neighborRank = neighbor->neighborRank();

  localIndex_array newNodesToSend;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > newSurfaceGhostsToSend;

  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > newElemsToSend;
  array1d< array1d< localIndex_array > > newElemsToSendData;

  arrayView1d< localIndex const > const & edgeGhostsToSend = edgeManager.getNeighborData( neighbor->neighborRank() ).ghostsToSend();
  arrayView1d< localIndex > const & parentIndex = nodeManager.getField< fields::parentEdgeIndex >();

  for( auto const ni : newObjects.newNodes )
  {
    forAll< serialPolicy >( edgeGhostsToSend.size(), [=, &newNodesToSend] ( localIndex const a )
    {
      if( edgeGhostsToSend[a] == parentIndex[ni] )
      {
        // a node is sent if the edge on which it was created had to be sent.
        newNodesToSend.emplace_back( ni );
      }
    } );
  }

  newElemsToSendData.resize( elemManager.numRegions() );
  newElemsToSend.resize( elemManager.numRegions() );
  //
  if( newObjects.newElements.size() > 0 )  // this map may be completely empty on some ranks.
  {
    elemManager.forElementSubRegionsComplete< EmbeddedSurfaceSubRegion >(
      [&]( localIndex const er, localIndex const esr, ElementRegionBase &, EmbeddedSurfaceSubRegion & subRegion )
    {

      OrderedVariableToManyElementRelation const & surfaceElementsToCells = subRegion.getToCellRelation();

      for( localIndex const & k : newObjects.newElements.at( {er, esr} ) )
      {
        localIndex const elemRegionIndex  = surfaceElementsToCells.m_toElementRegion[k][0];
        localIndex const elemSubRegionIndex = surfaceElementsToCells.m_toElementSubRegion[k][0];
        localIndex const elemIndex = surfaceElementsToCells.m_toElementIndex[k][0];

        CellElementRegion const & elemRegion = elemManager.getRegion< CellElementRegion >( elemRegionIndex );
        CellElementSubRegion const & elemSubRegion = elemRegion.getSubRegion< CellElementSubRegion >( elemSubRegionIndex );

        arrayView1d< localIndex const > const & elemGhostsToSend =
          elemSubRegion.getNeighborData( neighborRank ).ghostsToSend();

        forAll< serialPolicy >( elemGhostsToSend.size(), [=, &newSurfaceGhostsToSend] ( localIndex const a )
        {
          if( elemIndex == elemGhostsToSend[a] )
          {
            newSurfaceGhostsToSend[ {er, esr} ].insert( k );
            return;
          }
        } );
      }
    } );
  }

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newElemsToSendData[er].resize( elemRegion.numSubRegions() );
    newElemsToSend[er].resize( elemRegion.numSubRegions() );
    if( newSurfaceGhostsToSend.size() > 0 )  // this map may be completely empty on some ranks.
    {
      elemRegion.forElementSubRegionsIndex< EmbeddedSurfaceSubRegion >( [&]( localIndex const esr,
                                                                             EmbeddedSurfaceSubRegion & subRegion )
      {
        localIndex_array & surfaceElemGhostsToSend = subRegion.getNeighborData( neighborRank ).ghostsToSend();
        surfaceElemGhostsToSend.move( hostMemorySpace );

        for( localIndex const & k : newSurfaceGhostsToSend.at( {er, esr} ) )
        {
          newElemsToSendData[er][esr].emplace_back( k );
          surfaceElemGhostsToSend.emplace_back( k );
        }
        newElemsToSend[er][esr] = newElemsToSendData[er][esr];
      } );
    }
  }

  parallelDeviceEvents sizeEvents;
  int bufferSize = 0;

  bufferSize += nodeManager.packGlobalMapsSize( newNodesToSend, 0 );
  bufferSize += elemManager.packGlobalMapsSize( newElemsToSend );

  bufferSize += nodeManager.packUpDownMapsSize( newNodesToSend );
  bufferSize += elemManager.packUpDownMapsSize( newElemsToSend );

  bufferSize += nodeManager.packSize( newNodesToSend, 0, false, sizeEvents );
  bufferSize += elemManager.packSize( newElemsToSend );

  neighbor->resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor->sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  parallelDeviceEvents packEvents;
  int packedSize = 0;

  packedSize += nodeManager.packGlobalMaps( sendBufferPtr, newNodesToSend, 0 );
  packedSize += elemManager.packGlobalMaps( sendBufferPtr, newElemsToSend );

  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, newNodesToSend );
  packedSize += elemManager.packUpDownMaps( sendBufferPtr, newElemsToSend );

  packedSize += nodeManager.pack( sendBufferPtr, newNodesToSend, 0, false, packEvents );
  packedSize += elemManager.pack( sendBufferPtr, newElemsToSend );

  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );
}

void unpackNewToGhosts( NeighborCommunicator * const neighbor,
                        int commID,
                        MeshLevel & mesh )
{

  EmbeddedSurfaceNodeManager & nodeManager = mesh.getEmbSurfNodeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  buffer_type const & receiveBuffer = neighbor->receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  localIndex_array newGhostNodes;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > newGhostElems;
  array1d< array1d< localIndex_array > > newGhostElemsData;
  newGhostElems.resize( elemManager.numRegions() );
  newGhostElemsData.resize( elemManager.numRegions() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newGhostElemsData[er].resize( elemRegion.numSubRegions() );
    newGhostElems[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newGhostElems[er][esr].set( newGhostElemsData[er][esr] );
    }
  }

  parallelDeviceEvents events;

  nodeManager.unpackGlobalMaps( receiveBufferPtr, newGhostNodes, 0 );
  elemManager.unpackGlobalMaps( receiveBufferPtr, newGhostElems );

  nodeManager.unpackUpDownMaps( receiveBufferPtr, newGhostNodes, true, true );
  elemManager.unpackUpDownMaps( receiveBufferPtr, newGhostElems, true );

  nodeManager.unpack( receiveBufferPtr, newGhostNodes, 0, false, events );
  elemManager.unpack( receiveBufferPtr, newGhostElems );

  waitAllDeviceEvents( events );

  elemManager.forElementSubRegionsComplete< ElementSubRegionBase >(
    [&]( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    localIndex_array & elemGhostsToReceive = subRegion.getNeighborData( neighbor->neighborRank() ).ghostsToReceive();

    if( newGhostElemsData[er][esr].size() > 0 )
    {
      forAll< serialPolicy >( newGhostElemsData[er][esr].size(), [=, &elemGhostsToReceive] ( localIndex const k )
      {
        localIndex const & newElemIndex = newGhostElemsData[er][esr][k];
        elemGhostsToReceive.emplace_back( newElemIndex );
      } );
    }
  } );
}

void packFracturedToGhosts( NeighborCommunicator * const neighbor,
                            int commID,
                            MeshLevel & mesh,
                            string const fractureRegionName )
{
  ElementRegionManager & elemManager = mesh.getElemManager();

  int neighborRank = neighbor->neighborRank();

  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > elemsToSend;
  array1d< array1d< localIndex_array > > elemsToSendData;

  elemsToSendData.resize( elemManager.numRegions() );
  elemsToSend.resize( elemManager.numRegions() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    elemsToSendData[er].resize( elemRegion.numSubRegions() );
    elemsToSend[er].resize( elemRegion.numSubRegions() );

    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                       CellElementSubRegion & subRegion )
    {
      // we send all the ghosts
      arrayView1d< localIndex const >  const & elemsGhostsToSend = subRegion.getNeighborData( neighborRank ).ghostsToSend();
      elemsGhostsToSend.move( hostMemorySpace );
      forAll< serialPolicy >( elemsGhostsToSend.size(), [=, &elemsToSendData] ( localIndex const k )
      {
        elemsToSendData[er][esr].emplace_back( elemsGhostsToSend[k] );
      } );
      elemsToSend[er][esr] = elemsToSendData[er][esr];
    } );
  }

  parallelDeviceEvents sizeEvents;
  int bufferSize = 0;

  bufferSize += elemManager.packFracturedElementsSize( elemsToSend, fractureRegionName );

  neighbor->resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor->sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  parallelDeviceEvents packEvents;
  int packedSize = 0;

  packedSize += elemManager.packFracturedElements( sendBufferPtr, elemsToSend, fractureRegionName );

  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );
}

void unpackFracturedToGhosts( NeighborCommunicator * const neighbor,
                              int commID,
                              MeshLevel & mesh,
                              string const fractureRegionName )
{

  ElementRegionManager & elemManager = mesh.getElemManager();

  buffer_type const & receiveBuffer = neighbor->receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > ghostElems;
  array1d< array1d< localIndex_array > > ghostElemsData;
  ghostElems.resize( elemManager.numRegions() );
  ghostElemsData.resize( elemManager.numRegions() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    ghostElemsData[er].resize( elemRegion.numSubRegions() );
    ghostElems[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      ghostElems[er][esr].set( ghostElemsData[er][esr] );
    }
  }

  parallelDeviceEvents events;

  elemManager.unpackFracturedElements( receiveBufferPtr, ghostElems, fractureRegionName );

  waitAllDeviceEvents( events );
}

void synchronizeNewNodes( MeshLevel & mesh,
                          std::vector< NeighborCommunicator > & neighbors,
                          NewObjectLists & newObjects,
                          int const mpiCommOrder )
{
  //************************************************************************************************
  // We need to send over the new embedded surfaces and related objects for those whose parents are ghosts on neighbors.

  MPI_iCommData commData( CommunicationTools::getInstance().getCommID() );
  commData.resize( neighbors.size());
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewNodes( &neighbor,
                  commData.commID(),
                  mesh,
                  newObjects );

    neighbor.mpiISendReceiveBufferSizes( commData.commID(),
                                         commData.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOSX );

  }

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData.commID(),
                                     commData.mpiSendBufferRequest( neighborIndex ),
                                     commData.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOSX );
  }


  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;
    if( mpiCommOrder == 0 )
    {
      MpiWrapper::waitAny( commData.size(),
                           commData.mpiRecvBufferRequest(),
                           &neighborIndex,
                           commData.mpiRecvBufferStatus() );
    }
    else
    {
      MpiWrapper::wait( commData.mpiRecvBufferRequest() + count,
                        commData.mpiRecvBufferStatus() + count );
    }

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    unpackNewNodes( &neighbor, commData.commID(), mesh );
  }

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferSizeRequest(),
                       commData.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferRequest(),
                       commData.mpiSendBufferSizeStatus() );

}

void synchronizeNewSurfaces( MeshLevel & mesh,
                             std::vector< NeighborCommunicator > & neighbors,
                             NewObjectLists & newObjects,
                             int const mpiCommOrder )
{
  //************************************************************************************************
  // We need to send over the new embedded surfaces and related objects for those whose parents are ghosts on neighbors.

  MPI_iCommData commData( CommunicationTools::getInstance().getCommID() );
  commData.resize( neighbors.size());
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewObjectsToGhosts( &neighbor,
                            commData.commID(),
                            mesh,
                            newObjects );

    neighbor.mpiISendReceiveBufferSizes( commData.commID(),
                                         commData.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOSX );

  }

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData.commID(),
                                     commData.mpiSendBufferRequest( neighborIndex ),
                                     commData.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOSX );
  }


  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;
    if( mpiCommOrder == 0 )
    {
      MpiWrapper::waitAny( commData.size(),
                           commData.mpiRecvBufferRequest(),
                           &neighborIndex,
                           commData.mpiRecvBufferStatus() );
    }
    else
    {
      MpiWrapper::wait( commData.mpiRecvBufferRequest() + count,
                        commData.mpiRecvBufferStatus() + count );
    }

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    unpackNewToGhosts( &neighbor, commData.commID(), mesh );
  }

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferSizeRequest(),
                       commData.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferRequest(),
                       commData.mpiSendBufferSizeStatus() );
}

void synchronizeFracturedElements( MeshLevel & mesh,
                                   std::vector< NeighborCommunicator > & neighbors,
                                   string const fractureRegionName )
{
  MPI_iCommData commDataJunk( CommunicationTools::getInstance().getCommID() );
  MPI_iCommData commData( CommunicationTools::getInstance().getCommID() );
  commData.resize( neighbors.size());
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packFracturedToGhosts( &neighbor,
                           commData.commID(),
                           mesh,
                           fractureRegionName );

    neighbor.mpiISendReceiveBufferSizes( commData.commID(),
                                         commData.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOSX );

  }

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData.commID(),
                                     commData.mpiSendBufferRequest( neighborIndex ),
                                     commData.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOSX );
  }


  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;

    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    unpackFracturedToGhosts( &neighbor, commData.commID(), mesh, fractureRegionName );
  }

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferSizeRequest(),
                       commData.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferRequest(),
                       commData.mpiSendBufferSizeStatus() );
}



}  /* parallelSynchronizationHelpers */

using namespace parallelSynchronizationHelpers;

void sychronizeTopology( MeshLevel & mesh,
                         std::vector< NeighborCommunicator > & neighbors,
                         NewObjectLists & newObjects,
                         int const mpiCommOrder,
                         string const fractureRegionName )
{

  // Synchronize nodes
  synchronizeNewNodes( mesh,
                       neighbors,
                       newObjects,
                       mpiCommOrder );


  // Synchronize embedded Surfaces
  synchronizeNewSurfaces( mesh,
                          neighbors,
                          newObjects,
                          mpiCommOrder );

  synchronizeFracturedElements( mesh,
                                neighbors,
                                fractureRegionName );

}

} /* embeddedSurfacesParallelSynchronization */

} /* namespace geos */
