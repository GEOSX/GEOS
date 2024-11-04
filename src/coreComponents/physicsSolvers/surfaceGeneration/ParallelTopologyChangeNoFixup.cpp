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

/**
 * @file ParallelTopologyChange.cpp
 */

#include "ParallelTopologyChange.hpp"

#include "common/GeosxMacros.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/MeshFields.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"

#if PARALLEL_TOPOLOGY_CHANGE_METHOD==1
namespace geos
{

using namespace dataRepository;

namespace parallelTopologyChange
{


template< typename T >
void filterNonOwnedFromContainer( array1d< localIndex > & newList,
                                  T const & container,
                                  arrayView1d< localIndex const > const & ghostRank,
                                  integer const neighborRank )
{
  newList.resize( container.size());
  {
    localIndex a=0;
    for( auto index : container )
    {
      if( ghostRank[index] == neighborRank )
      {
        newList[a] = index;
        ++a;
      }
    }
    newList.resize( a );
  }
}

template< typename T >
void filterNonOwnedFromContainer( array1d< localIndex > & newList,
                                  T const & container,
                                  arrayView1d< localIndex const > const & parentIndices,
                                  arrayView1d< localIndex const > const & ghostRank,
                                  integer const neighborRank )
{
  newList.resize( container.size());
  {
    localIndex a=0;
    for( auto index : container )
    {
      localIndex const parentIndex = ObjectManagerBase::getParentRecursive( parentIndices, index );
      if( ghostRank[parentIndex] == neighborRank )
      {
        newList[a] = index;
        ++a;
      }
    }
    newList.resize( a );
  }
}

void filterNewObjectsForPackToGhosts( std::set< localIndex > const & objectList,
                                      arrayView1d< localIndex > const & parentIndices,
                                      localIndex_array & ghostsToSend,
                                      localIndex_array & objectsToSend )
{

  ghostsToSend.move( hostMemorySpace );
  //TODO this needs to be inverted since the ghostToSend list should be much longer....
  // and the objectList is a searchable set.
  for( auto const index : objectList )
  {
    localIndex const parentIndex = parentIndices[index];
    for( localIndex a=0; a<ghostsToSend.size(); ++a )
    {
      if( ghostsToSend[a]==parentIndex )
      {
        objectsToSend.emplace_back( index );
        ghostsToSend.emplace_back( index );
        break;
      }
    }
  }
}

void filterModObjectsForPackToGhosts( std::set< localIndex > const & objectList,
                                      localIndex_array const & ghostsToSend,
                                      localIndex_array & objectsToSend )
{
  ghostsToSend.move( hostMemorySpace );
  for( localIndex a=0; a<ghostsToSend.size(); ++a )
  {
    if( objectList.count( ghostsToSend[a] ) > 0 )
    {
      objectsToSend.emplace_back( ghostsToSend[a] );
    }
  }
}



//***** 1A *****//
void packNewAndModifiedObjectsToOwningRanks( NeighborCommunicator & neighbor,
                                             MeshLevel * const meshLevel,
                                             ModifiedObjectLists const & modifiedObjects,
                                             int const commID )
{
  int bufferSize = 0;

  NodeManager & nodeManager = meshLevel->getNodeManager();
  EdgeManager & edgeManager = meshLevel->getEdgeManager();
  FaceManager & faceManager = meshLevel->getFaceManager();
  ElementRegionManager & elemManager = meshLevel->getElemManager();

  arrayView1d< integer const > const & nodeGhostRank = nodeManager.ghostRank();
  arrayView1d< integer const > const & edgeGhostRank = edgeManager.ghostRank();
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

  arrayView1d< localIndex const > const & parentNodeIndices = nodeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex const > const & parentEdgeIndices = edgeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex const > const & parentFaceIndices = faceManager.getField< fields::parentIndex >();

  int const neighborRank = neighbor.neighborRank();

  array1d< localIndex > newNodePackListArray; filterNonOwnedFromContainer( newNodePackListArray, modifiedObjects.newNodes, parentNodeIndices, nodeGhostRank, neighborRank );
  array1d< localIndex > modNodePackListArray; filterNonOwnedFromContainer( modNodePackListArray, modifiedObjects.modifiedNodes, parentNodeIndices, nodeGhostRank, neighborRank );
  array1d< localIndex > newEdgePackListArray; filterNonOwnedFromContainer( newEdgePackListArray, modifiedObjects.newEdges, parentEdgeIndices, edgeGhostRank, neighborRank );
  array1d< localIndex > modEdgePackListArray; filterNonOwnedFromContainer( modEdgePackListArray, modifiedObjects.modifiedEdges, parentEdgeIndices, edgeGhostRank, neighborRank );
  array1d< localIndex > newFacePackListArray; filterNonOwnedFromContainer( newFacePackListArray, modifiedObjects.newFaces, parentFaceIndices, faceGhostRank, neighborRank );
  array1d< localIndex > modFacePackListArray; filterNonOwnedFromContainer( modFacePackListArray, modifiedObjects.modifiedFaces, parentFaceIndices, faceGhostRank, neighborRank );


  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > newElemPackList;
  array1d< array1d< localIndex_array > > newElemData;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modElemPackList;
  array1d< array1d< localIndex_array > > modElemData;
  newElemPackList.resize( elemManager.numRegions());
  newElemData.resize( elemManager.numRegions());
  modElemPackList.resize( elemManager.numRegions());
  modElemData.resize( elemManager.numRegions());
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newElemPackList[er].resize( elemRegion.numSubRegions() );
    newElemData[er].resize( elemRegion.numSubRegions() );
    modElemPackList[er].resize( elemRegion.numSubRegions() );
    modElemData[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
    {
      ElementSubRegionBase & subRegion = elemRegion.getSubRegion( esr );
      arrayView1d< integer > const & subRegionGhostRank = subRegion.ghostRank();
      if( modifiedObjects.modifiedElements.count( {er, esr} ) > 0 )
      {
        std::set< localIndex > const & elemList = modifiedObjects.modifiedElements.at( {er, esr} );
        filterNonOwnedFromContainer( modElemData[er][esr], elemList, subRegionGhostRank, neighborRank );
      }

      if( modifiedObjects.newElements.count( {er, esr} ) > 0 )
      {
        std::set< localIndex > const & elemList = modifiedObjects.newElements.at( {er, esr} );
//        std::cout<<"modifiedObjects.newElements{ "<<er<<", "<<esr<<" } ( "<<MpiWrapper::commRank()<<" -> "<<neighbor.neighborRank()<<"):
// ";
        // for( auto const & index : elemList )
        // {
        //   std::cout<<index<<", ";
        // }
        // std::cout<<std::endl;
        filterNonOwnedFromContainer( newElemData[er][esr], elemList, subRegionGhostRank, neighborRank );
      }

      newElemPackList[er][esr] = newElemData[er][esr];
      modElemPackList[er][esr].set( modElemData[er][esr] );
    }
  }

  // if we start packing sizing on device + async, poll for completion
  parallelDeviceEvents sizeEvents;
  bufferSize += nodeManager.packGlobalMapsSize( newNodePackListArray, 0 );
  bufferSize += edgeManager.packGlobalMapsSize( newEdgePackListArray, 0 );
  bufferSize += faceManager.packGlobalMapsSize( newFacePackListArray, 0 );
  bufferSize += elemManager.packGlobalMapsSize( newElemPackList );

  bufferSize += nodeManager.packParentChildMapsSize( newNodePackListArray );
  bufferSize += edgeManager.packParentChildMapsSize( newEdgePackListArray );
  bufferSize += faceManager.packParentChildMapsSize( newFacePackListArray );
  bufferSize += elemManager.packFaceElementToFaceSize( newElemPackList );

  bufferSize += nodeManager.packUpDownMapsSize( newNodePackListArray );
  bufferSize += edgeManager.packUpDownMapsSize( newEdgePackListArray );
  bufferSize += faceManager.packUpDownMapsSize( newFacePackListArray );
  bufferSize += elemManager.packUpDownMapsSize( newElemPackList );

  bufferSize += nodeManager.packSize( newNodePackListArray, 0, false, sizeEvents );
  bufferSize += edgeManager.packSize( newEdgePackListArray, 0, false, sizeEvents );
  bufferSize += faceManager.packSize( newFacePackListArray, 0, false, sizeEvents );
  bufferSize += elemManager.packSize( newElemPackList );

  bufferSize += nodeManager.packUpDownMapsSize( modNodePackListArray );
  bufferSize += edgeManager.packUpDownMapsSize( modEdgePackListArray );
  bufferSize += faceManager.packUpDownMapsSize( modFacePackListArray );
  bufferSize += elemManager.packUpDownMapsSize( modElemPackList );

  bufferSize += nodeManager.packParentChildMapsSize( modNodePackListArray );
  bufferSize += edgeManager.packParentChildMapsSize( modEdgePackListArray );
  bufferSize += faceManager.packParentChildMapsSize( modFacePackListArray );

  bufferSize += nodeManager.packSize( modNodePackListArray, 0, false, sizeEvents );
  bufferSize += edgeManager.packSize( modEdgePackListArray, 0, false, sizeEvents );
  bufferSize += faceManager.packSize( modFacePackListArray, 0, false, sizeEvents );

  waitAllDeviceEvents( sizeEvents );
  neighbor.resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor.sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  // empty event buffer
  int packedSize = 0;
  parallelDeviceEvents packEvents;

  packedSize += nodeManager.packGlobalMaps( sendBufferPtr, newNodePackListArray, 0 );
  packedSize += edgeManager.packGlobalMaps( sendBufferPtr, newEdgePackListArray, 0 );
  packedSize += faceManager.packGlobalMaps( sendBufferPtr, newFacePackListArray, 0 );
  packedSize += elemManager.packGlobalMaps( sendBufferPtr, newElemPackList );

  packedSize += nodeManager.packParentChildMaps( sendBufferPtr, newNodePackListArray );
  packedSize += edgeManager.packParentChildMaps( sendBufferPtr, newEdgePackListArray );
  packedSize += faceManager.packParentChildMaps( sendBufferPtr, newFacePackListArray );
  packedSize += elemManager.packFaceElementToFace( sendBufferPtr, newElemPackList );

//  std::cout<<"packedSize ( "<<MpiWrapper::commRank()<<" -> "<<neighbor.neighborRank()<<"): "<<packedSize<<std::endl;

  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, newNodePackListArray );
  packedSize += edgeManager.packUpDownMaps( sendBufferPtr, newEdgePackListArray );
  packedSize += faceManager.packUpDownMaps( sendBufferPtr, newFacePackListArray );
  packedSize += elemManager.packUpDownMaps( sendBufferPtr, newElemPackList );

  packedSize += nodeManager.pack( sendBufferPtr, newNodePackListArray, 0, false, packEvents );
  packedSize += edgeManager.pack( sendBufferPtr, newEdgePackListArray, 0, false, packEvents );
  packedSize += faceManager.pack( sendBufferPtr, newFacePackListArray, 0, false, packEvents );
  packedSize += elemManager.pack( sendBufferPtr, newElemPackList );

  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, modNodePackListArray );
  packedSize += edgeManager.packUpDownMaps( sendBufferPtr, modEdgePackListArray );
  packedSize += faceManager.packUpDownMaps( sendBufferPtr, modFacePackListArray );
  packedSize += elemManager.packUpDownMaps( sendBufferPtr, modElemPackList );

  packedSize += nodeManager.packParentChildMaps( sendBufferPtr, modNodePackListArray );
  packedSize += edgeManager.packParentChildMaps( sendBufferPtr, modEdgePackListArray );
  packedSize += faceManager.packParentChildMaps( sendBufferPtr, modFacePackListArray );

  packedSize += nodeManager.pack( sendBufferPtr, modNodePackListArray, 0, false, packEvents );
  packedSize += edgeManager.pack( sendBufferPtr, modEdgePackListArray, 0, false, packEvents );
  packedSize += faceManager.pack( sendBufferPtr, modFacePackListArray, 0, false, packEvents );

  // poll for pack completion here
  waitAllDeviceEvents( packEvents );
  GEOS_ERROR_IF( bufferSize != packedSize,
                 "Allocated Buffer Size ("<<bufferSize<<") is not equal to packed buffer size("<<packedSize<<")" );


}


//***** 1B *****//
localIndex unpackNewObjectsOnOwningRanks( NeighborCommunicator & neighbor,
                                          MeshLevel * const mesh,
                                          int const commID,
                                          ModifiedObjectLists & receivedObjects,
                                          TopologyChangeUnpackStepData & unpackStateData )
{
  GEOS_MARK_FUNCTION;

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  unpackStateData.init( neighbor.receiveBuffer( commID ), elemManager );
  buffer_unit_type const * & receiveBufferPtr = unpackStateData.m_bufferPtr;

  localIndex_array & newLocalNodes = unpackStateData.m_nodes;
  localIndex_array & newLocalEdges = unpackStateData.m_edges;
  localIndex_array & newLocalFaces = unpackStateData.m_faces;
  ElementRegionManager::ElementReferenceAccessor< array1d< localIndex > > & newLocalElements = unpackStateData.m_elements;
  array1d< array1d< localIndex_array > > & newLocalElementsData = unpackStateData.m_elementsData;

  newLocalNodes.resize( 0 );
  newLocalEdges.resize( 0 );
  newLocalFaces.resize( 0 );


  newLocalElements.resize( elemManager.numRegions());
  newLocalElementsData.resize( elemManager.numRegions());
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newLocalElements[er].resize( elemRegion.numSubRegions());
    newLocalElementsData[er].resize( elemRegion.numSubRegions());
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newLocalElementsData[er][esr].resize( 0 );
      newLocalElements[er][esr].set( newLocalElementsData[er][esr] );
    }
  }

  // if we move to device + async packing here, add polling of events or pass out
  buffer_type::size_type & unpackedSize = unpackStateData.m_size;
  unpackedSize += nodeManager.unpackGlobalMaps( receiveBufferPtr, newLocalNodes, 0 );
  unpackedSize += edgeManager.unpackGlobalMaps( receiveBufferPtr, newLocalEdges, 0 );
  unpackedSize += faceManager.unpackGlobalMaps( receiveBufferPtr, newLocalFaces, 0 );
  unpackedSize += elemManager.unpackGlobalMaps( receiveBufferPtr, newLocalElements );

  unpackedSize += nodeManager.unpackParentChildMaps( receiveBufferPtr, newLocalNodes );
  unpackedSize += edgeManager.unpackParentChildMaps( receiveBufferPtr, newLocalEdges );
  unpackedSize += faceManager.unpackParentChildMaps( receiveBufferPtr, newLocalFaces );
  unpackedSize += elemManager.unpackFaceElementToFace( receiveBufferPtr, newLocalElements, true );

  // std::cout<<"unpackedSize ( "<<neighbor.neighborRank()<<" -> "<<MpiWrapper::commRank()<<"): "<<unpackedSize<<std::endl;
  // std::cout<<" end of 1b receiveBufferPtr ("<<neighbor.neighborRank()<<" -> "<<MpiWrapper::commRank()<<" ) = "<<reinterpret_cast<void
  // const *>(receiveBufferPtr)<<std::endl;

  std::set< localIndex > & allNewNodes      = receivedObjects.newNodes;
  std::set< localIndex > & allNewEdges      = receivedObjects.newEdges;
  std::set< localIndex > & allNewFaces      = receivedObjects.newFaces;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > & allNewElements = receivedObjects.newElements;

  allNewNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
  allNewEdges.insert( newLocalEdges.begin(), newLocalEdges.end() );
  allNewFaces.insert( newLocalFaces.begin(), newLocalFaces.end() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
    {
      allNewElements[{er, esr}].insert( newLocalElements[er][esr].get().begin(),
                                        newLocalElements[er][esr].get().end() );
    }
  }

  return unpackedSize;
}


//***** 2a *****//
void packNewObjectsToGhosts( NeighborCommunicator & neighbor,
                             int commID,
                             MeshLevel * const mesh,
                             TopologyChangeStepData & packData,
                             ModifiedObjectLists & modifiedObjects )
{
  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  localIndex_array & newNodesToSend = packData.m_nodes;
  localIndex_array & newEdgesToSend = packData.m_edges;
  localIndex_array & newFacesToSend = packData.m_faces;
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > & newElemsToSend = packData.m_elementsView;
  array1d< array1d< localIndex_array > > & newElemsToSendData = packData.m_elementsData;


  localIndex_array & nodeGhostsToSend = nodeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  localIndex_array & edgeGhostsToSend = edgeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  localIndex_array & faceGhostsToSend = faceManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();

  arrayView1d< localIndex > const & nodalParentIndices = nodeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex > const & edgeParentIndices = edgeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex > const & faceParentIndices = faceManager.getField< fields::parentIndex >();

  filterNewObjectsForPackToGhosts( modifiedObjects.newNodes, nodalParentIndices, nodeGhostsToSend, newNodesToSend );
  filterNewObjectsForPackToGhosts( modifiedObjects.newEdges, edgeParentIndices, edgeGhostsToSend, newEdgesToSend );
  filterNewObjectsForPackToGhosts( modifiedObjects.newFaces, faceParentIndices, faceGhostsToSend, newFacesToSend );

  SortedArray< localIndex > faceGhostsToSendSet;
  for( localIndex const & kf : faceGhostsToSend )
  {
    faceGhostsToSendSet.insert( kf );
  }

  newElemsToSendData.resize( elemManager.numRegions() );
  newElemsToSend.resize( elemManager.numRegions() );
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newElemsToSendData[er].resize( elemRegion.numSubRegions() );
    newElemsToSend[er].resize( elemRegion.numSubRegions() );

    elemRegion.forElementSubRegionsIndex< FaceElementSubRegion >( [&]( localIndex const esr,
                                                                       FaceElementSubRegion & subRegion )
    {
      ArrayOfArraysView< localIndex const > const faceList = subRegion.faceList().toViewConst();
      localIndex_array & elemGhostsToSend = subRegion.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
      elemGhostsToSend.move( hostMemorySpace );
      for( localIndex const & k : modifiedObjects.newElements.at( {er, esr} ) )
      {
        if( faceGhostsToSendSet.count( faceList( k, 0 ) ) )
        {
          newElemsToSendData[er][esr].emplace_back( k );
          elemGhostsToSend.emplace_back( k );
        }
      }
      newElemsToSend[er][esr] = newElemsToSendData[er][esr];
    } );
  }

  int bufferSize = 0;

  bufferSize += nodeManager.packGlobalMapsSize( newNodesToSend, 0 );
  bufferSize += edgeManager.packGlobalMapsSize( newEdgesToSend, 0 );
  bufferSize += faceManager.packGlobalMapsSize( newFacesToSend, 0 );
  bufferSize += elemManager.packGlobalMapsSize( newElemsToSend );
  bufferSize += elemManager.packFaceElementToFaceSize( newElemsToSend );

  neighbor.resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor.sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = 0;

  packedSize += nodeManager.packGlobalMaps( sendBufferPtr, newNodesToSend, 0 );
  packedSize += edgeManager.packGlobalMaps( sendBufferPtr, newEdgesToSend, 0 );
  packedSize += faceManager.packGlobalMaps( sendBufferPtr, newFacesToSend, 0 );
  packedSize += elemManager.packGlobalMaps( sendBufferPtr, newElemsToSend );
  packedSize += elemManager.packFaceElementToFace( sendBufferPtr, newElemsToSend );

  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );
}


//***** 2b *****//
void unpackNewObjectsOnGhosts( NeighborCommunicator & neighbor,
                               int commID,
                               MeshLevel * const mesh,
                               ModifiedObjectLists & receivedObjects )
{

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  localIndex_array & nodeGhostsToRecv = nodeManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();
  localIndex_array & edgeGhostsToRecv = edgeManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();
  localIndex_array & faceGhostsToRecv = faceManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();

  buffer_type const & receiveBuffer = neighbor.receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  localIndex_array newGhostNodes;
  localIndex_array newGhostEdges;
  localIndex_array newGhostFaces;

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

  // if we move to device + async unoacking, poll these events for completion or pass out
  parallelDeviceEvents events;

  nodeManager.unpackGlobalMaps( receiveBufferPtr, newGhostNodes, 0 );
  edgeManager.unpackGlobalMaps( receiveBufferPtr, newGhostEdges, 0 );
  faceManager.unpackGlobalMaps( receiveBufferPtr, newGhostFaces, 0 );
  elemManager.unpackGlobalMaps( receiveBufferPtr, newGhostElems );
  elemManager.unpackFaceElementToFace( receiveBufferPtr, newGhostElems, true );

  waitAllDeviceEvents( events );

  if( newGhostNodes.size() > 0 )
  {
    nodeGhostsToRecv.move( hostMemorySpace );
    for( localIndex a=0; a<newGhostNodes.size(); ++a )
    {
      nodeGhostsToRecv.emplace_back( newGhostNodes[a] );
    }
  }

  if( newGhostEdges.size() > 0 )
  {
    edgeGhostsToRecv.move( hostMemorySpace );
    for( localIndex a=0; a<newGhostEdges.size(); ++a )
    {
      edgeGhostsToRecv.emplace_back( newGhostEdges[a] );
    }
  }

  if( newGhostFaces.size() > 0 )
  {
    faceGhostsToRecv.move( hostMemorySpace );
    for( localIndex a=0; a<newGhostFaces.size(); ++a )
    {
      faceGhostsToRecv.emplace_back( newGhostFaces[a] );
    }
  }

  elemManager.forElementSubRegionsComplete< ElementSubRegionBase >(
    [&]( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    localIndex_array & elemGhostsToReceive = subRegion.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();
    if( newGhostElemsData[er][esr].size() > 0 )
    {
      elemGhostsToReceive.move( hostMemorySpace );

      for( localIndex const & newElemIndex : newGhostElemsData[er][esr] )
      {
        elemGhostsToReceive.emplace_back( newElemIndex );
        receivedObjects.newElements[ { er, esr } ].insert( newElemIndex );
      }
    }
  } );

  receivedObjects.newNodes.insert( newGhostNodes.begin(), newGhostNodes.end() );
  receivedObjects.newEdges.insert( newGhostEdges.begin(), newGhostEdges.end() );
  receivedObjects.newFaces.insert( newGhostFaces.begin(), newGhostFaces.end() );
}



//***** 3a *****//
localIndex unpackNewAndModifiedObjectsDataOnOwningRanks( MeshLevel * const mesh,
                                                         ModifiedObjectLists & receivedObjects,
                                                         TopologyChangeUnpackStepData & unpackStateData )
{
  GEOS_MARK_FUNCTION;

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  buffer_unit_type const * & receiveBufferPtr = unpackStateData.m_bufferPtr;

  localIndex_array & newLocalNodes = unpackStateData.m_nodes;
  localIndex_array & newLocalEdges = unpackStateData.m_edges;
  localIndex_array & newLocalFaces = unpackStateData.m_faces;

  ElementRegionManager::ElementReferenceAccessor< array1d< localIndex > > & newLocalElements = unpackStateData.m_elements;

  localIndex_array modifiedLocalNodes;
  localIndex_array modifiedLocalEdges;
  localIndex_array modifiedLocalFaces;

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modifiedLocalElements;
  array1d< array1d< localIndex_array > > modifiedLocalElementsData;

  modifiedLocalElements.resize( elemManager.numRegions());
  modifiedLocalElementsData.resize( elemManager.numRegions());
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    modifiedLocalElements[er].resize( elemRegion.numSubRegions());
    modifiedLocalElementsData[er].resize( elemRegion.numSubRegions());
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      modifiedLocalElements[er][esr].set( modifiedLocalElementsData[er][esr] );
    }
  }

  // if we move to device + async packing here, add polling of events or pass out
  parallelDeviceEvents events;
  buffer_type::size_type & unpackedSize = unpackStateData.m_size;

  unpackedSize += nodeManager.unpackUpDownMaps( receiveBufferPtr, newLocalNodes, true, true );
  unpackedSize += edgeManager.unpackUpDownMaps( receiveBufferPtr, newLocalEdges, true, true );
  unpackedSize += faceManager.unpackUpDownMaps( receiveBufferPtr, newLocalFaces, true, true );
  unpackedSize += elemManager.unpackUpDownMaps( receiveBufferPtr, newLocalElements, true );

  unpackedSize += nodeManager.unpack( receiveBufferPtr, newLocalNodes, 0, false, events );
  unpackedSize += edgeManager.unpack( receiveBufferPtr, newLocalEdges, 0, false, events );
  unpackedSize += faceManager.unpack( receiveBufferPtr, newLocalFaces, 0, false, events );
  unpackedSize += elemManager.unpack( receiveBufferPtr, newLocalElements );

  unpackedSize += nodeManager.unpackUpDownMaps( receiveBufferPtr, modifiedLocalNodes, false, true );
  unpackedSize += edgeManager.unpackUpDownMaps( receiveBufferPtr, modifiedLocalEdges, false, true );
  unpackedSize += faceManager.unpackUpDownMaps( receiveBufferPtr, modifiedLocalFaces, false, true );
  unpackedSize += elemManager.unpackUpDownMaps( receiveBufferPtr, modifiedLocalElements, true );

  unpackedSize += nodeManager.unpackParentChildMaps( receiveBufferPtr, modifiedLocalNodes );
  unpackedSize += edgeManager.unpackParentChildMaps( receiveBufferPtr, modifiedLocalEdges );
  unpackedSize += faceManager.unpackParentChildMaps( receiveBufferPtr, modifiedLocalFaces );

  unpackedSize += nodeManager.unpack( receiveBufferPtr, modifiedLocalNodes, 0, false, events );
  unpackedSize += edgeManager.unpack( receiveBufferPtr, modifiedLocalEdges, 0, false, events );
  unpackedSize += faceManager.unpack( receiveBufferPtr, modifiedLocalFaces, 0, false, events );

  waitAllDeviceEvents( events );

  std::set< localIndex > & allModifiedNodes = receivedObjects.modifiedNodes;
  std::set< localIndex > & allModifiedEdges = receivedObjects.modifiedEdges;
  std::set< localIndex > & allModifiedFaces = receivedObjects.modifiedFaces;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > & allModifiedElements = receivedObjects.modifiedElements;

  allModifiedNodes.insert( modifiedLocalNodes.begin(), modifiedLocalNodes.end() );

  allModifiedEdges.insert( modifiedLocalEdges.begin(), modifiedLocalEdges.end() );

  allModifiedFaces.insert( modifiedLocalFaces.begin(), modifiedLocalFaces.end() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
    {
      allModifiedElements[{er, esr}].insert( modifiedLocalElements[er][esr].get().begin(),
                                             modifiedLocalElements[er][esr].get().end() );
    }
  }

  return unpackedSize;
}



//***** 3b *****//
void packNewModifiedObjectsToGhosts( NeighborCommunicator & neighbor,
                                     int commID,
                                     MeshLevel * const mesh,
                                     TopologyChangeStepData & packData,
                                     ModifiedObjectLists & receivedObjects )
{
  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  localIndex_array & newNodesToSend = packData.m_nodes;
  localIndex_array & newEdgesToSend = packData.m_edges;
  localIndex_array & newFacesToSend = packData.m_faces;
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > & newElemsToSend = packData.m_elementsView;
  // array1d< array1d< localIndex_array > > & newElemsToSendData = packData.m_elementsData;

  localIndex_array modNodesToSend;
  localIndex_array modEdgesToSend;
  localIndex_array modFacesToSend;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modElemsToSend;
  array1d< array1d< localIndex_array > > modElemsToSendData;

  localIndex_array & nodeGhostsToSend = nodeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  localIndex_array & edgeGhostsToSend = edgeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  localIndex_array & faceGhostsToSend = faceManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();

  filterModObjectsForPackToGhosts( receivedObjects.modifiedNodes, nodeGhostsToSend, modNodesToSend );
  filterModObjectsForPackToGhosts( receivedObjects.modifiedEdges, edgeGhostsToSend, modEdgesToSend );
  filterModObjectsForPackToGhosts( receivedObjects.modifiedFaces, faceGhostsToSend, modFacesToSend );

  // newElemsToSendData.resize( elemManager.numRegions() );
  // newElemsToSend.resize( elemManager.numRegions() );
  modElemsToSendData.resize( elemManager.numRegions() );
  modElemsToSend.resize( elemManager.numRegions() );
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    // newElemsToSendData[er].resize( elemRegion.numSubRegions() );
    // newElemsToSend[er].resize( elemRegion.numSubRegions() );
    modElemsToSendData[er].resize( elemRegion.numSubRegions() );
    modElemsToSend[er].resize( elemRegion.numSubRegions() );

    // elemRegion.forElementSubRegionsIndex< FaceElementSubRegion >( [&]( localIndex const esr,
    //                                                                    FaceElementSubRegion & subRegion )
    // {
    //   ArrayOfArraysView< localIndex const > const faceList = subRegion.faceList().toViewConst();
    //   localIndex_array & elemGhostsToSend = subRegion.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
    //   elemGhostsToSend.move( hostMemorySpace );
    //   for( localIndex const & k : receivedObjects.newElements.at( {er, esr} ) )
    //   {
    //     if( faceGhostsToSendSet.count( faceList( k, 0 ) ) )
    //     {
    //       newElemsToSendData[er][esr].emplace_back( k );
    //       elemGhostsToSend.emplace_back( k );
    //     }
    //   }
    //   newElemsToSend[er][esr] = newElemsToSendData[er][esr];
    // } );

    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >( [&]( localIndex const esr,
                                                                       ElementSubRegionBase const & subRegion )
    {
      modElemsToSend[er][esr].set( modElemsToSendData[er][esr] );
      arrayView1d< localIndex const > const & elemGhostsToSend = subRegion.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
      for( localIndex const ghostToSend : elemGhostsToSend )
      {
        if( receivedObjects.modifiedElements.at( { er, esr } ).count( ghostToSend ) > 0 )
        {
          modElemsToSendData[er][esr].emplace_back( ghostToSend );
        }
      }
    } );
  }

  parallelDeviceEvents sizeEvents;
  int bufferSize = 0;


  bufferSize += nodeManager.packUpDownMapsSize( newNodesToSend );
  bufferSize += edgeManager.packUpDownMapsSize( newEdgesToSend );
  bufferSize += faceManager.packUpDownMapsSize( newFacesToSend );
  bufferSize += elemManager.packUpDownMapsSize( newElemsToSend );

  bufferSize += nodeManager.packParentChildMapsSize( newNodesToSend );
  bufferSize += edgeManager.packParentChildMapsSize( newEdgesToSend );
  bufferSize += faceManager.packParentChildMapsSize( newFacesToSend );

  bufferSize += nodeManager.packSize( newNodesToSend, 0, false, sizeEvents );
  bufferSize += edgeManager.packSize( newEdgesToSend, 0, false, sizeEvents );
  bufferSize += faceManager.packSize( newFacesToSend, 0, false, sizeEvents );
  bufferSize += elemManager.packSize( newElemsToSend );

  bufferSize += nodeManager.packUpDownMapsSize( modNodesToSend );
  bufferSize += edgeManager.packUpDownMapsSize( modEdgesToSend );
  bufferSize += faceManager.packUpDownMapsSize( modFacesToSend );
  bufferSize += elemManager.packUpDownMapsSize( modElemsToSend );

  bufferSize += nodeManager.packParentChildMapsSize( modNodesToSend );
  bufferSize += edgeManager.packParentChildMapsSize( modEdgesToSend );
  bufferSize += faceManager.packParentChildMapsSize( modFacesToSend );

  waitAllDeviceEvents( sizeEvents );
  neighbor.resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor.sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  parallelDeviceEvents packEvents;
  int packedSize = 0;

  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, newNodesToSend );
  packedSize += edgeManager.packUpDownMaps( sendBufferPtr, newEdgesToSend );
  packedSize += faceManager.packUpDownMaps( sendBufferPtr, newFacesToSend );
  packedSize += elemManager.packUpDownMaps( sendBufferPtr, newElemsToSend );

  packedSize += nodeManager.packParentChildMaps( sendBufferPtr, newNodesToSend );
  packedSize += edgeManager.packParentChildMaps( sendBufferPtr, newEdgesToSend );
  packedSize += faceManager.packParentChildMaps( sendBufferPtr, newFacesToSend );

  packedSize += nodeManager.pack( sendBufferPtr, newNodesToSend, 0, false, packEvents );
  packedSize += edgeManager.pack( sendBufferPtr, newEdgesToSend, 0, false, packEvents );
  packedSize += faceManager.pack( sendBufferPtr, newFacesToSend, 0, false, packEvents );
  packedSize += elemManager.pack( sendBufferPtr, newElemsToSend );

  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, modNodesToSend );
  packedSize += edgeManager.packUpDownMaps( sendBufferPtr, modEdgesToSend );
  packedSize += faceManager.packUpDownMaps( sendBufferPtr, modFacesToSend );
  packedSize += elemManager.packUpDownMaps( sendBufferPtr, modElemsToSend );

  packedSize += nodeManager.packParentChildMaps( sendBufferPtr, modNodesToSend );
  packedSize += edgeManager.packParentChildMaps( sendBufferPtr, modEdgesToSend );
  packedSize += faceManager.packParentChildMaps( sendBufferPtr, modFacesToSend );

  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );

  waitAllDeviceEvents( packEvents );
}


//***** 3c *****
void unpackNewAndModifiedObjectsDataOnGhosts( NeighborCommunicator & neighbor,
                                              int commID,
                                              MeshLevel * const mesh,
                                              ModifiedObjectLists & receivedObjects )
{

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  buffer_type const & receiveBuffer = neighbor.receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  localIndex_array newGhostNodes;
  localIndex_array newGhostEdges;
  localIndex_array newGhostFaces;

  localIndex_array modGhostNodes;
  localIndex_array modGhostEdges;
  localIndex_array modGhostFaces;

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > newGhostElems;
  array1d< array1d< localIndex_array > > newGhostElemsData;
  newGhostElems.resize( elemManager.numRegions() );
  newGhostElemsData.resize( elemManager.numRegions() );
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modGhostElems;
  array1d< array1d< localIndex_array > > modGhostElemsData;
  modGhostElems.resize( elemManager.numRegions() );
  modGhostElemsData.resize( elemManager.numRegions() );
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newGhostElemsData[er].resize( elemRegion.numSubRegions() );
    newGhostElems[er].resize( elemRegion.numSubRegions() );
    modGhostElemsData[er].resize( elemRegion.numSubRegions() );
    modGhostElems[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newGhostElems[er][esr].set( newGhostElemsData[er][esr] );
      modGhostElems[er][esr].set( modGhostElemsData[er][esr] );
    }
  }

  // if we move to device + async unpacking, poll these events for completion or pass out
  parallelDeviceEvents events;
  nodeManager.unpackUpDownMaps( receiveBufferPtr, newGhostNodes, true, true );
  edgeManager.unpackUpDownMaps( receiveBufferPtr, newGhostEdges, true, true );
  faceManager.unpackUpDownMaps( receiveBufferPtr, newGhostFaces, true, true );
  elemManager.unpackUpDownMaps( receiveBufferPtr, newGhostElems, true );

  nodeManager.unpackParentChildMaps( receiveBufferPtr, newGhostNodes );
  edgeManager.unpackParentChildMaps( receiveBufferPtr, newGhostEdges );
  faceManager.unpackParentChildMaps( receiveBufferPtr, newGhostFaces );

  nodeManager.unpack( receiveBufferPtr, newGhostNodes, 0, false, events );
  edgeManager.unpack( receiveBufferPtr, newGhostEdges, 0, false, events );
  faceManager.unpack( receiveBufferPtr, newGhostFaces, 0, false, events );
  elemManager.unpack( receiveBufferPtr, newGhostElems );

  nodeManager.unpackUpDownMaps( receiveBufferPtr, modGhostNodes, false, true );
  edgeManager.unpackUpDownMaps( receiveBufferPtr, modGhostEdges, false, true );
  faceManager.unpackUpDownMaps( receiveBufferPtr, modGhostFaces, false, true );
  elemManager.unpackUpDownMaps( receiveBufferPtr, modGhostElems, true );

  nodeManager.unpackParentChildMaps( receiveBufferPtr, modGhostNodes );
  edgeManager.unpackParentChildMaps( receiveBufferPtr, modGhostEdges );
  faceManager.unpackParentChildMaps( receiveBufferPtr, modGhostFaces );

  waitAllDeviceEvents( events );


  elemManager.forElementSubRegionsComplete< ElementSubRegionBase >(
    [&]( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & )
  {
    receivedObjects.modifiedElements[ { er, esr } ].insert( modGhostElemsData[er][esr].begin(),
                                                            modGhostElemsData[er][esr].end() );
  } );

  receivedObjects.modifiedNodes.insert( modGhostNodes.begin(), modGhostNodes.end() );
  receivedObjects.modifiedEdges.insert( modGhostEdges.begin(), modGhostEdges.end() );
  receivedObjects.modifiedFaces.insert( modGhostFaces.begin(), modGhostFaces.end() );

}



void updateConnectorsToFaceElems( std::set< localIndex > const & newFaceElements,
                                  FaceElementSubRegion & faceElemSubRegion )
{
  ArrayOfArrays< localIndex > & connectorToElem = faceElemSubRegion.m_2dFaceTo2dElems;
  map< localIndex, localIndex > & edgesToConnectorEdges = faceElemSubRegion.m_edgesTo2dFaces;
  array1d< localIndex > & connectorEdgesToEdges = faceElemSubRegion.m_2dFaceToEdge;

  ArrayOfArraysView< localIndex const > const facesToEdges = faceElemSubRegion.edgeList().toViewConst();

  for( localIndex const & kfe : newFaceElements )
  {
    arraySlice1d< localIndex const > const faceToEdges = facesToEdges[kfe];
    for( localIndex ke=0; ke<faceToEdges.size(); ++ke )
    {
      localIndex const edgeIndex = faceToEdges[ke];

      auto connIter = edgesToConnectorEdges.find( edgeIndex );
      if( connIter==edgesToConnectorEdges.end() )
      {
        connectorToElem.appendArray( 0 );
        connectorEdgesToEdges.emplace_back( edgeIndex );
        edgesToConnectorEdges[edgeIndex] = connectorEdgesToEdges.size() - 1;
      }
      localIndex const connectorIndex = edgesToConnectorEdges.at( edgeIndex );

      localIndex const numExistingCells = connectorToElem.sizeOfArray( connectorIndex );
      bool cellExistsInMap = false;
      for( localIndex k=0; k<numExistingCells; ++k )
      {
        if( kfe == connectorToElem[connectorIndex][k] )
        {
          cellExistsInMap = true;
        }
      }
      if( !cellExistsInMap )
      {
        connectorToElem.resizeArray( connectorIndex, numExistingCells+1 );
        connectorToElem[connectorIndex][ numExistingCells ] = kfe;
        faceElemSubRegion.m_recalculateConnectionsFor2dFaces.insert( connectorIndex );
      }
    }
  }
}



void synchronizeTopologyChange( MeshLevel * const mesh,
                                std::vector< NeighborCommunicator > & neighbors,
                                ModifiedObjectLists & modifiedObjects,
                                ModifiedObjectLists & receivedObjects,
                                int mpiCommOrder )
{

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();


  /************************************************************************************************
   * The goal is to synchronize the changes from the rank that has topology changes to
   * ranks that have copies of the objects that were changed. In this "original" implementation, we
   * do this without map unpacking optimizations intended to reduce communications.
   *
   * Nomenclature is key to understanding the process:
   * - "New" objects are objects that have just been created on by the "active color rank (ACR)"
   * - "Modified" objects are objects that have been modified by the ACR.
   *
   * - ACR (active color rank) is the rank that has created the topology changes. Given the way we
   *   map the colors to ranks, the ACR are NOT neighbors...i.e. do not communicate with each other.
   * - OR (Owning rank/s) is the rank that owns the "new/modified" objects. This may or may not be
   *   the ACR.
   * - GR (Ghosted rank/s) is the rank that has a ghost copy of the "new/modified" object.
   *
   * note: object parents define the owning rank.
   * note: for any receive/unpack operation, the current rank is the rank performing the operation
   *       from each neighbor...i.e. the current rank is the OR and the GR.
   *
   * The sequence of steps are:
   * 1a) On the ACR, pack the new/modified objects that are not owned by the ACR and send them to
   *     their OR.
   * 1b) On the OR, unpack the new objects that are owned by the rank that has the changes. DO NOT
   *     unpack the maps as they will potentially contain indices that are not on the OR.
   *
   * At this point the OR has all the new objects that it owns...but not the maps or the fields.
   *
   * 2a) On the OR, pack the new objects that are owned by the rank and send them to the ranks
   *     where they are ghosted (GR). DO NOT PACK THE MAPS as they are incomplete.
   * 2b) On the GR, unpack the new objects.
   *
   * Now everyone has all the objects and we can pack/send/receive/unpack the maps.
   *
   * 3a) On the OR, unpack the map modification on owning ranks from 1b).
   *
   * Now the OR has the correct maps.
   *
   * 3b) On the OR, pack the map/field modification and send to the GR.
   * 3c) On the GR, unpack the map/field modifications.
   *
   ***********************************************************************************************/



  //***********************************************************************************************
  // 1a) On the ACR, pack the new/modified objects that are not owned by the ACR and send them to
  //     their OR.
  //***********************************************************************************************

//  std::cout<<"***** Step 1a *****"<<std::endl;;
  MPI_iCommData commData1;
  commData1.resize( neighbors.size() );
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewAndModifiedObjectsToOwningRanks( neighbor,
                                            mesh,
                                            modifiedObjects,
                                            commData1.commID() );

    neighbor.mpiISendReceiveBufferSizes( commData1.commID(),
                                         commData1.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData1.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOS );

  }

  // send/recv the buffers
  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData1.size(),
                         commData1.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData1.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData1.commID(),
                                     commData1.mpiSendBufferRequest( neighborIndex ),
                                     commData1.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOS );
  }

  MpiWrapper::barrier();
//  std::cout<<"***** Step 1b *****"<<std::endl;

  //***********************************************************************************************
  // 1b) On the OR, unpack the new objects that are owned by the rank that has the changes. DO NOT
  //     unpack the maps as they will potentially contain indices that are not on the OR.
  //***********************************************************************************************
  std::vector< TopologyChangeUnpackStepData > step1bUnpackData( neighbors.size() );
  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex = count;
    if( mpiCommOrder == 0 )
    {
      MpiWrapper::waitAny( commData1.size(),
                           commData1.mpiRecvBufferRequest(),
                           &neighborIndex,
                           commData1.mpiRecvBufferStatus() );
    }
    // Unpack buffers in set ordering for integration testing
    else
    {
      MpiWrapper::wait( commData1.mpiRecvBufferRequest() + count,
                        commData1.mpiRecvBufferStatus() + count );
    }

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    unpackNewObjectsOnOwningRanks( neighbor,
                                   mesh,
                                   commData1.commID(),
                                   receivedObjects,
                                   step1bUnpackData[neighborIndex] );

  }

  nodeManager.inheritGhostRankFromParent( receivedObjects.newNodes );
  edgeManager.inheritGhostRankFromParent( receivedObjects.newEdges );
  faceManager.inheritGhostRankFromParent( receivedObjects.newFaces );
  elemManager.forElementSubRegionsComplete< FaceElementSubRegion >( [&]( localIndex const er,
                                                                         localIndex const esr,
                                                                         ElementRegionBase &,
                                                                         FaceElementSubRegion & subRegion )
  {
    subRegion.inheritGhostRankFromParentFace( faceManager, receivedObjects.newElements[{er, esr}] );
  } );

  MpiWrapper::waitAll( commData1.size(),
                       commData1.mpiSendBufferSizeRequest(),
                       commData1.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData1.size(),
                       commData1.mpiSendBufferRequest(),
                       commData1.mpiSendBufferSizeStatus() );

  modifiedObjects.insert( receivedObjects );


  //************************************************************************************************
  // 2a) On the OR, pack the new objects that are owned by the rank and send them to the ranks
  //     where they are ghosted (GR). DO NOT PACK THE MAPS as they are incomplete.
  //************************************************************************************************

  MpiWrapper::barrier();
//  std::cout<<"***** Step 2a *****"<<std::endl;
  MpiWrapper::barrier();

  // a new MPI_iCommData object is created to avoid overwriting the previous one which isn't
  // finished unpacking
  MPI_iCommData commData2;
  commData2.resize( neighbors.size());
  std::vector< TopologyChangeStepData > step2and3PackData( neighbors.size() );

  // pack the new objects to send to ghost ranks
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewObjectsToGhosts( neighbor,
                            commData2.commID(),
                            mesh,
                            step2and3PackData[neighborIndex],
                            modifiedObjects );

    neighbor.mpiISendReceiveBufferSizes( commData2.commID(),
                                         commData2.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData2.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOS );
  }
  MpiWrapper::barrier();
//  std::cout<<"***** Step 2a - bp2 *****"<<std::endl;
  MpiWrapper::barrier();

  // send/recv buffer sizes
  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData2.size(),
                         commData2.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData2.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData2.commID(),
                                     commData2.mpiSendBufferRequest( neighborIndex ),
                                     commData2.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOS );
  }

  //************************************************************************************************
  // 2b) On the GR, unpack the new objects.
  //************************************************************************************************
  MpiWrapper::barrier();
//  std::cout<<"***** Step 2b *****"<<std::endl;
  MpiWrapper::barrier();

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;
    if( mpiCommOrder == 0 )
    {
      MpiWrapper::waitAny( commData2.size(),
                           commData2.mpiRecvBufferRequest(),
                           &neighborIndex,
                           commData2.mpiRecvBufferStatus() );
    }
    else
    {
      MpiWrapper::wait( commData2.mpiRecvBufferRequest() + count,
                        commData2.mpiRecvBufferStatus() + count );
    }
    unpackNewObjectsOnGhosts( neighbors[neighborIndex], commData2.commID(), mesh, receivedObjects );
  }

  MpiWrapper::waitAll( commData2.size(),
                       commData2.mpiSendBufferSizeRequest(),
                       commData2.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData2.size(),
                       commData2.mpiSendBufferRequest(),
                       commData2.mpiSendBufferSizeStatus() );

  //************************************************************************************************
  // 3a) On the OR, unpack the map modification on owning ranks from 1b).
  //************************************************************************************************
//  std::cout<<"***** Step 3a *****"<<std::endl;

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    unpackNewAndModifiedObjectsDataOnOwningRanks( mesh,
                                                  receivedObjects,
                                                  step1bUnpackData[count] );
  }

  modifiedObjects.insert( receivedObjects );



  //************************************************************************************************
  // 3b) On the OR, pack the map/field modification and send to the GR.
  //************************************************************************************************
//  std::cout<<"***** Step 3b *****"<<std::endl;

  // a new MPI_iCommData object is created...just because.
  MPI_iCommData commData3;
  commData3.resize( neighbors.size());

  // pack the new objects to send to ghost ranks
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewModifiedObjectsToGhosts( neighbor,
                                    commData3.commID(),
                                    mesh,
                                    step2and3PackData[neighborIndex],
                                    modifiedObjects );

    neighbor.mpiISendReceiveBufferSizes( commData3.commID(),
                                         commData3.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData3.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOS );
  }

//  std::cout<<"***** Step 3b - bp2 *****"<<std::endl;
  // send/recv buffer sizes
  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData3.size(),
                         commData3.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData3.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData3.commID(),
                                     commData3.mpiSendBufferRequest( neighborIndex ),
                                     commData3.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOS );
  }



  //************************************************************************************************
  // 3c) On the GR, unpack the map/field modifications.
  //************************************************************************************************

//  std::cout<<"***** Step 3c *****"<<std::endl;


  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;
    if( mpiCommOrder == 0 )
    {
      MpiWrapper::waitAny( commData3.size(),
                           commData3.mpiRecvBufferRequest(),
                           &neighborIndex,
                           commData3.mpiRecvBufferStatus() );
    }
    else
    {
      MpiWrapper::wait( commData3.mpiRecvBufferRequest() + count,
                        commData3.mpiRecvBufferStatus() + count );
    }
    unpackNewAndModifiedObjectsDataOnGhosts( neighbors[neighborIndex], commData3.commID(), mesh, receivedObjects );
  }



  MpiWrapper::waitAll( commData3.size(),
                       commData3.mpiSendBufferSizeRequest(),
                       commData3.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData3.size(),
                       commData3.mpiSendBufferRequest(),
                       commData3.mpiSendBufferSizeStatus() );

  modifiedObjects.insert( receivedObjects );


  // nodeManager.fixUpDownMaps( false );
  // edgeManager.fixUpDownMaps( false );
  // faceManager.fixUpDownMaps( false );


  // for( localIndex er = 0; er < elemManager.numRegions(); ++er )
  // {
  //   ElementRegionBase & elemRegion = elemManager.getRegion( er );
  //   for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
  //   {
  //     elemRegion.getSubRegion( esr ).fixUpDownMaps( false );
  //   }
  // }

  // for( int rank=0; rank<MpiWrapper::commSize(); ++rank )
  // {
  //   MpiWrapper::barrier();
  //   if( rank != MpiWrapper::commRank() )
  //   {
  //     for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  //     {
  //       NeighborCommunicator & neighbor = neighbors[neighborIndex];
  //       localIndex_array const & nodeGhostsToSend = nodeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  //       localIndex_array const & nodeGhostsToRecv = nodeManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();
  //       arrayView1d< globalIndex const > const & nodeGlobalIndices = nodeManager.localToGlobalMap();

  //       localIndex_array const & edgeGhostsToSend = edgeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  //       localIndex_array const & edgeGhostsToRecv = edgeManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();
  //       arrayView1d< globalIndex const > const & edgeGlobalIndices = edgeManager.localToGlobalMap();

  //       localIndex_array const & faceGhostsToSend = faceManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  //       localIndex_array const & faceGhostsToRecv = faceManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();
  //       arrayView1d< globalIndex const > const & faceGlobalIndices = faceManager.localToGlobalMap();

  //       std::cout<< "Rank: " << MpiWrapper::commRank() << " Neighbor: " << neighbor.neighborRank() << std::endl;

  //       std::cout<< "  nodeGhostsToSend: " ;
  //       for( localIndex const & k : nodeGhostsToSend )
  //       {
  //         std::cout << nodeGlobalIndices[k] << ", ";
  //       }
  //       std::cout<< std::endl;

  //       std::cout<< "  nodeGhostsToRecv: ";
  //       for( localIndex const & k : nodeGhostsToRecv )
  //       {
  //         std::cout << nodeGlobalIndices[k] << ", ";
  //       }
  //       std::cout<< std::endl;

  //       std::cout<< "  edgeGhostsToSend: ";
  //       for( localIndex const & k : edgeGhostsToSend )
  //       {
  //         std::cout << edgeGlobalIndices[k] << ", ";
  //       }
  //       std::cout<< std::endl;

  //       std::cout<< "  edgeGhostsToRecv: ";
  //       for( localIndex const & k : edgeGhostsToRecv )
  //       {
  //         std::cout << edgeGlobalIndices[k] << ", ";
  //       }
  //       std::cout<< std::endl;

  //       std::cout<< "  faceGhostsToSend: ";
  //       for( localIndex const & k : faceGhostsToSend )
  //       {
  //         std::cout << faceGlobalIndices[k] << ", ";
  //       }
  //       std::cout<< std::endl;

  //       std::cout<< "  faceGhostsToRecv: ";
  //       for( localIndex const & k : faceGhostsToRecv )
  //       {
  //         std::cout << faceGlobalIndices[k] << ", ";
  //       }
  //       std::cout<< std::endl;

  //     }
  //   }
  //   MpiWrapper::barrier();
  // }



  elemManager.forElementSubRegionsComplete< FaceElementSubRegion >( [&]( localIndex const er,
                                                                         localIndex const esr,
                                                                         ElementRegionBase const &,
                                                                         FaceElementSubRegion & subRegion )
  {
    updateConnectorsToFaceElems( receivedObjects.newElements.at( {er, esr} ), subRegion );
  } );


  std::set< localIndex > allTouchedNodes;
  allTouchedNodes.insert( modifiedObjects.newNodes.begin(), modifiedObjects.newNodes.end() );
  allTouchedNodes.insert( modifiedObjects.modifiedNodes.begin(), modifiedObjects.modifiedNodes.end() );
  nodeManager.depopulateUpMaps( allTouchedNodes,
                                edgeManager.nodeList(),
                                faceManager.nodeList().toViewConst(),
                                elemManager );

  std::set< localIndex > allTouchedEdges;
  allTouchedEdges.insert( modifiedObjects.newEdges.begin(), modifiedObjects.newEdges.end() );
  allTouchedEdges.insert( modifiedObjects.modifiedEdges.begin(), modifiedObjects.modifiedEdges.end() );
  edgeManager.depopulateUpMaps( allTouchedEdges,
                                faceManager.edgeList().toViewConst() );

  std::set< localIndex > allTouchedFaces;
  allTouchedFaces.insert( modifiedObjects.newFaces.begin(), modifiedObjects.newFaces.end() );
  allTouchedFaces.insert( modifiedObjects.modifiedFaces.begin(), modifiedObjects.modifiedFaces.end() );
  faceManager.depopulateUpMaps( allTouchedFaces, elemManager );

  nodeManager.enforceStateFieldConsistencyPostTopologyChange( modifiedObjects.modifiedNodes );
  edgeManager.enforceStateFieldConsistencyPostTopologyChange( modifiedObjects.modifiedEdges );
  faceManager.enforceStateFieldConsistencyPostTopologyChange( modifiedObjects.modifiedFaces );


}

}



} /* namespace geos */
#endif // PARALLEL_TOPOLOGY_CHANGE_METHOD==1
