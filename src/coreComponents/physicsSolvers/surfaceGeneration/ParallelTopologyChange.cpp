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

#if PARALLEL_TOPOLOGY_CHANGE_METHOD==0
namespace geos
{

using namespace dataRepository;

#define OWNING_RANK_UNPACK_METHOD 1
#define GHOST_RANK_UNPACK_METHOD 1

namespace
{

void print( ModifiedObjectLists const & objectList )
{
  for( int rank=0; rank<MpiWrapper::commSize(); ++rank )
  {
    MpiWrapper::barrier();
    if( rank==MpiWrapper::commRank() )
    {
      std::cout<<"rank "<<rank<<std::endl;

      std::cout<<"  newNodes: ";
      for( localIndex newNode : objectList.newNodes ) { std::cout<<newNode<<" "; }
      std::cout<<std::endl;

      std::cout<<"  newEdges: ";
      for( localIndex newEdge : objectList.newEdges ) { std::cout<<newEdge<<" "; }
      std::cout<<std::endl;

      std::cout<<"  newFaces: ";
      for( localIndex newFace : objectList.newFaces ) { std::cout<<newFace<<" "; }
      std::cout<<std::endl;

      std::cout<<"  newElements: "<<std::endl;
      for( auto newElementIter=objectList.newElements.begin(); newElementIter!=objectList.newElements.end(); ++newElementIter )
      {
        std::cout<<"             : ( "<<newElementIter->first.first<<", "<<newElementIter->first.second<<" ): ";
        for( auto elemIndex : newElementIter->second )
        {
          std::cout<<elemIndex<<", ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;

      std::cout<<"  modifiedNodes: ";
      for( localIndex modifiedNode : objectList.modifiedNodes ) { std::cout<<modifiedNode<<" "; }
      std::cout<<std::endl;

      std::cout<<"  modifiedEdges: ";
      for( localIndex modifiedEdge : objectList.modifiedEdges ) { std::cout<<modifiedEdge<<" "; }
      std::cout<<std::endl;

      std::cout<<"  modifiedFaces: ";
      for( localIndex modifiedFace : objectList.modifiedFaces ) { std::cout<<modifiedFace<<" "; }
      std::cout<<std::endl;

      std::cout<<"  modifiedElements: "<<std::endl;
      for( auto modifiedElementIter=objectList.modifiedElements.begin(); modifiedElementIter!=objectList.modifiedElements.end(); ++modifiedElementIter )
      {
        std::cout<<"             : ( "<<modifiedElementIter->first.first<<", "<<modifiedElementIter->first.second<<" ): ";
        for( auto elemIndex : modifiedElementIter->second )
        {
          std::cout<<elemIndex<<", ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;


    }

  }
}




void packNewAndModifiedObjectsToOwningRanks(  NeighborCommunicator & neighbor,
                                             MeshLevel * const meshLevel,
                                             ModifiedObjectLists const & modifiedObjects,
                                             int const commID )
{
  int bufferSize = 0;

  NodeManager & nodeManager = meshLevel->getNodeManager();
  EdgeManager & edgeManager = meshLevel->getEdgeManager();
  FaceManager & faceManager = meshLevel->getFaceManager();
  ElementRegionManager & elemManager = meshLevel->getElemManager();

  arrayView1d< integer > const & nodeGhostRank = nodeManager.ghostRank();
  arrayView1d< integer > const & edgeGhostRank = edgeManager.ghostRank();
  arrayView1d< integer > const & faceGhostRank = faceManager.ghostRank();

  arrayView1d< localIndex > const & parentNodeIndices = nodeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex > const & parentEdgeIndices = edgeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex > const & parentFaceIndices = faceManager.getField< fields::parentIndex >();

  int const neighborRank = neighbor.neighborRank();

  array1d< localIndex > newNodePackListArray( modifiedObjects.newNodes.size());
  {
    localIndex a=0;
    for( auto const index : modifiedObjects.newNodes )
    {
      localIndex const parentNodeIndex = ObjectManagerBase::getParentRecursive( parentNodeIndices, index );
      if( nodeGhostRank[parentNodeIndex] == neighborRank )
      {
        newNodePackListArray[a] = index;
        ++a;
      }
    }
    newNodePackListArray.resize( a );
  }
  array1d< localIndex > modNodePackListArray( modifiedObjects.modifiedNodes.size());
  {
    localIndex a=0;
    for( auto const index : modifiedObjects.modifiedNodes )
    {
      localIndex const parentNodeIndex = ObjectManagerBase::getParentRecursive( parentNodeIndices, index );
      if( nodeGhostRank[parentNodeIndex] == neighborRank )
      {
        modNodePackListArray[a] = index;
        ++a;
      }
    }
    modNodePackListArray.resize( a );
  }


  array1d< localIndex > newEdgePackListArray( modifiedObjects.newEdges.size());
  {
    localIndex a=0;
    for( auto const index : modifiedObjects.newEdges )
    {
      localIndex const parentIndex = ObjectManagerBase::getParentRecursive( parentEdgeIndices, index );
      if( edgeGhostRank[parentIndex] == neighborRank )
      {
        newEdgePackListArray[a] = index;
        ++a;
      }
    }
    newEdgePackListArray.resize( a );
  }
  array1d< localIndex > modEdgePackListArray( modifiedObjects.modifiedEdges.size());
  {
    localIndex a=0;
    for( auto const index : modifiedObjects.modifiedEdges )
    {
      localIndex const parentIndex = ObjectManagerBase::getParentRecursive( parentEdgeIndices, index );
      if( edgeGhostRank[parentIndex] == neighborRank )
      {
        modEdgePackListArray[a] = index;
        ++a;
      }
    }
    modEdgePackListArray.resize( a );
  }


  array1d< localIndex > newFacePackListArray( modifiedObjects.newFaces.size());
  {
    localIndex a=0;
    for( auto const index : modifiedObjects.newFaces )
    {
      localIndex const parentIndex = ObjectManagerBase::getParentRecursive( parentFaceIndices, index );
      if( faceGhostRank[parentIndex] == neighborRank )
      {
        newFacePackListArray[a] = index;
        ++a;
      }
    }
    newFacePackListArray.resize( a );
  }
  array1d< localIndex > modFacePackListArray( modifiedObjects.modifiedFaces.size());
  {
    localIndex a=0;
    for( auto const index : modifiedObjects.modifiedFaces )
    {
      localIndex const parentIndex = ObjectManagerBase::getParentRecursive( parentFaceIndices, index );
      if( faceGhostRank[parentIndex] == neighborRank )
      {
        modFacePackListArray[a] = index;
        ++a;
      }
    }
    modFacePackListArray.resize( a );
  }


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

        modElemData[er][esr].resize( elemList.size() );

        localIndex a=0;
        for( auto const index : elemList )
        {
          if( subRegionGhostRank[index] == neighborRank )
          {
            modElemData[er][esr][a] = index;
            ++a;
          }
        }
        modElemData[er][esr].resize( a );
      }

      if( modifiedObjects.newElements.count( {er, esr} ) > 0 )
      {
        std::set< localIndex > const & elemList = modifiedObjects.newElements.at( {er, esr} );

        newElemData[er][esr].resize( elemList.size() );

        localIndex a=0;
        for( auto const index : elemList )
        {
          if( subRegionGhostRank[index] == neighborRank )
          {
            newElemData[er][esr][a] = index;
            ++a;
          }
        }
        newElemData[er][esr].resize( a );
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

  bufferSize += nodeManager.packUpDownMapsSize( newNodePackListArray );
  bufferSize += edgeManager.packUpDownMapsSize( newEdgePackListArray );
  bufferSize += faceManager.packUpDownMapsSize( newFacePackListArray );
  bufferSize += elemManager.packUpDownMapsSize( newElemPackList );

  bufferSize += nodeManager.packParentChildMapsSize( newNodePackListArray );
  bufferSize += edgeManager.packParentChildMapsSize( newEdgePackListArray );
  bufferSize += faceManager.packParentChildMapsSize( newFacePackListArray );

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

  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, newNodePackListArray );
  packedSize += edgeManager.packUpDownMaps( sendBufferPtr, newEdgePackListArray );
  packedSize += faceManager.packUpDownMaps( sendBufferPtr, newFacePackListArray );
  packedSize += elemManager.packUpDownMaps( sendBufferPtr, newElemPackList );

  packedSize += nodeManager.packParentChildMaps( sendBufferPtr, newNodePackListArray );
  packedSize += edgeManager.packParentChildMaps( sendBufferPtr, newEdgePackListArray );
  packedSize += faceManager.packParentChildMaps( sendBufferPtr, newFacePackListArray );

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
#if OWNING_RANK_UNPACK_METHOD==0
localIndex unpackNewAndModifiedObjectsOnOwningRanks( NeighborCommunicator * const neighbor,
                                                     MeshLevel * const mesh,
                                                     int const commID,
//                                          array1d<array1d< std::set<localIndex> > > & allNewElements,
//                                          array1d<array1d< std::set<localIndex> > > & allModifiedElements,
                                                     ModifiedObjectLists & receivedObjects )
{
  GEOS_MARK_FUNCTION;

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();



  buffer_type const & receiveBuffer = neighbor->receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();
  localIndex_array newLocalNodes, modifiedLocalNodes;
  localIndex_array newLocalEdges, modifiedLocalEdges;
  localIndex_array newLocalFaces, modifiedLocalFaces;

  ElementRegionManager::ElementReferenceAccessor< array1d< localIndex > > newLocalElements;
  array1d< array1d< localIndex_array > > newLocalElementsData;

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modifiedLocalElements;
  array1d< array1d< localIndex_array > > modifiedLocalElementsData;

  newLocalElements.resize( elemManager.numRegions());
  newLocalElementsData.resize( elemManager.numRegions());
  modifiedLocalElements.resize( elemManager.numRegions());
  modifiedLocalElementsData.resize( elemManager.numRegions());
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newLocalElements[er].resize( elemRegion.numSubRegions());
    newLocalElementsData[er].resize( elemRegion.numSubRegions());
    modifiedLocalElements[er].resize( elemRegion.numSubRegions());
    modifiedLocalElementsData[er].resize( elemRegion.numSubRegions());
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newLocalElements[er][esr].set( newLocalElementsData[er][esr] );
      modifiedLocalElements[er][esr].set( modifiedLocalElementsData[er][esr] );
    }
  }

  // if we move to device + async packing here, add polling of events or pass out
  parallelDeviceEvents events;
  int unpackedSize = 0;
  unpackedSize += nodeManager.unpackGlobalMaps( receiveBufferPtr, newLocalNodes, 0 );
  unpackedSize += edgeManager.unpackGlobalMaps( receiveBufferPtr, newLocalEdges, 0 );
  unpackedSize += faceManager.unpackGlobalMaps( receiveBufferPtr, newLocalFaces, 0 );
  unpackedSize += elemManager.unpackGlobalMaps( receiveBufferPtr, newLocalElements );

  unpackedSize += nodeManager.unpackUpDownMaps( receiveBufferPtr, newLocalNodes, true, true );
  unpackedSize += edgeManager.unpackUpDownMaps( receiveBufferPtr, newLocalEdges, true, true );
  unpackedSize += faceManager.unpackUpDownMaps( receiveBufferPtr, newLocalFaces, true, true );
  unpackedSize += elemManager.unpackUpDownMaps( receiveBufferPtr, newLocalElements, true );

  unpackedSize += nodeManager.unpackParentChildMaps( receiveBufferPtr, newLocalNodes );
  unpackedSize += edgeManager.unpackParentChildMaps( receiveBufferPtr, newLocalEdges );
  unpackedSize += faceManager.unpackParentChildMaps( receiveBufferPtr, newLocalFaces );

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
//    unpackedSize += elemManager.Unpack( receiveBufferPtr, modifiedElements );

  waitAllDeviceEvents( events );

  std::set< localIndex > & allNewNodes      = receivedObjects.newNodes;
  std::set< localIndex > & allModifiedNodes = receivedObjects.modifiedNodes;
  std::set< localIndex > & allNewEdges      = receivedObjects.newEdges;
  std::set< localIndex > & allModifiedEdges = receivedObjects.modifiedEdges;
  std::set< localIndex > & allNewFaces      = receivedObjects.newFaces;
  std::set< localIndex > & allModifiedFaces = receivedObjects.modifiedFaces;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > & allNewElements = receivedObjects.newElements;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > & allModifiedElements = receivedObjects.modifiedElements;

  allNewNodes.insert( newLocalNodes.begin(), newLocalNodes.end() );
  allModifiedNodes.insert( modifiedLocalNodes.begin(), modifiedLocalNodes.end() );

  allNewEdges.insert( newLocalEdges.begin(), newLocalEdges.end() );
  allModifiedEdges.insert( modifiedLocalEdges.begin(), modifiedLocalEdges.end() );

  allNewFaces.insert( newLocalFaces.begin(), newLocalFaces.end() );
  allModifiedFaces.insert( modifiedLocalFaces.begin(), modifiedLocalFaces.end() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
    {
      allNewElements[{er, esr}].insert( newLocalElements[er][esr].get().begin(),
                                        newLocalElements[er][esr].get().end() );
      allModifiedElements[{er, esr}].insert( modifiedLocalElements[er][esr].get().begin(),
                                             modifiedLocalElements[er][esr].get().end() );
    }
  }

  return unpackedSize;
}
#else
localIndex unpackNewObjectsOnOwningRanks(  NeighborCommunicator & neighbor,
                                                     MeshLevel * const mesh,
                                                     int const commID,
                                                     ModifiedObjectLists & receivedObjects )
{
  GEOS_MARK_FUNCTION;

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  buffer_type const & receiveBuffer = neighbor.receiveBuffer( commID );
  neighbor.m_receiveBufferPtr = receiveBuffer.data();
  buffer_unit_type const * & receiveBufferPtr = neighbor.m_receiveBufferPtr;

  localIndex_array & newLocalNodes = neighbor.m_newLocalNodes;
  localIndex_array & newLocalEdges = neighbor.m_newLocalEdges;
  localIndex_array & newLocalFaces = neighbor.m_newLocalFaces;

  ElementRegionManager::ElementReferenceAccessor< array1d< localIndex > > & newLocalElements = neighbor.m_newLocalElements;
  array1d< array1d< localIndex_array > > & newLocalElementsData = neighbor.m_newLocalElementsData;

  newLocalNodes.resize(0);
  newLocalEdges.resize(0);
  newLocalFaces.resize(0);


  newLocalElements.resize( elemManager.numRegions());
  newLocalElementsData.resize( elemManager.numRegions());
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newLocalElements[er].resize( elemRegion.numSubRegions());
    newLocalElementsData[er].resize( elemRegion.numSubRegions());
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newLocalElementsData[er][esr].resize(0);
      newLocalElements[er][esr].set( newLocalElementsData[er][esr] );
    }
  }

  // if we move to device + async packing here, add polling of events or pass out
  parallelDeviceEvents events;
  int unpackedSize = 0;
  unpackedSize += nodeManager.unpackGlobalMaps( receiveBufferPtr, newLocalNodes, 0 );
  unpackedSize += edgeManager.unpackGlobalMaps( receiveBufferPtr, newLocalEdges, 0 );
  unpackedSize += faceManager.unpackGlobalMaps( receiveBufferPtr, newLocalFaces, 0 );
  unpackedSize += elemManager.unpackGlobalMaps( receiveBufferPtr, newLocalElements );

  waitAllDeviceEvents( events );

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

localIndex unpackNewAndModifiedObjectsDataOnOwningRanks(  NeighborCommunicator & neighbor,
                                                     MeshLevel * const mesh,
                                                     int const commID,
                                                     ModifiedObjectLists & receivedObjects )
{
  GEOS_MARK_FUNCTION;

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  buffer_unit_type const * & receiveBufferPtr = neighbor.m_receiveBufferPtr;

  localIndex_array & newLocalNodes = neighbor.m_newLocalNodes;
  localIndex_array & newLocalEdges = neighbor.m_newLocalEdges;
  localIndex_array & newLocalFaces = neighbor.m_newLocalFaces;

  ElementRegionManager::ElementReferenceAccessor< array1d< localIndex > > & newLocalElements = neighbor.m_newLocalElements;
  array1d< array1d< localIndex_array > > & newLocalElementsData = neighbor.m_newLocalElementsData;

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
  int unpackedSize = 0;

  unpackedSize += nodeManager.unpackUpDownMaps( receiveBufferPtr, newLocalNodes, true, true );
  unpackedSize += edgeManager.unpackUpDownMaps( receiveBufferPtr, newLocalEdges, true, true );
  unpackedSize += faceManager.unpackUpDownMaps( receiveBufferPtr, newLocalFaces, true, true );
  unpackedSize += elemManager.unpackUpDownMaps( receiveBufferPtr, newLocalElements, true );

  unpackedSize += nodeManager.unpackParentChildMaps( receiveBufferPtr, newLocalNodes );
  unpackedSize += edgeManager.unpackParentChildMaps( receiveBufferPtr, newLocalEdges );
  unpackedSize += faceManager.unpackParentChildMaps( receiveBufferPtr, newLocalFaces );

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

  std::set< localIndex > & allNewNodes      = receivedObjects.newNodes;
  std::set< localIndex > & allModifiedNodes = receivedObjects.modifiedNodes;
  std::set< localIndex > & allNewEdges      = receivedObjects.newEdges;
  std::set< localIndex > & allModifiedEdges = receivedObjects.modifiedEdges;
  std::set< localIndex > & allNewFaces      = receivedObjects.newFaces;
  std::set< localIndex > & allModifiedFaces = receivedObjects.modifiedFaces;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > & allNewElements = receivedObjects.newElements;
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
#endif

void FilterNewObjectsForPackToGhosts( std::set< localIndex > const & objectList,
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

void FilterModObjectsForPackToGhosts( std::set< localIndex > const & objectList,
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

void packNewModifiedObjectsToGhosts(  NeighborCommunicator & neighbor,
                                     int commID,
                                     MeshLevel * const mesh,
                                     ModifiedObjectLists & receivedObjects )
{
  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  localIndex_array newNodesToSend;
  localIndex_array newEdgesToSend;
  localIndex_array newFacesToSend;
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > newElemsToSend;
  array1d< array1d< localIndex_array > > newElemsToSendData;

  localIndex_array modNodesToSend;
  localIndex_array modEdgesToSend;
  localIndex_array modFacesToSend;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modElemsToSend;
  array1d< array1d< localIndex_array > > modElemsToSendData;

  localIndex_array & nodeGhostsToSend = nodeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  localIndex_array & edgeGhostsToSend = edgeManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
  localIndex_array & faceGhostsToSend = faceManager.getNeighborData( neighbor.neighborRank() ).ghostsToSend();

  arrayView1d< localIndex > const & nodalParentIndices = nodeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex > const & edgeParentIndices = edgeManager.getField< fields::parentIndex >();
  arrayView1d< localIndex > const & faceParentIndices = faceManager.getField< fields::parentIndex >();

  FilterNewObjectsForPackToGhosts( receivedObjects.newNodes, nodalParentIndices, nodeGhostsToSend, newNodesToSend );
  FilterModObjectsForPackToGhosts( receivedObjects.modifiedNodes, nodeGhostsToSend, modNodesToSend );

  FilterNewObjectsForPackToGhosts( receivedObjects.newEdges, edgeParentIndices, edgeGhostsToSend, newEdgesToSend );
  FilterModObjectsForPackToGhosts( receivedObjects.modifiedEdges, edgeGhostsToSend, modEdgesToSend );

  FilterNewObjectsForPackToGhosts( receivedObjects.newFaces, faceParentIndices, faceGhostsToSend, newFacesToSend );
  FilterModObjectsForPackToGhosts( receivedObjects.modifiedFaces, faceGhostsToSend, modFacesToSend );

  SortedArray< localIndex > faceGhostsToSendSet;
  for( localIndex const & kf : faceGhostsToSend )
  {
    faceGhostsToSendSet.insert( kf );
  }

  newElemsToSendData.resize( elemManager.numRegions() );
  newElemsToSend.resize( elemManager.numRegions() );
  modElemsToSendData.resize( elemManager.numRegions() );
  modElemsToSend.resize( elemManager.numRegions() );
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newElemsToSendData[er].resize( elemRegion.numSubRegions() );
    newElemsToSend[er].resize( elemRegion.numSubRegions() );
    modElemsToSendData[er].resize( elemRegion.numSubRegions() );
    modElemsToSend[er].resize( elemRegion.numSubRegions() );

    elemRegion.forElementSubRegionsIndex< FaceElementSubRegion >( [&]( localIndex const esr,
                                                                       FaceElementSubRegion & subRegion )
    {
      ArrayOfArraysView< localIndex const > const faceList = subRegion.faceList().toViewConst();
      localIndex_array & elemGhostsToSend = subRegion.getNeighborData( neighbor.neighborRank() ).ghostsToSend();
      elemGhostsToSend.move( hostMemorySpace );
      for( localIndex const & k : receivedObjects.newElements.at( {er, esr} ) )
      {
        if( faceGhostsToSendSet.count( faceList( k, 0 ) ) )
        {
          newElemsToSendData[er][esr].emplace_back( k );
          elemGhostsToSend.emplace_back( k );
        }
      }
      newElemsToSend[er][esr] = newElemsToSendData[er][esr];
    } );

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

  bufferSize += nodeManager.packGlobalMapsSize( newNodesToSend, 0 );
  bufferSize += edgeManager.packGlobalMapsSize( newEdgesToSend, 0 );
  bufferSize += faceManager.packGlobalMapsSize( newFacesToSend, 0 );
  bufferSize += elemManager.packGlobalMapsSize( newElemsToSend );

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

  packedSize += nodeManager.packGlobalMaps( sendBufferPtr, newNodesToSend, 0 );
  packedSize += edgeManager.packGlobalMaps( sendBufferPtr, newEdgesToSend, 0 );
  packedSize += faceManager.packGlobalMaps( sendBufferPtr, newFacesToSend, 0 );
  packedSize += elemManager.packGlobalMaps( sendBufferPtr, newElemsToSend );

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

#if GHOST_RANK_UNPACK_METHOD==0
void unpackNewModToGhosts( NeighborCommunicator * const neighbor,
                           int commID,
                           MeshLevel * const mesh,
                           ModifiedObjectLists & receivedObjects )
{

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  localIndex_array & nodeGhostsToRecv = nodeManager.getNeighborData( neighbor->neighborRank() ).ghostsToReceive();

  localIndex_array & edgeGhostsToRecv = edgeManager.getNeighborData( neighbor->neighborRank() ).ghostsToReceive();

  localIndex_array & faceGhostsToRecv = faceManager.getNeighborData( neighbor->neighborRank() ).ghostsToReceive();

  buffer_type const & receiveBuffer = neighbor->receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  localIndex_array newGhostNodes;
  localIndex_array newGhostEdges;
  localIndex_array newGhostFaces;

  localIndex_array modGhostNodes;
  localIndex_array modGhostEdges;
  localIndex_array modGhostFaces;

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > newGhostElems;
  array1d< array1d< localIndex_array > > newGhostElemsData;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modGhostElems;
  array1d< array1d< localIndex_array > > modGhostElemsData;
  newGhostElems.resize( elemManager.numRegions() );
  newGhostElemsData.resize( elemManager.numRegions() );
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

  // if we move to device + async unoacking, poll these events for completion or pass out
  parallelDeviceEvents events;

  nodeManager.unpackGlobalMaps( receiveBufferPtr, newGhostNodes, 0 );
  edgeManager.unpackGlobalMaps( receiveBufferPtr, newGhostEdges, 0 );
  faceManager.unpackGlobalMaps( receiveBufferPtr, newGhostFaces, 0 );
  elemManager.unpackGlobalMaps( receiveBufferPtr, newGhostElems );

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
    localIndex_array & elemGhostsToReceive = subRegion.getNeighborData( neighbor->neighborRank() ).ghostsToReceive();

    if( newGhostElemsData[er][esr].size() > 0 )
    {
      elemGhostsToReceive.move( hostMemorySpace );

      for( localIndex const & newElemIndex : newGhostElemsData[er][esr] )
      {
        elemGhostsToReceive.emplace_back( newElemIndex );
        receivedObjects.newElements[ { er, esr } ].insert( newElemIndex );
      }
    }

    receivedObjects.modifiedElements[ { er, esr } ].insert( modGhostElemsData[er][esr].begin(),
                                                            modGhostElemsData[er][esr].end() );
  } );

  receivedObjects.newNodes.insert( newGhostNodes.begin(), newGhostNodes.end() );
  receivedObjects.modifiedNodes.insert( modGhostNodes.begin(), modGhostNodes.end() );
  receivedObjects.newEdges.insert( newGhostEdges.begin(), newGhostEdges.end() );
  receivedObjects.modifiedEdges.insert( modGhostEdges.begin(), modGhostEdges.end() );
  receivedObjects.newFaces.insert( newGhostFaces.begin(), newGhostFaces.end() );
  receivedObjects.modifiedFaces.insert( modGhostFaces.begin(), modGhostFaces.end() );

}
#else
void unpackNewObjectsToGhosts(  NeighborCommunicator & neighbor,
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
  neighbor.m_receiveBufferPtr = receiveBuffer.data();
  buffer_unit_type const * & receiveBufferPtr = neighbor.m_receiveBufferPtr;

  localIndex_array & newGhostNodes = neighbor.m_newGhostNodes;
  localIndex_array & newGhostEdges = neighbor.m_newGhostEdges;
  localIndex_array & newGhostFaces = neighbor.m_newGhostFaces;

  newGhostNodes.resize(0);
  newGhostEdges.resize(0);
  newGhostFaces.resize(0);

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > & newGhostElems = neighbor.m_newGhostElems;
  array1d< array1d< localIndex_array > > & newGhostElemsData = neighbor.m_newGhostElemsData;
  newGhostElems.resize( elemManager.numRegions() );
  newGhostElemsData.resize( elemManager.numRegions() );
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newGhostElemsData[er].resize( elemRegion.numSubRegions() );
    newGhostElems[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newGhostElemsData[er][esr].resize(0);
      newGhostElems[er][esr].set( newGhostElemsData[er][esr] );
    }
  }

  // if we move to device + async unoacking, poll these events for completion or pass out
  parallelDeviceEvents events;

  nodeManager.unpackGlobalMaps( receiveBufferPtr, newGhostNodes, 0 );
  edgeManager.unpackGlobalMaps( receiveBufferPtr, newGhostEdges, 0 );
  faceManager.unpackGlobalMaps( receiveBufferPtr, newGhostFaces, 0 );
  elemManager.unpackGlobalMaps( receiveBufferPtr, newGhostElems );

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








void unpackNewAndModifiedObjectsDataToGhosts(  NeighborCommunicator & neighbor,
                           int commID,
                           MeshLevel * const mesh,
                           ModifiedObjectLists & receivedObjects )
{

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  localIndex_array & nodeGhostsToRecv = nodeManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();
\
  localIndex_array & edgeGhostsToRecv = edgeManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();

  localIndex_array & faceGhostsToRecv = faceManager.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();

  buffer_type const & receiveBuffer = neighbor.receiveBuffer( commID );
  buffer_unit_type const * & receiveBufferPtr = neighbor.m_receiveBufferPtr;

  localIndex_array & newGhostNodes = neighbor.m_newGhostNodes;
  localIndex_array & newGhostEdges = neighbor.m_newGhostEdges;
  localIndex_array & newGhostFaces = neighbor.m_newGhostFaces;

  localIndex_array modGhostNodes;
  localIndex_array modGhostEdges;
  localIndex_array modGhostFaces;

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > & newGhostElems = neighbor.m_newGhostElems;
  array1d< array1d< localIndex_array > > & newGhostElemsData = neighbor.m_newGhostElemsData;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > modGhostElems;
  array1d< array1d< localIndex_array > > modGhostElemsData;
  modGhostElems.resize( elemManager.numRegions() );
  modGhostElemsData.resize( elemManager.numRegions() );
  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    modGhostElemsData[er].resize( elemRegion.numSubRegions() );
    modGhostElems[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newGhostElems[er][esr].set( newGhostElemsData[er][esr] );
      modGhostElems[er][esr].set( modGhostElemsData[er][esr] );
    }
  }

  // if we move to device + async unoacking, poll these events for completion or pass out
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
    [&]( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    localIndex_array & elemGhostsToReceive = subRegion.getNeighborData( neighbor.neighborRank() ).ghostsToReceive();

    receivedObjects.modifiedElements[ { er, esr } ].insert( modGhostElemsData[er][esr].begin(),
                                                            modGhostElemsData[er][esr].end() );
  } );

  receivedObjects.modifiedNodes.insert( modGhostNodes.begin(), modGhostNodes.end() );
  receivedObjects.modifiedEdges.insert( modGhostEdges.begin(), modGhostEdges.end() );
  receivedObjects.modifiedFaces.insert( modGhostFaces.begin(), modGhostFaces.end() );

}
#endif







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

}

void parallelTopologyChange::synchronizeTopologyChange( MeshLevel * const mesh,
                                                        std::vector< NeighborCommunicator > & neighbors,
                                                        ModifiedObjectLists & modifiedObjects,
                                                        ModifiedObjectLists & receivedObjects,
                                                        int mpiCommOrder )
{

  NodeManager & nodeManager = mesh->getNodeManager();
  EdgeManager & edgeManager = mesh->getEdgeManager();
  FaceManager & faceManager = mesh->getFaceManager();
  ElementRegionManager & elemManager = mesh->getElemManager();

  //************************************************************************************************
  // 2) first we need to send over:
  //   a) the new objects to owning ranks. New objects are assumed to be
  //      owned by the rank that owned the parent.
  //   b) the modified objects to the owning ranks.

  // pack the buffers, and send the size of the buffers
  MPI_iCommData commData;
  commData.resize( neighbors.size() );
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewAndModifiedObjectsToOwningRanks( neighbor,
                                            mesh,
                                            modifiedObjects,
                                            commData.commID() );

    neighbor.mpiISendReceiveBufferSizes( commData.commID(),
                                         commData.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOS );

  }

  // send/recv the buffers
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
                                     MPI_COMM_GEOS );
  }

  // unpack the buffers and get lists of the new objects.
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
    // Unpack buffers in set ordering for integration testing
    else
    {
      MpiWrapper::wait( commData.mpiRecvBufferRequest() + count,
                        commData.mpiRecvBufferStatus() + count );
    }

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

#if OWNING_RANK_UNPACK_METHOD==0
    unpackNewAndModifiedObjectsOnOwningRanks( &neighbor,
                                              mesh,
                                              commData.commID(),
                                              receivedObjects );
#else
    unpackNewObjectsOnOwningRanks( neighbor,
                                              mesh,
                                              commData.commID(),
                                              receivedObjects );
#endif
  }

#if OWNING_RANK_UNPACK_METHOD==1
  for( NeighborCommunicator & neighbor : neighbors )
  {
    unpackNewAndModifiedObjectsDataOnOwningRanks( neighbor,
                                                  mesh,
                                                  commData.commID(),
                                                  receivedObjects );
  }
#endif

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

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferSizeRequest(),
                       commData.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData.size(),
                       commData.mpiSendBufferRequest(),
                       commData.mpiSendBufferSizeStatus() );

  modifiedObjects.insert( receivedObjects );


  //************************************************************************************************
  // 3) now we need to send over:
  //   a) the new objects whose parents are ghosts on neighbors.
  //   b) the modified objects whose parents are ghosts on neighbors.



  MPI_iCommData commData2;
  commData2.resize( neighbors.size());
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewModifiedObjectsToGhosts( neighbor,
                                    commData2.commID(),
                                    mesh,
                                    modifiedObjects );

    neighbor.mpiISendReceiveBufferSizes( commData2.commID(),
                                         commData2.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData2.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOS );

  }

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

    NeighborCommunicator & neighbor = neighbors[neighborIndex];


#if GHOST_RANK_UNPACK_METHOD==0
    unpackNewModToGhosts( &neighbor, commData2.commID(), mesh, receivedObjects );
#else
    unpackNewObjectsToGhosts( neighbor, commData2.commID(), mesh, receivedObjects );
#endif
  }

#if GHOST_RANK_UNPACK_METHOD==1
  for( NeighborCommunicator & neighbor : neighbors )
  {
    unpackNewAndModifiedObjectsDataToGhosts( neighbor, commData2.commID(), mesh, receivedObjects );
  }
#endif

  modifiedObjects.insert( receivedObjects );
  print( modifiedObjects );
  auto faceToEdge = faceManager.edgeList();
  for( int rank=0; rank<MpiWrapper::commSize(); ++rank )
  {
    MpiWrapper::barrier();
    if( rank==MpiWrapper::commRank() )
    {
      std::cout<<"rank "<<rank<<std::endl;
      std::cout<<faceToEdge<<std::endl;
    }
  }


  nodeManager.fixUpDownMaps( false );
  edgeManager.fixUpDownMaps( false );
  faceManager.fixUpDownMaps( false );


  for( localIndex er = 0; er < elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
    {
      elemRegion.getSubRegion( esr ).fixUpDownMaps( false );
    }
  }

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

  MpiWrapper::waitAll( commData2.size(),
                       commData2.mpiSendBufferSizeRequest(),
                       commData2.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( commData2.size(),
                       commData2.mpiSendBufferRequest(),
                       commData2.mpiSendBufferSizeStatus() );

}



} /* namespace geos */
#endif // PARALLEL_TOPOLOGY_CHANGE_METHOD==0
