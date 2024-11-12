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

#include "mesh/mpiCommunications/CommunicationTools.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <algorithm>

namespace geos
{


CommunicationTools * CommunicationTools::m_instance = nullptr;

using namespace dataRepository;

CommunicationTools::CommunicationTools()
{
  for( int i = 0; i < NeighborCommunicator::maxComm; ++i )
  {
    m_freeCommIDs.insert( i );
  }

  GEOS_ERROR_IF( m_instance != nullptr, "Only one CommunicationTools can exist at a time." );
  m_instance = this;
}

CommunicationTools::~CommunicationTools()
{
  GEOS_ERROR_IF( m_instance != this, "m_instance != this should not be possible." );
  m_instance = nullptr;
}

CommunicationTools & CommunicationTools::getInstance()
{
  GEOS_ERROR_IF( m_instance == nullptr,
                 "CommunicationTools has not been constructed, or is already been destructed." );
  return *m_instance;
}


void CommunicationTools::assignGlobalIndices( ObjectManagerBase & manager,
                                              NodeManager const & compositionManager,
                                              std::vector< NeighborCommunicator > & neighbors )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< integer > const & ghostRank = manager.ghostRank();
  ghostRank.setValues< serialPolicy >( -2 );

  int const commRank = MpiWrapper::commRank();

  localIndex const numberOfObjectsHere = manager.size();
  globalIndex const offset = MpiWrapper::prefixSum< globalIndex >( numberOfObjectsHere );

  arrayView1d< globalIndex > const localToGlobal = manager.localToGlobalMap();

  // set the global indices as if they were all local to this process
  for( localIndex a = 0; a < manager.size(); ++a )
  {
    localToGlobal[a] = offset + a;
  }

  // Get the relation to the composition object used that will be used to identify the main object.
  // For example, a face can be identified by its nodes.
  ArrayOfSets< globalIndex > const objectToCompositionObject =
    manager.extractMapFromObjectForAssignGlobalIndexNumbers( compositionManager );

  // Now arrange the data from objectToCompositionObject into a map "indexByFirstCompositionIndex",
  // such that the key is the lowest global index of the composition object that make up this object.
  // The value of the map is a pair, with the array being the remaining composition object global indices,
  // and the second being the global index of the object itself.
  map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, localIndex > > > indexByFirstCompositionIndex;

  localIndex bufferSize = 0;
  for( localIndex a = 0; a < objectToCompositionObject.size(); ++a )
  {
    arraySlice1d< globalIndex const > const nodeList = objectToCompositionObject[a];
    if( nodeList.size() > 0 )
    {
      // fill the array with the remaining composition object global indices
      std::vector< globalIndex > tempComp( nodeList.begin() + 1, nodeList.end() );

      // push the tempComp onto the map.
      indexByFirstCompositionIndex[nodeList[0]].emplace_back( std::make_pair( std::move( tempComp ), a ) );
      bufferSize += 2 + nodeList.size();
    }
  }

  array1d< globalIndex > objectToCompositionObjectSendBuffer;
  objectToCompositionObjectSendBuffer.reserve( bufferSize );

  // put the map into a buffer
  for( localIndex a = 0; a < objectToCompositionObject.size(); ++a )
  {
    arraySlice1d< globalIndex const > const nodeList = objectToCompositionObject[a];
    if( nodeList.size() > 0 )
    {
      objectToCompositionObjectSendBuffer.emplace_back( nodeList.size() );
      objectToCompositionObjectSendBuffer.emplace_back( localToGlobal[a] );
      objectToCompositionObjectSendBuffer.insert( objectToCompositionObjectSendBuffer.size(), nodeList.begin(), nodeList.end() );
    }
  }

  integer const numNeighbors = LvArray::integerConversion< integer >( neighbors.size() );

  MPI_iCommData commData;
  commData.resize( numNeighbors );

  array1d< int > receiveBufferSizes( numNeighbors );
  array1d< array1d< globalIndex > > receiveBuffers( numNeighbors );

  int const sendSize = LvArray::integerConversion< int >( objectToCompositionObjectSendBuffer.size() );

  for( integer neighborIndex = 0; neighborIndex < numNeighbors; ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    neighbor.mpiISendReceive( &sendSize,
                              1,
                              commData.mpiSendBufferSizeRequest( neighborIndex ),
                              &(receiveBufferSizes[neighborIndex]),
                              1,
                              commData.mpiRecvBufferSizeRequest( neighborIndex ),
                              commData.commID(),
                              MPI_COMM_GEOS );
  }

  for( std::size_t count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    receiveBuffers[neighborIndex].resize( receiveBufferSizes[neighborIndex] );
    neighbor.mpiISendReceive( objectToCompositionObjectSendBuffer.data(),
                              sendSize,
                              commData.mpiSendBufferRequest( neighborIndex ),
                              receiveBuffers[neighborIndex].data(),
                              receiveBufferSizes[neighborIndex],
                              commData.mpiRecvBufferRequest( neighborIndex ),
                              commData.commID(),
                              MPI_COMM_GEOS );
  }

  // unpack the data from neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager] to
  // the local arrays

  // object to receive the neighbor data
  // this baby is an Array (for each neighbor) of maps, with the key of lowest composition index, and a value
  // containing an array containing the std::pairs of the remaining composition indices, and the globalIndex of the
  // object.
  std::vector< map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, globalIndex > > > > neighborCompositionObjects( neighbors.size() );

  for( std::size_t count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    globalIndex const * recBuffer = receiveBuffers[neighborIndex].data();
    localIndex const recBufferSize = receiveBufferSizes[neighborIndex];
    globalIndex const * const endBuffer = recBuffer + recBufferSize;
    // iterate over data that was just received
    while( recBuffer < endBuffer )
    {
      // the first thing packed was the data size for a given object
      localIndex const dataSize = LvArray::integerConversion< localIndex >( *recBuffer++ );

      // the second thing packed was the globalIndex of that object
      globalIndex const neighborGlobalIndex = *recBuffer++;

      // the global indices of the composition objects were next. they are ordered, so the lowest one is first.
      globalIndex const firstCompositionIndex = *recBuffer++;

      // the remaining composition object indices.
      std::vector< globalIndex > temp( recBuffer, recBuffer + dataSize - 1 );
      recBuffer += dataSize - 1;

      // fill neighborCompositionObjects
      neighborCompositionObjects[neighborIndex][firstCompositionIndex].emplace_back( std::move( temp ), neighborGlobalIndex );
    }

    // Set iterators to the beginning of each indexByFirstCompositionIndex,
    // and neighborCompositionObjects[neighborNum].
    auto iter_local = indexByFirstCompositionIndex.begin();
    auto iter_neighbor = neighborCompositionObjects[neighborIndex].begin();

    // now we continue the while loop as long as both of our iterators are in range.
    while( iter_local != indexByFirstCompositionIndex.end() &&
           iter_neighbor != neighborCompositionObjects[neighborIndex].end() )
    {
      // check to see if the map keys (first composition index) are the same.
      if( iter_local->first == iter_neighbor->first )
      {
        // first we loop over all local composition arrays (objects with the matched key)
        for( auto const & localObj : iter_local->second )
        {
          // and loop over all the neighbor composition arrays (objects with the matched key)
          for( auto const & neighborObj : iter_neighbor->second )
          {
            // now compare the composition arrays
            if( localObj.first == neighborObj.first )
            {
              // they are equal, so we need to overwrite the global index for the object
              if( neighborObj.second < localToGlobal[localObj.second] )
              {
                if( neighbor.neighborRank() < commRank )
                {
                  localToGlobal[localObj.second] = neighborObj.second;
                  ghostRank[localObj.second] = neighbor.neighborRank();
                }
                else
                {
                  ghostRank[localObj.second] = -1;
                }
              }

              // we should break out of the iter_local2 loop since we aren't going to find another match.
              break;
            }
          }
        }
        ++iter_local;
        ++iter_neighbor;
      }
      else if( iter_local->first < iter_neighbor->first )
      {
        ++iter_local;
      }
      else if( iter_local->first > iter_neighbor->first )
      {
        ++iter_neighbor;
      }
    }
  }

  MpiWrapper::waitAll( neighbors.size(), commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
  MpiWrapper::waitAll( neighbors.size(), commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );


  manager.constructGlobalToLocalMap();

  manager.setMaxGlobalIndex();
}

void CommunicationTools::assignNewGlobalIndices( ObjectManagerBase & manager,
                                                 std::set< localIndex > const & indexList )
{
  globalIndex const glocalIndexOffset = MpiWrapper::prefixSum< globalIndex >( indexList.size(), MPI_COMM_GEOS );

  arrayView1d< globalIndex > const & localToGlobal = manager.localToGlobalMap();

  localIndex nIndicesAssigned = 0;
  for( localIndex const newLocalIndex : indexList )
  {
    GEOS_ERROR_IF( localToGlobal[newLocalIndex] != -1,
                   "Local object " << newLocalIndex << " should be new but already has a global index "
                                   << localToGlobal[newLocalIndex] );

    localToGlobal[newLocalIndex] = manager.maxGlobalIndex() + glocalIndexOffset + nIndicesAssigned + 1;
    manager.updateGlobalToLocalMap( newLocalIndex );

    nIndicesAssigned += 1;
  }

  manager.setMaxGlobalIndex();
}

void
CommunicationTools::assignNewGlobalIndices( ElementRegionManager & elementManager,
                                            std::map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems )
{
  localIndex numberOfNewObjectsHere = 0;
  for( auto const & iter : newElems )
  {
    numberOfNewObjectsHere += LvArray::integerConversion< localIndex >( iter.second.size() );
  }

  globalIndex const glocalIndexOffset = MpiWrapper::prefixSum< globalIndex >( numberOfNewObjectsHere, MPI_COMM_GEOS );

  localIndex nIndicesAssigned = 0;
  for( auto const & iter : newElems )
  {
    localIndex const er = iter.first.first;
    localIndex const esr = iter.first.second;
    std::set< localIndex > const & indexList = iter.second;

    ElementSubRegionBase & subRegion = elementManager.getRegion( er ).getSubRegion( esr );
    arrayView1d< globalIndex > const & localToGlobal = subRegion.localToGlobalMap();

    for( localIndex const newLocalIndex : indexList )
    {
      GEOS_ERROR_IF( localToGlobal[newLocalIndex] != -1,
                     "Local object " << newLocalIndex << " should be new but already has a global index "
                                     << localToGlobal[newLocalIndex] );

      localToGlobal[newLocalIndex] = elementManager.maxGlobalIndex() + glocalIndexOffset + nIndicesAssigned + 1;
      subRegion.updateGlobalToLocalMap( newLocalIndex );

      nIndicesAssigned += 1;
    }
  }

  elementManager.setMaxGlobalIndex();
}


/**
 * @brief Exchange some @p data with all the @p neighbors. The data received from the @p neighbors is the returned by the function.
 * @tparam DATA_PROVIDER Callable that takes neighbor index @p i and returns an <tt>array1d\< globalIndex \></tt> which will be sent.
 * Note that the index @p i is the index of the @p neighbor in the @p neighbors list, not its MPI rank.
 * @param commId MPI communicator Id.
 * @param neighbors List of all the concerned neighbor communicators.
 * @param data Provides the data to be sent to each neighbor.
 * @return The data received from all the @p neighbors. Data at index @p i coming from neighbor at index @p i in the list of @p neighbors.
 */
template< class DATA_PROVIDER >
array1d< array1d< globalIndex > > exchange( std::vector< NeighborCommunicator > & neighbors,
                                            DATA_PROVIDER const & data )
{
  MPI_iCommData commData;
  int commId = commData.commID();
  integer const numNeighbors = LvArray::integerConversion< integer >( neighbors.size() );
  commData.resize( numNeighbors );
  for( integer i = 0; i < numNeighbors; ++i )
  {
    neighbors[i].mpiISendReceiveSizes( data( i ),
                                       commData.mpiSendBufferSizeRequest( i ),
                                       commData.mpiRecvBufferSizeRequest( i ),
                                       commId,
                                       MPI_COMM_GEOS );
  }

  MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
  MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus() );

  array1d< array1d< globalIndex > > output( neighbors.size() );

  for( integer i = 0; i < numNeighbors; ++i )
  {
    neighbors[i].mpiISendReceiveData( data( i ),
                                      commData.mpiSendBufferRequest( i ),
                                      output[i],
                                      commData.mpiRecvBufferRequest( i ),
                                      commId,
                                      MPI_COMM_GEOS );
  }
  MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );
  MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus() );

  return output;
}


array1d< array1d< globalIndex > >
CommunicationTools::buildNeighborPartitionBoundaryObjects( ObjectManagerBase & manager,
                                                           std::vector< NeighborCommunicator > & allNeighbors )
{
  arrayView1d< integer > const & domainBoundaryIndicator = manager.getDomainBoundaryIndicator();
  array1d< globalIndex > const globalPartitionBoundaryObjectsIndices = manager.constructGlobalListOfBoundaryObjects();

  // The data functor will return the same information, whatever the rank `i` is.
  // This is because the `exchange` function is able to send different data to different neighbors.
  // Since in our precise case, we want to send the same information, we adapt the functor that way.
  auto const data = [&]( auto GEOS_UNUSED_PARAM( i ) ) -> array1d< globalIndex > const &
  {
    return std::cref( globalPartitionBoundaryObjectsIndices );
  };
  array1d< array1d< globalIndex > > const neighborPartitionBoundaryObjects = exchange( allNeighbors, data );

  integer const numNeighbors = LvArray::integerConversion< integer >( allNeighbors.size() );
  for( integer i = 0; i < numNeighbors; ++i )
  {
    NeighborCommunicator & neighbor = allNeighbors[i];
    array1d< localIndex > & matchedPartitionBoundaryObjects = manager.getNeighborData( neighbor.neighborRank() ).matchedPartitionBoundary();

    localIndex localCounter = 0;
    localIndex neighborCounter = 0;
    while( localCounter < globalPartitionBoundaryObjectsIndices.size() &&
           neighborCounter < neighborPartitionBoundaryObjects[i].size() )
    {
      globalIndex const & gi = globalPartitionBoundaryObjectsIndices[localCounter];
      globalIndex const & gni = neighborPartitionBoundaryObjects[i][neighborCounter];

      if( gi == gni )
      {
        localIndex const localMatchedIndex = manager.globalToLocalMap( gi );
        matchedPartitionBoundaryObjects.emplace_back( localMatchedIndex );
        domainBoundaryIndicator[localMatchedIndex] = 2;  // TODO check if useful.
        ++localCounter;
        ++neighborCounter;
      }
      else if( gi > gni )
      {
        ++neighborCounter;
      }
      else
      {
        ++localCounter;
      }
    }
  }

  return neighborPartitionBoundaryObjects;
}


void
CommunicationTools::
  findMatchedPartitionBoundaryObjects( ObjectManagerBase & objectManager,
                                       std::vector< NeighborCommunicator > & allNeighbors )
{
  GEOS_MARK_FUNCTION;
  this->buildNeighborPartitionBoundaryObjects( objectManager, allNeighbors );
}

/**
 * @brief During the ghosting process, some nodes may already be on multiple ranks in the same time.
 * But if we want to send those nodes, the rank with the lowest MPI rank really owns the node,
 * and therefore should be responsible for sending the node to the others.
 * @param mpiRankToNodes For each involved mpi rank, all the nodes that need to be sent.
 * @return The sanitized mapping: only the ranks with the minimal MPI rank will be in charged of sending the proper nodes.
 */
std::map< int, array1d< globalIndex > > reorganizeRequestedNodes( std::map< int, array1d< globalIndex > > const & mpiRankToNodes )
{
  class MinInt
  {
public:
    MinInt & operator=( int other )
    {
      m_value = std::min( other, m_value );
      return *this;
    }

    operator int() const
    {
      return m_value;
    }

private:
    int m_value = std::numeric_limits< int >::max();
  };

  std::map< int, array1d< globalIndex > > minMpiRankToNodes;  // Will be returned.

  std::map< globalIndex, MinInt > nodeToMpiRank;
  for( auto const & [mpiRank, nodes]: mpiRankToNodes )
  {
    for( globalIndex const & gi: nodes )
    {
      nodeToMpiRank[gi] = mpiRank;
    }
    minMpiRankToNodes[mpiRank];  // Explicitly allocate all the ranks. Important, otherwise mpi ranks without nodes will not be defined.
  }

  for( auto const & [node, mpiRank]: nodeToMpiRank )
  {
    minMpiRankToNodes[mpiRank].emplace_back( node );
  }

  return minMpiRankToNodes;
}


void CommunicationTools::findMatchedPartitionBoundaryNodes( NodeManager & nodeManager,
                                                            std::vector< NeighborCommunicator > & allNeighbors,
                                                            std::set< std::set< globalIndex > > const & collocatedNodesBuckets,
                                                            std::set< globalIndex > const & requestedNodes )
{
  GEOS_MARK_FUNCTION;
  auto const & g2l = nodeManager.globalToLocalMap();
  integer const numNeighbors = LvArray::integerConversion< integer >( allNeighbors.size() );

  std::map< int, array1d< globalIndex > > requestedMatchesMap;  // The key of the map is the MPI rank of the `neighbor`, not the index.

  {
    array1d< array1d< globalIndex > > const neighborBoundaryNodes = this->buildNeighborPartitionBoundaryObjects( nodeManager, allNeighbors );

    for( integer i = 0; i < numNeighbors; ++i )
    {
      array1d< localIndex > & matchedBoundaryNodes = nodeManager.getNeighborData( allNeighbors[i].neighborRank() ).matchedPartitionBoundary();
      array1d< localIndex > collocatedNodesToSend;  // The collocated nodes that match neighbor nodes and that we want to send.

      std::set< globalIndex > const neighborNodes( neighborBoundaryNodes[i].begin(), neighborBoundaryNodes[i].end() );
      // This loop searches for neighbor nodes that may be the collocated with some nodes on the current domain/rank.
      // If we find any match, then we tag the node on the current domain/rank in order to send it to the neighbor that will need it.
      for( std::set< globalIndex > const & collocatedNodes: collocatedNodesBuckets )
      {
        // First step is to find if there is any matching collocated node.
        std::vector< globalIndex > intersection;
        std::set_intersection( collocatedNodes.cbegin(), collocatedNodes.cend(),
                               neighborNodes.cbegin(), neighborNodes.cend(),
                               std::back_inserter( intersection ) );
        if( intersection.empty() )
        {
          continue;
        }

        std::vector< globalIndex > nodesMissingOnNeighbor;
        std::set_difference( collocatedNodes.cbegin(), collocatedNodes.cend(),
                             intersection.cbegin(), intersection.cend(),
                             std::back_inserter( nodesMissingOnNeighbor ) );
        for( globalIndex const & gi: nodesMissingOnNeighbor )
        {
          auto it = g2l.find( gi );
          if( it != g2l.cend() )
          {
            // We need to be sure that the lowest (owning) rank sends the node.
            collocatedNodesToSend.emplace_back( it->second );
          }
        }
      }

      for( localIndex const & li: collocatedNodesToSend )
      {
        if( std::find( matchedBoundaryNodes.begin(), matchedBoundaryNodes.end(), li ) == matchedBoundaryNodes.end() )
        {
          matchedBoundaryNodes.emplace_back( li );
        }
      }

      // Some ranks may also have some nodes that the current rank requires.
      // In that case, we store the information.
      // Later in this function, we'll send this request to the rank and it will send those to us.
      array1d< globalIndex > & requestedMatches = requestedMatchesMap[allNeighbors[i].neighborRank()];
      {
        std::vector< globalIndex > intersection;
        std::set_intersection( requestedNodes.cbegin(), requestedNodes.cend(),
                               neighborNodes.cbegin(), neighborNodes.cend(),
                               std::back_inserter( intersection ) );
        for( globalIndex const & gn: intersection )
        {
          requestedMatches.emplace_back( gn );
        }
      }
    }
  }

  // We send the requested nodes to the relevant ranks.
  auto const data = [req = reorganizeRequestedNodes( requestedMatchesMap ), &allNeighbors]( auto i ) -> array1d< globalIndex > const &
  {
    return req.at( allNeighbors[i].neighborRank() );
  };
  array1d< array1d< globalIndex > > const nodesRequestedByNeighbors = exchange( allNeighbors, data );

  // Then we store the requested nodes for each receiver.
  for( integer i = 0; i < numNeighbors; ++i )
  {
    std::set< globalIndex > nodesNotFound;
    array1d< localIndex > & matchedBoundaryNodes = nodeManager.getNeighborData( allNeighbors[i].neighborRank() ).matchedPartitionBoundary();
    for( globalIndex const & gi: nodesRequestedByNeighbors[i] )
    {
      auto const locIt = g2l.find( gi );
      if( locIt != g2l.cend() )
      {
        matchedBoundaryNodes.emplace_back( locIt->second );
      }
      else
      {
        nodesNotFound.insert( gi );
      }
    }
    GEOS_ERROR_IF( !nodesNotFound.empty(),
                   "Global nodes {" << stringutilities::join( nodesNotFound, ", " ) << "} requested by rank " << allNeighbors[i].neighborRank() << " were not found on this rank." );
  }
}

/**
 * @brief Check that the provided object's ghosts are consistent with the neighbors.
 * @param objectManager the owner of the objects to check.
 * @param neighbors list of all the neighbors.
 */
void verifyGhostingConsistency( ObjectManagerBase const & objectManager,
                                std::vector< NeighborCommunicator > const & neighbors )
{
  arrayView1d< integer const > const & ghostRank = objectManager.ghostRank();

  // Variable to track if an error has occurred.
  bool error = false;

  // For each neighbor make sure that the ghost rank is consistent with the send and receive lists.
  for( NeighborCommunicator const & neighbor : neighbors )
  {
    int const neighborRank = neighbor.neighborRank();
    NeighborData const & neighborData = objectManager.getNeighborData( neighborRank );

    arrayView1d< localIndex const > const & recvList = neighborData.ghostsToReceive();
    for( localIndex const recvIdx : recvList )
    {
      if( ghostRank[ recvIdx ] != neighborRank )
      {
        error = true;
        GEOS_LOG_RANK( "Receiving " << recvIdx << " from " << neighborRank <<
                       " but ghostRank[ " << recvIdx << " ] is " << ghostRank[ recvIdx ] );
      }
    }

    arrayView1d< localIndex const > const & sendList = neighborData.ghostsToSend();
    for( localIndex const sendIdx : sendList )
    {
      if( ghostRank[ sendIdx ] != -1 )
      {
        error = true;
        GEOS_LOG_RANK( "Sending " << sendIdx << " to " << neighborRank <<
                       " but ghostRank[ " << sendIdx << " ] is " << ghostRank[ sendIdx ] );
      }
    }

    arrayView1d< std::pair< globalIndex, int > const > const & nonLocalGhosts = neighborData.nonLocalGhosts();
    if( !nonLocalGhosts.empty() )
    {
      error = true;
      GEOS_LOG_RANK( "Expected to send 0 non local ghosts to rank " << neighborRank <<
                     " but sending " << nonLocalGhosts.size() );
    }
  }

  GEOS_ERROR_IF( error, "Encountered a ghosting inconsistency in " << objectManager.getName() );
}

/**
 * @brief Remove the given indices from the communication list.
 * @param indicesToAdd the local indices of objects to be removed.
 * @param commIndices the local indices of the existing objects to be communicated.
 */
void removeFromCommList( std::vector< localIndex > const & indicesToRemove, array1d< localIndex > & commIndices )
{
  localIndex * const itr = std::remove_if( commIndices.begin(), commIndices.end(), [&indicesToRemove]( localIndex const idx )
  {
    return std::find( indicesToRemove.begin(), indicesToRemove.end(), idx ) != indicesToRemove.end();
  } );

  localIndex const nRemoved = commIndices.end() - itr;
  GEOS_ERROR_IF_NE( nRemoved, localIndex( indicesToRemove.size() ) );
  commIndices.resize( commIndices.size() - nRemoved );
}

/**
 * @brief Fix up second neighbor ghosting issues by modifying the receive lists and ghost rank.
 * @param objectManager the owner of the objects to fix up.
 * @param neighbors array of neighbors.
 */
void fixReceiveLists( ObjectManagerBase & objectManager,
                      std::vector< NeighborCommunicator > const & neighbors )
{
  int nonLocalGhostsTag = 45;

  std::vector< MPI_Request > nonLocalGhostsRequests( neighbors.size(), MPI_REQUEST_NULL );

  /// For each neighbor send them the indices of their ghosts that they mistakenly believe are owned by this rank.
  for( std::size_t i = 0; i < neighbors.size(); ++i )
  {
    int const neighborRank = neighbors[ i ].neighborRank();

    MpiWrapper::iSend( objectManager.getNeighborData( neighborRank ).nonLocalGhosts().toView(),
                       neighborRank,
                       nonLocalGhostsTag,
                       MPI_COMM_GEOS,
                       &nonLocalGhostsRequests[ i ] );
  }

  for( NeighborCommunicator const & neighbor : neighbors )
  {
    int const neighborRank = neighbor.neighborRank();

    /// Receive the lists of ghosts we mistakenly though were owned by this neighbor.
    array1d< std::pair< globalIndex, int > > ghostsFromSecondNeighbor;
    MpiWrapper::recv( ghostsFromSecondNeighbor,
                      neighborRank,
                      nonLocalGhostsTag,
                      MPI_COMM_GEOS,
                      MPI_STATUS_IGNORE );

    /// Array of ghosts to fix.
    std::vector< localIndex > ghostsToFix;

    /// Map from owning MPI rank to an array of local objects we need to fix.
    std::unordered_map< int, std::vector< localIndex > > ghostsBySecondNeighbor;

    arrayView1d< integer > const & ghostRank = objectManager.ghostRank();

    /// Populate ghostsToFix and ghostsBySecondNeighbor while also updating ghostRank.
    for( std::pair< globalIndex, int > const & pair : ghostsFromSecondNeighbor )
    {
      localIndex const lid = objectManager.globalToLocalMap( pair.first );
      ghostsBySecondNeighbor[ pair.second ].emplace_back( lid );
      ghostsToFix.emplace_back( lid );
      ghostRank[ lid ] = pair.second;
    }

    /// Remove the ghosts to fix from the neighbor's receive list.
    removeFromCommList( ghostsToFix, objectManager.getNeighborData( neighborRank ).ghostsToReceive() );

    /// Iterate over the ranks that own the objects. For each rank add the new objects to the receive list.
    for( std::pair< int const, std::vector< localIndex > > const & pair : ghostsBySecondNeighbor )
    {
      array1d< localIndex > & trueOwnerRecvList = objectManager.getNeighborData( pair.first ).ghostsToReceive();
      trueOwnerRecvList.insert( trueOwnerRecvList.size(), pair.second.begin(), pair.second.end() );
    }
  }

  /// Wait on the initial send requests.
  MpiWrapper::waitAll( nonLocalGhostsRequests.size(), nonLocalGhostsRequests.data(), MPI_STATUSES_IGNORE );
}

/**
 * @brief Remove neighbors if there is no communication between them
 * @param nodeManager the NodeManager.
 * @param edgeManager the EdgeManager.
 * @param faceManager the FaceManager.
 * @param elemManager the ElementRegionManager.
 * @param neighbors the list of NeighborCommunicators, may be modified.
 */
void removeUnusedNeighbors( NodeManager & nodeManager,
                            EdgeManager & edgeManager,
                            FaceManager & faceManager,
                            ElementRegionManager & elemManager,
                            std::vector< NeighborCommunicator > & neighbors )
{
  for( std::size_t i = 0; i < neighbors.size(); )
  {
    int const neighborRank = neighbors[ i ].neighborRank();

    bool used = false;

    used = used || nodeManager.getNeighborData( neighborRank ).communicationExists();

    used = used || edgeManager.getNeighborData( neighborRank ).communicationExists();

    used = used || faceManager.getNeighborData( neighborRank ).communicationExists();

    elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
    {
      used = used || subRegion.getNeighborData( neighborRank ).communicationExists();
    } );

    if( used )
    {
      ++i;
    }
    else
    {
      nodeManager.removeNeighbor( neighborRank );
      edgeManager.removeNeighbor( neighborRank );
      faceManager.removeNeighbor( neighborRank );

      elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
      {
        subRegion.removeNeighbor( neighborRank );
      } );

      neighbors.erase( neighbors.begin() + i );
    }
  }
}

/**
 * @brief Complete each communication phase in order or with a wait any.
 * @param n the number of requests in each phase.
 * @param phases list of phases.
 * @param unorderedComms if true complete the communications of each phase in the order they are received.
 */
void waitOrderedOrWaitAll( int const n,
                           std::vector< std::tuple< MPI_Request *, MPI_Status *, std::function< MPI_Request ( int ) > > > const & phases,
                           bool const unorderedComms )
{
  if( unorderedComms )
  {
    MpiWrapper::activeWaitSomeCompletePhase( n, phases );
  }
  else
  {
    MpiWrapper::activeWaitOrderedCompletePhase( n, phases );
  }
}

void CommunicationTools::setupGhosts( MeshLevel & meshLevel,
                                      std::vector< NeighborCommunicator > & neighbors,
                                      bool const unorderedComms )
{
  GEOS_MARK_FUNCTION;
  MPI_iCommData commData;
  commData.resize( neighbors.size() );

  NodeManager & nodeManager = meshLevel.getNodeManager();
  EdgeManager & edgeManager = meshLevel.getEdgeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  auto sendGhosts = [&] ( int idx )
  {
    neighbors[idx].prepareAndSendGhosts( false,
                                         1,
                                         meshLevel,
                                         commData.commID(),
                                         commData.mpiRecvBufferSizeRequest( idx ),
                                         commData.mpiSendBufferSizeRequest( idx ),
                                         commData.mpiSendBufferRequest( idx ) );

    return commData.mpiRecvBufferSizeRequest( idx );
  };
  auto postRecv = [&] ( int idx )
  {
    neighbors[idx].postRecv( commData.commID(),
                             commData.mpiRecvBufferRequest( idx ) );
    return commData.mpiRecvBufferRequest( idx );
  };
  auto unpackGhosts = [&] ( int idx )
  {
    neighbors[idx].unpackGhosts( meshLevel, commData.commID() );
    return MPI_REQUEST_NULL;
  };

  waitOrderedOrWaitAll( neighbors.size(),
                        { std::make_tuple( static_cast< MPI_Request * >(nullptr), static_cast< MPI_Status * >(nullptr), sendGhosts ),
                          std::make_tuple( commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus(), postRecv ),
                          std::make_tuple( commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus(), unpackGhosts ) },
                        unorderedComms );

  // There are cases where the multiple waitOrderedOrWaitAll methods here will clash with
  // each other. This typically occurs at higher processor counts (>256) and large meshes
  // with ~14M elements. Adding the following waitAll methods ensures that the underlying
  // async MPI communication will not interfere with subsequent async communication calls.
  // The underlying problem is that for a given phase of async communication, the same
  // tag numbers are used. This will at least isolate the async calls from each other.
  MpiWrapper::waitAll( commData.size(), commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
  MpiWrapper::waitAll( commData.size(), commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );

  for( auto & neighbor : neighbors )
  {
    neighbor.unpackGhostsData( meshLevel, commData.commID() );
  }

  nodeManager.setReceiveLists();
  edgeManager.setReceiveLists();
  faceManager.setReceiveLists();

  auto sendSyncLists = [&] ( int idx )
  {
    MPI_Request & mpiSizeRecvRequest = commData.mpiRecvBufferSizeRequest( idx );
    MPI_Request & mpiSizeSendRequest = commData.mpiSendBufferSizeRequest( idx );
    MPI_Request & mpiSendRequest = commData.mpiSendBufferRequest( idx );
    neighbors[idx].prepareAndSendSyncLists( meshLevel,
                                            commData.commID(),
                                            mpiSizeRecvRequest,
                                            mpiSizeSendRequest,
                                            mpiSendRequest );
    return mpiSizeRecvRequest;
  };
  auto rebuildSyncLists = [&] ( int idx )
  {
    neighbors[idx].unpackAndRebuildSyncLists( meshLevel,
                                              commData.commID() );
    return MPI_REQUEST_NULL;
  };

  waitOrderedOrWaitAll( neighbors.size(),
                        { std::make_tuple( static_cast< MPI_Request * >(nullptr), static_cast< MPI_Status * >(nullptr), sendSyncLists ),
                          std::make_tuple( commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus(), postRecv ),
                          std::make_tuple( commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus(), rebuildSyncLists ) },
                        unorderedComms );

  // See above comments for the reason behind these waitAll commands
  // RE: isolate multiple async-wait calls
  MpiWrapper::waitAll( commData.size(), commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
  MpiWrapper::waitAll( commData.size(), commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );

  fixReceiveLists( nodeManager, neighbors );
  fixReceiveLists( edgeManager, neighbors );
  fixReceiveLists( faceManager, neighbors );

  waitOrderedOrWaitAll( neighbors.size(),
                        { std::make_tuple( static_cast< MPI_Request * >(nullptr), static_cast< MPI_Status * >(nullptr), sendSyncLists ),
                          std::make_tuple( commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus(), postRecv ),
                          std::make_tuple( commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus(), rebuildSyncLists ) },
                        unorderedComms );

  // See above comments for the reason behind these waitAll commands
  // RE: isolate multiple async-wait
  MpiWrapper::waitAll( commData.size(), commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
  MpiWrapper::waitAll( commData.size(), commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );

  // nodeManager.fixUpDownMaps( false );
  // verifyGhostingConsistency( nodeManager, neighbors );
  // edgeManager.fixUpDownMaps( false );
  // verifyGhostingConsistency( edgeManager, neighbors );
  // faceManager.fixUpDownMaps( false );
  // verifyGhostingConsistency( faceManager, neighbors );
  // elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
  // {
  //   subRegion.fixUpDownMaps( false );
  //   verifyGhostingConsistency( subRegion, neighbors );
  // } );
  // elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  // {
  //   subRegion.fixSecondaryMappings( nodeManager, edgeManager, faceManager, elemManager );
  // } );

  removeUnusedNeighbors( nodeManager, edgeManager, faceManager, elemManager, neighbors );

  nodeManager.compressRelationMaps();
  edgeManager.compressRelationMaps();
  faceManager.compressRelationMaps();
}

void CommunicationTools::synchronizePackSendRecvSizes( FieldIdentifiers const & fieldsToBeSync,
                                                       MeshLevel & mesh,
                                                       std::vector< NeighborCommunicator > & neighbors,
                                                       MPI_iCommData & icomm,
                                                       bool onDevice )
{
  GEOS_MARK_FUNCTION;
  icomm.setFieldsToBeSync( fieldsToBeSync );
  icomm.resize( neighbors.size() );

  parallelDeviceEvents events;
  for( std::size_t neighborIndex = 0; neighborIndex < neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    int const bufferSize = neighbor.packCommSizeForSync( fieldsToBeSync, mesh, icomm.commID(), onDevice, events );

    neighbor.mpiISendReceiveBufferSizes( icomm.commID(),
                                         icomm.mpiSendBufferSizeRequest( neighborIndex ),
                                         icomm.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOS );

    neighbor.resizeSendBuffer( icomm.commID(), bufferSize );
  }
  waitAllDeviceEvents( events );
}


void CommunicationTools::asyncPack( FieldIdentifiers const & fieldsToBeSync,
                                    MeshLevel & mesh,
                                    std::vector< NeighborCommunicator > & neighbors,
                                    MPI_iCommData & icomm,
                                    bool onDevice,
                                    parallelDeviceEvents & events )
{
  GEOS_MARK_FUNCTION;
  for( NeighborCommunicator & neighbor : neighbors )
  {
    neighbor.packCommBufferForSync( fieldsToBeSync, mesh, icomm.commID(), onDevice, events );
  }
}

void CommunicationTools::asyncSendRecv( std::vector< NeighborCommunicator > & neighbors,
                                        MPI_iCommData & icomm,
                                        bool onDevice,
                                        parallelDeviceEvents & events )
{
  GEOS_MARK_FUNCTION;
  if( onDevice )
  {
    waitAllDeviceEvents( events );
  }


  // could swap this to test and make this function call async as well, only launch the sends/recvs for
  // those we've already recv'd sizing for, go back to some usefule compute / launch some other compute, then
  // check this again
  for( std::size_t count = 0; count < neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( icomm.size(),
                         icomm.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         icomm.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( icomm.commID(),
                                     icomm.mpiSendBufferRequest( neighborIndex ),
                                     icomm.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOS );
  }
}

void CommunicationTools::synchronizePackSendRecv( FieldIdentifiers const & fieldsToBeSync,
                                                  MeshLevel & mesh,
                                                  std::vector< NeighborCommunicator > & neighbors,
                                                  MPI_iCommData & icomm,
                                                  bool onDevice )
{
  GEOS_MARK_FUNCTION;
  parallelDeviceEvents events;
  asyncPack( fieldsToBeSync, mesh, neighbors, icomm, onDevice, events );
  asyncSendRecv( neighbors, icomm, onDevice, events );
}


bool CommunicationTools::asyncUnpack( MeshLevel & mesh,
                                      std::vector< NeighborCommunicator > & neighbors,
                                      MPI_iCommData & icomm,
                                      bool onDevice,
                                      parallelDeviceEvents & events,
                                      MPI_Op op )
{
  GEOS_MARK_FUNCTION;

  int recvCount = 0;
  std::vector< int > neighborIndices;
  neighborIndices.reserve( icomm.size() );
  MpiWrapper::testSome( icomm.size(),
                        icomm.mpiRecvBufferRequest(),
                        &recvCount,
                        &neighborIndices[0],
                        icomm.mpiRecvBufferStatus() );

  for( int recvIdx = 0; recvIdx < recvCount; ++recvIdx )
  {
    NeighborCommunicator & neighbor = neighbors[ neighborIndices[ recvIdx ] ];
    neighbor.unpackBufferForSync( icomm.getFieldsToBeSync(), mesh, icomm.commID(), onDevice, events, op );
  }

  // we don't want to check if the request has completed,
  //  we want to check that we've processed the resulting buffer
  //  which means that we've tested the request and it has been
  //  deallocated and set to MPI_REQUEST_NULL
  int allDone = true;
  const MPI_Request * reqs = icomm.mpiRecvBufferRequest( );
  for( int idx = 0; idx < icomm.size(); ++idx )
  {
    if( reqs[ idx ] != MPI_REQUEST_NULL )
    {
      allDone = false;
      break;
    }
  }

  return allDone;
}

void CommunicationTools::finalizeUnpack( MeshLevel & mesh,
                                         std::vector< NeighborCommunicator > & neighbors,
                                         MPI_iCommData & icomm,
                                         bool onDevice,
                                         parallelDeviceEvents & events,
                                         MPI_Op op )
{
  GEOS_MARK_FUNCTION;

  // poll mpi for completion then wait 10 nanoseconds 6,000,000,000 times (60 sec timeout)
  GEOS_ASYNC_WAIT( 6000000000, 10, asyncUnpack( mesh, neighbors, icomm, onDevice, events, op ) );
  if( onDevice )
  {
    waitAllDeviceEvents( events );
  }

  MpiWrapper::waitAll( icomm.size(),
                       icomm.mpiSendBufferSizeRequest(),
                       icomm.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( icomm.size(),
                       icomm.mpiSendBufferRequest(),
                       icomm.mpiSendBufferStatus() );

}

void CommunicationTools::synchronizeUnpack( MeshLevel & mesh,
                                            std::vector< NeighborCommunicator > & neighbors,
                                            MPI_iCommData & icomm,
                                            bool onDevice )
{
  GEOS_MARK_FUNCTION;
  parallelDeviceEvents events;
  finalizeUnpack( mesh, neighbors, icomm, onDevice, events );
}

void CommunicationTools::synchronizeFields( FieldIdentifiers const & fieldsToBeSync,
                                            MeshLevel & mesh,
                                            std::vector< NeighborCommunicator > & neighbors,
                                            bool onDevice )
{
  MPI_iCommData icomm;
  icomm.resize( neighbors.size() );
  synchronizePackSendRecvSizes( fieldsToBeSync, mesh, neighbors, icomm, onDevice );
  synchronizePackSendRecv( fieldsToBeSync, mesh, neighbors, icomm, onDevice );
  synchronizeUnpack( mesh, neighbors, icomm, onDevice );
}

} /* namespace geos */
