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
 * @file EmbeddedSurfaceNodeManager.cpp
 */

#include "EmbeddedSurfaceNodeManager.hpp"

#include "common/TimingMacros.hpp"
#include "common/MpiWrapper.hpp"
#include "mesh/BufferOps.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ToElementRelation.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"

namespace geos
{

using namespace dataRepository;

EmbeddedSurfaceNodeManager::EmbeddedSurfaceNodeManager( string const & name,
                                                        Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition );
  this->registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  this->registerWrapper( viewKeyStruct::elementRegionListString(), &elementRegionList() );
  this->registerWrapper( viewKeyStruct::elementSubRegionListString(), &elementSubRegionList() );
  this->registerWrapper( viewKeyStruct::elementListString(), &elementList() );
  this->registerWrapper( viewKeyStruct::parentEdgeGlobalIndexString(), &m_parentEdgeGlobalIndex );

  excludeWrappersFromPacking( { viewKeyStruct::edgeListString(),
                                viewKeyStruct::elementRegionListString(),
                                viewKeyStruct::elementSubRegionListString(),
                                viewKeyStruct::elementListString() } );
}


void EmbeddedSurfaceNodeManager::resize( localIndex const newSize )
{
  m_toEdgesRelation.resize( newSize, 2 * edgeMapOverallocation() );
  m_toElements.m_toElementRegion.resize( newSize, 2 * elemMapOverallocation() );
  m_toElements.m_toElementSubRegion.resize( newSize, 2 * elemMapOverallocation() );
  m_toElements.m_toElementIndex.resize( newSize, 2 * elemMapOverallocation() );
  ObjectManagerBase::resize( newSize );
}


void EmbeddedSurfaceNodeManager::setEdgeMaps( EdgeManager const & embSurfEdgeManager )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< localIndex const > const edgeToNodeMap = embSurfEdgeManager.nodeList();

  ArrayOfArrays< localIndex > nodeToEdges =
    meshMapUtilities::transposeIndexMap< parallelHostPolicy >( edgeToNodeMap,
                                                               size(),
                                                               edgeMapOverallocation() );

  m_toEdgesRelation.assimilate< parallelHostPolicy >( std::move( nodeToEdges ),
                                                      LvArray::sortedArrayManipulation::UNSORTED_NO_DUPLICATES );
  m_toEdgesRelation.setRelatedObject( embSurfEdgeManager );
}

void EmbeddedSurfaceNodeManager::setElementMaps( ElementRegionManager const & elementRegionManager )
{
  GEOS_MARK_FUNCTION;

  ArrayOfArrays< localIndex > & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays< localIndex > & toElementSubRegionList = m_toElements.m_toElementSubRegion;
  ArrayOfArrays< localIndex > & toElementList = m_toElements.m_toElementIndex;
  localIndex const numNodes = size();

  // The number of elements attached to the each node.
  array1d< localIndex > elemsPerNode( numNodes );

  // The total number of elements, the sum of elemsPerNode.
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeElems = 0;

  elementRegionManager.
    forElementSubRegions< EmbeddedSurfaceSubRegion >( [&elemsPerNode, &totalNodeElems]( EmbeddedSurfaceSubRegion const & subRegion )
  {
    EmbeddedSurfaceSubRegion::NodeMapType const & elemToNodeMap = subRegion.nodeList();
    forAll< parallelHostPolicy >( subRegion.size(), [&elemsPerNode, totalNodeElems, &elemToNodeMap ] ( localIndex const k )
    {
      localIndex const numElementNodes = elemToNodeMap.sizeOfArray( k );
      totalNodeElems += numElementNodes;
      for( localIndex a = 0; a < numElementNodes; ++a )
      {
        localIndex const nodeIndex = elemToNodeMap( k, a );
        RAJA::atomicInc< parallelHostAtomic >( &elemsPerNode[ nodeIndex ] );
      }
    } );
  } );

  // Resize the node to elem map.
  toElementRegionList.resize( 0 );
  toElementSubRegionList.resize( 0 );
  toElementList.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  toElementRegionList.reserve( entriesToReserve );
  toElementSubRegionList.reserve( entriesToReserve );
  toElementList.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeElems.get() + numNodes * elemMapOverallocation() * ( 1 + 2 * overAllocationFactor );
  toElementRegionList.reserveValues( valuesToReserve );
  toElementSubRegionList.reserveValues( valuesToReserve );
  toElementList.reserveValues( valuesToReserve );

  // Append an array for each node with capacity to hold the appropriate number of elements plus some wiggle room.
  for( localIndex nodeIndex = 0; nodeIndex < numNodes; ++nodeIndex )
  {
    toElementRegionList.appendArray( 0 );
    toElementSubRegionList.appendArray( 0 );
    toElementList.appendArray( 0 );

    toElementRegionList.setCapacityOfArray( nodeIndex, elemsPerNode[ nodeIndex ] + elemMapOverallocation() );
    toElementSubRegionList.setCapacityOfArray( nodeIndex, elemsPerNode[ nodeIndex ] + elemMapOverallocation() );
    toElementList.setCapacityOfArray( nodeIndex, elemsPerNode[ nodeIndex ] + elemMapOverallocation() );
  }

  // Populate the element maps.
  // Note that this can't be done in parallel because the three element lists must be in the same order.
  // If this becomes a bottleneck create a temporary ArrayOfArrays of tuples and insert into that first then copy over.
  elementRegionManager.
    forElementSubRegionsComplete< EmbeddedSurfaceSubRegion >( [&toElementRegionList, &toElementSubRegionList, &toElementList]
                                                                ( localIndex const er, localIndex const esr, ElementRegionBase const &,
                                                                EmbeddedSurfaceSubRegion const & subRegion )
  {
    EmbeddedSurfaceSubRegion::NodeMapType const & elemToNodeMap = subRegion.nodeList();
    for( localIndex k = 0; k < subRegion.size(); ++k )
    {
      for( localIndex a=0; a<elemToNodeMap.sizeOfArray( k ); ++a )
      {
        localIndex const nodeIndex = elemToNodeMap( k, a );
        toElementRegionList.emplaceBack( nodeIndex, er );
        toElementSubRegionList.emplaceBack( nodeIndex, esr );
        toElementList.emplaceBack( nodeIndex, k );
      }
    }
  } );

  this->m_toElements.setElementRegionManager( elementRegionManager );
}


void EmbeddedSurfaceNodeManager::compressRelationMaps()
{
  m_toEdgesRelation.compress();
  m_toElements.m_toElementRegion.compress();
  m_toElements.m_toElementSubRegion.compress();
  m_toElements.m_toElementIndex.compress();
}


void EmbeddedSurfaceNodeManager::appendNode( arraySlice1d< real64 const > const & pointCoord,
                                             integer const & pointGhostRank )
{
  localIndex nodeIndex =  this->size();
  this->resize( nodeIndex + 1 );
  LvArray::tensorOps::copy< 3 >( m_referencePosition[nodeIndex], pointCoord );
  m_ghostRank[ nodeIndex ] = pointGhostRank;
}


localIndex EmbeddedSurfaceNodeManager::packNewNodesGlobalMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packNewNodesGlobalMapsImpl< false >( junk, packList );
}

localIndex EmbeddedSurfaceNodeManager::packNewNodesGlobalMaps( buffer_unit_type * & buffer,
                                                               arrayView1d< localIndex const > const & packList ) const
{
  return packNewNodesGlobalMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex EmbeddedSurfaceNodeManager::packNewNodesGlobalMapsImpl( buffer_unit_type * & buffer,
                                                                   arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer, this->getName() );

  // this doesn't link without the string()...no idea why.
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::localToGlobalMapString() ) );

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, rank );

  localIndex const numPackedIndices = packList.size();
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, numPackedIndices );

  if( numPackedIndices > 0 )
  {
    // We pack 3 things:
    // 1. the global indices
    globalIndex_array globalIndices;
    globalIndices.resize( numPackedIndices );
    // 2. the ghostRank
    array1d< integer > ghostRanks;
    ghostRanks.resize( numPackedIndices );
    for( localIndex a=0; a<numPackedIndices; ++a )
    {
      globalIndices[a] = this->m_localToGlobalMap[packList[a]];
      ghostRanks[a]= this->m_ghostRank[ packList[a] ];
    }
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, globalIndices );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, ghostRanks );

    // 3. the referencePosition
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePosition = this->referencePosition();
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::referencePositionString() ) );
    packedSize += bufferOps::PackByIndex< DO_PACKING >( buffer, referencePosition, packList );
  }

  return packedSize;
}


localIndex EmbeddedSurfaceNodeManager::unpackNewNodesGlobalMaps( buffer_unit_type const * & buffer,
                                                                 localIndex_array & packList )
{
  GEOS_MARK_FUNCTION;

  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOS_ERROR_IF( groupName != this->getName(), "EmbeddedSurfaceNodeManager::unpackGlobalMaps(): group names do not match" );

  string localToGlobalString;
  unpackedSize += bufferOps::Unpack( buffer, localToGlobalString );
  GEOS_ERROR_IF( localToGlobalString != viewKeyStruct::localToGlobalMapString(), "ObjectManagerBase::unpack(): label incorrect" );

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  int sendingRank;
  unpackedSize += bufferOps::Unpack( buffer, sendingRank );

  localIndex numUnpackedIndices;
  unpackedSize += bufferOps::Unpack( buffer, numUnpackedIndices );

  if( numUnpackedIndices > 0 )
  {
    localIndex_array unpackedLocalIndices;
    unpackedLocalIndices.resize( numUnpackedIndices );

    globalIndex_array globalIndices;
    unpackedSize += bufferOps::Unpack( buffer, globalIndices );

    array1d< integer > ghostRankOnSendingRank;
    unpackedSize += bufferOps::Unpack( buffer, ghostRankOnSendingRank );

    // Unpack referencePosition
    localIndex_array indicesOnBuffer( numUnpackedIndices );
    for( localIndex i = 0; i < numUnpackedIndices; i++ )
    {
      indicesOnBuffer[i] = i;
    }
    array2d< real64, nodes::REFERENCE_POSITION_PERM > referencePositionData( 0, 3 );
    referencePositionData.resize( numUnpackedIndices );
    arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & referencePosition = referencePositionData.toView();
    string referencePositionString;
    unpackedSize += bufferOps::Unpack( buffer, referencePositionString );
    GEOS_ERROR_IF( referencePositionString != viewKeyStruct::referencePositionString(), "EmbeddedSurfaceNodeManager::unpackGlobalMaps(): label incorrect" );
    unpackedSize += bufferOps::UnpackByIndex( buffer, referencePosition, indicesOnBuffer );

    localIndex numNewIndices = 0;
    globalIndex_array newGlobalIndices;
    newGlobalIndices.reserve( numUnpackedIndices );
    array1d< integer > ghostRank;
    ghostRank.reserve( numUnpackedIndices );
    array1d< integer > originalIndex;
    originalIndex.reserve( numUnpackedIndices );
    localIndex const oldSize = this->size();
    for( localIndex a = 0; a < numUnpackedIndices; ++a )
    {
      // check to see if the node already exists by comparing the coordinates by checking for the global
      localIndex nodeIndexOnThisRank = -1;
      real64 nodeCoord[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( referencePosition[a] );
      nodeExistsOnThisRank( nodeCoord, nodeIndexOnThisRank );
      if( nodeIndexOnThisRank > -1 )
      {
        if( m_ghostRank[nodeIndexOnThisRank] == sendingRank )
        {
          // object already exists on this rank and it's a ghost
          unpackedLocalIndices( a ) = nodeIndexOnThisRank;

          globalIndex const globalIndexOnthisRank = m_localToGlobalMap[ nodeIndexOnThisRank ];

          m_localToGlobalMap[ nodeIndexOnThisRank ] = globalIndices( a );

          //modify global to local map
          m_globalToLocalMap.erase( globalIndexOnthisRank );
          m_globalToLocalMap.insert( {globalIndices( a ), nodeIndexOnThisRank} );

        }
      }
      else
      {
        // object does not exist on this domain
        const localIndex newLocalIndex = oldSize + numNewIndices;

        // add the global index of the new object to the globalToLocal map
        m_globalToLocalMap[ globalIndices[a] ] = newLocalIndex;

        unpackedLocalIndices( a ) = newLocalIndex;

        newGlobalIndices.emplace_back( globalIndices[a] );

        if( ghostRankOnSendingRank[a] == rank )
        {
          ghostRank.emplace_back( -1 );
          originalIndex.emplace_back( a );
        }
        else
        {
          ghostRank.emplace_back( sendingRank );
          originalIndex.emplace_back( -1 );
        }

        ++numNewIndices;

        GEOS_ERROR_IF( packList.size() != 0,
                       "EmbeddedSurfaceNodeManager::unpackGlobalMaps(): packList specified, "
                       "but a new globalIndex is unpacked" );
      }
    }

    ///

    // figure out new size of object container, and resize it
    const localIndex newSize = oldSize + numNewIndices;
    this->resize( newSize );

    // add the new indices to the maps.
    for( int a=0; a<numNewIndices; ++a )
    {
      localIndex const b = oldSize + a;
      m_localToGlobalMap[b] = newGlobalIndices( a );
      m_ghostRank[b] = ghostRank( a );
      if( originalIndex[a] > -1 && m_ghostRank[b] == -1 )
      {
        LvArray::tensorOps::copy< 3 >( m_referencePosition[b], referencePosition[originalIndex[a]] );
      }
    }


    packList = unpackedLocalIndices;
  }

  return unpackedSize;
}


localIndex EmbeddedSurfaceNodeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}


localIndex EmbeddedSurfaceNodeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                                       arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}


template< bool DO_PACKING >
localIndex EmbeddedSurfaceNodeManager::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                           arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elementListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               this->m_toElements,
                                               packList,
                                               m_toElements.getElementRegionManager() );
  return packedSize;
}


localIndex EmbeddedSurfaceNodeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                         array1d< localIndex > & packList,
                                                         bool const overwriteUpMaps,
                                                         bool const )
{
  localIndex unPackedSize = 0;

  string temp;

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ERROR_IF( temp != viewKeyStruct::elementListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     this->m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void EmbeddedSurfaceNodeManager::nodeExistsOnThisRank( real64 const (&nodeCoord)[3],
                                                       localIndex & nodeIndexOnThisRank ) const
{
  for( localIndex a=0; a<size(); a++ )
  {
    real64 const tolerance = 100.0 * std::numeric_limits< real64 >::epsilon();
    real64 distance[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_referencePosition[a] );
    LvArray::tensorOps::subtract< 3 >( distance, nodeCoord );
    if( LvArray::tensorOps::l2Norm< 3 >( distance ) < tolerance )
    {
      nodeIndexOnThisRank = a;
    }
  }
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, EmbeddedSurfaceNodeManager, string const &, Group * const )

}
