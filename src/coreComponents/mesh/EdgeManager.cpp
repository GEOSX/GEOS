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
 * @file EdgeManager.cpp
 */

#include "EdgeManager.hpp"

#include "BufferOps.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"
#include "common/TimingMacros.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "mesh/generators/CellBlockUtilities.hpp"

namespace geos
{
using namespace dataRepository;

EdgeManager::EdgeManager( string const & name,
                          Group * const parent ):
  ObjectManagerBase( name, parent )
{
  registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation );
  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation ).
    setSizedFromParent( 0 );

  m_toNodesRelation.resize( 0, 2 );

  excludeWrappersFromPacking( { viewKeyStruct::nodeListString(),
                                viewKeyStruct::faceListString(),
                                viewKeyStruct::elementRegionListString(),
                                viewKeyStruct::elementSubRegionListString(),
                                viewKeyStruct::elementListString() } );
}

void EdgeManager::resize( localIndex const newSize )
{
  m_toFacesRelation.resize( newSize, 2 * faceMapOverallocation() );
  ObjectManagerBase::resize( newSize );
}

void EdgeManager::buildSets( NodeManager const & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // Make sets from node sets.
  auto const & nodeSets = nodeManager.sets().wrappers();
  for( int i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    createSet( setName );
  }

  // Then loop over them in parallel.
  forAll< parallelHostPolicy >( nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    SortedArrayView< localIndex const > const targetSet = nodeManager.sets().getReference< SortedArray< localIndex > >( setName ).toViewConst();
    constructSetFromSetAndMap( targetSet, m_toNodesRelation, setName );
  } );
}

void EdgeManager::buildEdges( localIndex const numNodes,
                              ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                              ArrayOfArrays< localIndex > & faceToEdgeMap )
{
  ArrayOfArrays< localIndex > edgeToFace;
  localIndex const numEdges = buildEdgeMaps( numNodes,
                                             faceToNodeMap,
                                             faceToEdgeMap,
                                             edgeToFace,
                                             m_toNodesRelation );

  m_toFacesRelation.base().assimilate< parallelHostPolicy >( std::move( edgeToFace ),
                                                             LvArray::sortedArrayManipulation::UNSORTED_NO_DUPLICATES );
  resize( numEdges );
}

void EdgeManager::setGeometricalRelations( CellBlockManagerABC const & cellBlockManager, bool isBaseMeshLevel )
{
  GEOS_MARK_FUNCTION;
  if( isBaseMeshLevel )
  {
    resize( cellBlockManager.numEdges() );
  }

  m_toNodesRelation.base() = cellBlockManager.getEdgeToNodes();
  m_toFacesRelation.base().assimilate< parallelHostPolicy >( cellBlockManager.getEdgeToFaces(),
                                                             LvArray::sortedArrayManipulation::UNSORTED_NO_DUPLICATES );
}

void EdgeManager::setupRelatedObjectsInRelations( NodeManager const & nodeManager,
                                                  FaceManager const & faceManager )
{
  m_toNodesRelation.setRelatedObject( nodeManager );
  m_toFacesRelation.setRelatedObject( faceManager );
}

void EdgeManager::setDomainBoundaryObjects( FaceManager const & faceManager )
{
  // get the "isDomainBoundary" field from the faceManager. This should have been set already!
  arrayView1d< integer const > const isFaceOnDomainBoundary = faceManager.getDomainBoundaryIndicator();

  // get the "isDomainBoundary" field from for *this, and set it to zero
  arrayView1d< integer > const isEdgeOnDomainBoundary = this->getDomainBoundaryIndicator();
  isEdgeOnDomainBoundary.zero();

  ArrayOfArraysView< localIndex const > const faceToEdges = faceManager.edgeList().toViewConst();

  forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const faceIndex )
  {
    if( isFaceOnDomainBoundary[faceIndex] == 1 )
    {
      for( localIndex const edgeIndex : faceToEdges[faceIndex] )
      {
        isEdgeOnDomainBoundary[edgeIndex] = 1;
      }
    }
  } );
}

bool EdgeManager::hasNode( const localIndex edgeIndex, const localIndex nodeIndex ) const
{
  return m_toNodesRelation( edgeIndex, 0 ) == nodeIndex || m_toNodesRelation( edgeIndex, 1 ) == nodeIndex;
}

void EdgeManager::setIsExternal( FaceManager const & faceManager )
{
  // get the "isExternal" field from the faceManager...This should have been
  // set already!
  arrayView1d< integer const > const & isExternalFace = faceManager.isExternal();

  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  // get the "isExternal" field from for *this, and set it to zero
  m_isExternal.zero();

  // loop through all faces
  for( localIndex kf=0; kf<faceManager.size(); ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isExternalFace[kf] == 1 )
    {
      // loop over all nodes connected to face, and set isNodeDomainBoundary
      localIndex const numEdges = faceToEdges.sizeOfArray( kf );
      for( localIndex a = 0; a < numEdges; ++a )
      {
        m_isExternal[ faceToEdges( kf, a ) ] = 1;
      }
    }
  }
}

ArrayOfSets< globalIndex >
EdgeManager::extractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & nodeManager )
{
  GEOS_MARK_FUNCTION;

  localIndex const numEdges = size();

  arrayView2d< localIndex const > const edgeNodes = this->nodeList();
  arrayView1d< integer const > const isDomainBoundary = this->getDomainBoundaryIndicator();
  arrayView1d< globalIndex const > const nodeLocalToGlobal = nodeManager.localToGlobalMap();

  ArrayOfSets< globalIndex > globalEdgeNodes( numEdges, 2 );

  forAll< parallelHostPolicy >( numEdges, [globalEdgeNodes = globalEdgeNodes.toView(),
                                           isDomainBoundary, edgeNodes, nodeLocalToGlobal]( localIndex const edgeIndex )
  {
    if( isDomainBoundary( edgeIndex ) )
    {
      for( localIndex const nodeIndex : edgeNodes[edgeIndex] )
      {
        globalEdgeNodes.insertIntoSet( edgeIndex, nodeLocalToGlobal[nodeIndex] );
      }
    }
  } );

  return globalEdgeNodes;
}

localIndex EdgeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex EdgeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex EdgeManager::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                            arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();

  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toNodesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToNodes,
                                               packList,
                                               localToGlobal,
                                               nodeLocalToGlobal );


  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toFacesRelation.base().toArrayOfArraysView(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  return packedSize;
}



localIndex EdgeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  GEOS_MARK_FUNCTION;

  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void EdgeManager::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation.base(),
                                    m_toFacesRelation.relatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}

void EdgeManager::compressRelationMaps()
{
  m_toFacesRelation.compress();
}

void EdgeManager::depopulateUpMaps( std::set< localIndex > const & receivedEdges,
                                    ArrayOfArraysView< localIndex const > const & facesToEdges )
{
  ObjectManagerBase::cleanUpMap( receivedEdges, m_toFacesRelation.toView(), facesToEdges );
}


} /// namespace geos
