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
 * @file FaceManager.cpp
 */

#include "FaceManager.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/Logger.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "mesh/BufferOps.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/MeshFields.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"
#include "utilities/ComputationalGeometry.hpp"
#include "CellElementRegion.hpp"

namespace geos
{
using namespace dataRepository;

FaceManager::FaceManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  this->registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation );
  this->registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );

  this->registerWrapper( viewKeyStruct::elementRegionListString(), &m_toElements.m_toElementRegion ).
    setApplyDefaultValue( -1 );

  this->registerWrapper( viewKeyStruct::elementSubRegionListString(), &m_toElements.m_toElementSubRegion ).
    setApplyDefaultValue( -1 );

  // Do we really want this to be resized and accessed by anyone?
  this->registerWrapper( viewKeyStruct::elementListString(), &m_toElements.m_toElementIndex ).
    setApplyDefaultValue( -1 );

  this->registerWrapper( viewKeyStruct::faceAreaString(), &m_faceArea );

  this->registerWrapper( viewKeyStruct::faceCenterString(), &m_faceCenter ).
    reference().resizeDimension< 1 >( 3 );

  this->registerWrapper( viewKeyStruct::faceNormalString(), &m_faceNormal ).
    reference().resizeDimension< 1 >( 3 );

  excludeWrappersFromPacking( { viewKeyStruct::nodeListString(),
                                viewKeyStruct::edgeListString(),
                                viewKeyStruct::elementRegionListString(),
                                viewKeyStruct::elementSubRegionListString(),
                                viewKeyStruct::elementListString() } );

  m_toElements.resize( 0, 2 );

}

void FaceManager::resize( localIndex const newSize )
{
  m_toNodesRelation.resize( newSize, 2 * nodeMapOverallocation() );
  m_toEdgesRelation.resize( newSize, 2 * edgeMapOverallocation() );
  ObjectManagerBase::resize( newSize );
}

void FaceManager::buildSets( NodeManager const & nodeManager )
{
  GEOS_MARK_FUNCTION;

  // First create the sets
  auto const & nodeSets = nodeManager.sets().wrappers();
  for( auto const & setWrapper : nodeSets )
  {
    createSet( setWrapper.second->getName() );
  }

  // Then loop over them in parallel and fill them in.
  forAll< parallelHostPolicy >( nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    SortedArrayView< localIndex const > const targetSet =
      nodeManager.sets().getReference< SortedArray< localIndex > >( setName ).toViewConst();
    constructSetFromSetAndMap( targetSet, m_toNodesRelation.toViewConst(), setName );
  } );
}

void FaceManager::setDomainBoundaryObjects( ElementRegionManager const & elemRegionManager )
{
  arrayView1d< integer > const isFaceOnDomainBoundary = getDomainBoundaryIndicator();
  isFaceOnDomainBoundary.zero();

  arrayView2d< localIndex const > const toElementRegion = m_toElements.m_toElementRegion.toViewConst();

  forAll< parallelHostPolicy >( size(), [=]( localIndex const kf )
  {
    if( toElementRegion( kf, 1 ) == -1 )
    {
      isFaceOnDomainBoundary( kf ) = 1;
    }
  } );

  // We want to tag as boundary faces all the faces that touch a surface element (mainly a fracture),
  // if this element only has one unique neighbor.
  auto const f = [&]( SurfaceElementRegion const & region )
  {
    if( region.subRegionType() != SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      return;
    }

    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    ArrayOfArraysView< localIndex const > const elem2dToFaces = subRegion.faceList().toViewConst();
    for( int ei = 0; ei < elem2dToFaces.size(); ++ei )
    {
      if( elem2dToFaces.sizeOfArray( ei ) == 2 )
      {
        continue;
      }

      for( localIndex const & face: elem2dToFaces[ei] )
      {
        isFaceOnDomainBoundary[face] = 1;
      }
    }
  };
  elemRegionManager.forElementRegions< SurfaceElementRegion >( f );
}

void FaceManager::setGeometricalRelations( CellBlockManagerABC const & cellBlockManager,
                                           ElementRegionManager const & elemRegionManager,
                                           NodeManager const & nodeManager,
                                           bool isBaseMeshLevel )
{
  GEOS_MARK_FUNCTION;

  if( isBaseMeshLevel )
  {
    resize( cellBlockManager.numFaces() );
  }

  m_toNodesRelation.base() = cellBlockManager.getFaceToNodes();
  m_toEdgesRelation.base() = cellBlockManager.getFaceToEdges();

  ToCellRelation< array2d< localIndex > > const toCellBlock = cellBlockManager.getFaceToElements();
  array2d< localIndex > const blockToSubRegion = elemRegionManager.getCellBlockToSubRegionMap( cellBlockManager );
  meshMapUtilities::transformCellBlockToRegionMap< parallelHostPolicy >( blockToSubRegion.toViewConst(),
                                                                         toCellBlock,
                                                                         m_toElements );

  // Since the mappings of the current FaceManager instance were only filled based on the {face -> elements} mapping,
  // they do not take into account any connection between the 3d elements and the 2d elements of fracture meshes.
  // The following function adds those connections between the 3d and 2d elements.
  auto const connect2dElems = [&]( localIndex er,
                                   SurfaceElementRegion const & region )
  {
    if( region.subRegionType() != SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      return;
    }

    constexpr char err[] = "Internal error when trying to connect matrix mapping and fracture mapping. Face {} seems wrongly connected.";

    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    int const esr = 0;  // Since there's only on unique subregion, the index is always 0.
    // The fracture subregion knows the faces it's connected to.
    // And since a 2d element is connected to a given face, and since a face can only have 2 neighbors,
    // then the second neighbor of the face is bound to be undefined (i.e. -1).
    ArrayOfArraysView< localIndex const > const & elem2dToFaces = subRegion.faceList().toViewConst();
    for( localIndex ei = 0; ei < elem2dToFaces.size(); ++ei )
    {
      for( localIndex const & face: elem2dToFaces[ei] )
      {
        GEOS_ERROR_IF_EQ_MSG( m_toElements.m_toElementRegion( face, 0 ), -1, GEOS_FMT( err, face ) );
        GEOS_ERROR_IF_EQ_MSG( m_toElements.m_toElementSubRegion( face, 0 ), -1, GEOS_FMT( err, face ) );
        GEOS_ERROR_IF_EQ_MSG( m_toElements.m_toElementIndex( face, 0 ), -1, GEOS_FMT( err, face ) );

        GEOS_ERROR_IF_NE_MSG( m_toElements.m_toElementRegion( face, 1 ), -1, GEOS_FMT( err, face ) );
        GEOS_ERROR_IF_NE_MSG( m_toElements.m_toElementSubRegion( face, 1 ), -1, GEOS_FMT( err, face ) );
        GEOS_ERROR_IF_NE_MSG( m_toElements.m_toElementIndex( face, 1 ), -1, GEOS_FMT( err, face ) );

        m_toElements.m_toElementRegion( face, 1 ) = er;
        m_toElements.m_toElementSubRegion( face, 1 ) = esr;
        m_toElements.m_toElementIndex( face, 1 ) = ei;
      }
    }
  };
  // Connecting all the 3d elements (information is already in the m_toElements mappings) and all the 2d elements.
  elemRegionManager.forElementRegionsComplete< SurfaceElementRegion >( connect2dElems );

  if( isBaseMeshLevel )
  {
    computeGeometry( nodeManager );
  }
}

void FaceManager::setupRelatedObjectsInRelations( NodeManager const & nodeManager,
                                                  EdgeManager const & edgeManager,
                                                  ElementRegionManager const & elemRegionManager )
{
  m_toNodesRelation.setRelatedObject( nodeManager );
  m_toEdgesRelation.setRelatedObject( edgeManager );

  m_toElements.setElementRegionManager( elemRegionManager );
}

void FaceManager::computeGeometry( NodeManager const & nodeManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  forAll< parallelHostPolicy >( this->size(), [&]( localIndex const faceIndex )
  {
    m_faceArea[ faceIndex ] = computationalGeometry::centroid_3DPolygon( m_toNodesRelation[ faceIndex ],
                                                                         X,
                                                                         m_faceCenter[ faceIndex ],
                                                                         m_faceNormal[ faceIndex ] );

  } );
}

void FaceManager::setIsExternal()
{
  arrayView1d< integer const > const isDomainBoundary = this->getDomainBoundaryIndicator();

  m_isExternal.zero();
  for( localIndex k=0; k<size(); ++k )
  {
    if( isDomainBoundary[k]==1 )
    {
      m_isExternal[k] = 1;
    }
  }
}

void FaceManager::sortAllFaceNodes( NodeManager const & nodeManager,
                                    ElementRegionManager const & elemManager )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< localIndex const > const facesToElementRegions = elementRegionList();
  arrayView2d< localIndex const > const facesToElementSubRegions = elementSubRegionList();
  arrayView2d< localIndex const > const facesToElements = elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();

  ArrayOfArraysView< localIndex > const facesToNodes = nodeList().toView();

  elemManager.forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >( [&] ( auto const & subRegion )
  {
    subRegion.calculateElementCenters( X );
  } );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > elemCenter =
    elemManager.constructArrayViewAccessor< real64, 2 >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

  forAll< parallelHostPolicy >( size(), [=, elemCenter = elemCenter.toNestedViewConst()]( localIndex const faceIndex )
  {
    // The face should be connected to at least one element.
    if( facesToElements( faceIndex, 0 ) < 0 && facesToElements( faceIndex, 1 ) < 0 )
    {
      GEOS_ERROR( getDataContext() << ": Face " << faceIndex <<
                  " is not connected to any cell.\n"
                  "You might have:\n"
                  "- an invalid mesh,\n"
                  "- not enough CellElementRegions to describe your input mesh (all regions, simulated or not, must be listed),\n"
                  "- forgotten one cell type in an existing \"" << CellElementRegion::viewKeyStruct::sourceCellBlockNamesString() << "\'." );
    }

    // Take the first defined face-to-(elt/region/sub region) to sorting direction.
    localIndex const iElemLoc = facesToElements( faceIndex, 0 ) >= 0 ? 0 : 1;

    localIndex const er = facesToElementRegions( faceIndex, iElemLoc );
    localIndex const esr = facesToElementSubRegions( faceIndex, iElemLoc );
    localIndex const ei = facesToElements( faceIndex, iElemLoc );

    if( er < 0 || esr < 0 || ei < 0 )
    {
      GEOS_ERROR( GEOS_FMT( "{0}: Face {1} is connected to an invalid element ({2}/{3}/{4}).",
                            getDataContext().toString(), faceIndex, er, esr, ei ) );
    }

    try
    {
      sortFaceNodes( X, elemCenter[er][esr][ei], facesToNodes[faceIndex] );
    } catch( std::runtime_error const & e )
    {
      throw std::runtime_error( getDataContext().toString() + ": " + e.what() );
    }
  } );
}

void FaceManager::sortFaceNodes( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                 arraySlice1d< real64 const > const elementCenter,
                                 Span< localIndex > const faceNodes )
{
  localIndex const numFaceNodes = LvArray::integerConversion< localIndex >( faceNodes.size() );
  GEOS_THROW_IF_GT_MSG( numFaceNodes, MAX_FACE_NODES, "The number of maximum nodes allocated per cell face has been reached.", std::runtime_error );

  localIndex const firstNodeIndex = faceNodes[0];

  // get face center (average vertex location)
  real64 fc[3] = { 0 };
  for( localIndex n = 0; n < numFaceNodes; ++n )
  {
    LvArray::tensorOps::add< 3 >( fc, X[faceNodes[n]] );
  }
  LvArray::tensorOps::scale< 3 >( fc, 1.0 / numFaceNodes );

  // Approximate face normal direction (unscaled)

  if( numFaceNodes == 2 )  //2D only.
  {
    real64 ex[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[faceNodes[1]] );
    LvArray::tensorOps::subtract< 3 >( ex, X[faceNodes[0]] );

    real64 ey[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( elementCenter );
    LvArray::tensorOps::subtract< 3 >( ey, fc );

    real64 ez[3];
    LvArray::tensorOps::crossProduct( ez, ex, ey );

    // The element should be on the right hand side of the vector from node 0 to node 1.
    // This ensure that the normal vector of an external face points to outside the element.
    if( ez[2] > 0 )
    {
      localIndex itemp = faceNodes[0];
      faceNodes[0] = faceNodes[1];
      faceNodes[1] = itemp;
    }
  }
  else
  {
    real64 ez[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
    LvArray::tensorOps::subtract< 3 >( ez, elementCenter );

    // Approximate in-plane axis
    real64 ex[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[faceNodes[0]] );
    LvArray::tensorOps::subtract< 3 >( ex, fc );
    LvArray::tensorOps::normalize< 3 >( ex );

    real64 ey[3];
    LvArray::tensorOps::crossProduct( ey, ez, ex );
    LvArray::tensorOps::normalize< 3 >( ey );

    std::pair< real64, localIndex > thetaOrder[MAX_FACE_NODES];

    // Sort nodes counterclockwise around face center
    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      real64 v[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[faceNodes[n]] );
      LvArray::tensorOps::subtract< 3 >( v, fc );
      thetaOrder[n] = std::make_pair( atan2( LvArray::tensorOps::AiBi< 3 >( v, ey ), LvArray::tensorOps::AiBi< 3 >( v, ex ) ), faceNodes[n] );
    }

    std::sort( thetaOrder, thetaOrder + numFaceNodes );
    // Reorder nodes on face
    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      faceNodes[n] = thetaOrder[n].second;
    }

    localIndex tempFaceNodes[MAX_FACE_NODES];

    localIndex firstIndexIndex = 0;
    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      tempFaceNodes[n] = thetaOrder[n].second;
      if( tempFaceNodes[n] == firstNodeIndex )
      {
        firstIndexIndex = n;
      }
    }

    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      localIndex const index = (firstIndexIndex + n) % numFaceNodes;
      faceNodes[n] = tempFaceNodes[index];
    }
  }
}

ArrayOfSets< globalIndex >
FaceManager::extractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & nodeManager )
{
  GEOS_MARK_FUNCTION;

  localIndex const numFaces = size();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = this->nodeList().toViewConst();
  arrayView1d< integer const > const isDomainBoundary = this->getDomainBoundaryIndicator();
  arrayView1d< globalIndex const > const nodeLocalToGlobal = nodeManager.localToGlobalMap();

  array1d< localIndex > nodeCounts( numFaces );
  forAll< parallelHostPolicy >( numFaces, [nodeCounts = nodeCounts.toView(),
                                           faceToNodeMap, isDomainBoundary]( localIndex const & faceIndex )
  {
    if( isDomainBoundary( faceIndex ) )
    {
      nodeCounts[ faceIndex ] = faceToNodeMap.sizeOfArray( faceIndex );
    }
  } );

  ArrayOfSets< globalIndex > globalFaceNodes;
  globalFaceNodes.resizeFromCapacities< parallelHostPolicy >( nodeCounts.size(), nodeCounts.data() );

  forAll< parallelHostPolicy >( numFaces, [globalFaceNodes = globalFaceNodes.toView(),
                                           faceToNodeMap, isDomainBoundary, nodeLocalToGlobal]( localIndex const & faceIndex )
  {
    if( isDomainBoundary( faceIndex ) )
    {
      for( localIndex const nodeIndex : faceToNodeMap[faceIndex] )
      {
        globalFaceNodes.insertIntoSet( faceIndex, nodeLocalToGlobal[nodeIndex] );
      }
    }
  } );

  return globalFaceNodes;
}


localIndex FaceManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex FaceManager::packUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex FaceManager::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                            arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::nodeListString() ) );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toNodesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToNodes,
                                               packList,
                                               this->localToGlobalMap(),
                                               m_toNodesRelation.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toEdgesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToEdges,
                                               packList,
                                               this->localToGlobalMap(),
                                               m_toEdgesRelation.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elementListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               this->m_toElements,
                                               packList,
                                               m_toElements.getElementRegionManager() );

  return packedSize;
}

localIndex FaceManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
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

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::elementListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void FaceManager::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );
}

void FaceManager::compressRelationMaps()
{
  m_toNodesRelation.compress();
  m_toEdgesRelation.compress();
}

void FaceManager::enforceStateFieldConsistencyPostTopologyChange( std::set< localIndex > const & targetIndices )
{
  arrayView1d< localIndex const > const childFaceIndices = getField< fields::childIndex >();

  ObjectManagerBase::enforceStateFieldConsistencyPostTopologyChange ( targetIndices );

  for( localIndex const targetIndex : targetIndices )
  {
    localIndex const childIndex = childFaceIndices[targetIndex];
    if( childIndex != -1 )
    {
      LvArray::tensorOps::scaledCopy< 3 >( m_faceNormal[ targetIndex ], m_faceNormal[ childIndex ], -1 );
    }
  }
}

void FaceManager::depopulateUpMaps( std::set< localIndex > const & receivedFaces,
                                    ElementRegionManager const & elemRegionManager )
{
  for( auto const & receivedFaceIdx: receivedFaces )
  {
    for( localIndex k = 0; k < m_toElements.m_toElementRegion.size( 1 ); ++k )
    {
      localIndex const elemRegionIdx    = m_toElements.m_toElementRegion[receivedFaceIdx][k];
      localIndex const elemSubRegionIdx = m_toElements.m_toElementSubRegion[receivedFaceIdx][k];
      localIndex const elemIdx          = m_toElements.m_toElementIndex[receivedFaceIdx][k];

      if( elemRegionIdx != -1 && elemSubRegionIdx != -1 && elemIdx != -1 )
      {
        CellElementSubRegion const & subRegion = elemRegionManager.getRegion( elemRegionIdx ).getSubRegion< CellElementSubRegion >( elemSubRegionIdx );
        array2d< localIndex > const & downMap = subRegion.faceList();
        bool hasTargetIndex = false;

        for( localIndex a = 0; a < downMap.size( 1 ); ++a )
        {
          localIndex const compositeLocalIdx = downMap[elemIdx][a];
          if( compositeLocalIdx == receivedFaceIdx )
          {
            hasTargetIndex = true;
          }
        }
        if( !hasTargetIndex )
        {
          m_toElements.m_toElementRegion[receivedFaceIdx][k] = -1;
          m_toElements.m_toElementSubRegion[receivedFaceIdx][k] = -1;
          m_toElements.m_toElementIndex[receivedFaceIdx][k] = -1;
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceManager, string const &, Group * const )

}
