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
 * @file EmbeddedSurfaceSubRegion.cpp
 */

#include "EmbeddedSurfaceSubRegion.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "common/MpiWrapper.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "mesh/MeshFields.hpp"


#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{
using namespace dataRepository;


void surfaceWithGhostNodes::insert( globalIndex const & edgeIndex )
{
  parentEdgeIndex.push_back( edgeIndex );
  numOfNodes += 1;
}

EmbeddedSurfaceSubRegion::EmbeddedSurfaceSubRegion( string const & name,
                                                    dataRepository::Group * const parent ):
  SurfaceElementSubRegion( name, parent ),
  m_numOfJumpEnrichments( 3 ),
  m_connectivityIndex(),
  m_parentPlaneName()
{
  m_elementType = ElementType::Polygon;

  registerWrapper( viewKeyStruct::elementCenterString(), &m_elementCenter ).
    setDescription( "The center of each EmbeddedSurface element." );

  registerWrapper( viewKeyStruct::elementVolumeString(), &m_elementVolume ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_0 ).
    setDescription( "The volume of each EmbeddedSurface element." );

  registerWrapper( viewKeyStruct::connectivityIndexString(), &m_connectivityIndex ).
    setApplyDefaultValue( 1 ).
    setDescription( "Connectivity index of each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::surfaceElementToParentPlaneString(), &m_parentPlaneName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setDescription( "A map of surface element to the parent fracture name" );

  m_normalVector.resizeDimension< 1 >( 3 );
  m_tangentVector1.resizeDimension< 1 >( 3 );
  m_tangentVector2.resizeDimension< 1 >( 3 );
  m_2dElemToElems.resize( 0, 1 );
}

void EmbeddedSurfaceSubRegion::calculateElementGeometricQuantities( NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                                                                    FaceManager const & GEOS_UNUSED_PARAM( facemanager ) )
{
  // loop over the elements
  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
  } );
}

void EmbeddedSurfaceSubRegion::calculateElementGeometricQuantities( arrayView2d< real64 const > const intersectionPoints,
                                                                    localIndex const k )
{
  for( localIndex p = 0; p < intersectionPoints.size( 0 ); p++ )
  {
    LvArray::tensorOps::add< 3 >( m_elementCenter[k], intersectionPoints[ p ] );
  }

  // update area
  m_elementArea[ k ] = computationalGeometry::ComputeSurfaceArea( intersectionPoints, m_normalVector[k] );

  LvArray::tensorOps::scale< 3 >( m_elementCenter[ k ], 1.0 / intersectionPoints.size( 0 ) );

  // update volume
  m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
}

bool EmbeddedSurfaceSubRegion::addNewEmbeddedSurface( localIndex const cellIndex,
                                                      localIndex const regionIndex,
                                                      localIndex const subRegionIndex,
                                                      NodeManager const & nodeManager,
                                                      EmbeddedSurfaceNodeManager & embSurfNodeManager,
                                                      EdgeManager const & edgeManager,
                                                      FixedOneToManyRelation const & cellToEdges,
                                                      PlanarGeometricObject const * fracture )
{
  /* The goal is to add an embeddedSurfaceElem if it is contained within the PlanarGeometricObject
   *
   * A. Identify whether the cell falls within the bounded planar object or not
   *
   * we loop over each edge:
   *
   *   1. check if it is cut by the planar object using the Dot product between the distance of each node
   *   from the origin and the normal vector.
   *   2. If an edge is cut by the planar object we look for the intersection between a line and a plane.
   *   3. Once we have the intersection we check if it falls inside the bounded planar object.
   *
   * Only elements for which all intersection points fall within the fracture plane limits will be added.
   * If the fracture does not cut through the entire element we will just chop it (it's a discretization error).
   *
   * B. Once we know the element has to be added we compute it's geometric properties:
   * - Surface Area: this is trivial given the intersection points as long as they are ordered either CW or CCW.
   * - centre:
   * - Volume:
   */

  bool addEmbeddedElem = true;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager.referencePosition();
  arrayView2d< localIndex const > const edgeToNodes = edgeManager.nodeList();
  arrayView1d< integer const > const edgeGhostRank = edgeManager.ghostRank();
  arrayView1d< globalIndex const > const & edgeLocalToGlobal = edgeManager.localToGlobalMap();

  R1Tensor origin        = fracture->getCenter();
  R1Tensor normalVector  = fracture->getNormal();
  localIndex edgeIndex;
  real64 lineDir[3], dist[3], point[3], distance[3], prodScalarProd;

  array2d< real64 > intersectionPoints( 0, 3 );
  array1d< integer > pointGhostRank;
  array1d< localIndex > pointParentIndex;

  localIndex numPoints = 0;
  for( localIndex ke = 0; ke < cellToEdges.size( 1 ); ke++ )
  {
    edgeIndex = cellToEdges[cellIndex][ke];
    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][0]] );
    LvArray::tensorOps::subtract< 3 >( dist, origin );
    prodScalarProd = LvArray::tensorOps::AiBi< 3 >( dist, normalVector );
    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][1]] );
    LvArray::tensorOps::subtract< 3 >( dist, origin );
    prodScalarProd *= LvArray::tensorOps::AiBi< 3 >( dist, normalVector );

    // check if the plane intersects the edge
    if( prodScalarProd < 0 )
    {
      LvArray::tensorOps::copy< 3 >( lineDir, nodesCoord[edgeToNodes[edgeIndex][0]] );
      LvArray::tensorOps::subtract< 3 >( lineDir, nodesCoord[edgeToNodes[edgeIndex][1]] );
      LvArray::tensorOps::normalize< 3 >( lineDir );
      //find the intersection point
      computationalGeometry::LinePlaneIntersection( lineDir,
                                                    nodesCoord[edgeToNodes[edgeIndex][0]],
                                                    normalVector,
                                                    origin,
                                                    point );

      // Check if the point is inside the fracture (bounded plane)
      if( !(fracture->isCoordInObject( point )) )
      {
        addEmbeddedElem = false;
      }
      intersectionPoints.resizeDimension< 0 >( numPoints+1 );
      pointGhostRank.resize( numPoints+1 );
      pointParentIndex.resize( numPoints+1 );
      pointGhostRank[numPoints] = edgeGhostRank[edgeIndex];
      pointParentIndex[numPoints] = edgeIndex;
      LvArray::tensorOps::copy< 3 >( intersectionPoints[numPoints], point );
      numPoints++;
    }
  } //end of edge loop


  if( addEmbeddedElem && intersectionPoints.size( 0 ) > 0 )
  {
    // resize
    localIndex surfaceIndex = this->size();
    this->resize( surfaceIndex + 1 );

    // Reorder the points CCW and then add the point to the list in the nodeManager if it is a new one.
    array1d< int > originalIndices = computationalGeometry::orderPointsCCW( intersectionPoints, normalVector );

    // Get location of embedded surfaces nodes.
    array2d< real64, nodes::REFERENCE_POSITION_PERM > & embSurfNodesPos =
      embSurfNodeManager.referencePosition();

    // fill out elemNodes array with the previously found intersection points
    // add new nodes to embSurfNodes
    bool isNew;
    localIndex nodeIndex;
    array1d< localIndex > elemNodes( intersectionPoints.size( 0 ) );

    for( localIndex j=0; j < intersectionPoints.size( 0 ); j++ )
    {
      isNew = true;
      for( localIndex h=0; h < embSurfNodeManager.size(); h++ )
      {
        LvArray::tensorOps::copy< 3 >( distance, intersectionPoints[ j ] );
        LvArray::tensorOps::subtract< 3 >( distance, embSurfNodesPos[ h ] );
        if( LvArray::tensorOps::l2Norm< 3 >( distance ) < 1e-9 )
        {
          isNew = false;
          nodeIndex = h;
          break;
        }
      }
      if( isNew )
      {
        // Add the point to the node Manager
        globalIndex parentEdgeID = edgeLocalToGlobal[ pointParentIndex[ originalIndices[ j ] ] ];
        nodeIndex = embSurfNodeManager.size();

        embSurfNodeManager.appendNode( intersectionPoints[ j ],
                                       pointGhostRank[ originalIndices[ j ] ] );

        arrayView1d< localIndex > const & parentIndex =
          embSurfNodeManager.getField< fields::parentEdgeIndex >();

        parentIndex[nodeIndex] = pointParentIndex[ originalIndices[ j ] ];

        array1d< globalIndex > & parentEdgeGlobalIndex = embSurfNodeManager.getParentEdgeGlobalIndex();
        parentEdgeGlobalIndex[nodeIndex] = parentEdgeID;
      }
      elemNodes[ j ] =  nodeIndex;
    }

    m_toNodesRelation.resizeArray( surfaceIndex, elemNodes.size() );
    for( localIndex inode = 0; inode < elemNodes.size(); inode++ )
    {
      m_toNodesRelation( surfaceIndex, inode ) = elemNodes[ inode ];
    }

    // For now 2d elements are always only connected to a single 3d element.
    if( m_2dElemToElems.m_toElementIndex[surfaceIndex].size() > 0 )
    {
      m_2dElemToElems.m_toElementIndex[surfaceIndex][0] = cellIndex;
      m_2dElemToElems.m_toElementSubRegion[surfaceIndex][0] = subRegionIndex;
      m_2dElemToElems.m_toElementRegion[surfaceIndex][0] = regionIndex;
    }
    else
    {
      m_2dElemToElems.m_toElementIndex.emplaceBack( surfaceIndex, cellIndex );
      m_2dElemToElems.m_toElementSubRegion.emplaceBack( surfaceIndex, subRegionIndex );
      m_2dElemToElems.m_toElementRegion.emplaceBack( surfaceIndex, regionIndex );
    }

    m_parentPlaneName[ surfaceIndex ] = fracture->getName();
    LvArray::tensorOps::copy< 3 >( m_normalVector[ surfaceIndex ], normalVector );
    LvArray::tensorOps::copy< 3 >( m_tangentVector1[ surfaceIndex ], fracture->getWidthVector());
    LvArray::tensorOps::copy< 3 >( m_tangentVector2[ surfaceIndex ], fracture->getLengthVector());
    this->calculateElementGeometricQuantities( intersectionPoints.toViewConst(), this->size()-1 );
  }
  return addEmbeddedElem;
}
array1d< localIndex > EmbeddedSurfaceSubRegion::getEdfmNodeParentEdgeIndex( ArrayOfArraysView< real64 > const & elemNodesLocations,
                                                                            ArrayOfArraysView< localIndex > const & elemToNodes,
                                                                            ToCellRelation< ArrayOfArrays< localIndex > > const & elemTo3dElem,
                                                                            FixedOneToManyRelation const & cellToEdges,
                                                                            arrayView2d< localIndex const > const & edgeToNodes,
                                                                            arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord )
{
  // first we build a node index to parent cell index map (we only need one)
  array1d< int > nodeToElem( elemNodesLocations.size());
  for( int i = 0; i<elemNodesLocations.size(); i++ )
  {
    nodeToElem[i] = -1;
  }
  // Assign each node to the first element that contains it
  for( int elementIndex = 0; elementIndex < elemToNodes.size(); ++elementIndex )
  {
    for( int node : elemToNodes[elementIndex] )
    {
      if( nodeToElem[node] == -1 )
      {
        nodeToElem[node] = elementIndex;
      }
    }
  }

  array1d< localIndex > parentEdgeIndex( elemNodesLocations.size());
  for( localIndex i = 0; i < elemNodesLocations.size(); ++i )
  {

    auto elementIndex =nodeToElem[i];
    auto cell3dIndex = elemTo3dElem.toCellIndex[elementIndex][0];
    double elemMinDistance = std::numeric_limits< double >::max();
    localIndex targetEdgeIndex = -1;

    R1Tensor nodexyz;
    LvArray::tensorOps::copy< 3 >( nodexyz, elemNodesLocations[i] );
    for( localIndex ke = 0; ke < cellToEdges.size( 1 ); ke++ )
    {
      R1Tensor v1;
      R1Tensor v2;
      R1Tensor crossProd;
      localIndex edgeIndex = cellToEdges[cell3dIndex][ke];
      LvArray::tensorOps::copy< 3 >( v1, nodesCoord[edgeToNodes[edgeIndex][0]] );
      LvArray::tensorOps::subtract< 3 >( v1, nodesCoord[edgeToNodes[edgeIndex][1]] );

      LvArray::tensorOps::copy< 3 >( v2, nodesCoord[edgeToNodes[edgeIndex][1]] );
      LvArray::tensorOps::subtract< 3 >( v2, nodexyz );

      LvArray::tensorOps::crossProduct( crossProd, v1, v2 );
      auto numeratorNorm = LvArray::tensorOps::l2Norm< 3 >( crossProd );
      auto denomNorm = LvArray::tensorOps::l2Norm< 3 >( v1 );
      auto dist = numeratorNorm / denomNorm;
      if( dist < elemMinDistance )
      {
        elemMinDistance = dist;
        targetEdgeIndex = edgeIndex;
      }
    }

    parentEdgeIndex[i] = targetEdgeIndex;
  }
  return parentEdgeIndex;
}

bool EmbeddedSurfaceSubRegion::addAllEmbeddedSurfaces(
  localIndex const regionIndex,
  localIndex const subRegionIndex,
  NodeManager const & nodeManager,
  EmbeddedSurfaceNodeManager & embSurfNodeManager,
  EdgeManager const & edgeManager,
  FixedOneToManyRelation const & cellToEdges,
  EmbeddedSurfaceBlockABC const & embeddedSurfaceBlock )
{

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager.referencePosition();
  arrayView2d< localIndex const > const edgeToNodes = edgeManager.nodeList();
  arrayView1d< integer const > const edgeGhostRank = edgeManager.ghostRank();
  arrayView1d< globalIndex const > const & edgeLocalToGlobal = edgeManager.localToGlobalMap();

  localIndex const numElems = embeddedSurfaceBlock.numEmbeddedSurfElem();

  resize( embeddedSurfaceBlock.numEmbeddedSurfElem());

  ToCellRelation< ArrayOfArrays< localIndex > > elemTo3dElem = embeddedSurfaceBlock.getEmbeddedSurfElemTo3dElem();
  ArrayOfArrays< localIndex > elemToNodes = embeddedSurfaceBlock.getEmbeddedSurfElemToNodes();
  ArrayOfArrays< real64 > elemNodesLocations = embeddedSurfaceBlock.getEmbeddedSurfElemNodesCoords();
  ArrayOfArrays< real64 > elemNormalVectors = embeddedSurfaceBlock.getEmbeddedSurfElemNormalVectors();
  ArrayOfArrays< real64 > elemTangentialVectors1 = embeddedSurfaceBlock.getEmbeddedSurfElemTangentialWidthVectors();
  ArrayOfArrays< real64 > elemTangentialVectors2 = embeddedSurfaceBlock.getEmbeddedSurfElemTangentialLengthVectors();

  localIndex numEdfmNodes = elemNodesLocations.size();
  array1d< localIndex > allPointsParentIndex =  getEdfmNodeParentEdgeIndex( elemNodesLocations.toView(), elemToNodes.toView(), elemTo3dElem, cellToEdges, edgeToNodes, nodesCoord );

  // EDFM nodes get their gost rank and global parent index
  // from their parent 3d element local edge index.
  array1d< integer > allPointsGhostRank( numEdfmNodes );
  array1d< integer > allPointsGlobalParenIndex( numEdfmNodes );
  for( auto i = 0; i < numEdfmNodes; ++i )
  {
    integer pointGhostRank = edgeGhostRank[allPointsParentIndex[i]];
    localIndex pointLocalParentEdgeIndex = allPointsParentIndex[i];
    globalIndex pointGlobalParentEdgeIndex = edgeLocalToGlobal[allPointsParentIndex[i]];

    embSurfNodeManager.appendNode( elemNodesLocations[i].toSliceConst(), pointGhostRank );

    arrayView1d< localIndex > const & parentIndex = embSurfNodeManager.getField< fields::parentEdgeIndex >();

    parentIndex[i] = pointLocalParentEdgeIndex;

    array1d< globalIndex > & parentEdgeGlobalIndex = embSurfNodeManager.getParentEdgeGlobalIndex();
    parentEdgeGlobalIndex[i] = pointGlobalParentEdgeIndex;
  }

  for( auto i = 0; i < numElems; ++i )
  {
    m_toNodesRelation.resizeArray( i, 4 );
    localIndex elem3dIndex = elemTo3dElem.toCellIndex[i][0];
    //localIndex elemRegionIndex = elemTo3dElem.toBlockIndex[i][0];
    // region and subregion will be filled later.
    m_2dElemToElems.m_toElementIndex.emplaceBack( i, elem3dIndex );
    m_2dElemToElems.m_toElementSubRegion.emplaceBack( i, subRegionIndex );
    m_2dElemToElems.m_toElementRegion.emplaceBack( i, regionIndex );

    // we only accept quads
    array2d< real64 > elemNodeCoords( 4, 3 );
    for( localIndex inode = 0; inode < 4; inode++ )
    {
      auto nodeGlobalIndex = elemToNodes[i][inode];
      m_toNodesRelation( i, inode ) = nodeGlobalIndex;
      LvArray::tensorOps::copy< 3 >( elemNodeCoords[inode], elemNodesLocations[nodeGlobalIndex] );
    }

    m_parentPlaneName[i] = "elem_" + std::to_string( i );
    LvArray::tensorOps::copy< 3 >( m_normalVector[i], elemNormalVectors[i] );
    LvArray::tensorOps::copy< 3 >( m_tangentVector1[i], elemTangentialVectors1[i] );
    LvArray::tensorOps::copy< 3 >( m_tangentVector2[i], elemTangentialVectors2[i] );
    this->calculateElementGeometricQuantities( elemNodeCoords.toViewConst(), i );
  }

  return true;

}

void EmbeddedSurfaceSubRegion::inheritGhostRank( array1d< array1d< arrayView1d< integer const > > > const & cellGhostRank )
{
  arrayView1d< integer > const & ghostRank = this->ghostRank();
  for( localIndex k=0; k < size(); ++k )
  {
    localIndex regionIndex    = m_2dElemToElems.m_toElementRegion[k][0];
    localIndex subRegionIndex = m_2dElemToElems.m_toElementSubRegion[k][0];
    localIndex cellIndex      = m_2dElemToElems.m_toElementIndex[k][0];

    ghostRank[k] = cellGhostRank[regionIndex][subRegionIndex][cellIndex];
  }
}

void EmbeddedSurfaceSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getEmbSurfNodeManager() );
}

localIndex EmbeddedSurfaceSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex EmbeddedSurfaceSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex EmbeddedSurfaceSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                         arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               nodeList().base().toViewConst(),
                                               m_unmappedGlobalIndicesInToNodes,
                                               packList,
                                               localToGlobal,
                                               nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::surfaceElementsToCellRegionsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               this->m_2dElemToElems,
                                               packList,
                                               m_2dElemToElems.getElementRegionManager() );

  return packedSize;
}


localIndex EmbeddedSurfaceSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                       localIndex_array & packList,
                                                       bool const overwriteUpMaps,
                                                       bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
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

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::surfaceElementsToCellRegionsString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToElems,
                                     packList.toViewConst(),
                                     m_2dElemToElems.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

} /* namespace geos */
