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

#include "WellElementSubRegion.hpp"

#include "mesh/MeshLevel.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "common/MpiWrapper.hpp"
#include "LvArray/src/output.hpp"

#include <unordered_set>

namespace geos
{

WellElementSubRegion::WellElementSubRegion( string const & name, Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_wellControlsName( "" ),
  m_toNodesRelation(),
  m_topWellElementIndex( -1 ),
  m_perforationData( groupKeyStruct::perforationDataString(), this ),
  m_topRank( -1 ),
  m_searchDepth( 10 )
{
  m_elementType = ElementType::Line;

  registerWrapper( viewKeyStruct::wellControlsString(), &m_wellControlsName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef );

  registerWrapper( viewKeyStruct::wellNodeListString(), &m_toNodesRelation );
  registerWrapper( viewKeyStruct::nextWellElementIndexString(), &m_nextWellElementIndex );
  registerWrapper( viewKeyStruct::nextWellElementIndexGlobalString(), &m_nextWellElementIndexGlobal );
  registerWrapper( viewKeyStruct::topWellElementIndexString(), &m_topWellElementIndex );
  registerWrapper( viewKeyStruct::topRankString(), &m_topRank );
  registerWrapper( viewKeyStruct::radiusString(), &m_radius );

  registerGroup( groupKeyStruct::perforationDataString(), &m_perforationData );

  excludeWrappersFromPacking( { viewKeyStruct::nodeListString() } );

  m_numNodesPerElement = 2;
  m_numFacesPerElement = 0;
  m_toNodesRelation.resizeDimension< 1 >( this->numNodesPerElement() );
}

void WellElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
}

namespace
{

/**
 * @brief Now that the well elements are assigned, collect the nodes and tag the boundary nodes between ranks
          The function WellElementSubRegion::AssignUnownedElements must have been called before this function
 * @param[in] lineBlock the LineBlockABC containing the global well topology
 * @param[in] localElems set of local well elems. At this point all the well elems have been assigned
 * @param[out] localNodes set of local well nodes (includes boundary nodes)
 * @param[out] boundaryNodes set of local well nodes that are at the boundary between this rank
               and another rank
 */
void collectLocalAndBoundaryNodes( LineBlockABC const & lineBlock,
                                   SortedArray< globalIndex >      const & localElems,
                                   SortedArray< globalIndex > & localNodes,
                                   SortedArray< globalIndex > & boundaryNodes )
{
  // get the well connectivity
  arrayView1d< globalIndex const >                      const & nextElemIdGlobal  = lineBlock.getNextElemIndex();
  arrayView1d< arrayView1d< globalIndex const > const > const & prevElemIdsGlobal = lineBlock.getPrevElemIndices();
  arrayView2d< globalIndex const >                      const & elemToNodesGlobal = lineBlock.getElemToNodesMap();

  // loop over the local elements and collect the local and boundary nodes
  for( globalIndex currGlobal : localElems )
  {

    // if the element is local, its two nodes are also local
    globalIndex const inodeTopGlobal    = elemToNodesGlobal[currGlobal][LineBlockABC::NodeLocation::TOP];
    globalIndex const inodeBottomGlobal = elemToNodesGlobal[currGlobal][LineBlockABC::NodeLocation::BOTTOM];
    localNodes.insert( inodeTopGlobal );
    localNodes.insert( inodeBottomGlobal );

    localIndex const nextGlobal =
      LvArray::integerConversion< localIndex >( nextElemIdGlobal[ LvArray::integerConversion< localIndex >( currGlobal ) ] );

    // if the next well elem is not local, add the node in between curr and next to boundaryNodes
    if( nextGlobal >= 0 && !localElems.contains( nextGlobal ))
    {
      boundaryNodes.insert( inodeTopGlobal );
    }

    // if the prev well elem is not local, add the node in between curr and prev to boundaryNodes (relevant for
    // branches)
    for( localIndex iwelem = 0; iwelem < prevElemIdsGlobal[currGlobal].size(); ++iwelem )
    {
      globalIndex const prevGlobal = prevElemIdsGlobal[currGlobal][iwelem];
      if( prevGlobal >= 0 && !localElems.contains( prevGlobal ))
      {
        boundaryNodes.insert( inodeBottomGlobal );
      }
    }
  }
}

/**
 * @brief Collect the nodes of reservoir element ei
 * @param[in] subRegion the subRegion of reservoir element ei
 * @param[in] ei the index of the reservoir element
 * @param[inout] nodes the nodes that have already been visited
 */
template< typename SUBREGION_TYPE >
void collectElementNodes( SUBREGION_TYPE const & subRegion,
                          localIndex ei,
                          SortedArray< localIndex > & nodes )
{
  // get all the nodes belonging to this element
  for( localIndex a = 0; a < subRegion.numNodesPerElement(); ++a )
  {
    localIndex const inode = subRegion.nodeList( ei, a );

    // if not already visited, store the newly found node
    if( !nodes.contains( inode ))
    {
      nodes.insert( inode );
    }
  }
}

template< typename SUBREGION_TYPE >
bool isPointInsideElement( SUBREGION_TYPE const & GEOS_UNUSED_PARAM( subRegion ),
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & GEOS_UNUSED_PARAM( referencePosition ),
                           localIndex const & GEOS_UNUSED_PARAM( eiLocal ),
                           ArrayOfArraysView< localIndex const > const & GEOS_UNUSED_PARAM( facesToNodes ),
                           real64 const (&GEOS_UNUSED_PARAM( location ))[3] )
{
  // only CellElementSubRegion is currently supported
  return false;
}

bool isPointInsideElement( CellElementSubRegion const & subRegion,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePosition,
                           localIndex const & eiLocal,
                           ArrayOfArraysView< localIndex const > const & facesToNodes,
                           real64 const (&location)[3] )
{
  arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
  arrayView2d< real64 const > const elemCenters = subRegion.getElementCenter();
  real64 const elemCenter[3] = { elemCenters[eiLocal][0],
                                 elemCenters[eiLocal][1],
                                 elemCenters[eiLocal][2] };
  return computationalGeometry::isPointInsidePolyhedron( referencePosition,
                                                         elemsToFaces[eiLocal],
                                                         facesToNodes,
                                                         elemCenter,
                                                         location );
}

// Define a hash function
template< typename POINT_TYPE >
struct PointHash
{
  std::size_t operator()( POINT_TYPE const point ) const
  {
    std::size_t h1 = std::hash< real64 >()( point[0] );
    std::size_t h2 = std::hash< real64 >()( point[1] );
    std::size_t h3 = std::hash< real64 >()( point[2] );
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

// Define equality operator
template< typename POINT_TYPE >
struct PointsEqual
{
  bool operator()( POINT_TYPE const & p1, POINT_TYPE const & p2 ) const
  {
    return (std::abs( p1[0] - p2[0] ) < 1e-10) && (std::abs( p1[1] - p2[1] ) < 1e-10) && (std::abs( p1[2] - p2[2] ) < 1e-10);
  }
};

bool isPointInsideElement( SurfaceElementSubRegion const & subRegion,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePosition,
                           localIndex const & eiLocal,
                           ArrayOfArraysView< localIndex const > const & GEOS_UNUSED_PARAM( facesToNodes ),
                           real64 const (&location)[3] )
{
  typedef std::array< real64, 3 > Point3d;

  // collect element nodes
  integer const nV = subRegion.numNodesPerElement( eiLocal );
  SurfaceElementSubRegion::NodeMapType const & nodeList = subRegion.nodeList();
  std::vector< Point3d > polygon( nV );
  for( integer i = 0; i < nV; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      polygon[i][j] = referencePosition[nodeList( eiLocal, i )][j];
    }
  }

  // remove duplicates
  std::unordered_set< Point3d, PointHash< Point3d >, PointsEqual< Point3d > >
  unique_points( polygon.begin(), polygon.end());
  polygon.clear();
  std::copy( unique_points.begin(), unique_points.end(), std::back_inserter( polygon ));

  return computationalGeometry::isPointInPolygon3d( polygon, polygon.size(), location );
}

/**
 * @brief Search the reservoir elements that can be accessed from the set "nodes".
          Stop if a reservoir element containing the perforation is found.
          If not, enlarge the set "nodes"
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[inout] nodes the nodes that have already been visited
 * @param[inout] elements the reservoir elements that have already been visited
 * @param[inout] erMatched the region index of the reservoir element that contains "location", if any
 * @param[inout] esrMatched the subregion index of the reservoir element that contains "location", if any
 * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
 * @param[inout] giMatched the element global index of the reservoir element that contains "location", if any
 */
bool visitNeighborElements( MeshLevel const & mesh,
                            real64 const (&location)[3],
                            SortedArray< localIndex > & nodes,
                            SortedArray< globalIndex > & elements,
                            localIndex & erMatched,
                            localIndex & esrMatched,
                            localIndex & eiMatched,
                            globalIndex & giMatched )
{
  ElementRegionManager const & elemManager = mesh.getElemManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  ArrayOfArraysView< localIndex const > const & toElementRegionList    = nodeManager.elementRegionList();
  ArrayOfArraysView< localIndex const > const & toElementSubRegionList = nodeManager.elementSubRegionList();
  ArrayOfArraysView< localIndex const > const & toElementList          = nodeManager.elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition =
    nodeManager.referencePosition().toViewConst();

  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

  bool matched = false;

  // In this function, we loop over all the nodes that we have collected so far.
  // For each node, we grab all the reservoir elements that contain the node
  //     For reservoir element that has not been visited yet, we check if it contains "location"
  //           If the reservoir element contains "location" we can stop the search
  //           If the reservoir element does not contain "location", we collect its nodes

  // we will enlarge the set of nodes in the loop below
  // to do this we have to create a new set, "currNodes"
  // that contains only the nodes that have already been visited
  // the newly added nodes will be added to "nodes"
  SortedArray< localIndex > currNodes = nodes;
  giMatched = -1;
  // for all the nodes already visited
  for( localIndex currNode : currNodes )
  {
    // collect the elements that have not been visited yet
    for( localIndex b=0; b<toElementRegionList.sizeOfArray( currNode ); ++b )
    {
      localIndex const er      = toElementRegionList[currNode][b];
      localIndex const esr     = toElementSubRegionList[currNode][b];
      localIndex const eiLocal = toElementList[currNode][b];
      CellElementRegion const & region = elemManager.getRegion< CellElementRegion >( er );
      CellElementSubRegion const & subRegion = region.getSubRegion< CellElementSubRegion >( esr );
      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
      arrayView2d< real64 const > const elemCenters = subRegion.getElementCenter();

      globalIndex const eiGlobal = subRegion.localToGlobalMap()[eiLocal];

      // if this element has not been visited yet, save it
      if( !elements.contains( eiGlobal ))
      {
        elements.insert( eiGlobal );

        real64 const elemCenter[3] = { elemCenters[eiLocal][0],
                                       elemCenters[eiLocal][1],
                                       elemCenters[eiLocal][2] };

        // perform the test to see if the point is in this reservoir element
        // if the point is in the resevoir element, save the indices and stop the search
        if( computationalGeometry::isPointInsidePolyhedron( referencePosition,
                                                            elemsToFaces[eiLocal],
                                                            facesToNodes,
                                                            elemCenter,
                                                            location ) )
        {
          erMatched  = er;
          esrMatched = esr;
          eiMatched  = eiLocal;
          giMatched  = eiGlobal;
          matched    = true;
          break;
        }
        // otherwise add the nodes of this element to the set of new nodes to visit
        else
        {
          collectElementNodes( subRegion, eiLocal, nodes );
        }
      }
    }

    if( matched )
    {
      break;
    }
  }

  // if not matched, insert the new nodes
  return matched;
}

/**
 * @brief Search the reservoir elements that can be accessed from the set "nodes".
          Stop if a reservoir element containing the perforation is found.
          If not, enlarge the set "nodes"
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[inout] nodes the nodes that have already been visited
 * @param[inout] elements the reservoir elements that have already been visited
 * @param[in] targetRegionIndex the target region index for the reservoir element
 * @param[in] targetSubRegionIndex the target subregion index for the reservoir element
 * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
 * @param[inout] giMatched the element global index of the reservoir element that contains "location", if any
 */
template< typename SUBREGION_TYPE >
bool visitNeighborElements( MeshLevel const & mesh,
                            real64 const (&location)[3],
                            SortedArray< localIndex > & nodes,
                            SortedArray< globalIndex > & elements,
                            localIndex const & targetRegionIndex,
                            localIndex const & targetSubRegionIndex,
                            localIndex & eiMatched,
                            globalIndex & giMatched )
{
  ElementRegionManager const & elemManager = mesh.getElemManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  ArrayOfArraysView< localIndex const > const & toElementRegionList    = nodeManager.elementRegionList();
  ArrayOfArraysView< localIndex const > const & toElementSubRegionList = nodeManager.elementSubRegionList();
  ArrayOfArraysView< localIndex const > const & toElementList          = nodeManager.elementList();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition =
    nodeManager.referencePosition().toViewConst();

  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

  bool matched = false;

  // In this function, we loop over all the nodes that we have collected so far.
  // For each node, we grab all the reservoir elements that contain the node
  //     For reservoir element that has not been visited yet, we check if it contains "location"
  //           If the reservoir element contains "location" we can stop the search
  //           If the reservoir element does not contain "location", we collect its nodes

  // we will enlarge the set of nodes in the loop below
  // to do this we have to create a new set, "currNodes"
  // that contains only the nodes that have already been visited
  // the newly added nodes will be added to "nodes"
  SortedArray< localIndex > currNodes = nodes;
  giMatched = -1;
  // for all the nodes already visited
  for( localIndex currNode : currNodes )
  {
    // collect the elements that have not been visited yet
    for( localIndex b = 0; b < toElementRegionList.sizeOfArray( currNode ); ++b )
    {
      localIndex const er      = toElementRegionList[currNode][b];
      localIndex const esr     = toElementSubRegionList[currNode][b];
      localIndex const eiLocal = toElementList[currNode][b];

      if( er != targetRegionIndex || esr != targetSubRegionIndex )
        continue;

      ElementRegionBase const & region = elemManager.getRegion< ElementRegionBase >( er );
      SUBREGION_TYPE const & subRegion = region.getSubRegion< SUBREGION_TYPE >( esr );
      globalIndex const eiGlobal = subRegion.localToGlobalMap()[eiLocal];

      // if this element has not been visited yet, save it
      if( !elements.contains( eiGlobal ))
      {
        elements.insert( eiGlobal );

        // perform the test to see if the point is in this reservoir element
        // if the point is in the resevoir element, save the indices and stop the search
        if( isPointInsideElement( subRegion, referencePosition, eiLocal, facesToNodes, location ) )
        {
          eiMatched = eiLocal;
          giMatched = eiGlobal;
          matched   = true;
          break;
        }
        // otherwise add the nodes of this element to the set of new nodes to visit
        else
        {
          collectElementNodes( subRegion, eiLocal, nodes );
        }
      }
    }

    if( matched )
    {
      break;
    }
  }

  // if not matched, insert the new nodes
  return matched;
}

/**
 * @brief Search for the reservoir element that is the *closest* from the center of well element.
          Note that this reservoir element does not necessarily contain the center of the well element.
          This "init" reservoir element will be used in SearchLocalElements to find the reservoir element that
          contains the well element.
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[inout] erInit the region index of the reservoir element from which we start the search
 * @param[inout] esrInit the subregion index of the reservoir element from which we start the search
 * @param[inout] eiInit the element index of the reservoir element from which we start the search
 */
void initializeLocalSearch( MeshLevel const & mesh,
                            real64 const (&location)[3],
                            localIndex & erInit,
                            localIndex & esrInit,
                            localIndex & eiInit )
{
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
  resElemCenter = mesh.getElemManager().constructViewAccessor< array2d< real64 >,
                                                               arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );
  // to initialize the local search for the reservoir element that contains "location",
  // we find the reservoir element that minimizes the distance from "location" to the reservoir element center
  auto ret = minLocOverElemsInMesh( mesh, [&] ( localIndex const er,
                                                localIndex const esr,
                                                localIndex const ei )
  {
    real64 v[3] = { location[0], location[1], location[2] };
    LvArray::tensorOps::subtract< 3 >( v, resElemCenter[er][esr][ei] );
    return LvArray::tensorOps::l2Norm< 3 >( v );
  } );

  // save the region, subregion and index of the reservoir element
  // note that this reservoir element does not necessarily contains "location"
  erInit  = std::get< 0 >( ret.second );
  esrInit = std::get< 1 >( ret.second );
  eiInit  = std::get< 2 >( ret.second );
}

/**
 * @brief Search for the reservoir element that is the *closest* from the center of well element.
          Note that this reservoir element does not necessarily contain the center of the well element.
          This "init" reservoir element will be used in SearchLocalElements to find the reservoir element that
          contains the well element.
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[in] targetRegionIndex the region index of the reservoir element from which we start the search
 * @param[in] targetSubRegionIndex the subregion index of the reservoir element from which we start the search
 * @param[inout] eiInit the element index of the reservoir element from which we start the search
 */
void initializeLocalSearch( MeshLevel const & mesh,
                            real64 const (&location)[3],
                            localIndex const & targetRegionIndex,
                            localIndex const & targetSubRegionIndex,
                            localIndex & eiInit )
{
  ElementSubRegionBase const & subRegion = mesh.getElemManager().getRegion( targetRegionIndex ).getSubRegion( targetSubRegionIndex );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
  resElemCenter = mesh.getElemManager().constructViewAccessor< array2d< real64 >,
                                                               arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );
  // to initialize the local search for the reservoir element that contains "location",
  // we find the reservoir element that minimizes the distance from "location" to the reservoir element center
  auto ret = minLocOverElemsInSubRegion( subRegion, [&] ( localIndex const ei )
  {
    real64 v[3] = { location[0], location[1], location[2] };
    LvArray::tensorOps::subtract< 3 >( v, resElemCenter[targetRegionIndex][targetSubRegionIndex][ei] );
    auto dist = LvArray::tensorOps::l2Norm< 3 >( v );
    return dist;
  } );

  // save the index of the reservoir element
  // note that this reservoir element does not necessarily contains "location"
  eiInit  = ret.second;
}

/**
 * @brief Search for the reservoir element that contains the well element.
          To do that, loop over the reservoir elements that are in the neighborhood of (erInit,esrInit,eiInit)
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[in] erInit the region index of the reservoir element from which we start the search
 * @param[in] esrInit the subregion index of the reservoir element from which we start the search
 * @param[in] eiInit the element index of the reservoir element from which we start the search
 * @param[inout] erMatched the region index of the reservoir element that contains "location", if any
 * @param[inout] esrMatched the subregion index of the reservoir element that contains "location", if any
 * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
 * @param[inout] giMatched the element global index of the reservoir element that contains "location", if any
 */
bool searchLocalElements( MeshLevel const & mesh,
                          real64 const (&location)[3],
                          localIndex const & searchDepth,
                          localIndex const & erInit,
                          localIndex const & esrInit,
                          localIndex const & eiInit,
                          localIndex & erMatched,
                          localIndex & esrMatched,
                          localIndex & eiMatched,
                          globalIndex & giMatched )
{
  // search locally, starting from the location of the previous perforation
  // the assumption here is that perforations have been entered in order of depth
  bool resElemFound = false;

  CellElementRegion const & region = mesh.getElemManager().getRegion< CellElementRegion >( erInit );
  CellElementSubRegion const & subRegion = region.getSubRegion< CellElementSubRegion >( esrInit );

  SortedArray< localIndex >  nodes;
  SortedArray< globalIndex > elements;

  // here is how the search is done:
  //   1 - We check if "location" is within the "init" reservoir element defined by (erInit,esrMatched,eiMatched)
  //   2 - If yes, stop
  //     - If not, a) collect the nodes of the reservoir element defined by (erInit,esrMatched,eiMatched)
  //               b) use these nodes to grab the neighbors of (erInit,esrMatched,eiMatched)
  //               c) check if "location" is within the neighbors. If not, grab the neighbors of the neighbors, and so
  // on...

  // collect the nodes of the current element
  // they will be used to access the neighbors and check if they contain the perforation
  collectElementNodes( subRegion, eiInit, nodes );
  // if no match is found, enlarge the neighborhood m_searchDepth'th times
  for( localIndex d = 0; d < searchDepth; ++d )
  {
    localIndex nNodes = nodes.size();

    // search the reservoir elements that can be accessed from the set "nodes"
    // stop if a reservoir element containing the perforation is found
    // if not, enlarge the set "nodes"
    resElemFound = visitNeighborElements( mesh, location, nodes, elements,
                                          erMatched, esrMatched, eiMatched, giMatched );
    if( resElemFound || nNodes == nodes.size())
    {
      break;
    }
  }
  return resElemFound;
}

/**
 * @brief Search for the reservoir element that contains the well element.
          To do that, loop over the reservoir elements that are in the neighborhood of (erInit,esrInit,eiInit)
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[in] targetRegionIndex the target region index for the reservoir element
 * @param[inout] esrMatched the subregion index of the reservoir element that contains "location", if any
 * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
 * @param[inout] giMatched the element global index of the reservoir element that contains "location", if any
 */
bool searchLocalElements( MeshLevel const & mesh,
                          real64 const (&location)[3],
                          localIndex const & searchDepth,
                          localIndex const & targetRegionIndex,
                          localIndex & esrMatched,
                          localIndex & eiMatched,
                          globalIndex & giMatched )
{
  ElementRegionBase const & region = mesh.getElemManager().getRegion< ElementRegionBase >( targetRegionIndex );

  bool resElemFound = false;
  for( localIndex esr = 0; esr < region.numSubRegions(); ++esr )
  {
    ElementSubRegionBase const & subRegionBase = region.getSubRegion( esr );
    region.applyLambdaToContainer< CellElementSubRegion, SurfaceElementSubRegion >( subRegionBase, [&]( auto const & subRegion )
    {
      GEOS_LOG( GEOS_FMT( "  searching well connections with region/subregion: {}/{}", region.getName(), subRegion.getName() ) );

      // first, we search for the reservoir element that is the *closest* from the center of well element
      // note that this reservoir element does not necessarily contain the center of the well element
      // this "init" reservoir element will be used later to find the reservoir element that
      // contains the well element
      localIndex eiInit     = -1;
      initializeLocalSearch( mesh, location, targetRegionIndex, esr, eiInit );

      // loop over the reservoir elements that are in the neighborhood of (esrInit,eiInit)
      // search locally, starting from the location of the previous perforation
      // the assumption here is that perforations have been entered in order of depth

      SortedArray< localIndex >  nodes;
      SortedArray< globalIndex > elements;

      // here is how the search is done:
      //   1 - We check if "location" is within the "init" reservoir element defined by (erInit,esrMatched,eiMatched)
      //   2 - If yes, stop
      //     - If not, a) collect the nodes of the reservoir element defined by (erInit,esrMatched,eiMatched)
      //               b) use these nodes to grab the neighbors of (erInit,esrMatched,eiMatched)
      //               c) check if "location" is within the neighbors. If not, grab the neighbors of the neighbors, and so
      // on...

      // collect the nodes of the current element
      // they will be used to access the neighbors and check if they contain the perforation
      collectElementNodes( subRegion, eiInit, nodes );

      // if no match is found, enlarge the neighborhood m_searchDepth'th times
      for( localIndex d = 0; d < searchDepth; ++d )
      {
        localIndex nNodes = nodes.size();

        // search the reservoir elements that can be accessed from the set "nodes"
        // stop if a reservoir element containing the perforation is found
        // if not, enlarge the set "nodes"

        resElemFound =
          visitNeighborElements< TYPEOFREF( subRegion ) >( mesh, location, nodes, elements, targetRegionIndex, esr, eiMatched, giMatched );

        if( resElemFound || nNodes == nodes.size())
        {
          esrMatched = esr;
          GEOS_LOG( GEOS_FMT( "    found {}/{}/{}", region.getName(), subRegion.getName(), eiMatched ) );
          // TODO learn how to exit forElementSubRegionsIndex
          break;
        }
      }
    } );

    if( resElemFound )
    {
      break;
    }
  }

  return resElemFound;
}

}

void WellElementSubRegion::generate( MeshLevel & mesh,
                                     LineBlockABC const & lineBlock,
                                     arrayView1d< integer > & elemStatusGlobal,
                                     globalIndex nodeOffsetGlobal,
                                     globalIndex elemOffsetGlobal )
{

  map< integer, SortedArray< globalIndex > > elemSetsByStatus;

  // convert elemStatus list into sets of indices
  for( localIndex iwelemGlobal = 0; iwelemGlobal < elemStatusGlobal.size(); ++iwelemGlobal )
  {
    elemSetsByStatus[elemStatusGlobal[iwelemGlobal]].insert( iwelemGlobal );
  }

  // initialize the sets using the classification of well elems
  // localElems will be enlarged once boundary elements ownership is determined
  SortedArray< globalIndex > & localElems   = elemSetsByStatus[WellElemStatus::LOCAL];
  SortedArray< globalIndex > & sharedElems  = elemSetsByStatus[WellElemStatus::SHARED];
  SortedArray< globalIndex > & unownedElems = elemSetsByStatus[WellElemStatus::UNOWNED];

  // here we make sure that there are no shared elements
  // this is enforced in the LineBlockABC that currently merges two perforations
  // if they belong to the same well element. This is a temporary solution.
  // TODO: split the well elements that contain multiple perforations, so that no element is shared
  GEOS_THROW_IF( sharedElems.size() > 0,
                 "Well " << lineBlock.getDataContext() << " contains shared well elements",
                 InputError );

  // In Steps 1 and 2 we determine the local objects on this rank (elems and nodes)
  // Once this is done, in Steps 3, 4, and 5, we update the nodeManager and wellElementSubRegion (size, maps)


  // 1) First assign the unowned elements to a rank
  // this is done in two steps

  // 1.a) First assign unowned elements in the reservoir based on location
  //      ie., if the center of the well element falls in the domain owned by rank k
  //      then the well element is assigned to rank k
  assignUnownedElementsInReservoir( mesh,
                                    lineBlock,
                                    unownedElems,
                                    localElems,
                                    elemStatusGlobal );
  // 1.b) Then we check that all the well elements have been assigned (and assigned once)
  //      This is needed because if the center of the well element falls on the boundary of
  //      a reservoir element, the assignment algorithm of 1.a) can assign the same well element
  //      to two ranks, or to no rank at all (which will break the solver).
  //      In this function we also check that the resulting well partitioning is valid, that is,
  //      we make sure that if two ranks are neighbors in the well, that are also neighbors in the
  //      reservoir mesh
  checkPartitioningValidity( lineBlock,
                             localElems,
                             elemStatusGlobal );

  SortedArray< globalIndex > localNodes;
  SortedArray< globalIndex > boundaryNodes;

  // 2) collect the local nodes and tag the boundary nodes using element info
  // now that all the elements have been assigned, we collected the local nodes
  // and tag the boundary nodes (i.e., the nodes in contact with both local and remote elems)
  collectLocalAndBoundaryNodes( lineBlock,
                                localElems,
                                localNodes,
                                boundaryNodes );

  // 3) size update in the nodeManager
  // this is necessary to later use the node matching procedure
  // to place ghosts in DomainPartition::SetupCommunications
  updateNodeManagerSize( mesh,
                         lineBlock,
                         localNodes,
                         boundaryNodes,
                         nodeOffsetGlobal );

  // 4) resize the well element subregion
  // and construct local to global, global to local, maps, etc
  constructSubRegionLocalElementMaps( mesh,
                                      lineBlock,
                                      localElems,
                                      nodeOffsetGlobal,
                                      elemOffsetGlobal );

  // 5) node-to-elem map update in the nodeManager
  // This map will be used by MeshLevel::GenerateAdjacencyLists
  // this assumes that the elemToNodes maps has been filled at Step 5)
  updateNodeManagerNodeToElementMap( mesh );

  // Store local to global index mapping
  integer n_localElems = localElems.size();
  m_globalWellElementIndex.resize( n_localElems );
  for( integer i=0; i<n_localElems; i++ )
  {
    m_globalWellElementIndex[i] = localElems[i];
  }

}


void WellElementSubRegion::assignUnownedElementsInReservoir( MeshLevel & mesh,
                                                             LineBlockABC const & lineBlock,
                                                             SortedArray< globalIndex > const & unownedElems,
                                                             SortedArray< globalIndex > & localElems,
                                                             arrayView1d< integer > & elemStatusGlobal ) const
{
  // get the well and reservoir element coordinates
  arrayView2d< real64 const > const & wellElemCoordsGlobal = lineBlock.getElemCoords();

  // assign the well elements based on location wrt the reservoir elements
  // if the center of the well element falls in the domain owned by rank k
  // then the well element is assigned to rank k
  for( globalIndex currGlobal : unownedElems )
  {
    real64 const location[3] = { wellElemCoordsGlobal[currGlobal][0],
                                 wellElemCoordsGlobal[currGlobal][1],
                                 wellElemCoordsGlobal[currGlobal][2] };

    // this will contain the indices of the reservoir element
    // in which the center of the well element is located
    localIndex erMatched  = -1;
    localIndex esrMatched = -1;
    localIndex eiMatched  = -1;

    // this will contain the indices of the reservoir element
    // from which we are going to start the search
    localIndex erInit     = -1;
    localIndex esrInit    = -1;
    localIndex eiInit     = -1;
    globalIndex giMatched   = -1;

    // Step 1: first, we search for the reservoir element that is the *closest* from the center of well element
    //         note that this reservoir element does not necessarily contain the center of the well element
    //         this "init" reservoir element will be used in SearchLocalElements to find the reservoir element that
    //         contains the well element
    initializeLocalSearch( mesh, location,
                           erInit, esrInit, eiInit );

    // Step 2: then, search for the reservoir element that contains the well element
    //         to do that, we loop over the reservoir elements that are in the neighborhood of (erInit,esrInit,eiInit)
    bool resElemFound = searchLocalElements( mesh, location, m_searchDepth,
                                             erInit, esrInit, eiInit,
                                             erMatched, esrMatched, eiMatched, giMatched );

    // if the element was found
    if( resElemFound )
    {
      // the well element is in the reservoir element (erMatched,esrMatched,eiMatched), so tag it as local
      localElems.insert( currGlobal );
      elemStatusGlobal[currGlobal] = WellElemStatus::LOCAL;
    }
  }
}


void WellElementSubRegion::checkPartitioningValidity( LineBlockABC const & lineBlock,
                                                      SortedArray< globalIndex > & localElems,
                                                      arrayView1d< integer > & elemStatusGlobal ) const
{
  arrayView1d< arrayView1d< globalIndex const > const > const & prevElemIdsGlobal = lineBlock.getPrevElemIndices();

  // we are going to make sure that the partitioning is good,
  // well element per well element, starting from the bottom of the well
  for( globalIndex iwelemGlobal = lineBlock.numElements()-1; iwelemGlobal >= 0; --iwelemGlobal )
  {

    // communicate the status of this element
    array1d< integer > thisElemStatusGlobal;
    MpiWrapper::allGather( elemStatusGlobal[iwelemGlobal],
                           thisElemStatusGlobal );
    // group the ranks by well element status
    map< integer, SortedArray< globalIndex > > rankSetsByStatus;
    for( globalIndex irank = 0; irank < thisElemStatusGlobal.size(); ++irank )
    {
      rankSetsByStatus[thisElemStatusGlobal[irank]].insert( irank );
    }
    globalIndex const numLocalRanks = rankSetsByStatus[WellElemStatus::LOCAL].size();

    // in this case, this element has not been assigned
    //    => we assign it to the rank that owns
    //       the well element below iwelemGlobal (prevGlobal, already assigned and checked)
    if( numLocalRanks == 0 )
    {
      globalIndex const numBranches = prevElemIdsGlobal[iwelemGlobal].size();
      globalIndex const prevGlobal  = prevElemIdsGlobal[iwelemGlobal][numBranches-1];

      GEOS_THROW_IF( prevGlobal <= iwelemGlobal || prevGlobal < 0,
                     "The structure of well " << lineBlock.getDataContext() << " is invalid. " <<
                     " The main reason for this error is that there may be no perforation" <<
                     " in the bottom well element of the well, which is required to have" <<
                     " a well-posed problem.",
                     InputError );

      if( elemStatusGlobal[prevGlobal] == WellElemStatus::LOCAL )
      {
        localElems.insert( iwelemGlobal );
        elemStatusGlobal[iwelemGlobal] = WellElemStatus::LOCAL;
      }
      else
      {
        elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
      }
    }
    // in this case, everything is fine,
    // we just update the elemStatusGlobal array for all ranks
    else if( numLocalRanks == 1 )
    {

      for( globalIndex iownerRank : rankSetsByStatus[WellElemStatus::LOCAL] )
      {
        if( MpiWrapper::commRank( MPI_COMM_GEOS ) != iownerRank )
        {
          elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
        }
      }

    }
    // in this case, this element has been assigned to more than rank
    //    => the smallest rank keeps it
    else // (numLocalRanks > 1)
    {

      localIndex rankCount = 0;
      for( globalIndex iownerRank : rankSetsByStatus[WellElemStatus::LOCAL] )
      {
        if( rankCount == 0 )
        {
          // update the elemStatusGlobal array for all ranks
          if( MpiWrapper::commRank( MPI_COMM_GEOS ) != iownerRank )
          {
            elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
          }
        }
        else // (rankCount > 0)
        {
          // remove the duplicate elements
          if( MpiWrapper::commRank( MPI_COMM_GEOS ) == iownerRank )
          {
            localElems.remove( iwelemGlobal );
          }
        }
        rankCount++;
      }

    }

    // TODO: check neighbor rank
  }
}

void WellElementSubRegion::updateNodeManagerSize( MeshLevel & mesh,
                                                  LineBlockABC const & lineBlock,
                                                  SortedArray< globalIndex >      const & localNodes,
                                                  SortedArray< globalIndex >      const & boundaryNodes,
                                                  globalIndex nodeOffsetGlobal )
{

  // get the node manager to compute the total number of mesh nodes
  NodeManager & nodeManager    = mesh.getNodeManager();
  localIndex const numWellNodesLocal = localNodes.size();
  localIndex const oldNumNodesLocal  = nodeManager.size();

  // resize nodeManager to account for the new well nodes and update the properties
  nodeManager.resize( oldNumNodesLocal + numWellNodesLocal );

  arrayView1d< integer > const & isDomainBoundary = nodeManager.getDomainBoundaryIndicator();

  arrayView1d< globalIndex > const & nodeLocalToGlobal = nodeManager.localToGlobalMap();

  arrayView2d< real64 const > const & nodeCoordsGlobal = lineBlock.getNodeCoords();

  // local *well* index
  localIndex iwellNodeLocal = 0;
  // loop over global *well* indices

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  for( globalIndex iwellNodeGlobal : localNodes )
  {
    // local *nodeManager* index
    localIndex const inodeLocal = oldNumNodesLocal + iwellNodeLocal;

    // update node manager maps and position
    nodeLocalToGlobal[inodeLocal]  = nodeOffsetGlobal + iwellNodeGlobal; // global *nodeManager* index
    LvArray::tensorOps::copy< 3 >( X[inodeLocal], nodeCoordsGlobal[ iwellNodeGlobal ] );

    // mark the boundary nodes for ghosting in DomainPartition::SetupCommunications
    if( boundaryNodes.contains( iwellNodeGlobal ) )
    {
      isDomainBoundary[inodeLocal] = 1;
    }

    iwellNodeLocal++;
  }

  // now with update the relevant node indices in nodeManager.globalToLocalMap
  // this is to avoid a call to nodeManager.ConstructGlobalToLocalMap everytime we add a well
  for( iwellNodeLocal = 0; iwellNodeLocal < numWellNodesLocal; ++iwellNodeLocal )
  {
    // local *nodeManager* index
    localIndex const inodeLocal = oldNumNodesLocal + iwellNodeLocal;
    nodeManager.updateGlobalToLocalMap( inodeLocal );
  }
}

void WellElementSubRegion::constructSubRegionLocalElementMaps( MeshLevel & mesh,
                                                               LineBlockABC const & lineBlock,
                                                               SortedArray< globalIndex > const & localElems,
                                                               globalIndex nodeOffsetGlobal,
                                                               globalIndex elemOffsetGlobal )
{
  // get the well geometry
  arrayView1d< globalIndex const > const & nextElemIdGlobal  = lineBlock.getNextElemIndex();
  arrayView2d< real64 const >      const & elemCoordsGlobal  = lineBlock.getElemCoords();
  arrayView2d< globalIndex const > const & elemToNodesGlobal = lineBlock.getElemToNodesMap();
  arrayView1d< real64 const >      const & elemVolumeGlobal  = lineBlock.getElemVolume();

  NodeManager const & nodeManager = mesh.getNodeManager();

  resize( localElems.size() );

  // create local elem numbering

  // local well elem ordering
  localIndex iwelemLocal = 0;
  // loop over global well elem indices
  for( globalIndex iwelemGlobal : localElems )
  {
    // create a global *elemManager* index
    m_localToGlobalMap[iwelemLocal++] = elemOffsetGlobal + iwelemGlobal;
  }
  constructGlobalToLocalMap();

  // recreate local wellbore tree by connecting locally relevant elems
  for( iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal )
  {
    globalIndex const ielemGlobal      = m_localToGlobalMap[iwelemLocal];     // global index in elemManager ordering
    globalIndex const iwelemGlobal     = ielemGlobal - elemOffsetGlobal;      // global index in well ordering
    globalIndex const iwelemNextGlobal = nextElemIdGlobal[iwelemGlobal];      // global index in well ordering
    globalIndex const ielemNextGlobal  = elemOffsetGlobal + iwelemNextGlobal; // global index in elemManager ordering

    if( iwelemNextGlobal < 0 )
    {
      m_nextWellElementIndexGlobal[iwelemLocal] = -1; // wellhead
      m_nextWellElementIndex[iwelemLocal]       = -1;
      m_topWellElementIndex = iwelemLocal;
    }
    else
    {
      m_nextWellElementIndexGlobal[iwelemLocal] = ielemNextGlobal; // wellhead

      if( globalToLocalMap().count( ielemNextGlobal ) > 0 )
      {
        m_nextWellElementIndex[iwelemLocal] = globalToLocalMap( ielemNextGlobal );
      }
      else
      {
        m_nextWellElementIndex[iwelemLocal] = -2; // remote elem
      }
    }

    LvArray::tensorOps::copy< 3 >( m_elementCenter[ iwelemLocal ], elemCoordsGlobal[ iwelemGlobal ] );

    m_elementVolume[iwelemLocal] = elemVolumeGlobal[iwelemGlobal];
    m_radius[iwelemLocal] = lineBlock.getElementRadius();

    // update local well elem to node map (note: nodes are in nodeManager ordering)

    // first get the global node indices in nodeManager ordering
    globalIndex const inodeTopGlobal    = nodeOffsetGlobal + elemToNodesGlobal[iwelemGlobal][LineBlockABC::NodeLocation::TOP];
    globalIndex const inodeBottomGlobal = nodeOffsetGlobal + elemToNodesGlobal[iwelemGlobal][LineBlockABC::NodeLocation::BOTTOM];

    // then get the local node indices in nodeManager ordering
    localIndex const inodeTopLocal    = nodeManager.globalToLocalMap( inodeTopGlobal );
    localIndex const inodeBottomLocal = nodeManager.globalToLocalMap( inodeBottomGlobal );

    m_toNodesRelation[iwelemLocal][LineBlockABC::NodeLocation::TOP]    = inodeTopLocal;
    m_toNodesRelation[iwelemLocal][LineBlockABC::NodeLocation::BOTTOM] = inodeBottomLocal;
  }

}

void WellElementSubRegion::updateNodeManagerNodeToElementMap( MeshLevel & mesh )
{
  ElementRegionManager const & elemManager = mesh.getElemManager();
  NodeManager & nodeManager = mesh.getNodeManager();

  // at this point, NodeManager::SetElementMaps has already been called for the mesh nodes
  // we have to update the following maps for the well nodes
  ArrayOfArrays< localIndex > & toElementRegionList    = nodeManager.elementRegionList();
  ArrayOfArrays< localIndex > & toElementSubRegionList = nodeManager.elementSubRegionList();
  ArrayOfArrays< localIndex > & toElementList          = nodeManager.elementList();

  // we get the region and subregion indices in the elemManager
  WellElementRegion const & elemRegion = dynamicCast< WellElementRegion & >( this->getParent().getParent() );
  string const & elemRegionName = elemRegion.getName();

  localIndex const iregion    = elemManager.getRegions().getIndex( elemRegionName );
  localIndex const isubRegion = elemRegion.getSubRegions().getSubGroupIndex( getName() );

  // for each (new) well element
  for( localIndex iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal )
  {
    for( localIndex a=0; a < numNodesPerElement(); ++a )
    {
      // get the local node index (in nodeManager ordering) using the elem-to-nodes maps constructed above
      localIndex const inodeLocal = m_toNodesRelation[iwelemLocal][a];

      // update the reverse map from well node to well element
      // this is needed to generate the adjacency list in communication setup phase
      toElementRegionList.emplaceBack( inodeLocal, iregion );
      toElementSubRegionList.emplaceBack( inodeLocal, isubRegion );
      toElementList.emplaceBack( inodeLocal, iwelemLocal );
    }
  }

  setupRelatedObjectsInRelations( mesh );
}

void WellElementSubRegion::connectPerforationsToMeshElements( MeshLevel & mesh,
                                                              LineBlockABC const & lineBlock )
{
  arrayView2d< real64 const > const perfCoordsGlobal = lineBlock.getPerfCoords();
  arrayView1d< real64 const > const perfWellTransmissibilityGlobal = lineBlock.getPerfTransmissibility();
  arrayView1d< real64 const > const perfWellSkinFactorGlobal = lineBlock.getPerfSkinFactor();
  arrayView1d< string const > const perfTargetTegionGlobal = lineBlock.getPerfTargetRegion();

  m_perforationData.resize( perfCoordsGlobal.size( 0 ) );
  localIndex iperfLocal = 0;

  arrayView2d< real64 > const perfLocation = m_perforationData.getLocation();

  ElementRegionManager const & elemManager = mesh.getElemManager();

  // loop over all the perforations
  for( globalIndex iperfGlobal = 0; iperfGlobal < perfCoordsGlobal.size( 0 ); ++iperfGlobal )
  {
    real64 const location[3] = { perfCoordsGlobal[iperfGlobal][0],
                                 perfCoordsGlobal[iperfGlobal][1],
                                 perfCoordsGlobal[iperfGlobal][2] };
    GEOS_LOG( GEOS_FMT( "{}: perforation {} location = ({}, {}, {})", lineBlock.getName(), iperfGlobal,
                        location[0], location[1], location[2] ) );

    localIndex erStart = -1, erEnd = -1;

    localIndex const targetRegionIndex = elemManager.getRegions().getIndex( perfTargetTegionGlobal[iperfGlobal] );
    if( targetRegionIndex >= 0 )
    {
      erStart = targetRegionIndex;
      erEnd = erStart + 1;
    }
    else // default is all regions
    {
      erStart = 0;
      erEnd = elemManager.numRegions();
    }

    // for each perforation, we have to find the reservoir element that contains the perforation
    for( localIndex er = erStart; er < erEnd; er++ )
    {
      // search for the reservoir element that contains the well element
      localIndex esrMatched = -1;
      localIndex eiMatched  = -1;
      globalIndex giMatched  = -1;
      bool resElemFound = searchLocalElements( mesh, location, m_searchDepth, er, esrMatched, eiMatched, giMatched );

      // if the element was found
      if( resElemFound )
      {
        // set the indices for the matched reservoir element
        m_perforationData.getMeshElements().m_toElementRegion[iperfLocal] = er;
        m_perforationData.getMeshElements().m_toElementSubRegion[iperfLocal] = esrMatched;
        m_perforationData.getMeshElements().m_toElementIndex[iperfLocal] = eiMatched;
        m_perforationData.getReservoirElementGlobalIndex()[iperfLocal] = giMatched;

        // construct the local wellTransmissibility and location maps
        m_perforationData.getWellTransmissibility()[iperfLocal] = perfWellTransmissibilityGlobal[iperfGlobal];
        m_perforationData.getWellSkinFactor()[iperfLocal] = perfWellSkinFactorGlobal[iperfGlobal];
        LvArray::tensorOps::copy< 3 >( perfLocation[iperfLocal], location );

        // increment the local to global map
        m_perforationData.localToGlobalMap()[iperfLocal++] = iperfGlobal;

        break;
      }
    }
  }

  // set the size based on the number of perforations matched with local reservoir elements
  m_perforationData.resize( iperfLocal );
  m_perforationData.constructGlobalToLocalMap();
}

void WellElementSubRegion::reconstructLocalConnectivity()
{
  // here we reconstruct the array m_nextWellElementIndexGlobal
  // this is needed after the addition of ghost well elements

  for( localIndex iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal )
  {
    globalIndex const nextGlobal = m_nextWellElementIndexGlobal[iwelemLocal];

    if( nextGlobal < 0 )  // well head
    {
      m_nextWellElementIndex[iwelemLocal] = -1;
      m_topWellElementIndex = iwelemLocal; // reset this is case top element was added as ghost
    }
    else if( globalToLocalMap().count( nextGlobal ) == 0 )  // next is remote
    {
      m_nextWellElementIndex[iwelemLocal] = -2;
    }
    else // local
    {
      m_nextWellElementIndex[iwelemLocal] = this->globalToLocalMap( nextGlobal );
    }
  }
}


bool WellElementSubRegion::isLocallyOwned() const
{
  return m_topRank == MpiWrapper::commRank( MPI_COMM_GEOS );
}


localIndex WellElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex WellElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex WellElementSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > const nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
  return bufferOps::Pack< DO_PACKING >( buffer,
                                        nodeList().base().toViewConst(),
                                        m_unmappedGlobalIndicesInNodelist,
                                        packList,
                                        localToGlobal.toSliceConst(),
                                        nodeLocalToGlobal );
}

localIndex WellElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const GEOS_UNUSED_PARAM( overwriteUpMaps ),
                                                   bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  return bufferOps::Unpack( buffer,
                            nodeList().base().toView(),
                            packList,
                            m_unmappedGlobalIndicesInNodelist,
                            this->globalToLocalMap(),
                            nodeList().relatedObjectGlobalToLocal() );
}

void WellElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );
}

}
