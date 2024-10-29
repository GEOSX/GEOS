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

#include "WellGeneratorBase.hpp"
#include "LineBlockABC.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/Perforation.hpp"

#include "LvArray/src/genericTensorOps.hpp"
#include "fileIO/Table/TableLayout.hpp"
#include "fileIO/Table/TableData.hpp"
#include "fileIO/Table/TableFormatter.hpp"
#include "common/Format.hpp"
namespace geos
{
using namespace dataRepository;

WellGeneratorBase::WellGeneratorBase( string const & name, Group * const parent ):
  WellGeneratorABC( name, parent )
  , m_numPerforations( 0 )
  , m_numElemsPerSegment( 0 )
  , m_minSegmentLength( 1e-2 )
  , m_minElemLength( 1e-3 )
  , m_radius( 0 )
  , m_wellRegionName( "" )
  , m_wellControlsName( "" )
  , m_numElems( 0 )
  , m_numNodesPerElem( 2 )
  , m_numNodes( 0 )
  , m_nDims( 3 )
  , m_polylineHeadNodeId( -1 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::radiusString(), &m_radius ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Radius of the well [m]" );

  registerWrapper( viewKeyStruct::numElementsPerSegmentString(), &m_numElemsPerSegment ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of well elements per polyline segment" );

  registerWrapper( viewKeyStruct::minSegmentLengthString(), &m_minSegmentLength ).
    setApplyDefaultValue( 1e-2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Minimum length of a well segment [m]" );

  registerWrapper( viewKeyStruct::minElementLengthString(), &m_minElemLength ).
    setApplyDefaultValue( 1e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Minimum length of a well element, computed as (segment length / number of elements per segment ) [m]" );

  registerWrapper( viewKeyStruct::wellRegionNameString(), &m_wellRegionName ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Name of the well element region" );

  registerWrapper( viewKeyStruct::wellControlsNameString(), &m_wellControlsName ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Name of the set of constraints associated with this well" );
}

Group * WellGeneratorBase::createChild( string const & childKey, string const & childName )
{
  if( childKey == viewKeyStruct::perforationString() )
  {
    ++m_numPerforations;

    // keep track of the perforations that have been added
    m_perforationList.emplace_back( childName );
    GEOS_LOG_RANK_0( "Adding Well attribute: " << childKey << ", " << childName );
    return &registerGroup< Perforation >( childName );
  }
  else
  {
    GEOS_THROW( "Unrecognized node: " << childKey, InputError );
  }
  return nullptr;
}

void WellGeneratorBase::expandObjectCatalogs()
{
  createChild( viewKeyStruct::perforationString(), viewKeyStruct::perforationString() );
}

WellGeneratorBase::CatalogInterface::CatalogType & WellGeneratorBase::getCatalog()
{
  static WellGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void WellGeneratorBase::generateWellGeometry( )
{
  fillPolylineDataStructure();

  // count the number of well elements to create
  m_numElems = m_numElemsPerSegment * m_segmentToPolyNodeMap.size( 0 );
  m_numNodes = m_numElems + 1;

  // resize the well element, node, and perforation arrays
  m_elemCenterCoords.resize( m_numElems, 3 );
  m_nextElemId.resize( m_numElems );
  m_prevElemId.resize( m_numElems );
  m_elemToNodesMap.resizeDimension< 0 >( m_numElems );
  m_elemToNodesMap.resizeDimension< 1 >( m_numNodesPerElem );
  m_elemVolume.resize( m_numElems );

  m_nodeDistFromHead.resize( m_numNodes );
  m_nodeCoords.resize( m_numNodes, 3 );

  m_perfCoords.resize( m_numPerforations, 3 );
  m_perfDistFromHead.resize( m_numPerforations );
  m_perfTransmissibility.resize( m_numPerforations );
  m_perfSkinFactor.resize( m_numPerforations );
  m_perfElemId.resize( m_numPerforations );

  // construct a reverse map from the polyline nodes to the segments
  constructPolylineNodeToSegmentMap();

  // detect the head polyline node based on depth
  findPolylineHeadNodeIndex();

  // compute the location and distance from top of the well elements
  discretizePolyline();

  // map the perforations to the well elements
  connectPerforationsToWellElements();

  // make sure that the perforation locations are valid
  checkPerforationLocationsValidity();

  if( getLogLevel() >= 1 && MpiWrapper::commRank() == 0 )
  {
    logInternalWell();
    logPerforationTable();
  }

}

void WellGeneratorBase::postInputInitialization()
{
  GEOS_THROW_IF( m_radius <= 0,
                 "Invalid " << viewKeyStruct::radiusString() << " in well " << getName(),
                 InputError );

  GEOS_THROW_IF( m_wellRegionName.empty(),
                 "Invalid well region name in well " << getName(),
                 InputError );

  GEOS_THROW_IF( m_wellControlsName.empty(),
                 "Invalid well constraint name in well " << getName(),
                 InputError );
}

void WellGeneratorBase::constructPolylineNodeToSegmentMap()
{
  m_polyNodeToSegmentMap.resize( m_polyNodeCoords.size( 0 ) );

  bool foundSmallElem = false;

  // loop over the segments
  for( globalIndex iseg = 0; iseg < m_segmentToPolyNodeMap.size( 0 ); ++iseg )
  {
    globalIndex const ipolyNode_a = m_segmentToPolyNodeMap[iseg][0];
    globalIndex const ipolyNode_b = m_segmentToPolyNodeMap[iseg][1];

    real64 vSeg[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_polyNodeCoords[ipolyNode_a] );
    LvArray::tensorOps::subtract< 3 >( vSeg, m_polyNodeCoords[ipolyNode_b] );

    real64 const segmentLength = LvArray::tensorOps::l2Norm< 3 >( vSeg );

    // various checks and warnings on the segment and element length
    GEOS_THROW_IF( segmentLength < m_minSegmentLength,
                   "Error in the topology of well '" << getName() <<
                   "': we detected a polyline segment measuring less than " << m_minSegmentLength << "m. \n" <<
                   "You can change the minimum segment length using the field " << viewKeyStruct::minSegmentLengthString(),
                   InputError );

    GEOS_THROW_IF( m_polyNodeCoords[ipolyNode_a][2] < m_polyNodeCoords[ipolyNode_b][2],
                   "Error in the topology of well '" << getName() <<
                   "': in the polyline, each segment must be going down. \n" <<
                   "This is not the case between polyline nodes " << m_polyNodeCoords[ipolyNode_a] << " and " << m_polyNodeCoords[ipolyNode_b],
                   InputError );

    if( segmentLength / m_numElemsPerSegment < m_minElemLength )
    {
      foundSmallElem = true;
    }

    // map the polyline node ids to the polyline segment ids
    m_polyNodeToSegmentMap[ipolyNode_a].insert( iseg );
    m_polyNodeToSegmentMap[ipolyNode_b].insert( iseg );
  }

  if( foundSmallElem )
  {
    GEOS_LOG_RANK_0( "\nWarning: the chosen number of well elements per polyline segment (" << m_numElemsPerSegment <<
                     ") leads to well elements measuring less than " << m_minElemLength << "m in the topology of well '" << getName() << "'.\n" <<
                     "The simulation can proceed like that, but small well elements may cause numerical issues, so it is something to keep an eye on.\n" <<
                     "You can get rid of this message by changing the field " << viewKeyStruct::minElementLengthString() );
  }
}

void WellGeneratorBase::findPolylineHeadNodeIndex()
{
  // we assume that the first *polyline segment* in the XML input is the well head segment
  globalIndex const wellHeadSegId = 0;

  // get the corresponding node indices
  globalIndex const ipolyNode_a = m_segmentToPolyNodeMap[wellHeadSegId][0];
  globalIndex const ipolyNode_b = m_segmentToPolyNodeMap[wellHeadSegId][1];

  // we determine which node is the head node based on depth
  // therefore here we throw an error if the well head segment is horizontal
  GEOS_THROW_IF( !(m_polyNodeCoords[ipolyNode_a][2] < m_polyNodeCoords[ipolyNode_b][2])
                 && !(m_polyNodeCoords[ipolyNode_a][2] > m_polyNodeCoords[ipolyNode_b][2]),
                 "The head polyline segment cannot be horizontal in well '" << getName()
                                                                            << "' since we use depth to determine which of its nodes is to head node of the well.\n"
                                                                            << "If you are trying to set up a horizontal well, please simply add a non-horizontal segment at the top of the well,"
                                                                            << " and this error will go away",
                 InputError );

  // detect the top node, assuming z oriented upwards
  m_polylineHeadNodeId =
    ( m_polyNodeCoords[ipolyNode_a][2] > m_polyNodeCoords[ipolyNode_b][2] )
    ? ipolyNode_a
    : ipolyNode_b;

  real64 const headZcoord = m_polyNodeCoords[m_polylineHeadNodeId][2];
  for( globalIndex inode = 0; inode < m_polyNodeCoords.size( 0 ); ++inode )
  {
    if( inode == m_polylineHeadNodeId )
    {
      continue;
    }
    real64 const currentZcoord = m_polyNodeCoords[inode][2];

    GEOS_THROW_IF( !(currentZcoord < headZcoord),
                   "Error in the topology of well '" << getName()
                                                     << "' since we found a well node that is above the head node",
                   InputError );
  }
}

void WellGeneratorBase::discretizePolyline()
{
  // initialize well elements and node ids
  globalIndex ipolyNodeTop  = m_polylineHeadNodeId;
  globalIndex iwelemCurrent = 0;
  globalIndex isegCurrent   = 0;

  // set the location of the first well node and distance from well head
  LvArray::tensorOps::copy< 3 >( m_nodeCoords[0], m_polyNodeCoords[ipolyNodeTop] );
  m_nodeDistFromHead[0] = 0.0;

  // note: this part of the code does not support well branching
  // TODO: check that there is only one branch
  // TODO: read wells with branching (already supported elsewhere in the code)

  // go through the well from top to bottom
  for( globalIndex is = 0; is < m_segmentToPolyNodeMap.size( 0 ); ++is )
  {

    GEOS_THROW_IF( isegCurrent == -1,
                   "Invalid segmentToNode map in well " << getName(),
                   InputError );

    globalIndex const ipolyNodeBottom = ( ipolyNodeTop == m_segmentToPolyNodeMap[isegCurrent][0] )
                                      ? m_segmentToPolyNodeMap[isegCurrent][1]
                                      : m_segmentToPolyNodeMap[isegCurrent][0];

    real64 vPoly[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_polyNodeCoords[ipolyNodeBottom] );
    LvArray::tensorOps::subtract< 3 >( vPoly, m_polyNodeCoords[ipolyNodeTop] );

    // add the well elements and well nodes corresponding to this polyline segment
    for( localIndex iw = 0; iw < m_numElemsPerSegment; ++iw )
    {

      // 1) set the element location, connectivity
      real64 const scaleCenter = (iw + 0.5) / static_cast< real64 >(m_numElemsPerSegment);
      LvArray::tensorOps::copy< 3 >( m_elemCenterCoords[iwelemCurrent], vPoly );
      LvArray::tensorOps::scale< 3 >( m_elemCenterCoords[iwelemCurrent], scaleCenter );
      LvArray::tensorOps::add< 3 >( m_elemCenterCoords[iwelemCurrent], m_polyNodeCoords[ipolyNodeTop] );

      GEOS_THROW_IF( iwelemCurrent >= m_numElems,
                     "Invalid well topology in well " << getName(),
                     InputError );

      globalIndex const iwelemTop    = iwelemCurrent - 1;
      globalIndex const iwelemBottom = iwelemCurrent + 1;
      m_nextElemId[iwelemCurrent] = iwelemTop;
      m_prevElemId[iwelemCurrent].resize( 1 );
      m_prevElemId[iwelemCurrent][0] = (iwelemBottom < m_numElems)
                                     ? iwelemBottom
                                     : -1;

      // 2) set the node properties
      globalIndex const iwellNodeTop    = iwelemCurrent;
      globalIndex const iwellNodeBottom = iwelemCurrent+1;
      m_elemToNodesMap[iwelemCurrent][LineBlockABC::NodeLocation::TOP]    = iwellNodeTop;
      m_elemToNodesMap[iwelemCurrent][LineBlockABC::NodeLocation::BOTTOM] = iwellNodeBottom;

      real64 const scaleBottom = (iw + 1.0) / m_numElemsPerSegment;
      LvArray::tensorOps::copy< 3 >( m_nodeCoords[iwellNodeBottom], vPoly );
      LvArray::tensorOps::scale< 3 >( m_nodeCoords[iwellNodeBottom], scaleBottom );
      LvArray::tensorOps::add< 3 >( m_nodeCoords[iwellNodeBottom], m_polyNodeCoords[ipolyNodeTop] );

      // 3) increment the distance from the well head to the bottom of current element
      real64 vWellElem[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_nodeCoords[iwellNodeBottom] );
      LvArray::tensorOps::subtract< 3 >( vWellElem, m_nodeCoords[iwellNodeTop] );
      m_nodeDistFromHead[iwellNodeBottom]  = LvArray::tensorOps::l2Norm< 3 >( vWellElem );
      m_nodeDistFromHead[iwellNodeBottom] += m_nodeDistFromHead[iwellNodeTop];

      // 4) set element volume
      m_elemVolume[iwelemCurrent] = LvArray::tensorOps::l2Norm< 3 >( vWellElem ) * M_PI * m_radius * m_radius;

      // 4) increment the element counter
      ++iwelemCurrent;
    }

    // then consider the next polyline segment
    ipolyNodeTop = ipolyNodeBottom;
    isegCurrent  = getNextSegmentIndex( isegCurrent, ipolyNodeTop );
  }

  // set the previous index for the bottom segment
  m_prevElemId[iwelemCurrent-1].resize( 1 );
  m_prevElemId[iwelemCurrent-1][0] = -1;

}


void WellGeneratorBase::connectPerforationsToWellElements()
{

  // assign a well element to each perforation
  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {

    // get the perforation and its properties
    Perforation const & perf = this->getGroup< Perforation >( m_perforationList[iperf] );
    m_perfDistFromHead[iperf]  = perf.getDistanceFromWellHead();
    m_perfTransmissibility[iperf] = perf.getWellTransmissibility();
    m_perfSkinFactor[iperf] = perf.getWellSkinFactor();

    // search in all the elements of this well between head and bottom
    globalIndex iwelemTop    = 0;
    globalIndex iwelemBottom = m_numElems - 1;

    // check the validity of the perforation before starting
    real64 const wellLength = m_nodeDistFromHead[m_elemToNodesMap[iwelemBottom][LineBlockABC::NodeLocation::BOTTOM]];

    GEOS_THROW_IF( m_perfDistFromHead[iperf] > wellLength,
                   "Distance from perforation " << perf.getName() << " to head is larger than well polyline length for well " << getName() << "\n \n"
                                                << "Here is how the \"distanceFromHead\" keyword is used in the definition of the perforation location: \n"
                                                << "We start from the well head (top of the well) and we measure the linear distance along the well polyline as we go down the well.\n"
                                                << "When we reach the distanceFromHead specified by the user, we place a perforation on the well at this location of the polyline, and connect it to the reservoir element that contains this perforation",
                   InputError );

    // start binary search
    const globalIndex maxNumSteps = m_numElems + 1;
    globalIndex currentNumSteps = 0;
    while( iwelemTop < iwelemBottom )
    {
      globalIndex iwelemMid =
        static_cast< globalIndex >(floor( static_cast< real64 >(iwelemTop + iwelemBottom) / 2.0 ));
      real64 const headToBottomDist =
        m_nodeDistFromHead[m_elemToNodesMap[iwelemMid][LineBlockABC::NodeLocation::BOTTOM]];

      if( headToBottomDist < m_perfDistFromHead[iperf] )
      {
        iwelemTop = iwelemMid + 1;
      }
      else
      {
        iwelemBottom = iwelemMid;
      }

      GEOS_THROW_IF( currentNumSteps > maxNumSteps,
                     "Perforation " << perf.getName() << " cannot be mapped to a well element in well " << getName(),
                     InputError );

      currentNumSteps++;
    }

    // set the index of the matched element
    globalIndex iwelemMatched = iwelemTop;
    GEOS_THROW_IF( iwelemMatched >= m_numElems,
                   "Invalid topology in well " << getName(),
                   InputError );

    m_perfElemId[iperf] = iwelemMatched;

    // compute the physical location of the perforation
    globalIndex const inodeTop    = m_elemToNodesMap[iwelemMatched][LineBlockABC::NodeLocation::TOP];
    globalIndex const inodeBottom = m_elemToNodesMap[iwelemMatched][LineBlockABC::NodeLocation::BOTTOM];
    real64 const elemLength       = m_nodeDistFromHead[inodeBottom] - m_nodeDistFromHead[inodeTop];
    real64 const topToPerfDist    = m_perfDistFromHead[iperf] - m_nodeDistFromHead[inodeTop];

    GEOS_THROW_IF( (elemLength <= 0) || (topToPerfDist < 0),
                   "Invalid topology in well " << getName(),
                   InputError );

    LvArray::tensorOps::copy< 3 >( m_perfCoords[iperf], m_nodeCoords[inodeBottom] );
    LvArray::tensorOps::subtract< 3 >( m_perfCoords[iperf], m_nodeCoords[inodeTop] );
    LvArray::tensorOps::scale< 3 >( m_perfCoords[iperf], topToPerfDist / elemLength );
    LvArray::tensorOps::add< 3 >( m_perfCoords[iperf], m_nodeCoords[inodeTop] );

  }
}


globalIndex WellGeneratorBase::getNextSegmentIndex( globalIndex topSegId,
                                                    globalIndex currentPolyNodeId ) const
{
  globalIndex nextSegId = -1;

  // get the index of the two segments sharing this node
  for( auto is : m_polyNodeToSegmentMap[currentPolyNodeId] )
  {
    // if this is not equal to the current segId, we found the next segId
    if( is != topSegId )
    {
      nextSegId = is;
      break;
    }
  }

  return nextSegId;
}

void WellGeneratorBase::checkPerforationLocationsValidity()
{
  array1d< array1d< localIndex > > elemToPerfMap;
  elemToPerfMap.resize( m_numElems );

  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {
    elemToPerfMap[m_perfElemId[iperf]].emplace_back( iperf );
  }

  // merge perforations to make sure that no well element is shared between two MPI domains
  // TODO: instead of merging perforations, split the well elements and do not change the physical location of the
  // perforation
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  if( mpiSize > 1 )
  {
    mergePerforations( elemToPerfMap );
  }

  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    // check that there is always a perforation in the last well element (otherwise, the problem is not well posed)
    for( localIndex iwelemPrev = 0; iwelemPrev < m_prevElemId[iwelem].size(); ++iwelemPrev )
    {
      GEOS_THROW_IF( m_prevElemId[iwelem][iwelemPrev] == -1 && elemToPerfMap[iwelem].size() == 0,
                     "The bottom element of well " << getName() << " does not have a perforation. "
                                                   << "This is needed to have a well-posed problem. \n\n"
                                                   << "Here are the two possible ways to solve this problem: \n\n"
                                                   << "1) Adding a perforation located close to the bottom of the well. "
                                                   << "To do that, compute the total length of the well polyline (by summing the length of the well segments defined by the keywords \"polylineNodeCoords\" and \"polylineSegmentConn\") "
                                                   << "and place a perforation whose \"distanceFromHead\" is slightly smaller than this total length. \n \n"
                                                   << "2) Shorten  the well polyline. "
                                                   << "To do that, reduce the length of the well polyline by shortening the segments defined by the keywords \"polylineNodeCoords\" and \"polylineSegmentConn\", or by removing a segment.",
                     InputError );
    }
  }
}

void WellGeneratorBase::mergePerforations( array1d< array1d< localIndex > > const & elemToPerfMap )
{

  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    // collect the indices of the elems with more that one perforation
    if( elemToPerfMap[iwelem].size() > 1 )
    {
      // find the perforation with the largest Peaceman index and keep its location
      globalIndex iperfMaxTransmissibility = elemToPerfMap[iwelem][0];
      real64 maxTransmissibility = m_perfTransmissibility[iperfMaxTransmissibility];
      for( localIndex ip = 1; ip < elemToPerfMap[iwelem].size(); ++ip )
      {
        if( m_perfTransmissibility[elemToPerfMap[iwelem][ip]] > maxTransmissibility )
        {
          iperfMaxTransmissibility = elemToPerfMap[iwelem][ip];
          maxTransmissibility = m_perfTransmissibility[iperfMaxTransmissibility];
        }
      }

      // assign the coordinates of the perf with the largest trans to the other perfs on this elem
      for( localIndex ip = 0; ip < elemToPerfMap[iwelem].size(); ++ip )
      {
        if( elemToPerfMap[iwelem][ip] == iperfMaxTransmissibility )
        {
          continue;
        }

        GEOS_LOG_RANK_0( "\n \nThe GEOSX wells currently have the following limitation in parallel: \n"
                         << "We cannot allow an element of the well mesh to have two or more perforations associated with it. \n"
                         << "So, in the present simulation, perforation #" << elemToPerfMap[iwelem][ip]
                         << " of well " << getName()
                         << " is moved from " << m_perfCoords[elemToPerfMap[iwelem][ip]]
                         << " to " << m_perfCoords[iperfMaxTransmissibility]
                         << " to make sure that no element of the well mesh has two perforations associated with it. \n"
                         << "To circumvent this issue, please increase the value of \"numElementsPerSegment\" that controls the number of (uniformly distributed) well mesh elements per segment of the well polyline. \n"
                         << "Our recommendation is to choose \"numElementsPerSegment\" such that each well mesh element has at most one perforation. \n\n" );
        LvArray::tensorOps::copy< 3 >( m_perfCoords[elemToPerfMap[iwelem][ip]], m_perfCoords[iperfMaxTransmissibility] );
      }
    }
  }
}

void WellGeneratorBase::logInternalWell() const
{
  TableData tableWellData;
  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    std::optional< globalIndex > nextElement;
    std::optional< globalIndex > prevElement;

    if( m_nextElemId[iwelem] >= 0 )
    {
      nextElement =  m_nextElemId[iwelem];
    }

    if( m_prevElemId[iwelem][0] >= 0 )
    {
      prevElement =  m_prevElemId[iwelem][0];
    }

    tableWellData.addRow( iwelem,
                          m_elemCenterCoords[iwelem][0],
                          m_elemCenterCoords[iwelem][1],
                          m_elemCenterCoords[iwelem][2],
                          prevElement,
                          nextElement );
  }

  string const wellTitle = GEOS_FMT( "Well '{}' Element Table", getName() );
  TableLayout const tableWellLayout = TableLayout( {
      TableLayout::ColumnParam{"Element no.", TableLayout::Alignment::right},
      TableLayout::ColumnParam{"CoordX", TableLayout::Alignment::right},
      TableLayout::ColumnParam{"CoordY", TableLayout::Alignment::right},
      TableLayout::ColumnParam{"CoordZ", TableLayout::Alignment::right},
      TableLayout::ColumnParam{"Prev\nElement", TableLayout::Alignment::right},
      TableLayout::ColumnParam{"Next\nElement", TableLayout::Alignment::right},
    }, wellTitle );

  TableTextFormatter const tableFormatter( tableWellLayout );
  GEOS_LOG_RANK_0( tableFormatter.toString( tableWellData ));
}

void WellGeneratorBase::logPerforationTable() const
{
  TableData tablePerfoData;
  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {
    tablePerfoData.addRow( iperf, m_perfCoords[iperf], m_perfElemId[iperf] );
  }

  TableLayout const tableLayoutPerfo ( {"Perforation no.", "Coordinates", "connected to"},
                                       GEOS_FMT( "Well '{}' Perforation Table", getName() ) );
  TableTextFormatter const tablePerfoLog( tableLayoutPerfo );
  GEOS_LOG_RANK_0( tablePerfoLog.toString( tablePerfoData ));
}

}
