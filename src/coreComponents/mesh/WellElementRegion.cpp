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
 * @file WellElementRegion.cpp
 */

#include "WellElementRegion.hpp"

#include "common/MpiWrapper.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/generators/InternalWellGenerator.hpp"

namespace geos
{
using namespace dataRepository;

WellElementRegion::WellElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent ),
  m_subRegionName( name+"UniqueSubRegion" ),
  m_wellControlsName( "" ),
  m_wellGeneratorName( "" )
{
  registerWrapper( viewKeyStruct::wellControlsString(), &m_wellControlsName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef );

  registerWrapper( viewKeyStruct::wellGeneratorString(), &m_wellGeneratorName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef );

  this->getGroup( viewKeyStruct::elementSubRegions() ).registerGroup< WellElementSubRegion >( m_subRegionName );

}

WellElementRegion::~WellElementRegion()
{}

void WellElementRegion::generateWell( MeshLevel & mesh,
                                      LineBlockABC const & lineBlock,
                                      globalIndex nodeOffsetGlobal,
                                      globalIndex elemOffsetGlobal )
{
  // get the (unique) subregion
  WellElementSubRegion &
  subRegion = this->getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                .getGroup< WellElementSubRegion >( m_subRegionName );

  subRegion.setWellControlsName( m_wellControlsName );

  PerforationData * const perforationData = subRegion.getPerforationData();
  perforationData->setNumPerforationsGlobal( lineBlock.numPerforations() );

  globalIndex const numElemsGlobal        = lineBlock.numElements();
  globalIndex const numPerforationsGlobal = lineBlock.numPerforations();

  // 1) select the local perforations based on connectivity to the local reservoir elements
  subRegion.connectPerforationsToMeshElements( mesh, lineBlock );

  globalIndex const matchedPerforations = MpiWrapper::sum( perforationData->size() );
  GEOS_THROW_IF( matchedPerforations != numPerforationsGlobal,
                 "Invalid mapping perforation-to-element in "<<
                 InternalWellGenerator::catalogName() << " " << getWellGeneratorName() << "." <<
                 " This happens when GEOSX cannot match a perforation with a reservoir element." <<
                 " There are two common reasons for this error:\n" <<
                 " 1- The most common reason for this error is that a perforation is on a section of " <<
                 " the well polyline located outside the domain.\n" <<
                 " 2- This error can also happen if a perforation falls on a mesh face or a mesh vertex." <<
                 " Please try to move the perforation slightly (to the interior of the perforated cell) to see if it fixes the problem.",
                 InputError );

  // 2) classify well elements based on connectivity to local mesh partition
  array1d< integer > elemStatusGlobal;
  elemStatusGlobal.resizeDefault( numElemsGlobal, WellElementSubRegion::WellElemStatus::UNOWNED );

  arrayView1d< globalIndex const > const & perfElemIdGlobal = lineBlock.getPerfElemIndex();

  for( localIndex iperfGlobal = 0; iperfGlobal < numPerforationsGlobal; ++iperfGlobal )
  {
    globalIndex const iwelemGlobal = perfElemIdGlobal[iperfGlobal];

    if( perforationData->globalToLocalMap().count( iperfGlobal ) > 0 )
    {
      elemStatusGlobal[iwelemGlobal] |= WellElementSubRegion::WellElemStatus::LOCAL;
    }
    else
    {
      elemStatusGlobal[iwelemGlobal] |= WellElementSubRegion::WellElemStatus::REMOTE;
    }
  }


  // 3) select the local well elements and mark boundary nodes (for ghosting)
  subRegion.generate( mesh,
                      lineBlock,
                      elemStatusGlobal,
                      nodeOffsetGlobal,
                      elemOffsetGlobal );


  // 4) find out which rank is the owner of the top segment
  localIndex const topElemIdLocal = subRegion.getTopWellElementIndex();

  array1d< localIndex > allRankTopElem;
  MpiWrapper::allGather( topElemIdLocal, allRankTopElem );
  int topRank = -1;
  for( int irank = 0; irank < allRankTopElem.size(); ++irank )
  {
    if( allRankTopElem[irank] >= 0 )
    {
      GEOS_ASSERT( topRank < 0 );
      topRank = irank;
    }
  }
  GEOS_ASSERT( topRank >= 0 );
  subRegion.setTopRank( topRank );


  // 5) construct the local perforation to well element map
  perforationData->connectToWellElements( lineBlock,
                                          subRegion.globalToLocalMap(),
                                          elemOffsetGlobal );

}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, WellElementRegion, string const &, Group * const )

} /* namespace geos */
