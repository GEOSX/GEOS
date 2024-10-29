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

#include "CellElementRegion.hpp"
#include "CellElementSubRegion.hpp"
#include "mesh/generators/CellBlockABC.hpp"

namespace geos
{
using namespace dataRepository;

CellElementRegion::CellElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceCellBlockNamesString(), &m_cellBlockNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::coarseningRatioString(), &m_coarseningRatio ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 );
}

CellElementRegion::~CellElementRegion()
{}

void CellElementRegion::generateMesh( Group const & cellBlocks )
{
  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );

  for( string const & cellBlockName : this->m_cellBlockNames )
  {
    CellElementSubRegion & subRegion = elementSubRegions.registerGroup< CellElementSubRegion >( cellBlockName );
    CellBlockABC const & source = cellBlocks.getGroup< CellBlockABC >( subRegion.getName() );
    subRegion.copyFromCellBlock( source );
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geos */
