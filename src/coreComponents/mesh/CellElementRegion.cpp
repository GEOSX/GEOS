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

#include "CellElementRegion.hpp"
#include "CellElementSubRegion.hpp"
#include "mesh/generators/CellBlockABC.hpp"

#include <fnmatch.h>

namespace geos
{
using namespace dataRepository;

CellElementRegion::CellElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  std::vector< string > elementNames;
  for( int i = 0; i < numElementTypes(); ++i )
  {
    if( getElementDim( (ElementType)i ) == 3 )
    {
      elementNames.push_back( getElementTypeName( (ElementType)i ) );
    }
  }

  registerWrapper( viewKeyStruct::sourceCellBlockNamesString(), &m_cellBlockNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( GEOS_FMT( "List of the desired cell-blocks qualifiers to contain in this {}. Qualifiers can be either cell-block "
                              "names, region attribute values, or fnmatch pattern."
                              "The form of loaded cell-block names is of \"regionAttribute_elementType\", so \"1_tetrahedra\" "
                              " contains the tetrahedric elements for which the regionAttribute is 1.\n"
                              "The element types are: {}.",
                              catalogName(), stringutilities::join( elementNames, ", " ) ) );

  registerWrapper( viewKeyStruct::coarseningRatioString(), &m_coarseningRatio ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 );
}

CellElementRegion::~CellElementRegion()
{}

void CellElementRegion::generateMesh( Group const & cellBlocks )
{
  GEOS_THROW_IF( m_cellBlockNames.empty(),
                 GEOS_FMT( "{}: No cellBlock selected in this region.",
                           getDataContext() ),
                 InputError );
  Group & subRegions = this->getGroup( viewKeyStruct::elementSubRegions() );
  for( string const & cbName : m_cellBlockNames )
  {
    CellBlockABC const * cellBlock = cellBlocks.getGroupPointer< CellBlockABC >( cbName );
    GEOS_THROW_IF( cellBlock == nullptr,
                   GEOS_FMT( "{}: No cellBlock named '{}' found.\nAvailable cellBlock list: {{ {} }}\nNo CellElementRegionSelector has been used to verify the cellBlock selection.",
                             getDataContext(), cbName, stringutilities::join( m_cellBlockNames, ", " ) ),
                   InputError );

    // subRegion name must be the same as the cell-block (so we can match them and reference them in errors).
    CellElementSubRegion & subRegion = subRegions.registerGroup< CellElementSubRegion >( cbName );
    subRegion.copyFromCellBlock( *cellBlock );
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geos */
