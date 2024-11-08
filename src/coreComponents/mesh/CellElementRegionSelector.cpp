/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "mesh/CellElementRegionSelector.hpp"

#include <fnmatch.h>
#include <optional>


namespace geos
{
using namespace dataRepository;
using ViewKeys = CellElementRegion::viewKeyStruct;


CellElementRegionSelector::CellElementRegionSelector(
  Group const & cellBlocks,
  std::map< integer, std::set< string > > const & regionsCellBlocks )
{
  // The owners lists need to be initialized so we will be able to verify later that it is not empty.

  cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
  {
    string const name = cellBlock.getName();
    m_cellBlocksOwners.emplace( name, std::vector< CellElementRegion const * >() );
  } );

  for( auto const & regionCellBlocks : regionsCellBlocks )
  {
    string const regionAttributeStr = std::to_string( regionCellBlocks.first );
    m_regionAttributesCellBlocks.emplace( regionAttributeStr, regionCellBlocks.second );
    m_regionAttributesOwners.emplace( regionAttributeStr, std::vector< CellElementRegion const * >() );
  }
}


std::set< string >
CellElementRegionSelector::getMatchingCellblocks( CellElementRegion const & region,
                                                  string_view matchPattern ) const
{
  std::set< string > matchedCellBlocks;
  bool matching = false;
  for( auto const & [cellBlockName, owners] : m_cellBlocksOwners )
  {
    // if the pattern matches the tested cellBlock name
    if( fnmatch( matchPattern.data(), cellBlockName.c_str(), 0 ) == 0 )
    {
      matching = true;
      matchedCellBlocks.emplace( cellBlockName );
    }
  }

  GEOS_THROW_IF( !matching,
                 GEOS_FMT( "{}: No cellBlock name is satisfying the qualifier '{}'.\n"
                           "Available cellBlock list: {{ {} }}\nAvailable region attribute list: {{ {} }}",
                           region.getWrapperDataContext( ViewKeys::sourceCellBlockNamesString() ),
                           matchPattern,
                           stringutilities::joinLambda( m_regionAttributesOwners, ", ",
                                                        []( auto pair ) { return pair->first; } ),
                           stringutilities::joinLambda( m_cellBlocksOwners, ", ",
                                                        []( auto pair ) { return pair->first; } ) ),
                 InputError );
  return matchedCellBlocks;
}

void
CellElementRegionSelector::verifyRequestedCellBlocks( CellElementRegion const & region,
                                                      std::set< string > const & cellBlockNames ) const
{
  for( string const & requestedCellBlockName : cellBlockNames )
  {
    // if cell block does not exist in the mesh
    GEOS_THROW_IF( m_cellBlocksOwners.count( requestedCellBlockName ) == 0,
                   GEOS_FMT( "{}: No cellBlock named '{}'.\nAvailable cellBlock list: {{ {} }}",
                             region.getWrapperDataContext( ViewKeys::sourceCellBlockNamesString() ),
                             requestedCellBlockName,
                             stringutilities::joinLambda( m_cellBlocksOwners, ", ",
                                                          []( auto pair ) { return pair->first; } ) ),
                   InputError );
  }
}


void
CellElementRegionSelector::registerRegionSelection( CellElementRegion const & region,
                                                    std::set< string > const & cellBlockNames,
                                                    std::set< string > const & attributeValues )
{
  for( string attributeValue : attributeValues )
  {
    m_regionAttributesOwners[attributeValue].push_back( &region );
  }

  for( string const & cellBlockName : cellBlockNames )
  {
    m_cellBlocksOwners[cellBlockName].push_back( &region );
  }
}


std::set< string > CellElementRegionSelector::buildCellBlocksSelection( CellElementRegion const & region )
{
  std::set< string > cellBlocksSelection;
  std::set< string > regionAttributeSelection;

  auto const qualifiers = region.getCellBlockNames();
  for( string const & qualifier : qualifiers )
  {
    auto const regionCellBlocks = m_regionAttributesCellBlocks.find( qualifier );
    if( regionCellBlocks != m_regionAttributesCellBlocks.end() )
    { // if the qualifier is a region attribute value, let's select it
      regionAttributeSelection.emplace( regionCellBlocks->first );
      for( string const & cellBlock : regionCellBlocks->second )
      {
        cellBlocksSelection.emplace( cellBlock );
      }
    }
    else
    { // the qualifier is a match pattern, or a simple cellblock name, let's select all matching cell-blocks
      std::set< string > const matchedCellBlocks = getMatchingCellblocks( region, qualifier );
      cellBlocksSelection.insert( matchedCellBlocks.begin(), matchedCellBlocks.end() );
    }
  }

  verifyRequestedCellBlocks( region, cellBlocksSelection );
  registerRegionSelection( region, cellBlocksSelection, regionAttributeSelection );

  return cellBlocksSelection;
}

void CellElementRegionSelector::checkSelectionConsistency() const
{
  auto const getRegionStr = []( auto regionPtrIterator ) -> string {
    return GEOS_FMT( "- {}", (*regionPtrIterator)->getDataContext() );
  };

  auto const checkOwnerCount = [&]( string_view qualifierType,
                                    auto const & qualifiersOwners,
                                    auto & orphanList ) {
    // Search of never or multiple selected attribute values
    std::vector< string > multipleRefsErrors;
    for( auto const & [qualifier, owningRegions] : qualifiersOwners )
    {
      if( owningRegions.size() == 0 )
      {
        orphanList.insert( qualifier );
      }
      else if( owningRegions.size() > 1 )
      {
        multipleRefsErrors.push_back(
          GEOS_FMT( "The {} '{}' has been referenced in multiple {}:\n{}",
                    qualifierType, qualifier, CellElementRegion::catalogName(),
                    stringutilities::joinLambda( owningRegions, '\n', getRegionStr ) ) );
      }
    }
    GEOS_THROW_IF( !multipleRefsErrors.empty(), stringutilities::join( multipleRefsErrors, "\n\n" ), InputError );
  };

  std::set< string > orphanRegionAttributes;
  std::set< string > orphanCellBlockNames;
  checkOwnerCount( "region attribute", m_regionAttributesOwners, orphanRegionAttributes );
  checkOwnerCount( "cell-block", m_cellBlocksOwners, orphanCellBlockNames );

  if( !orphanCellBlockNames.empty() )
  {
    std::ostringstream oss;
    if( !orphanRegionAttributes.empty() )
    {
      oss << GEOS_FMT( "The region attributes {{ {} }} has not been referenced in any {}.\n",
                       stringutilities::join( orphanRegionAttributes, ", " ),
                       CellElementRegion::catalogName() );
    }
    oss << GEOS_FMT( "The following cell-blocks has not been referenced in any region: {{ {} }}.\n",
                     stringutilities::join( orphanCellBlockNames, ", " ) );
    oss << GEOS_FMT( "Please add it in an existing {} (through the '{}' attribute), or consider creating a new one to describe your model.",
                     CellElementRegion::catalogName(), ViewKeys::sourceCellBlockNamesString() );
    GEOS_THROW( oss.str(), InputError );
  }
}


} /* namespace geos */
