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

#include "MeshGeneratorBase.hpp"
#include "mesh/generators/CellBlockManager.hpp"

namespace geos
{
using namespace dataRepository;

MeshGeneratorBase::MeshGeneratorBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

Group * MeshGeneratorBase::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding Mesh attribute: " << childKey << ", " << childName );
  std::unique_ptr< WellGeneratorBase > wellGen = WellGeneratorBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< WellGeneratorBase >( childName, std::move( wellGen ) );
}

void MeshGeneratorBase::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from WellGeneratorBase here
  for( auto & catalogIter: WellGeneratorBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

MeshGeneratorBase::CatalogInterface::CatalogType & MeshGeneratorBase::getCatalog()
{
  static MeshGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void MeshGeneratorBase::generateMesh( Group & parent, array1d< int > const & partition )
{
  CellBlockManager & cellBlockManager = parent.registerGroup< CellBlockManager >( keys::cellManager );
  fillCellBlockManager( cellBlockManager, partition );
  this->attachWellInfo( cellBlockManager );
}

void MeshGeneratorBase::attachWellInfo( CellBlockManager & cellBlockManager )
{
  forSubGroups< WellGeneratorBase >( [&]( WellGeneratorBase & wellGen ) {
    wellGen.generateWellGeometry( );
    LineBlock & lb = cellBlockManager.registerLineBlock( wellGen.getWellRegionName() );
    lb.setNumElements( wellGen.numElements() );
    lb.setElemCoords( wellGen.getElemCoords() );
    lb.setNextElemIndex( wellGen.getNextElemIndex() );
    lb.setPrevElemIndices( wellGen.getPrevElemIndices() );
    lb.setElemToNodesMap( wellGen.getElemToNodesMap() );
    lb.setElemVolume( wellGen.getElemVolume() );
    lb.setElementRadius( wellGen.getElementRadius() );
    lb.setNumNodes( wellGen.numNodes() );
    lb.setNodeCoords( wellGen.getNodeCoords() );
    lb.setNumPerforations( wellGen.numPerforations() );
    lb.setPerfCoords( wellGen.getPerfCoords() );
    lb.setPerfTransmissibility( wellGen.getPerfTransmissibility() );
    lb.setPerfSkinFactor( wellGen.getPerfSkinFactor() );
    lb.setPerfElemIndex( wellGen.getPerfElemIndex() );
    lb.setWellControlsName( wellGen.getWellControlsName() );
    lb.setWellGeneratorName( wellGen.getName() );
  } );
}
}
