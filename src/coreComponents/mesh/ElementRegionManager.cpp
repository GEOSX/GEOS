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

#include <map>
#include <vector>

#include "ElementRegionManager.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "SurfaceElementRegion.hpp"
#include "FaceManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"
#include "schema/schemaUtilities.hpp"
#include "mesh/generators/LineBlockABC.hpp"
#include "mesh/CellElementRegionSelector.hpp"

namespace geos
{
using namespace dataRepository;

ElementRegionManager::ElementRegionManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  this->registerGroup< Group >( ElementRegionManager::groupKeyStruct::elementRegionsGroup() );
}

ElementRegionManager::~ElementRegionManager()
{
  // TODO Auto-generated destructor stub
}

void ElementRegionManager::resize( integer_array const & numElements,
                                   string_array const & regionNames,
                                   string_array const & GEOS_UNUSED_PARAM( elementTypes ) )
{
  localIndex const n_regions = LvArray::integerConversion< localIndex >( regionNames.size());
  for( localIndex reg=0; reg<n_regions; ++reg )
  {
    this->getRegion( reg ).resize( numElements[reg] );
  }
}

void ElementRegionManager::setMaxGlobalIndex()
{
  forElementSubRegions< ElementSubRegionBase >( [this] ( ElementSubRegionBase const & subRegion )
  {
    m_localMaxGlobalIndex = std::max( m_localMaxGlobalIndex, subRegion.maxGlobalIndex() );
  } );

  m_maxGlobalIndex = MpiWrapper::max( m_localMaxGlobalIndex, MPI_COMM_GEOS );
}



Group * ElementRegionManager::createChild( string const & childKey, string const & childName )
{
  GEOS_ERROR_IF( !(CatalogInterface::hasKeyName( childKey )),
                 "KeyName ("<<childKey<<") not found in ObjectManager::Catalog" );
  GEOS_LOG_RANK_0( "Adding Object " << childKey<<" named "<< childName<<" from ObjectManager::Catalog." );

  Group & elementRegions = this->getGroup( ElementRegionManager::groupKeyStruct::elementRegionsGroup() );
  return &elementRegions.registerGroup( childName,
                                        CatalogInterface::factory( childKey, childName, &elementRegions ) );
}

void ElementRegionManager::expandObjectCatalogs()
{
  ObjectManagerBase::CatalogInterface::CatalogType const & catalog = ObjectManagerBase::getCatalog();
  for( ObjectManagerBase::CatalogInterface::CatalogType::const_iterator iter = catalog.begin();
       iter!=catalog.end();
       ++iter )
  {
    string const key = iter->first;
    if( key.find( "ElementRegion" ) != string::npos )
    {
      this->createChild( key, key );
    }
  }
}


void ElementRegionManager::setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                                xmlWrapper::xmlNode schemaParent,
                                                integer documentationType )
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child( "xsd:choice" );
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child( "xsd:choice" );
    targetChoiceNode.append_attribute( "minOccurs" ) = "0";
    targetChoiceNode.append_attribute( "maxOccurs" ) = "unbounded";
  }

  std::set< string > names;
  this->forElementRegions( [&]( ElementRegionBase & elementRegion )
  {
    names.insert( elementRegion.getName() );
  } );

  for( string const & name: names )
  {
    schemaUtilities::SchemaConstruction( getRegion( name ), schemaRoot, targetChoiceNode, documentationType );
  }
}


void ElementRegionManager::generateMesh( CellBlockManagerABC const & cellBlockManager )
{
  { // cellBlocks loading
    Group const & cellBlocks = cellBlockManager.getCellBlocks();
    CellElementRegionSelector cellBlockSelector{ cellBlocks,
                                                 cellBlockManager.getRegionAttributesCellBlocks() };
    this->forElementRegions< CellElementRegion >( [&]( CellElementRegion & elemRegion )
    {
      elemRegion.setCellBlockNames( cellBlockSelector.buildCellBlocksSelection( elemRegion ) );
      elemRegion.generateMesh( cellBlocks );
    } );
    // selecting all cellblocks is mandatory
    cellBlockSelector.checkSelectionConsistency();
  }


  this->forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion & elemRegion )
  {
    elemRegion.generateMesh( cellBlockManager.getFaceBlocks() );
  } );

  // Some mappings of the surfaces subregions point to elements in other subregions and regions.
  // For the moment, those mappings only point to cell block indices.
  // The following makes use of cell block to subregions mappings to finalize the surfaces information.
  // It also creates the mappings concerning the regions.
  array2d< localIndex > const blockToSubRegion = this->getCellBlockToSubRegionMap( cellBlockManager );
  this->forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion & elemRegion )
  {
    SurfaceElementSubRegion & subRegion = elemRegion.getUniqueSubRegion< SurfaceElementSubRegion >();
    // While indicated as containing element subregion information,
    // `relation` currently contains cell block information
    // that will be transformed into element subregion information.
    // This is why we copy the information into a temporary,
    // which frees space for the final information (of same size).
    // OrderedVariableToManyElementRelation & relation = subRegion.getToCellRelation();
    // ToCellRelation< ArrayOfArrays< localIndex > > const tmp( relation.m_toElementSubRegion,
    //                                                          relation.m_toElementIndex );
    // meshMapUtilities::transformCellBlockToRegionMap< parallelHostPolicy >( blockToSubRegion.toViewConst(),
    //                                                                        tmp,
    //                                                                        relation );
  } );
}

void ElementRegionManager::generateWells( CellBlockManagerABC const & cellBlockManager,
                                          MeshLevel & meshLevel )
{
  NodeManager & nodeManager = meshLevel.getNodeManager();

  // get the offsets to construct local-to-global maps for well nodes and elements
  nodeManager.setMaxGlobalIndex();
  globalIndex const nodeOffsetGlobal = nodeManager.maxGlobalIndex() + 1;
  localIndex const elemOffsetLocal  = this->getNumberOfElements();
  globalIndex const elemOffsetGlobal = MpiWrapper::sum( elemOffsetLocal );

  globalIndex wellElemCount = 0;
  globalIndex wellNodeCount = 0;

  // construct the wells one by one
  forElementRegions< WellElementRegion >( [&]( WellElementRegion & wellRegion )
  {

    // get the global well geometry from the well generator
    LineBlockABC const & lineBlock = cellBlockManager.getLineBlock( wellRegion.getName() );
    wellRegion.setWellGeneratorName( lineBlock.getWellGeneratorName() );
    wellRegion.setWellControlsName( lineBlock.getWellControlsName() );
    // generate the local data (well elements, nodes, perforations) on this well
    // note: each MPI rank knows the global info on the entire well (constructed earlier in InternalWellGenerator)
    // so we only need node and element offsets to construct the local-to-global maps in each wellElemSubRegion
    wellRegion.generateWell( meshLevel, lineBlock, nodeOffsetGlobal + wellNodeCount, elemOffsetGlobal + wellElemCount );

    // increment counters with global number of nodes and elements
    wellElemCount += lineBlock.numElements();
    wellNodeCount += lineBlock.numNodes();

    string const & subRegionName = wellRegion.getSubRegionName();
    WellElementSubRegion &
    subRegion = wellRegion.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                  .getGroup< WellElementSubRegion >( subRegionName );

    globalIndex const numWellElemsGlobal = MpiWrapper::sum( subRegion.size() );

    GEOS_ERROR_IF( numWellElemsGlobal != lineBlock.numElements(),
                   "Invalid partitioning in well " << lineBlock.getDataContext() <<
                   ", subregion " << subRegion.getDataContext() );

  } );

  // communicate to rebuild global node info since we modified global ordering
  nodeManager.setMaxGlobalIndex();
}

void ElementRegionManager::buildSets( NodeManager const & nodeManager )
{
  GEOS_MARK_FUNCTION;

  dataRepository::Group const & nodeSets = nodeManager.sets();

  map< string, array1d< bool > > nodeInSet; // map to contain indicator of whether a node is in a set.
  string_array setNames; // just a holder for the names of the sets

  // loop over all wrappers and fill the nodeIndSet arrays for each set
  for( auto & wrapper: nodeSets.wrappers() )
  {
    string const & name = wrapper.second->getName();
    nodeInSet[name].resize( nodeManager.size() );
    nodeInSet[name].setValues< serialPolicy >( false );

    if( nodeSets.hasWrapper( name ) )
    {
      setNames.emplace_back( name );
      SortedArrayView< localIndex const > const & set = nodeSets.getReference< SortedArray< localIndex > >( name );
      for( localIndex const a: set )
      {
        nodeInSet[name][a] = true;
      }
    }
  }

  this->forElementSubRegions(
    [&]( auto & subRegion ) -> void
  {
    dataRepository::Group & elementSets = subRegion.sets();

    auto const & elemToNodeMap = subRegion.nodeList();

    for( string const & setName: setNames )
    {
      arrayView1d< bool const > const nodeInCurSet = nodeInSet[setName];

      SortedArray< localIndex > & targetSet = elementSets.registerWrapper< SortedArray< localIndex > >( setName ).reference();
      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        localIndex const numNodes = subRegion.numNodesPerElement( k );

        localIndex elementInSet = true;
        for( localIndex i = 0; i < numNodes; ++i )
        {
          if( !nodeInCurSet( elemToNodeMap[k][i] ) )
          {
            elementInSet = false;
            break;
          }
        }

        if( elementInSet )
        {
          targetSet.insert( k );
        }
      }
    }
  } );
}

int ElementRegionManager::packSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packImpl< false >( junk, packList );
}

int ElementRegionManager::pack( buffer_unit_type * & buffer,
                                ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return packImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
int ElementRegionManager::packImpl( buffer_unit_type * & buffer,
                                    ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  int packedSize = 0;

//  packedSize += Group::Pack( buffer, wrapperNames, {}, 0, 0);

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, this->getName() );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, numRegions() );

  parallelDeviceEvents events;
  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.getName() );

    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.numSubRegions() );

    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const elemList = packList[kReg][esr];
      if( DO_PACKING )
      {
        packedSize += subRegion.pack( buffer, elemList, 0, false, events );
      }
      else
      {
        packedSize += subRegion.packSize( elemList, 0, false, events );
      }
    } );
  }

  waitAllDeviceEvents( events );
  return packedSize;
}


int ElementRegionManager::unpack( buffer_unit_type const * & buffer,
                                  ElementViewAccessor< arrayView1d< localIndex > > & packList )
{
  return unpackImpl( buffer, packList );
}

int ElementRegionManager::unpack( buffer_unit_type const * & buffer,
                                  ElementReferenceAccessor< array1d< localIndex > > & packList )
{
  return unpackImpl( buffer, packList );
}

template< typename T >
int ElementRegionManager::unpackImpl( buffer_unit_type const * & buffer,
                                      T & packList )
{
  int unpackedSize = 0;

  string name;
  unpackedSize += bufferOps::Unpack( buffer, name );

  GEOS_ERROR_IF( name != this->getName(), "Unpacked name (" << name << ") does not equal object name (" << this->getName() << ")" );

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  parallelDeviceEvents events;
  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase & elemRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG??
      arrayView1d< localIndex > & elemList = packList[kReg][esr];

      unpackedSize += subRegion.unpack( buffer, elemList, 0, false, events );
    } );
  }

  waitAllDeviceEvents( events );
  return unpackedSize;
}

int ElementRegionManager::packGlobalMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packGlobalMapsImpl< false >( junk, packList );
}

int ElementRegionManager::packGlobalMaps( buffer_unit_type * & buffer,
                                          ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return packGlobalMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
int ElementRegionManager::packGlobalMapsImpl( buffer_unit_type * & buffer,
                                              ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.getName() );

    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.numSubRegions() );
    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const elemList = packList[kReg][esr];
      if( DO_PACKING )
      {
        packedSize += subRegion.packGlobalMaps( buffer, elemList, 0 );
      }
      else
      {
        packedSize += subRegion.packGlobalMapsSize( elemList, 0 );
      }
    } );
  }

  return packedSize;
}


int
ElementRegionManager::unpackGlobalMaps( buffer_unit_type const * & buffer,
                                        ElementViewAccessor< ReferenceWrapper< localIndex_array > > & packList )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  packList.resize( numRegionsRead );
  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase & elemRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    packList[kReg].resize( numSubRegionsRead );
    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & elemList = packList[kReg][esr].get();

      unpackedSize += subRegion.unpackGlobalMaps( buffer, elemList, 0 );
    } );
  }

  return unpackedSize;
}



int ElementRegionManager::packUpDownMapsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}
int ElementRegionManager::packUpDownMapsSize( ElementReferenceAccessor< array1d< localIndex > > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

int ElementRegionManager::packUpDownMaps( buffer_unit_type * & buffer,
                                          ElementViewAccessor< arrayView1d< localIndex > > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}
int ElementRegionManager::packUpDownMaps( buffer_unit_type * & buffer,
                                          ElementReferenceAccessor< array1d< localIndex > > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING, typename T >
int ElementRegionManager::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                              T const & packList ) const
{
  int packedSize = 0;

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.getName() );

    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.numSubRegions() );
    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, subRegion.getName() );

      arrayView1d< localIndex > const elemList = packList[kReg][esr];
      if( DO_PACKING )
      {
        packedSize += subRegion.packUpDownMaps( buffer, elemList );
      }
      else
      {
        packedSize += subRegion.packUpDownMapsSize( elemList );
      }
    } );
  }

  return packedSize;
}

//template int
//ElementRegionManager::
//PackUpDownMapsImpl<true>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;
//template int
//ElementRegionManager::
//PackUpDownMapsImpl<false>( buffer_unit_type * & buffer,
//                             ElementViewAccessor<arrayView1d<localIndex>> const & packList ) const;

int
ElementRegionManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                        ElementReferenceAccessor< localIndex_array > & packList,
                                        bool const overwriteMap )
{
  int unpackedSize = 0;

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase & elemRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const kSubReg, ElementSubRegionBase & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & elemList = packList[kReg][kSubReg];
      unpackedSize += subRegion.unpackUpDownMaps( buffer, elemList, false, overwriteMap );
    } );
  }

  return unpackedSize;
}

int ElementRegionManager::packFracturedElementsSize( ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                                     string const fractureRegionName ) const
{
  buffer_unit_type * junk = nullptr;
  return packFracturedElementsImpl< false >( junk, packList, fractureRegionName );
}

int ElementRegionManager::packFracturedElements( buffer_unit_type * & buffer,
                                                 ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                                 string const fractureRegionName ) const
{
  return packFracturedElementsImpl< true >( buffer, packList, fractureRegionName );
}

template< bool DO_PACKING >
int ElementRegionManager::packFracturedElementsImpl( buffer_unit_type * & buffer,
                                                     ElementViewAccessor< arrayView1d< localIndex > > const & packList,
                                                     string const fractureRegionName ) const
{
  int packedSize = 0;

  SurfaceElementRegion const & embeddedSurfaceRegion =
    this->getRegion< SurfaceElementRegion >( fractureRegionName );
  EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
    embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  arrayView1d< globalIndex const > const embeddedSurfacesLocalToGlobal =
    embeddedSurfaceSubRegion.localToGlobalMap();

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, numRegions() );

  for( typename dataRepository::indexType kReg=0; kReg<numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = getRegion( kReg );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.getName() );

    packedSize += bufferOps::Pack< DO_PACKING >( buffer, elemRegion.numSubRegions() );
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >(
      [&]( localIndex const esr, CellElementSubRegion const & subRegion )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, subRegion.getName() );

      arrayView1d< localIndex const > const elemList = packList[kReg][esr];
      if( DO_PACKING )
      {
        packedSize += subRegion.packFracturedElements( buffer, elemList, embeddedSurfacesLocalToGlobal );
      }
      else
      {
        packedSize += subRegion.packFracturedElementsSize( elemList, embeddedSurfacesLocalToGlobal );
      }
    } );
  }

  return packedSize;
}

int
ElementRegionManager::unpackFracturedElements( buffer_unit_type const * & buffer,
                                               ElementReferenceAccessor< localIndex_array > & packList,
                                               string const fractureRegionName )
{
  int unpackedSize = 0;

  SurfaceElementRegion & embeddedSurfaceRegion =
    this->getRegion< SurfaceElementRegion >( fractureRegionName );
  EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion =
    embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  unordered_map< globalIndex, localIndex > embeddedSurfacesGlobalToLocal =
    embeddedSurfaceSubRegion.globalToLocalMap();

  localIndex numRegionsRead;
  unpackedSize += bufferOps::Unpack( buffer, numRegionsRead );

  for( localIndex kReg=0; kReg<numRegionsRead; ++kReg )
  {
    string regionName;
    unpackedSize += bufferOps::Unpack( buffer, regionName );

    ElementRegionBase & elemRegion = getRegion( regionName );

    localIndex numSubRegionsRead;
    unpackedSize += bufferOps::Unpack( buffer, numSubRegionsRead );
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >(
      [&]( localIndex const kSubReg, CellElementSubRegion & subRegion )
    {
      string subRegionName;
      unpackedSize += bufferOps::Unpack( buffer, subRegionName );

      /// THIS IS WRONG
      localIndex_array & elemList = packList[kReg][kSubReg];
      unpackedSize += subRegion.unpackFracturedElements( buffer, elemList, embeddedSurfacesGlobalToLocal );
    } );
  }

  return unpackedSize;
}

array2d< localIndex >
ElementRegionManager::getCellBlockToSubRegionMap( CellBlockManagerABC const & cellBlockManager ) const
{
  array2d< localIndex > blockMap( cellBlockManager.getCellBlocks().numSubGroups(), 2 );
  blockMap.setValues< serialPolicy >( -1 );

  Group::subGroupMap const & cellBlocks = cellBlockManager.getCellBlocks().getSubGroups();

  forElementSubRegionsComplete< CellElementSubRegion >( [blockMap = blockMap.toView(),
                                                         &cellBlocks]( localIndex const er,
                                                                       localIndex const esr,
                                                                       ElementRegionBase const & region,
                                                                       CellElementSubRegion const & subRegion )
  {
    GEOS_UNUSED_VAR( region ); // unused if geos_error_if is nulld
    localIndex const blockIndex = cellBlocks.getIndex( subRegion.getName() );
    GEOS_ERROR_IF( blockIndex == Group::subGroupMap::KeyIndex::invalid_index,
                   GEOS_FMT( "{}, subregion {}: Cell block not found at index {}.",
                             region.getDataContext().toString(), subRegion.getName(), blockIndex ) );
    GEOS_ERROR_IF( blockMap( blockIndex, 1 ) != -1,
                   GEOS_FMT( "{}, subregion {}: Cell block at index {} is mapped to more than one subregion.",
                             region.getDataContext().toString(), subRegion.getName(), blockIndex ) );

    blockMap( blockIndex, 0 ) = er;
    blockMap( blockIndex, 1 ) = esr;
  } );

  return blockMap;
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, ElementRegionManager, string const &, Group * const )
}
