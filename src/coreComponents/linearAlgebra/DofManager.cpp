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
 * @file DofManager.cpp
 */

#include "DofManager.hpp"

#include "common/FieldSpecificationOps.hpp"
#include "common/TypeDispatch.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/ReverseCutHillMcKeeOrdering.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/NodeManager.hpp"

#include "DofManagerHelpers.hpp"

#include <numeric>
#include <functional>

namespace geos
{

using namespace dataRepository;

DofManager::DofManager( string name )
  : m_name( std::move( name ) )
{}

void DofManager::clear()
{
  // deallocate index arrays from the mesh
  for( FieldDescription const & field : m_fields )
  {
    removeIndexArray( field );
  }

  // delete internal data
  m_fields.clear();
  m_coupling.clear();
  m_reordered = false;
}

void DofManager::setDomain( DomainPartition & domain )
{
  clear();
  m_domain = &domain;
}

localIndex DofManager::getFieldIndex( string const & name ) const
{
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  GEOS_ASSERT_MSG( it != m_fields.end(), "DofManager: field does not exist: " << name );
  return std::distance( m_fields.begin(), it );
}

bool DofManager::fieldExists( string const & name ) const
{
  auto const it = std::find_if( m_fields.begin(), m_fields.end(),
                                [&]( FieldDescription const & f ) { return f.name == name; } );
  return it != m_fields.end();
}

string const & DofManager::getKey( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].key;
}

globalIndex DofManager::numGlobalDofs( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].numGlobalDof;
}

globalIndex DofManager::numGlobalDofs() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), globalIndex( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.numGlobalDof; } );
}

localIndex DofManager::numLocalDofs( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].numLocalDof;
}

localIndex DofManager::numLocalDofs() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), localIndex( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.numLocalDof; } );
}

globalIndex DofManager::rankOffset( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].rankOffset;
}

globalIndex DofManager::rankOffset() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), globalIndex( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.rankOffset; } );
}

integer DofManager::numComponents( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].numComponents;
}

integer DofManager::numComponents() const
{
  return std::accumulate( m_fields.begin(), m_fields.end(), integer( 0 ),
                          []( globalIndex const n, FieldDescription const & f ) { return n + f.numComponents; } );
}

array1d< integer > DofManager::numComponentsPerField() const
{
  array1d< integer > ret( m_fields.size() );
  std::transform( m_fields.begin(), m_fields.end(), ret.begin(), std::mem_fn( &FieldDescription::numComponents ) );
  return ret;
}

FieldLocation DofManager::location( string const & fieldName ) const
{
  return m_fields[getFieldIndex( fieldName )].location;
}

globalIndex DofManager::globalOffset( string const & fieldName ) const
{
  GEOS_ASSERT_MSG( m_reordered, "Global offset not available until after reorderByRank() has been called." );
  return m_fields[getFieldIndex( fieldName )].globalOffset;
}

namespace
{

template< typename FUNC >
void forMeshSupport( std::vector< DofManager::FieldSupport > const & support,
                     DomainPartition & domain,
                     FUNC && func )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel & meshLevel = meshBody.getMeshLevel( regions.meshLevelName );
    func( meshBody, meshLevel, regions.regionNames );
  }
}

template< typename FUNC >
void forMeshSupport( std::vector< DofManager::FieldSupport > const & support,
                     DomainPartition const & domain,
                     FUNC && func )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & meshLevel = meshBody.getMeshLevel( regions.meshLevelName );
    func( meshBody, meshLevel, regions.regionNames );
  }
}

void fillTrivialPermutation( arrayView1d< localIndex > const permutation )
{
  forAll< parallelHostPolicy >( permutation.size(), [&]( localIndex const i )
  {
    permutation[i] = i;
  } );
}

} // namespace

array1d< localIndex > DofManager::computePermutation( FieldDescription & field )
{
  localIndex const fieldIndex = getFieldIndex( field.name );

  // step 1: save the number of components, and then set it to 1 temporarily
  //         do not forget to restore at the end
  //         we set the number of components to 1 to compute the reordering on a smaller matrix

  localIndex const numComps = field.numComponents;
  CompMask const globallyCoupledComps = field.globallyCoupledComponents;
  field.numComponents = 1;
  field.globallyCoupledComponents = CompMask( 1, true );

  // step 2: compute field dimensions (local dofs, global dofs, etc)
  //         this is needed to make sure that the sparsity pattern function work properly
  //         in particular, this function defines the rankOffset (computed with number of components = 1)

  computeFieldDimensions( fieldIndex );

  // the number of local dofs is available at this point, we allocate space for the permutation
  array1d< localIndex > permutation( numLocalDofs( field.name ) );

  // if no reordering is requesting, we just return the identity permutation
  if( field.reorderingType == LocalReorderingType::None )
  {
    fillTrivialPermutation( permutation );
  }
  else
  {
    fillTrivialPermutation( permutation );
    computePermutation( field, permutation );
  }

  // reset the number of components
  field.numComponents = numComps;
  field.globallyCoupledComponents = globallyCoupledComps;
  // reset the offsets, since they will be recomputed with the proper number of components
  // this is important to get the right reordering for multiphysics problems
  field.numLocalDof  = 0;
  field.rankOffset   = 0;
  field.numGlobalDof = 0;
  field.globalOffset = 0;
  field.blockOffset = 0;

  return permutation;
}

void DofManager::computePermutation( FieldDescription const & field,
                                     arrayView1d< localIndex > const permutation )
{
  localIndex const fieldIndex = getFieldIndex( field.name );

  // step 3: allocate and fill the dofNumber array
  //         this is needed to have dofNumber computed with numComps = 1
  //         note that the dofNumbers will be recomputed once we have obtained the permutation

  createIndexArray( field, permutation.toViewConst() );

  // step 4: compute the local sparsity pattern for this field

  SparsityPattern< globalIndex > pattern;
  array1d< localIndex > rowSizes( numLocalDofs( field.name ) );
  // in a first pass, count the row lengths to allocate enough space
  countRowLengthsOneBlock( rowSizes, fieldIndex, fieldIndex );

  // resize the sparsity pattern now that we know the row sizes
  pattern.resizeFromRowCapacities< parallelHostPolicy >( numLocalDofs( field.name ),
                                                         numGlobalDofs( field.name ),
                                                         rowSizes.data() );

  // compute the sparsity pattern
  setSparsityPatternOneBlock( pattern.toView(), fieldIndex, fieldIndex );

  // step 5: call the reordering function
  //         the goal of this step is to fill the permutation array

  localIndex const * const offsets = pattern.getOffsets();
  globalIndex const * const columns = pattern.getColumns();
  array1d< localIndex > reversePermutation( permutation.size() );

  if( field.reorderingType == LocalReorderingType::ReverseCutHillMcKee )
  {
    reverseCutHillMcKeeOrdering::
      computePermutation( offsets, columns, rankOffset(), reversePermutation );
  }
  else
  {
    GEOS_ERROR( "This local ordering type is not supported yet" );
  }

  forAll< parallelHostPolicy >( permutation.size(), [&]( localIndex const i )
  {
    permutation[reversePermutation[i]] = i;
  } );
}

void DofManager::createIndexArray( FieldDescription const & field,
                                   arrayView1d< localIndex const > const permutation )
{
  LocationSwitch( field.location, [&]( auto const loc )
  {
    localIndex index = 0;

    forMeshSupport( field.support, *m_domain, [&]( MeshBody &, MeshLevel & mesh, auto const & regions )
    {
      FieldLocation constexpr LOC = decltype(loc)::value;
      using helper = ArrayHelper< globalIndex, LOC >;

      FieldIdentifiers fieldsToBeSync;
      if( field.location == FieldLocation::Elem )
      {
        fieldsToBeSync.addElementFields( { field.key }, regions );
      }
      else
      {
        fieldsToBeSync.addFields( field.location, { field.key } );
      }

      // register index array
      helper::template create<>( mesh, field.key, field.docstring, regions );
      typename helper::Accessor indexArray = helper::get( mesh, field.key );

      // populate index array using a sequential counter
      forMeshLocation< LOC, false, serialPolicy >( mesh, regions, [&]( auto const locIdx )
      {
        helper::reference( indexArray, locIdx ) = field.rankOffset + field.numComponents * permutation[index++];
      } );

      // synchronize across ranks
      CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, m_domain->getNeighbors(), false );
    } );
  } );
}

void DofManager::removeIndexArray( FieldDescription const & field )
{
  LocationSwitch( field.location, [&]( auto const loc )
  {
    forMeshSupport( field.support, *m_domain, [&]( MeshBody &, MeshLevel & mesh, auto const & regions )
    {
      FieldLocation constexpr LOC = decltype(loc)::value;
      ArrayHelper< globalIndex, LOC >::template remove<>( mesh, field.key, regions );
    } );
  } );
}

void DofManager::computeFieldDimensions( localIndex const fieldIndex )
{
  FieldDescription & field = m_fields[fieldIndex];

  // determine number of local support points
  localIndex numLocalSupport = 0;
  forMeshSupport( field.support, *m_domain, [&]( MeshBody &, MeshLevel & mesh, auto const & regions )
  {
    numLocalSupport += countMeshObjects< false >( field.location, mesh, regions );
  } );
  field.numLocalDof = field.numComponents * numLocalSupport;

  // gather dof counts across ranks
  field.rankOffset = MpiWrapper::prefixSum< globalIndex >( field.numLocalDof );
  field.numGlobalDof = field.rankOffset + field.numLocalDof;
  MpiWrapper::broadcast( field.numGlobalDof, MpiWrapper::commSize() - 1 );

  // determine field's offsets
  if( fieldIndex > 0 )
  {
    FieldDescription const & prev = m_fields[fieldIndex - 1];
    field.blockOffset = prev.blockOffset + prev.numGlobalDof;
  }
  field.globalOffset = field.rankOffset; // actual value computed in reorderByRank()
}

namespace
{

std::set< string >
processFieldRegionList( MeshLevel const & mesh,
                        std::set< string > inputList )
{
  std::set< string > regions( std::move( inputList ) );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  if( regions.empty() )
  {
    elemManager.forElementRegions( [&]( ElementRegionBase const & region )
    {
      regions.insert( region.getName() );
    } );
  }
  else
  {
    for( string const & regionName : regions )
    {
      // Check existence, discard return value
      elemManager.getRegion( regionName );
    }
  }

  return regions;
}

std::vector< DofManager::FieldSupport >
processFieldSupportList( DomainPartition const & domain,
                         std::vector< DofManager::FieldSupport > const & inputList )
{
  std::vector< DofManager::FieldSupport > result;
  std::set< std::pair< string, string > > processedMeshLevels;
  for( DofManager::FieldSupport const & r : inputList )
  {
    GEOS_ERROR_IF( processedMeshLevels.count( { r.meshBodyName, r.meshLevelName } ) > 0,
                   GEOS_FMT( "Duplicate mesh: {}/{}", r.meshBodyName, r.meshLevelName ) );
    processedMeshLevels.insert( { r.meshBodyName, r.meshLevelName } );

    MeshBody const & body = domain.getMeshBody( r.meshBodyName );
    MeshLevel const & mesh = body.getMeshLevel( r.meshLevelName );
    result.push_back( { r.meshBodyName, r.meshLevelName, processFieldRegionList( mesh, r.regionNames ) } );
  }
  return result;
}

void addNewSupports( std::vector< DofManager::FieldSupport > const & inputSupport,
                     std::vector< DofManager::FieldSupport > & fieldSupport )
{
  for( auto const & newRegion : inputSupport )
  {
    bool added = false;
    for( auto & support : fieldSupport )
    {
      added = support.add( newRegion );
      if( added )
        break;
    }
    if( !added )
    {
      fieldSupport.push_back( newRegion );
    }
  }
}

} // namespace

void DofManager::addField( string const & fieldName,
                           FieldLocation const location,
                           integer const components,
                           std::vector< FieldSupport > const & regions )
{
  GEOS_ASSERT_MSG( m_domain != nullptr, "Domain has not been set" );
  GEOS_ERROR_IF( m_reordered, "Cannot add fields after reorderByRank() has been called." );
  GEOS_ERROR_IF_GT_MSG( components, MAX_COMP, "Number of components limit exceeded" );

  std::vector< FieldSupport > processedSupports = processFieldSupportList( *m_domain, regions );

  if( !fieldExists( fieldName ))
  {
    // populate basic info
    m_fields.emplace_back();
    FieldDescription & field = m_fields.back();
    field.name = fieldName;
    field.location = location;
    field.numComponents = components;
    field.globallyCoupledComponents = CompMask( components, true ); // everything globally coupled
    field.key = m_name + '_' + fieldName + "_dofIndex";
    field.docstring = fieldName + " DoF indices";
    // advanced processing
    field.support = processedSupports;
  }
  else
  {
    addNewSupports( processedSupports, m_fields[getFieldIndex( fieldName )].support );
  }
}

void DofManager::addField( string const & fieldName,
                           FieldLocation const location,
                           integer const components,
                           map< std::pair< string, string >, array1d< string > > const & regions )
{
  // Convert input into internal format
  std::vector< FieldSupport > support;
  for( auto const & p : regions )
  {
    MeshBody const & meshBody = m_domain->getMeshBody( p.first.first );
    MeshLevel const & mesh = meshBody.getMeshLevel( p.first.second ).getShallowParent();
    std::set< string > regionNames( p.second.begin(), p.second.end() );
    support.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
  }
  addField( fieldName, location, components, support );
}

void DofManager::setLocalReorderingType( string const & fieldName,
                                         LocalReorderingType const reorderingType )
{
  FieldDescription & field = m_fields[getFieldIndex( fieldName )];
  field.reorderingType = reorderingType;
}

void DofManager::disableGlobalCouplingForEquation( string const & fieldName,
                                                   integer const c )
{
  FieldDescription & field = m_fields[getFieldIndex( fieldName )];
  field.globallyCoupledComponents.unset( c ); // this equation will not interact with neighbors
}

void DofManager::disableGlobalCouplingForEquations( string const & fieldName,
                                                    arrayView1d< integer const > const components )
{
  FieldDescription & field = m_fields[getFieldIndex( fieldName )];
  for( integer c = 0; c < components.size(); ++c )
  {
    field.globallyCoupledComponents.unset( c ); // this equation will not interact with neighbors
  }
}

namespace
{

std::set< string >
processCouplingRegionList( std::set< string > inputList,
                           std::set< string > const & rowFieldRegions,
                           string const & rowFieldName,
                           std::set< string > const & colFieldRegions,
                           string const & colFieldName )
{
  std::set< string > regions( std::move( inputList ) );
  if( regions.empty() )
  {
    // Populate with regions common between two fields
    std::set_intersection( rowFieldRegions.begin(), rowFieldRegions.end(),
                           colFieldRegions.begin(), colFieldRegions.end(),
                           std::inserter( regions, regions.begin() ) );
  }
  else
  {
    // Check that both fields exist on all regions in the list
    auto const checkSupport = [&regions]( std::set< string > const & fieldRegions, string const & fieldName )
    {
      GEOS_UNUSED_VAR( fieldName ); // unused if geos_error_if is nulld
      // Both regions lists are sorted at this point
      GEOS_ERROR_IF( !std::includes( fieldRegions.begin(), fieldRegions.end(), regions.begin(), regions.end() ),
                     GEOS_FMT( "Coupling domain is not a subset of {}'s support:\nCoupling: {}\nField: {}",
                               fieldName, stringutilities::join( regions ), stringutilities::join( fieldRegions ) ) );
    };
    checkSupport( rowFieldRegions, rowFieldName );
    checkSupport( colFieldRegions, colFieldName );
  }

  return regions;
}

/// Comparison function that ignores region lists and only looks at mesh body/level names
template< typename OP >
struct RegionComp
{
  bool operator()( DofManager::FieldSupport const & l, DofManager::FieldSupport const & r ) const
  {
    return OP{} ( std::tie( l.meshBodyName, l.meshLevelName ), std::tie( r.meshBodyName, r.meshLevelName ) );
  }
};

std::vector< DofManager::FieldSupport >
processCouplingRegionList( std::vector< DofManager::FieldSupport > inputList,
                           std::vector< DofManager::FieldSupport > const & rowFieldRegions,
                           string const & rowFieldName,
                           std::vector< DofManager::FieldSupport > const & colFieldRegions,
                           string const & colFieldName )
{
  std::vector< DofManager::FieldSupport > regions( std::move( inputList ) );

  if( regions.empty() )
  {
    // First, compute a common set of mesh body/level pairs
    std::set_intersection( rowFieldRegions.begin(), rowFieldRegions.end(),
                           colFieldRegions.begin(), colFieldRegions.end(),
                           std::back_inserter( regions ), RegionComp< std::less<> >{} );
    // Now, compute intersections of region lists
    for( DofManager::FieldSupport & r : regions )
    {
      // Find the body/level pair in the column field list (unsorted range)
      auto const comp = [&r]( auto const & c ){ return RegionComp< std::equal_to<> >{} ( r, c ); };
      DofManager::FieldSupport const & colRegions = *std::find_if( colFieldRegions.begin(), colFieldRegions.end(), comp );
      // Intersect row field regions (already copied into result) with col field regions (found above)
      r.regionNames = processCouplingRegionList( {}, r.regionNames, rowFieldName, colRegions.regionNames, colFieldName );
    }
  }
  else
  {
    // Check that each input entry is included in both row and col field supports
    auto const checkSupport = [&regions]( std::vector< DofManager::FieldSupport > const & fieldRegions, string const & fieldName )
    {
      GEOS_UNUSED_VAR( fieldName ); // unused if geos_error_if is nulled
      for( DofManager::FieldSupport const & r : regions )
      {
        auto const comp = [&r]( auto const & f ){ return RegionComp< std::equal_to<> >{} ( r, f ); };
        auto const it = std::find_if( fieldRegions.begin(), fieldRegions.end(), comp );
        GEOS_ERROR_IF( it == fieldRegions.end(),
                       GEOS_FMT( "Mesh {}/{} not found in support of field {}", r.meshBodyName, r.meshLevelName, fieldName ) );

        GEOS_ERROR_IF( !std::includes( it->regionNames.begin(), it->regionNames.end(), r.regionNames.begin(), r.regionNames.end() ),
                       GEOS_FMT( "Coupling domain is not a subset of {}'s support:\nCoupling: {}\nField: {}",
                                 fieldName, stringutilities::join( r.regionNames ), stringutilities::join( it->regionNames ) ) );
      }
    };
    checkSupport( rowFieldRegions, rowFieldName );
    checkSupport( colFieldRegions, colFieldName );
  }

  return regions;
}

} // namespace

void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              Connector const connectivity,
                              std::vector< FieldSupport > const & supports,
                              bool const symmetric )
{
  GEOS_ASSERT_MSG( m_domain != nullptr, "Domain has not been set" );
  localIndex const rowFieldIndex = getFieldIndex( rowFieldName );
  localIndex const colFieldIndex = getFieldIndex( colFieldName );

  // Check if already defined
  std::vector< FieldSupport > processSupportList = processCouplingRegionList( supports,
                                                                              m_fields[rowFieldIndex].support,
                                                                              rowFieldName,
                                                                              m_fields[colFieldIndex].support,
                                                                              colFieldName );

  if( m_coupling.count( {rowFieldIndex, colFieldIndex} ) == 0 )
  {
    CouplingDescription & coupling = m_coupling[ { rowFieldIndex, colFieldIndex } ];
    coupling.connector = connectivity;
    if( connectivity == Connector::None && rowFieldIndex == colFieldIndex )
    {
      coupling.connector = static_cast< Connector >( m_fields[rowFieldIndex].location );
    }

    // Make a list of supports on which coupling is defined
    coupling.support = processSupportList;

    // Set connectivity with active symmetry flag
    if( symmetric && colFieldIndex != rowFieldIndex )
    {
      m_coupling.insert( { { colFieldIndex, rowFieldIndex }, coupling } );
    }
  }
  else
  {
    CouplingDescription & coupling = m_coupling[ { rowFieldIndex, colFieldIndex } ];
    addNewSupports( processSupportList, m_coupling[ { rowFieldIndex, colFieldIndex } ].support );

    // Set connectivity with active symmetry flag
    if( symmetric && colFieldIndex != rowFieldIndex )
    {
      m_coupling.insert( { { colFieldIndex, rowFieldIndex }, coupling } );
    }
  }
}

void DofManager::addCoupling( string const & fieldName,
                              FluxApproximationBase const & stencils )
{
  localIndex const fieldIndex = getFieldIndex( fieldName );
  FieldDescription const & field = m_fields[fieldIndex];

  GEOS_ERROR_IF( field.location != FieldLocation::Elem, "Field must be supported on elements in order to use stencil sparsity" );

  CouplingDescription & coupling = m_coupling[ {fieldIndex, fieldIndex} ];
  coupling.connector = Connector::Stencil;
  coupling.support = field.support;
  coupling.stencils = &stencils;
}

void DofManager::addCoupling( string const & rowFieldName,
                              string const & colFieldName,
                              DofManager::Connector connectivity,
                              map< std::pair< string, string >, array1d< string > > const & supports,
                              bool symmetric )
{
  // Convert input into internal format
  std::vector< FieldSupport > support;
  for( auto const & p : supports )
  {
    MeshBody const & meshBody = m_domain->getMeshBody( p.first.first );
    MeshLevel const & mesh = meshBody.getMeshLevel( p.first.second ).getShallowParent();
    std::set< string > regionNames( p.second.begin(), p.second.end() );
    support.push_back( { meshBody.getName(), mesh.getName(), std::move( regionNames ) } );
  }
  addCoupling( rowFieldName, colFieldName, connectivity, support, symmetric );
}

namespace
{

/**
 * @brief Implements a mesh loop that populates a local connector-to-location pattern.
 * @tparam LOC type of mesh locations
 * @tparam CONN type of mesh connector
 * @tparam SUBREGIONTYPES types of subregions to loop over
 *
 * The algorithm used requires that a connector-to-location map is present and populated in
 * the mesh data structure. We loop over connectors first, then visit adjacent locations from
 * each connector and thus populate sparsity pattern row-by-row. The numbering of connectors
 * is implicit, i.e. implied by the looping order, thus we don't have to compute/store it.
 *
 * The sparsity pattern thus constructed should only be contracted with another sparsity pattern that
 * was built with the same type of connectors and list of region names and types, otherwise connector
 * numbering will be incompatible between the two.
 */
template< FieldLocation LOC, FieldLocation CONN, typename ... SUBREGIONTYPES >
struct ConnLocPatternBuilder
{
  static void build( MeshLevel const & mesh,
                     string const & key,
                     DofManager::CompMask const & mask,
                     std::set< string > const & regions,
                     localIndex const rowOffset,
                     SparsityPattern< globalIndex > & connLocPattern )
  {
    using helper = ArrayHelper< globalIndex const, LOC >;
    typename helper::Accessor dofIndexArray = helper::get( mesh, key );
    auto connIdxPrev = MeshHelper< CONN >::invalid_local_index;

    localIndex connectorCount = -1;

    forMeshLocation< CONN, LOC, true, serialPolicy, SUBREGIONTYPES... >( mesh, regions,
                                                                         [&]( auto const & connIdx,
                                                                              auto const & locIdx,
                                                                              localIndex const GEOS_UNUSED_PARAM( k ) )
    {
      if( connIdx != connIdxPrev )
      {
        ++connectorCount;
        connIdxPrev = connIdx;
      }

      globalIndex const dofOffset = helper::value( dofIndexArray, locIdx );
      if( dofOffset >= 0 )
      {
        for( integer const c : mask )
        {
          connLocPattern.insertNonZero( rowOffset + connectorCount, dofOffset + c );
        }
      }
    } );
  }
};

/**
 * @brief A specialization of ConnLocPatternBuilder for elements connected through edges.
 * @tparam SUBREGIONTYPES types of subregions to loop over
 *
 * This is required because edge-to-element map does not exist in the mesh. Therefore, we have
 * to invert the loop order and go through elements first. This requires us to create a unique
 * connector numbering for edges in the regions of interest first - which comes with a memory
 * and runtime overhead, so we only use this method when absolutely have to, i.e. in this case.
 */
template< typename ... SUBREGIONTYPES >
struct ConnLocPatternBuilder< FieldLocation::Elem, FieldLocation::Edge, SUBREGIONTYPES... >
{
  static void build( MeshLevel const & mesh,
                     string const & key,
                     DofManager::CompMask const & mask,
                     std::set< string > const & regions,
                     localIndex const rowOffset,
                     SparsityPattern< globalIndex > & connLocPattern )
  {
    FieldLocation constexpr ELEM  = FieldLocation::Elem;
    FieldLocation constexpr EDGE = FieldLocation::Edge;

    EdgeManager const & edgeManager = mesh.getEdgeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofIndex =
      elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( key );

    array1d< localIndex > edgeConnectorIndex( edgeManager.size() );
    edgeConnectorIndex.setValues< serialPolicy >( -1 );

    localIndex edgeCount = 0;
    forMeshLocation< EDGE, true, serialPolicy, SUBREGIONTYPES... >( mesh, regions, [&] ( localIndex const edgeIdx )
    {
      edgeConnectorIndex[edgeIdx] = edgeCount++;
    } );

    forMeshLocation< ELEM, EDGE, true, serialPolicy, SUBREGIONTYPES... >( mesh, regions,
                                                                          [&]( auto const & elemIdx,
                                                                               auto const & edgeIdx,
                                                                               localIndex const GEOS_UNUSED_PARAM( k ) )
    {
      globalIndex const dofOffset = dofIndex[elemIdx[0]][elemIdx[1]][elemIdx[2]];
      if( dofOffset >= 0 )
      {
        for( integer const c : mask )
        {
          connLocPattern.insertNonZero( rowOffset + edgeConnectorIndex[edgeIdx], dofOffset + c );
        }
      }
    } );
  }
};

template< FieldLocation LOC, FieldLocation CONN >
void makeConnLocPattern( MeshLevel const & mesh,
                         string const & key,
                         localIndex const numComp,
                         DofManager::CompMask const & mask,
                         globalIndex const numGlobalDof,
                         std::set< string > const & regions,
                         SparsityPattern< globalIndex > & connLocPattern )
{
  using Loc = FieldLocation;

  // 1. Count number of connectors in order to size the sparsity pattern
  localIndex const numConnectors = countMeshObjects< CONN, true >( mesh, regions );
  localIndex numConnectorsTotal = numConnectors;

  // TPFA+fracture special case
  if( LOC == Loc::Elem && CONN == Loc::Face )
  {
    numConnectorsTotal += countMeshObjects< Loc::Edge, true, FaceElementSubRegion >( mesh, regions );
  }

  // 2. Resize the pattern appropriately
  localIndex const numEntriesPerRow = MeshIncidence< CONN, LOC >::max * numComp;
  connLocPattern.resize( numConnectorsTotal,
                         numGlobalDof,
                         std::min( LvArray::integerConversion< globalIndex >( numEntriesPerRow ), numGlobalDof ) );

  // 3. Populate the local CL pattern
  ConnLocPatternBuilder< LOC, CONN >::build( mesh, key, mask, regions, 0, connLocPattern );

  // TPFA+fracture special case
  if( LOC == Loc::Elem && CONN == Loc::Face )
  {
    ConnLocPatternBuilder< Loc::Elem, Loc::Edge, FaceElementSubRegion >::build( mesh, key, mask, regions, numConnectors, connLocPattern );
  }
}

} // namespace

void DofManager::setSparsityPatternFromStencil( SparsityPatternView< globalIndex > const & pattern,
                                                localIndex const fieldIndex ) const
{
  FieldDescription const & field = m_fields[fieldIndex];
  CouplingDescription const & coupling = m_coupling.at( {fieldIndex, fieldIndex} );
  localIndex const numComp = field.numComponents;
  globalIndex const rankDofOffset = rankOffset();
  CompMask const & globallyCoupledComps = field.globallyCoupledComponents;

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const dofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( field.key );

    array1d< globalIndex > rowDofIndices;
    array1d< globalIndex > colDofIndices;

    // 1. Assemble diagonal and off-diagonal blocks for elements in stencil
    coupling.stencils->forAllStencils( mesh, [&]( auto const & stencil )
    {
      using StenciType = typename std::decay< decltype( stencil ) >::type;
      constexpr localIndex maxNumFluxElems = StenciType::maxNumPointsInFlux;
      constexpr localIndex maxStencilSize = StenciType::maxStencilSize;

      typename StenciType::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      rowDofIndices.reserve( maxNumFluxElems );
      colDofIndices.reserve( maxStencilSize * numComp );

      forAll< serialPolicy >( stencil.size(), [&]( localIndex const iconn )
      {
        // This weirdness is because of fracture stencils, which don't have separate
        // getters for num flux elems vs stencil size... it won't work for MPFA though
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        localIndex const stencilSize = numFluxElems;

        rowDofIndices.resize( numFluxElems );
        for( localIndex i = 0; i < numFluxElems; ++i )
        {
          rowDofIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        }

        colDofIndices.resize( stencilSize * numComp );
        for( localIndex i = 0; i < stencilSize; ++i )
        {
          for( localIndex c = 0; c < numComp; ++c )
          {
            colDofIndices[i * numComp + c] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + c;
          }
        }

        std::sort( colDofIndices.begin(), colDofIndices.end() );

        for( localIndex i = 0; i < numFluxElems; ++i )
        {
          localIndex const localDofNumber = rowDofIndices[i] - rankDofOffset;
          if( localDofNumber >= 0 && localDofNumber < pattern.numRows() )
          {
            for( integer const c : globallyCoupledComps ) // add a non-zero for globally coupled components only
            {
              pattern.insertNonZeros( localDofNumber + c, colDofIndices.begin(), colDofIndices.end() );
            }
          }
        }
      } );
    } );

    // 2. Insert diagonal blocks, in case there are elements not included in stencil
    // (e.g. a single fracture element not connected to any other)
    auto dofNumberView = dofNumber.toNestedViewConst();
    colDofIndices.resize( numComp );
    forMeshLocation< FieldLocation::Elem, false, parallelHostPolicy >( mesh, regions, [=]( auto const & elemIdx )
    {
      globalIndex const elemDof = dofNumberView[elemIdx[0]][elemIdx[1]][elemIdx[2]];
      for( localIndex c = 0; c < numComp; ++c )
      {
        colDofIndices[c] = elemDof + c;
      }
      for( localIndex c = 0; c < numComp; ++c )
      {
        pattern.insertNonZeros( elemDof - rankDofOffset + c, colDofIndices.begin(), colDofIndices.end() );
      }
    } );
  } );
}

void DofManager::setSparsityPatternOneBlock( SparsityPatternView< globalIndex > const & pattern,
                                             localIndex const rowFieldIndex,
                                             localIndex const colFieldIndex ) const
{
  GEOS_ASSERT( rowFieldIndex >= 0 );
  GEOS_ASSERT( colFieldIndex >= 0 );

  FieldDescription const & rowField = m_fields[rowFieldIndex];
  FieldDescription const & colField = m_fields[colFieldIndex];

  if( m_coupling.count( {rowFieldIndex, colFieldIndex} ) == 0 )
  {
    return;
  }
  CouplingDescription const & coupling = m_coupling.at( {rowFieldIndex, colFieldIndex} );

  // Special treatment for stencil-based sparsity
  if( rowFieldIndex == colFieldIndex && coupling.connector == Connector::Stencil )
  {
    setSparsityPatternFromStencil( pattern, rowFieldIndex );
    return;
  }

  forMeshSupport( coupling.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    SparsityPattern< globalIndex > connLocRow( 0, 0, 0 ), connLocCol( 0, 0, 0 );

    LocationSwitch( rowField.location, static_cast< FieldLocation >( coupling.connector ),
                    [&]( auto const locType, auto const connType )
    {
      FieldLocation constexpr LOC  = decltype(locType)::value;
      FieldLocation constexpr CONN = decltype(connType)::value;

      makeConnLocPattern< LOC, CONN >( mesh,
                                       rowField.key,
                                       rowField.numComponents,
                                       rowField.globallyCoupledComponents, // select only globally coupled components
                                       rowField.numGlobalDof,
                                       regions,
                                       connLocRow );
    } );

    // we create a temporary mask including all components
    // if the colFieldIndex is equal to the rowFieldIndex and all components are in the mask,
    // we don't need to recomponent the connLocPattern
    CompMask const allComponentsMask( colField.numComponents, true );
    bool const areAllComponentsCoupled = colField.globallyCoupledComponents.begin() == allComponentsMask.begin();

    if( colFieldIndex == rowFieldIndex && areAllComponentsCoupled )
    {
      connLocCol = connLocRow; // TODO avoid copying
    }
    else
    {
      LocationSwitch( colField.location, static_cast< FieldLocation >( coupling.connector ),
                      [&]( auto const locType, auto const connType )
      {
        FieldLocation constexpr LOC = decltype(locType)::value;
        FieldLocation constexpr CONN = decltype(connType)::value;

        // note that below we activate all components using allComponentsMask and not colField.globallyCoupledComponents
        // for the following reasons:
        // - we want the **equation** in LOC to be decoupled from its neighbors (think volume balance constraint)
        // - however, we want the corresponding **dof** to always be coupled (think last component density)

        makeConnLocPattern< LOC, CONN >( mesh,
                                         colField.key,
                                         colField.numComponents,
                                         allComponentsMask, // select all components
                                         colField.numGlobalDof,
                                         regions,
                                         connLocCol );
      } );
    }
    GEOS_ASSERT_EQ( connLocRow.numRows(), connLocCol.numRows() );

    // First, we perform assembly/multiply patterns (at this step, locally coupled equations are excluded by construction)
    globalIndex const globalDofOffset = rankOffset();
    forAll< serialPolicy >( connLocRow.numRows(), [&]( localIndex const irow )
    {
      arraySlice1d< globalIndex const > const dofIndicesRow = connLocRow.getColumns( irow );
      arraySlice1d< globalIndex const > const dofIndicesCol = connLocCol.getColumns( irow );
      for( globalIndex const globalRow : dofIndicesRow )
      {
        localIndex const localRow = globalRow - globalDofOffset;
        if( localRow >= 0 && localRow < pattern.numRows() )
        {
          pattern.insertNonZeros( localRow, dofIndicesCol.begin(), dofIndicesCol.end() );
        }
      }
    } );

    // Finally, we add the diagonal terms of the locally coupled equations to the sparsity pattern
    // Note: this is only needed if some components in the row field are not globally coupled,
    // because the diagonal terms of the locally coupled components have not been added at the previous step
    if( rowFieldIndex == colFieldIndex && !areAllComponentsCoupled )
    {
      LocationSwitch( rowField.location, [&]( auto const loc )
      {
        FieldLocation constexpr LOC = decltype(loc)::value;
        using ArrayHelper = ArrayHelper< globalIndex const, LOC >;

        typename ArrayHelper::Accessor dofIndexArray = ArrayHelper::get( mesh, rowField.key );

        CompMask locallyCoupledComponents = rowField.globallyCoupledComponents;
        locallyCoupledComponents.invert();

        array1d< globalIndex > colDofIndices( rowField.numComponents );
        globalIndex const rankDofOffset = rankOffset();
        forMeshLocation< LOC, false, serialPolicy >( mesh, regions, [&]( auto const locIdx )
        {
          globalIndex const dofNumber = ArrayHelper::value( dofIndexArray, locIdx );
          for( localIndex c = 0; c < rowField.numComponents; ++c )
          {
            colDofIndices[c] = dofNumber + c;
          }
          // add a non-zero for locally coupled components since diagonal terms for globally coupled components have been added already
          for( integer const c : locallyCoupledComponents )
          {
            pattern.insertNonZeros( dofNumber - rankDofOffset + c, colDofIndices.begin(), colDofIndices.end() );
          }
        } );
      } );
    }
  } );
}

void DofManager::countRowLengthsFromStencil( arrayView1d< localIndex > const & rowLengths,
                                             localIndex const fieldIndex ) const
{
  FieldDescription const & field = m_fields[fieldIndex];
  CouplingDescription const & coupling = m_coupling.at( {fieldIndex, fieldIndex} );
  GEOS_ASSERT( coupling.connector == Connector::Stencil );
  GEOS_ASSERT( coupling.stencils != nullptr );

  localIndex const numComp = field.numComponents;
  globalIndex const rankDofOffset = rankOffset();
  CompMask const & globallyCoupledComps = field.globallyCoupledComponents;

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const dofNumberAccessor =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( field.key );

    // 1. Count row contributions from stencil
    coupling.stencils->forAllStencils( mesh, [&]( auto const & stencil )
    {
      using StenciType = typename std::decay< decltype( stencil ) >::type;
      typename StenciType::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename StenciType::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      forAll< parallelHostPolicy >( stencil.size(), [&]( localIndex const iconn )
      {
        // This weirdness is because of fracture stencils, which don't have separate
        // getters for num flux elems vs stencil size... it won't work for MPFA though
        localIndex const stencilSize = stencil.stencilSize( iconn );
        localIndex const numFluxElems = stencilSize;

        for( localIndex i = 0; i < numFluxElems; ++i )
        {
          localIndex const localDofNumber = dofNumberAccessor[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] - rankDofOffset;
          if( localDofNumber >= 0 && localDofNumber < rowLengths.size() )
          {
            for( integer const c : globallyCoupledComps ) // add a non-zero for globally coupled components only
            {
              RAJA::atomicAdd( parallelHostAtomic{}, &rowLengths[localDofNumber + c], ( stencilSize - 1 ) * numComp );
            }
          }
        }
      } );
    } );

    // 2. Add diagonal contributions to account for elements not in stencil
    mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex const, ElementSubRegionBase const & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( field.key );
      forAll< parallelHostPolicy >( subRegion.size(), [&]( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          localIndex const localDofNumber = dofNumber[ei] - rankDofOffset;
          for( localIndex c = 0; c < numComp; ++c )
          {
            rowLengths[localDofNumber + c] += numComp;
          }
        }
      } );
    } );
  } );
}

void DofManager::countRowLengthsOneBlock( arrayView1d< localIndex > const & rowLengths,
                                          localIndex const rowFieldIndex,
                                          localIndex const colFieldIndex ) const
{
  GEOS_ASSERT( rowFieldIndex >= 0 );
  GEOS_ASSERT( colFieldIndex >= 0 );

  FieldDescription const & rowField = m_fields[rowFieldIndex];
  FieldDescription const & colField = m_fields[colFieldIndex];

  if( m_coupling.count( {rowFieldIndex, colFieldIndex} ) == 0 )
  {
    return;
  }
  CouplingDescription const & coupling = m_coupling.at( {rowFieldIndex, colFieldIndex} );

  if( rowFieldIndex == colFieldIndex && coupling.connector == Connector::Stencil )
  {
    countRowLengthsFromStencil( rowLengths, rowFieldIndex );
    return;
  }

  forMeshSupport( coupling.support, *m_domain, [&]( MeshBody const &, MeshLevel const & mesh, auto const & regions )
  {
    SparsityPattern< globalIndex > connLocRow( 0, 0, 0 ), connLocCol( 0, 0, 0 );

    LocationSwitch( rowField.location, static_cast< FieldLocation >( coupling.connector ),
                    [&]( auto const locType, auto const connType )
    {
      FieldLocation constexpr LOC  = decltype(locType)::value;
      FieldLocation constexpr CONN = decltype(connType)::value;

      makeConnLocPattern< LOC, CONN >( mesh,
                                       rowField.key,
                                       rowField.numComponents,
                                       rowField.globallyCoupledComponents, // select only globally coupled components
                                       rowField.numGlobalDof,
                                       regions,
                                       connLocRow );
    } );

    // we create a temporary mask including all components
    // if the colFieldIndex is equal to the rowFieldIndex and all components are in the mask,
    // we don't need to recomponent the connLocPattern
    CompMask const allComponentsMask( colField.numComponents, true );
    bool const areAllComponentsCoupled = colField.globallyCoupledComponents.begin() == allComponentsMask.begin();

    if( colFieldIndex == rowFieldIndex && areAllComponentsCoupled )
    {
      connLocCol = connLocRow;
    }
    else
    {
      LocationSwitch( colField.location, static_cast< FieldLocation >( coupling.connector ),
                      [&]( auto const locType, auto const connType )
      {
        FieldLocation constexpr LOC = decltype(locType)::value;
        FieldLocation constexpr CONN = decltype(connType)::value;

        // note that below we activate all components using allComponentsMask and not colField.globallyCoupledComponents
        // for the following reasons:
        // - we want the **equation** in LOC to be decoupled from its neighbors (think volume balance constraint)
        // - however, we want the corresponding dof to always be coupled (think last component density)

        makeConnLocPattern< LOC, CONN >( mesh,
                                         colField.key,
                                         colField.numComponents,
                                         allComponentsMask, // select all components
                                         colField.numGlobalDof,
                                         regions,
                                         connLocCol );
      } );
    }

    GEOS_ASSERT_EQ( connLocRow.numRows(), connLocCol.numRows() );

    // Estimate an upper bound on row length by adding adjacent entries on a connector
    globalIndex const rankDofOffset = rankOffset();
    forAll< parallelHostPolicy >( connLocRow.numRows(), [&]( localIndex const iconn )
    {
      arraySlice1d< globalIndex const > const dofIndicesRow = connLocRow.getColumns( iconn );
      for( globalIndex const globalRow : dofIndicesRow )
      {
        localIndex const localRow = globalRow - rankDofOffset;
        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          RAJA::atomicAdd( parallelHostAtomic{}, &rowLengths[localRow], connLocCol.numNonZeros( iconn ) );
        }
      }
    } );

    // Finally, we add the diagonal terms of the locally coupled equations to the sparsity pattern
    // Note: this is only needed if some components in the row field are not globally coupled,
    // because the diagonal terms of the locally coupled components have not been added at the previous step
    if( rowFieldIndex == colFieldIndex && rowField.globallyCoupledComponents.begin() != allComponentsMask.begin() )
    {
      LocationSwitch( rowField.location, [&]( auto const loc )
      {
        FieldLocation constexpr LOC = decltype(loc)::value;
        using ArrayHelper = ArrayHelper< globalIndex const, LOC >;

        typename ArrayHelper::Accessor dofIndexArray = ArrayHelper::get( mesh, rowField.key );

        CompMask locallyCoupledComponents = rowField.globallyCoupledComponents;
        locallyCoupledComponents.invert();

        forMeshLocation< LOC, false, serialPolicy >( mesh, regions, [&]( auto const locIdx )
        {
          globalIndex const dofNumber = ArrayHelper::value( dofIndexArray, locIdx );
          localIndex const localDofNumber = dofNumber - rankDofOffset;
          if( localDofNumber >= 0 && localDofNumber < rowLengths.size() )
          {
            // add a non-zero for locally coupled components since diagonal terms for globally coupled components have been added already
            for( integer const c : locallyCoupledComponents )
            {
              RAJA::atomicAdd( parallelHostAtomic{}, &rowLengths[localDofNumber + c], rowField.numComponents );
            }
          }
        } );
      } );
    }
  } );
}

// Create the sparsity pattern (location-location). Low level interface
void DofManager::setSparsityPattern( SparsityPattern< globalIndex > & pattern ) const
{
  GEOS_ERROR_IF( !m_reordered, "Cannot set monolithic sparsity pattern before reorderByRank() has been called." );

  localIndex const numLocalRows = numLocalDofs();
  localIndex const numFields = LvArray::integerConversion< localIndex >( m_fields.size() );

  // Step 1. Do a dry run of sparsity construction to get the total number of nonzeros in each row
  array1d< localIndex > rowSizes( numLocalRows );
  for( localIndex blockRow = 0; blockRow < numFields; ++blockRow )
  {
    for( localIndex blockCol = 0; blockCol < numFields; ++blockCol )
    {
      countRowLengthsOneBlock( rowSizes, blockRow, blockCol );
    }
  }

  // Step 2. Allocate enough capacity for all nonzero entries in each row
  pattern.resizeFromRowCapacities< parallelHostPolicy >( numLocalRows, numGlobalDofs(), rowSizes.data() );

  // Step 3. Fill the sparsity block-by-block
  for( localIndex blockRow = 0; blockRow < numFields; ++blockRow )
  {
    for( localIndex blockCol = 0; blockCol < numFields; ++blockCol )
    {
      setSparsityPatternOneBlock( pattern.toView(), blockRow, blockCol );
    }
  }

  // Step 4. Compress to remove unused space between rows
  pattern.compress();
}

namespace
{

template< typename FIELD_OP, typename POLICY, typename FIELD_VIEW >
void vectorToFieldKernel( ObjectManagerBase & GEOS_UNUSED_PARAM( manager ),
                          arrayView1d< real64 const > const & localVector,
                          FIELD_VIEW const & field,
                          arrayView1d< globalIndex const > const & dofNumber,
                          arrayView1d< integer const > const & ghostRank,
                          real64 const scalingFactor,
                          localIndex const dofOffset,
                          DofManager::CompMask const mask )
{
  forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    if( ghostRank[i] < 0 && dofNumber[i] >= 0 )
    {
      localIndex const lid = dofNumber[i] - dofOffset;
      GEOS_ASSERT_GE( lid, 0 );

      integer fieldComp = 0;
      for( integer const vecComp : mask )
      {
        FIELD_OP::template SpecifyFieldValue( field,
                                              i,
                                              fieldComp++,
                                              scalingFactor * localVector[lid + vecComp] );
      }
    }
  } );
}

template< typename FIELD_OP, typename POLICY, typename FIELD_VIEW >
void vectorToFieldKernel( ObjectManagerBase & manager,
                          arrayView1d< real64 const > const & localVector,
                          FIELD_VIEW const & field,
                          arrayView1d< globalIndex const > const & dofNumber,
                          arrayView1d< integer const > const & ghostRank,
                          string const & scalingFactorName,
                          localIndex const dofOffset,
                          DofManager::CompMask const mask )
{
  arrayView1d< real64 const > const & scalingFactor = manager.getReference< array1d< real64 > >( scalingFactorName );

  forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    if( ghostRank[i] < 0 && dofNumber[i] >= 0 )
    {
      localIndex const lid = dofNumber[i] - dofOffset;
      GEOS_ASSERT_GE( lid, 0 );

      integer fieldComp = 0;
      for( integer const vecComp : mask )
      {
        FIELD_OP::template SpecifyFieldValue( field,
                                              i,
                                              fieldComp++,
                                              scalingFactor[i] * localVector[lid + vecComp] );
      }
    }
  } );
}

template< typename FIELD_OP, typename POLICY, typename SCALING_FACTOR_TYPE >
void vectorToFieldImpl( arrayView1d< real64 const > const & localVector,
                        ObjectManagerBase & manager,
                        string const & dofKey,
                        string const & fieldName,
                        SCALING_FACTOR_TYPE const & scalingFactor,
                        localIndex const dofOffset,
                        DofManager::CompMask const mask )
{
  arrayView1d< globalIndex const > const dofNumber = manager.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const ghostRank = manager.ghostRank();

  WrapperBase & wrapper = manager.getWrapperBase( fieldName );

  // Restrict primary solution fields to 1-2D real arrays,
  // because applying component index is not well defined for 3D and higher
  using FieldTypes = types::ListofTypeList< types::ArrayTypes< types::RealTypes, types::DimsUpTo< 2 > > >;
  types::dispatch( FieldTypes{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    vectorToFieldKernel< FIELD_OP, POLICY >( manager,
                                             localVector,
                                             wrapperT.reference().toView(),
                                             dofNumber,
                                             ghostRank,
                                             scalingFactor,
                                             dofOffset,
                                             mask );
  }, wrapper );
}

template< typename FIELD_OP, typename POLICY, typename FIELD_VIEW >
void fieldToVectorKernel( arrayView1d< real64 > const & localVector,
                          FIELD_VIEW const & field,
                          arrayView1d< globalIndex const > const & dofNumber,
                          arrayView1d< integer const > const & ghostRank,
                          real64 const GEOS_UNUSED_PARAM( scalingFactor ),
                          localIndex const dofOffset,
                          DofManager::CompMask const mask )
{
  forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    if( ghostRank[i] < 0 && dofNumber[i] >= 0 )
    {
      localIndex const lid = dofNumber[i] - dofOffset;
      GEOS_ASSERT_GE( lid, 0 );

      integer fieldComp = 0;
      for( integer const vecComp : mask )
      {
        FIELD_OP::template readFieldValue( field,
                                           i,
                                           fieldComp++,
                                           localVector[lid + vecComp] );
      }
    }
  } );
}

template< typename FIELD_OP, typename POLICY, typename SCALING_FACTOR_TYPE >
void fieldToVectorImpl( arrayView1d< real64 > const & localVector,
                        ObjectManagerBase const & manager,
                        string const & dofKey,
                        string const & fieldName,
                        SCALING_FACTOR_TYPE const & scalingFactor,
                        localIndex const dofOffset,
                        DofManager::CompMask const mask )
{
  arrayView1d< globalIndex const > const & dofNumber = manager.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = manager.ghostRank();

  WrapperBase const & wrapper = manager.getWrapperBase( fieldName );

  // Restrict primary solution fields to 1-2D real arrays,
  // because applying component index is not well defined for 3D and higher
  using FieldTypes = types::ListofTypeList< types::ArrayTypes< types::RealTypes, types::DimsUpTo< 2 > > >;
  types::dispatch( FieldTypes{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    Wrapper< ArrayType > const & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    fieldToVectorKernel< FIELD_OP, POLICY >( localVector,
                                             wrapperT.reference(),
                                             dofNumber,
                                             ghostRank,
                                             scalingFactor,
                                             dofOffset,
                                             mask );
  }, wrapper );
}

} // namespace

template< typename FIELD_OP, typename POLICY, typename SCALING_FACTOR_TYPE >
void DofManager::vectorToField( arrayView1d< real64 const > const & localVector,
                                string const & srcFieldName,
                                string const & dstFieldName,
                                SCALING_FACTOR_TYPE const & scalingFactor,
                                CompMask mask ) const
{
  FieldDescription const & field = m_fields[getFieldIndex( srcFieldName )];
  mask.setNumComp( field.numComponents );

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel & mesh, auto const & regions )
  {
    if( field.location == FieldLocation::Elem )
    {
      mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
      {
        vectorToFieldImpl< FIELD_OP, POLICY >( localVector,
                                               subRegion,
                                               field.key,
                                               dstFieldName,
                                               scalingFactor,
                                               rankOffset(),
                                               mask );
      } );
    }
    else
    {
      vectorToFieldImpl< FIELD_OP, POLICY >( localVector,
                                             getObjectManager( field.location, mesh ),
                                             field.key,
                                             dstFieldName,
                                             scalingFactor,
                                             rankOffset(),
                                             mask );
    }
  } );
}

void DofManager::copyVectorToField( arrayView1d< real64 const > const & localVector,
                                    string const & srcFieldName,
                                    string const & dstFieldName,
                                    real64 const scalingFactor,
                                    CompMask const mask ) const
{
  vectorToField< FieldSpecificationEqual, parallelDevicePolicy<> >( localVector,
                                                                    srcFieldName,
                                                                    dstFieldName,
                                                                    scalingFactor,
                                                                    mask );
}

template< typename SCALING_FACTOR_TYPE >
void DofManager::addVectorToField( arrayView1d< real64 const > const & localVector,
                                   string const & srcFieldName,
                                   string const & dstFieldName,
                                   SCALING_FACTOR_TYPE const & scalingFactor,
                                   CompMask const mask ) const
{
  vectorToField< FieldSpecificationAdd, parallelDevicePolicy<> >( localVector,
                                                                  srcFieldName,
                                                                  dstFieldName,
                                                                  scalingFactor,
                                                                  mask );
}
template void DofManager::addVectorToField< real64 >( arrayView1d< real64 const > const & localVector,
                                                      string const & srcFieldName,
                                                      string const & dstFieldName,
                                                      real64 const & scalingFactor,
                                                      CompMask const mask ) const;
template void DofManager::addVectorToField< char const * >( arrayView1d< real64 const > const & localVector,
                                                            string const & srcFieldName,
                                                            string const & dstFieldName,
                                                            char const * const & scalingFactor,
                                                            CompMask const mask ) const;

template< typename FIELD_OP, typename POLICY >
void DofManager::fieldToVector( arrayView1d< real64 > const & localVector,
                                string const & srcFieldName,
                                string const & dstFieldName,
                                real64 const scalingFactor,
                                CompMask mask ) const
{
  FieldDescription const & field = m_fields[getFieldIndex( srcFieldName )];
  mask.setNumComp( field.numComponents );

  forMeshSupport( field.support, *m_domain, [&]( MeshBody const &, MeshLevel & mesh, auto const & regions )
  {
    if( field.location == FieldLocation::Elem )
    {
      mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex const,
                                                                ElementSubRegionBase const & subRegion )
      {
        fieldToVectorImpl< FIELD_OP, POLICY >( localVector,
                                               subRegion,
                                               field.key,
                                               dstFieldName,
                                               scalingFactor,
                                               rankOffset(),
                                               mask );
      } );
    }
    else
    {
      fieldToVectorImpl< FIELD_OP, POLICY >( localVector,
                                             getObjectManager( field.location, mesh ),
                                             field.key,
                                             dstFieldName,
                                             scalingFactor,
                                             rankOffset(),
                                             mask );
    }
  } );
}

void DofManager::copyFieldToVector( arrayView1d< real64 > const & localVector,
                                    string const & srcFieldName,
                                    string const & dstFieldName,
                                    real64 const scalingFactor,
                                    CompMask const mask ) const
{
  fieldToVector< FieldSpecificationEqual, parallelDevicePolicy<> >( localVector,
                                                                    srcFieldName,
                                                                    dstFieldName,
                                                                    scalingFactor,
                                                                    mask );
}

void DofManager::addFieldToVector( arrayView1d< real64 > const & localVector,
                                   string const & srcFieldName,
                                   string const & dstFieldName,
                                   real64 const scalingFactor,
                                   CompMask const mask ) const
{
  fieldToVector< FieldSpecificationAdd, parallelDevicePolicy<> >( localVector,
                                                                  srcFieldName,
                                                                  dstFieldName,
                                                                  scalingFactor,
                                                                  mask );
}

void DofManager::reorderByRank()
{
  GEOS_LAI_ASSERT( !m_reordered );

  std::map< string, array1d< localIndex > > permutations;

  // First loop: compute the local permutation
  for( FieldDescription & field : m_fields )
  {
    // compute local permutation of dofs, if needed
    permutations[ field.name ] = computePermutation( field );
  }

  // Second loop: compute the dof number array
  for( FieldDescription & field : m_fields )
  {
    // compute field dimensions (local dofs, global dofs, etc)
    computeFieldDimensions( static_cast< localIndex >( getFieldIndex( field.name ) ) );
    // allocate and fill index array
    createIndexArray( field, permutations.at( field.name ).toViewConst() );
  }

  // update field offsets to account for renumbering
  globalIndex dofOffset = rankOffset();
  for( FieldDescription & field : m_fields )
  {
    field.globalOffset = dofOffset;
    dofOffset += field.numLocalDof;
  }

  // This is a map with a key that is the pair of strings specifying the
  // ( MeshBody name, MeshLevel name), and a value that is another map with a
  // key that indicates the name of the object that contains the field to be
  // synced, and a value that contans the name of the field to be synced.
  std::map< std::pair< string, string >, FieldIdentifiers > fieldsToBeSync;

  // adjust index arrays for owned locations
  for( FieldDescription const & field : m_fields )
  {
    globalIndex const adjustment = field.globalOffset - field.rankOffset;

    forMeshSupport( field.support, *m_domain, [&]( MeshBody const & body, MeshLevel & mesh, auto const & regions )
    {
      LocationSwitch( field.location, [&]( auto const loc )
      {
        FieldLocation constexpr LOC = decltype(loc)::value;
        using ArrayHelper = ArrayHelper< globalIndex, LOC >;

        typename ArrayHelper::Accessor indexArray = ArrayHelper::get( mesh, field.key );

        forMeshLocation< LOC, false, parallelHostPolicy >( mesh, regions, [&]( auto const locIdx )
        {
          ArrayHelper::reference( indexArray, locIdx ) += adjustment;
        } );

        if( field.location == FieldLocation::Elem )
        {
          fieldsToBeSync[{ body.getName(), mesh.getName() }].addElementFields( {field.key}, regions );

        }
        else
        {
          fieldsToBeSync[{ body.getName(), mesh.getName() }].addFields( field.location, {field.key} );
        }

      } );
    } );
  }

  for( auto const & meshFieldPair : fieldsToBeSync )
  {
    MeshLevel & mesh = m_domain->getMeshBody( meshFieldPair.first.first ).getMeshLevel( meshFieldPair.first.second );
    CommunicationTools::getInstance().synchronizeFields( meshFieldPair.second, mesh, m_domain->getNeighbors(), false );
  }

  m_reordered = true;
}

std::vector< DofManager::SubComponent >
DofManager::filterDofs( std::vector< SubComponent > const & excluded ) const
{
  std::vector< DofManager::SubComponent > result;
  for( const auto & field : m_fields )
  {
    auto const it = std::find_if( excluded.begin(), excluded.end(),
                                  [&]( SubComponent const & sc ) { return sc.fieldName == field.name; } );
    CompMask const includedComp = ( it != excluded.end() ) ? ~it->mask : ~CompMask( field.numComponents );
    if( !includedComp.empty() )
    {
      result.push_back( { field.name, includedComp } );
    }
  }

  return result;
}

void DofManager::setupFrom( DofManager const & source,
                            std::vector< SubComponent > const & selection )
{
  clear();
  for( FieldDescription const & field : source.m_fields )
  {
    auto const it = std::find_if( selection.begin(), selection.end(),
                                  [&]( SubComponent const & sc ) { return sc.fieldName == field.name; } );
    if( it != selection.end() )
    {
      SubComponent const & sc = *it;
      GEOS_ASSERT_GE( field.numComponents, sc.mask.size() );
      addField( field.name, field.location, sc.mask.size(), field.support );
    }
  }
  for( auto const & entry : source.m_coupling )
  {
    string const & rowFieldName = source.m_fields[entry.first.first].name;
    string const & colFieldName = source.m_fields[entry.first.second].name;
    if( fieldExists( rowFieldName ) && fieldExists( colFieldName ) )
    {
      localIndex const rowFieldIndex = getFieldIndex( rowFieldName );
      localIndex const colFieldIndex = getFieldIndex( colFieldName );
      m_coupling.insert( { {rowFieldIndex, colFieldIndex}, entry.second } );
    }
  }
  reorderByRank();
}

template< typename MATRIX >
void DofManager::makeRestrictor( std::vector< SubComponent > const & selection,
                                 MPI_Comm const & comm,
                                 bool const transpose,
                                 MATRIX & restrictor ) const
{
  GEOS_ERROR_IF( !m_reordered, "Cannot make restrictors before reorderByRank() has been called." );

  // 1. Populate selected fields and compute some basic dimensions
  // array1d< FieldDescription > fieldsSelected( selection.size() );
  std::vector< FieldDescription > fieldsSelected( selection.size() );

  for( std::size_t k = 0; k < fieldsSelected.size(); ++k )
  {
    SubComponent const & dof = selection[k];
    FieldDescription const & fieldOld = m_fields[getFieldIndex( dof.fieldName )];
    FieldDescription & fieldNew = fieldsSelected[k];

    fieldNew.name = selection[k].fieldName;
    fieldNew.numComponents = dof.mask.size();
    fieldNew.numLocalDof = fieldOld.numLocalDof / fieldOld.numComponents * fieldNew.numComponents;
    fieldNew.numGlobalDof = fieldOld.numGlobalDof / fieldOld.numComponents * fieldNew.numComponents;
    fieldNew.rankOffset = fieldOld.rankOffset / fieldOld.numComponents * fieldNew.numComponents;
  }

  // 2. Compute remaining offsets (we mostly just need globalOffset, but it depends on others)

  globalIndex blockOffset = 0;
  for( auto & field : fieldsSelected )
  {
    field.blockOffset = blockOffset;
    blockOffset += field.numGlobalDof;
  }

  globalIndex globalOffset = std::accumulate( fieldsSelected.begin(), fieldsSelected.end(), globalIndex( 0 ),
                                              []( localIndex const n, FieldDescription const & f )
  { return n + f.rankOffset; } );

  for( auto & field : fieldsSelected )
  {
    field.globalOffset = globalOffset;
    globalOffset += field.numLocalDof;
  }

  // 3. Build the restrictor field by field

  localIndex const numLocalDofSelected = std::accumulate( fieldsSelected.begin(), fieldsSelected.end(), localIndex( 0 ),
                                                          []( localIndex const n, FieldDescription const & f )
  { return n + f.numLocalDof; } );

  localIndex const rowSize = transpose ? numLocalDofs() : numLocalDofSelected;
  localIndex const colSize = transpose ? numLocalDofSelected : numLocalDofs();

  restrictor.createWithLocalSize( rowSize, colSize, 1, comm );
  restrictor.open();

  array1d< globalIndex > rows;
  array1d< globalIndex > cols;
  array1d< real64 > values;

  for( std::size_t k = 0; k < fieldsSelected.size(); ++k )
  {
    FieldDescription const & fieldNew = fieldsSelected[k];
    FieldDescription const & fieldOld = m_fields[getFieldIndex( fieldNew.name )];

    FieldDescription const & fieldRow = transpose ? fieldOld : fieldNew;
    FieldDescription const & fieldCol = transpose ? fieldNew : fieldOld;

    CompMask const mask = selection[k].mask;
    localIndex const numLocalNodes = fieldNew.numLocalDof / fieldNew.numComponents;


    rows.resize( numLocalNodes*mask.size() );
    cols.resize( numLocalNodes*mask.size() );
    values.resize( numLocalNodes*mask.size() );

    for( localIndex i = 0; i < numLocalNodes; ++i )
    {
      integer newComp = 0;
      for( integer const oldComp : mask )
      {
        globalIndex const row = fieldRow.globalOffset + i * fieldRow.numComponents + ( transpose ? oldComp : newComp );
        globalIndex const col = fieldCol.globalOffset + i * fieldCol.numComponents + ( transpose ? newComp : oldComp );

        rows( i*mask.size() + newComp ) = row;
        cols( i*mask.size() + newComp ) = col;
        values[ i*mask.size() + newComp ] = 1.0;
        ++newComp;
      }
    }

    restrictor.insert( rows.toViewConst(),
                       cols.toViewConst(),
                       values.toViewConst() );


  }
  restrictor.close();
}


void DofManager::printFieldInfo( std::ostream & os ) const
{
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
  {
    localIndex const numFields = LvArray::integerConversion< localIndex >( m_fields.size() );

    os << "Fields:" << std::endl;
    os << " # | " << std::setw( 20 ) << "name" << " | " << "comp" << " | " << "N global DOF" << std::endl;
    os << "---+----------------------+------+-------------" << std::endl;
    for( localIndex i = 0; i < numFields; ++i )
    {
      FieldDescription const & f = m_fields[i];
      os << ' ' << i << " | "
         << std::setw( 20 ) << f.name << " | "
         << std::setw( 4 ) << f.numComponents << " | "
         << std::setw( 12 ) << f.numGlobalDof << std::endl;
    }
    os << "---+----------------------+------+-------------" << std::endl;

    os << std::endl << "Connectivity:" << std::endl;
    for( localIndex i = 0; i < numFields; ++i )
    {
      for( localIndex j = 0; j < numFields; ++j )
      {
        if( m_coupling.count( {i, j} ) == 0 )
        {
          continue;
        }
        switch( m_coupling.at( {i, j} ).connector )
        {
          case Connector::Elem:
          {
            os << " E ";
            break;
          }
          case Connector::Face:
          {
            os << " F ";
            break;
          }
          case Connector::Edge:
          {
            os << " G ";
            break;
          }
          case Connector::Node:
          {
            os << " N ";
            break;
          }
          case Connector::None:
          {
            os << "   ";
            break;
          }
          case Connector::Stencil:
          {
            os << " S ";
            break;
          }
          default:
          {
            GEOS_ERROR( "Invalid connector type: " << static_cast< int >( m_coupling.at( {i, j} ).connector ) );
          }
        }
        if( j < numFields - 1 )
        {
          os << "|";
        }
      }
      os << std::endl;
      if( i < numFields - 1 )
      {
        for( localIndex j = 0; j < numFields - 1; ++j )
        {
          os << "---|";
        }
        os << "---";
      }
      os << std::endl;
    }
  }
}

#define MAKE_DOFMANAGER_METHOD_INST( LAI ) \
  template void DofManager::makeRestrictor( std::vector< SubComponent > const & selection, \
                                            MPI_Comm const & comm, \
                                            bool const transpose, \
                                            LAI::ParallelMatrix & restrictor ) const;

#ifdef GEOS_USE_TRILINOS
MAKE_DOFMANAGER_METHOD_INST( TrilinosInterface )
#endif

#ifdef GEOS_USE_HYPRE
MAKE_DOFMANAGER_METHOD_INST( HypreInterface )
#endif

#ifdef GEOS_USE_PETSC
MAKE_DOFMANAGER_METHOD_INST( PetscInterface )
#endif

} // namespace geos
