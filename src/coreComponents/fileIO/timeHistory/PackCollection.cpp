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

#include "PackCollection.hpp"

namespace geos
{
PackCollection::PackCollection ( string const & name, Group * parent )
  : HistoryCollectionBase( name, parent )
  , m_setsIndices( )
  , m_objectPath( )
  , m_fieldName( )
  , m_setNames( )
  , m_setChanged( true )
  , m_onlyOnSetChange( 0 )
  , m_initialized( false )
{
  registerWrapper( PackCollection::viewKeysStruct::objectPathString(), &m_objectPath ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the object from which to retrieve field values." );

  registerWrapper( PackCollection::viewKeysStruct::fieldNameString(), &m_fieldName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the (packable) field associated with the specified object to retrieve data from" );

  registerWrapper( PackCollection::viewKeysStruct::setNamesString(), &m_setNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The set(s) for which to retrieve data." );

  registerWrapper( PackCollection::viewKeysStruct::onlyOnSetChangeString(), &m_onlyOnSetChange ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "Whether or not to only collect when the collected sets of indices change in any way." );

  registerWrapper( PackCollection::viewKeysStruct::disableCoordCollectionString(), &m_disableCoordCollection ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "Whether or not to create coordinate meta-collectors if collected objects are mesh objects." );
}

void PackCollection::initializePostSubGroups( )
{
  if( !m_initialized )
  {
    DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
    m_collectionCount = collectAll() ? 1 : m_setNames.size();
    // determine whether we're collecting from a mesh object manager
    Group const * const targetObject = this->getTargetObject( domain, m_objectPath );
    ObjectManagerBase const * const objectManagerTarget = dynamic_cast< ObjectManagerBase const * >( targetObject );
    m_targetIsMeshObject = objectManagerTarget != nullptr;
    // update sets after we know whether to filter ghost indices ( m_targetIsMeshObject )
    updateSetsIndices( domain );
    HistoryCollectionBase::initializePostSubGroups();
    // build any persistent meta collectors that are required
    buildMetaDataCollectors();
    for( auto & metaCollector : m_metaDataCollectors )
    {
      // coord meta collectors should have m_disableCoordCollection == true to avoid
      //  infinite recursive init calls here
      metaCollector->initializePostSubGroups();
    }
    m_initialized = true;
  }
}

HistoryMetadata PackCollection::getMetaData( DomainPartition const & domain, localIndex collectionIdx ) const
{
  Group const * targetObject = this->getTargetObject( domain, m_objectPath );
  WrapperBase const & targetField = targetObject->getWrapperBase( m_fieldName );

  if( collectAll() )
  {
    // If we collect all the data, then we have a unique field: "all".
    // So it's safe to grab index 0.
    localIndex const packCount = m_setsIndices[0].size() == 0 ? targetField.size() : m_setsIndices[0].size();
    return targetField.getHistoryMetadata( packCount );
  }
  else
  {
    GEOS_ERROR_IF( collectionIdx < 0 || collectionIdx >= m_setNames.size(), "Invalid collection index specified." );
    localIndex collectionSize = m_setsIndices[collectionIdx].size();
    if( ( m_onlyOnSetChange != 0 ) && ( !m_setChanged ) ) // if we're only collecting when the set changes but the set hasn't changed
    {
      collectionSize = 0;
    }
    HistoryMetadata metadata = targetField.getHistoryMetadata( collectionSize );
    metadata.setName( metadata.getName() + " " + m_setNames[collectionIdx] );
    return metadata;
  }
}

/**
 * @brief Extracts the wrapper names that exist in the @p sets group of @p group.
 * @param setNames The candidate wrapper names.
 * @param omb The  @p ObjectManagerBase instance holding the @p sets from which we'll extract the @p setNames.
 * @return The set names, guarantied to be unique.
 * @note This function simply discards requested sets which are not available. No warning or error is provided...
 */
std::vector< string > getExistingWrapperNames( arrayView1d< string const > setNames, ObjectManagerBase const * omb )
{
  // Extract available wrapper names from `omb`.
  std::set< string > available;
  omb->sets().forWrappers( [&]( WrapperBase const & w ) { available.insert( w.getName() ); } );

  // Convert input `setNames` to a `std::set` to allow set intersection computation.
  std::set< string > const requested( setNames.begin(), setNames.end() );

  // Compute the intersection of requested and available.
  std::vector< string > intersection;
  std::set_intersection( requested.cbegin(), requested.cend(), available.cbegin(), available.cend(), std::back_inserter( intersection ) );

  return intersection;
}


// if we add a check for the setsizes changing we can use that data to determine whether the metadata collector
//  should only collect on size changes, otherwise not collect anything
void PackCollection::updateSetsIndices( DomainPartition const & domain )
{
  // In the current function, depending on the context, the target can be considered as a `Group` or as an `ObjectManagerBase`.
  // This convenience lambda function is a shortcut to get lighter syntax while not changing the whole function.
  auto asOMB = []( Group const * grp ) -> ObjectManagerBase const *
  {
    ObjectManagerBase const * omb = dynamicCast< ObjectManagerBase const * >( grp );
    GEOS_ERROR_IF( omb == nullptr, "Group " << grp->getName() << " could not be converted to an ObjectManagerBase during the `PackCollection` process." );
    return omb;
  };

  Group const * targetGrp = this->getTargetObject( domain, m_objectPath );
  try
  {
    WrapperBase const & targetField = targetGrp->getWrapperBase( m_fieldName );
    GEOS_ERROR_IF( !targetField.isPackable( false ), "The object targeted for collection must be packable!" );
  }
  catch( std::exception const & e )
  {
    throw InputError( e, getWrapperDataContext( viewKeysStruct::fieldNameString() ).toString() +
                      ": Target not found !\n" );
  }

  // If no set or "all" is specified we retrieve the entire field.
  // If sets are specified we retrieve the field only from those sets.
  bool const collectAll = this->collectAll();

  // In the wake of previous trick about `m_setNames`, another small trick not to compute the real set names if they are not needed.
  // This is questionable but lets me define `setNames` as `const` variable.
  // Note that the third operator will be evaluated iff `collectAll` is `false` (C++ paragraph 6.5.15).
  // So the `asOMB` function will not be called inappropriately and kill the simulation.
  std::vector< string > const setNames = collectAll ? std::vector< string >{} : getExistingWrapperNames( m_setNames.toViewConst(), asOMB( targetGrp ) );

  std::size_t const numSets = collectAll ? 1 : setNames.size();
  m_setsIndices.resize( numSets );
  // `oldSetSizes` will help us check if the sets have changed.
  std::vector< localIndex > oldSetSizes( numSets );
  for( std::size_t setIdx = 0; setIdx < numSets; ++setIdx )
  {
    oldSetSizes[setIdx] = m_setsIndices[setIdx].size();
  }

  if( collectAll )
  {
    // Here we only have one "all" set.
    array1d< localIndex > & setIndices = m_setsIndices.front();
    setIndices.resize( targetGrp->size() );
    for( localIndex i = 0; i < targetGrp->size(); ++i )
    {
      setIndices[i] = i;
    }
  }
  else
  {
    ObjectManagerBase const * targetOMB = asOMB( targetGrp );
    for( std::size_t setIdx = 0; setIdx < numSets; ++setIdx )
    {
      string const & setName = setNames[setIdx];
      array1d< localIndex > & setIndices = m_setsIndices[setIdx];

      SortedArrayView< localIndex const > const & set = targetOMB->getSet( setName );
      localIndex const setSize = set.size();
      setIndices.resize( setSize );
      if( setSize > 0 )
      {
        for( localIndex idx = 0; idx < setSize; ++idx )
        {
          setIndices[idx] = set[idx];
        }
      }
    }
  }

  // filter out the ghost indices immediately when we update the index sets
  if( m_targetIsMeshObject )
  {
    arrayView1d< integer const > const ghostRank = asOMB( targetGrp )->ghostRank();
    for( std::size_t setIdx = 0; setIdx < numSets; ++setIdx )
    {
      array1d< localIndex > & setIndices = m_setsIndices[ setIdx ];
      array1d< localIndex > ownedIndices( setIndices.size() );

      localIndex ownIdx = 0;
      for( localIndex idx = 0; idx < setIndices.size(); ++idx )
      {
        if( ghostRank[ setIndices[ idx ] ] < 0 )
        {
          ownedIndices[ ownIdx ] = setIndices[ idx ];
          ++ownIdx;
        }
      }

      localIndex const newSetSize = ownIdx;

      if( oldSetSizes[ setIdx ] != newSetSize )
      {
        m_setChanged = true;
      }
      if( !m_setChanged )  // if the size hasn't changed check the individual values
      {
        for( localIndex idx = 0; idx < ownedIndices.size(); ++idx )
        {
          if( setIndices[ idx ] != ownedIndices[ idx ] )
          {
            m_setChanged = true;
            break;
          }
        }
      }
      setIndices.resize( newSetSize );
      for( localIndex idx = 0; idx < newSetSize; ++idx )
      {
        setIndices[ idx ] = ownedIndices[ idx ];
      }
    }
  }
}

localIndex PackCollection::numMetaDataCollectors() const
{
  return m_metaDataCollectors.size();
}

void PackCollection::buildMetaDataCollectors()
{
  if( !m_disableCoordCollection )
  {
    char const * coordField = nullptr;
    if( m_objectPath.find( "nodeManager" ) != string::npos )
    {
      coordField = NodeManager::viewKeyStruct::referencePositionString();
    }
    else if( m_objectPath.find( "edgeManager" ) != string::npos )
    {
      GEOS_ERROR( "Edge coordinate data collection is unimplemented." );
    }
    else if( m_objectPath.find( "faceManager" ) != string::npos )
    {
      coordField = FaceManager::viewKeyStruct::faceCenterString();
    }
    else if( m_objectPath.find( "ElementRegions" ) != string::npos )
    {
      coordField = ElementSubRegionBase::viewKeyStruct::elementCenterString();
    }

    // "metaCollector" is a dummy name that should not appear in the results.
    string metaName( "metaCollector" );
    std::unique_ptr< PackCollection > collector = std::make_unique< PackCollection >( metaName, this );
    collector->m_objectPath = m_objectPath;
    collector->m_fieldName = coordField ? string( coordField ) : m_fieldName;
    collector->m_setNames = m_setNames;
    collector->m_onlyOnSetChange = true;
    // don't recursively keep creating metaDataCollectors
    collector->disableCoordCollection();
    m_metaDataCollectors.push_back( std::move( collector ) );
  }
}

void PackCollection::collect( DomainPartition const & domain,
                              localIndex const collectionIdx,
                              buffer_unit_type * & buffer )
{
  GEOS_MARK_FUNCTION;
  GEOS_ERROR_IF( collectionIdx < 0 || collectionIdx >= numCollectors(), "Attempting to collection from an invalid collection index!" );
  Group const * targetObject = this->getTargetObject( domain, m_objectPath );
  WrapperBase const & targetField = targetObject->getWrapperBase( m_fieldName );
  // If we have any indices to collect, and we're either collecting every time or we're only collecting
  // when the set changes and the set has changed.
  parallelDeviceEvents events;
  if( m_setsIndices[collectionIdx].size() > 0 )
  {
    if( ( ( m_onlyOnSetChange != 0 ) && m_setChanged ) || ( m_onlyOnSetChange == 0 ) )
    {
      targetField.packByIndex< true >( buffer, m_setsIndices[collectionIdx], false, true, events );
    }
  }
  // If we're not collecting from a set of indices, we're collecting the entire object.
  else if( !m_targetIsMeshObject && collectAll() )
  {
    targetField.pack< true >( buffer, false, true, events );
  }
  m_setChanged = false;
  GEOS_ASYNC_WAIT( 6000000000, 10, testAllDeviceEvents( events ) );
}

bool PackCollection::collectAll() const
{
  // The current pattern is that an empty `m_setNames` means that we have to collect all the sets.
  // But is `m_setNames` contains the "all" keyword, then we must collect all the sets too.
  // Otherwise, we rely on the content of `m_setNames`.
  return m_setNames.empty() or std::find( m_setNames.begin(), m_setNames.end(), "all" ) != m_setNames.end();
}

REGISTER_CATALOG_ENTRY( TaskBase, PackCollection, string const &, Group * const )
}
