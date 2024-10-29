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

#include "HistoryCollectionBase.hpp"

namespace geos
{

void HistoryCollectionBase::initializePostSubGroups()
{
  m_bufferProviders.resize( m_collectionCount );
}

localIndex HistoryCollectionBase::numCollectors() const
{
  return m_collectionCount;
}

void HistoryCollectionBase::registerBufferProvider( localIndex collectionIdx,
                                                    BufferProvider bufferProvider )
{
  GEOS_ERROR_IF( collectionIdx < 0 || collectionIdx >= this->numCollectors(), "Invalid collection index specified." );
  m_bufferProviders[collectionIdx] = bufferProvider;
}

HistoryCollection & HistoryCollectionBase::getMetaDataCollector( localIndex metaIdx )
{
  GEOS_ASSERT_MSG( metaIdx >= 0 && metaIdx < numMetaDataCollectors(), "Requesting nonexistent meta collector index." );
  return *m_metaDataCollectors[ metaIdx ].get( );
}

HistoryMetadata HistoryCollectionBase::getTimeMetaData() const
{
  return HistoryMetadata( this->getTargetName() + " Time", 1, std::type_index( typeid( real64 ) ) );
}

void HistoryCollectionBase::registerTimeBufferProvider( TimeBufferProvider timeBufferProvider )
{
  m_timeBufferProvider = timeBufferProvider;
}

bool HistoryCollectionBase::execute( real64 const time_n,
                                     real64 const dt,
                                     integer const cycleNumber,
                                     integer const eventCounter,
                                     real64 const eventProgress,
                                     DomainPartition & domain )
{
  for( localIndex collectionIdx = 0; collectionIdx < numCollectors(); ++collectionIdx )
  {
    // std::function defines the == and =! comparable against nullptr_t to check the
    // function pointer is actually assigned (an error would be thrown on the call attempt even so)
    GEOS_ERROR_IF( m_bufferProviders[collectionIdx] == nullptr,
                   "History collection buffer retrieval function is unassigned, did you declare a related TimeHistoryOutput event?" );
    // using GEOS_ERROR_IF_EQ caused type issues since the values are used in streams
    this->updateSetsIndices( domain );
    HistoryMetadata hmd = this->getMetaData( domain, collectionIdx );
    buffer_unit_type * buffer = m_bufferProviders[collectionIdx]( hmd.size( 0 ) );
    collect( domain, collectionIdx, buffer );
  }
  for( auto & metaCollector: m_metaDataCollectors )
  {
    metaCollector->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  int rank = MpiWrapper::commRank();
  if( rank == 0 && m_timeBufferProvider )
  {
    buffer_unit_type * timeBuffer = m_timeBufferProvider();
    memcpy( timeBuffer, &time_n, sizeof( time_n ) );
  }
  return false;
}

dataRepository::Group const * HistoryCollectionBase::getTargetObject( DomainPartition const & domain, string const & objectPath ) const
{
  try
  {
    // absolute objectPaths can be used to target anything in the data repo that is packable
    if( objectPath[0] == '/' )
    {
      return &Group::getGroupByPath( objectPath );
    }
    else // relative objectPaths use relative lookup identical to fieldSpecification to make xml input spec easier
    {
      std::vector< string > targetTokens = stringutilities::tokenize( objectPath, "/" );
      localIndex targetTokenLength = LvArray::integerConversion< localIndex >( targetTokens.size() );

      dataRepository::Group const * targetGroup = nullptr;
      int const numMeshBodies = domain.getMeshBodies().numSubGroups();


      if( numMeshBodies==1 )
      {
        string const singleMeshBodyName = domain.getMeshBody( 0 ).getName();
        if( targetTokens[0] != singleMeshBodyName )
        {
          ++targetTokenLength;
          targetTokens.insert( targetTokens.begin(), singleMeshBodyName );
        }
      }
      else
      {
        bool bodyFound = false;
        domain.forMeshBodies( [&]( MeshBody const & meshBody )
        {
          if( meshBody.getName()==targetTokens[0] )
          {
            bodyFound=true;
          }
        } );

        GEOS_THROW_IF( !bodyFound,
                       GEOS_FMT( "MeshBody ({}) is specified, but not found.",
                                 targetTokens[0] ),
                       std::domain_error );
      }

      string const meshBodyName = targetTokens[0];
      MeshBody const & meshBody = domain.getMeshBody( meshBodyName );

      // set mesh level in path
      localIndex const numMeshLevels = meshBody.getMeshLevels().numSubGroups();
      if( numMeshLevels==1 )
      {
        string const singleMeshLevelName = meshBody.getMeshLevels().getGroup< MeshLevel >( 0 ).getName();
        if( targetTokens[1] != singleMeshLevelName )
        {
          ++targetTokenLength;
          targetTokens.insert( targetTokens.begin()+1, singleMeshLevelName );
        }
        else
        {
          bool levelFound = false;
          meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
          {
            if( meshLevel.getName()==targetTokens[1] )
            {
              levelFound=true;
            }
          } );

          GEOS_THROW_IF( !levelFound,
                         GEOS_FMT( "MeshLevel ({}) is specified, but not found.",
                                   targetTokens[1] ),
                         std::domain_error );
        }
      }
      else if( !meshBody.getMeshLevels().hasGroup< MeshLevel >( targetTokens[1] ) )
      {
        string const baseMeshLevelName = MeshBody::groupStructKeys::baseDiscretizationString();
        ++targetTokenLength;
        targetTokens.insert( targetTokens.begin()+1, baseMeshLevelName );
      }

      string meshLevelName = targetTokens[1];
      MeshLevel const & meshLevel = meshBody.getMeshLevel( meshLevelName );
      targetGroup = &meshLevel;


      if( targetTokens[2]== MeshLevel::groupStructKeys::elemManagerString() )
      {
        ElementRegionManager const & elemRegionManager = meshLevel.getElemManager();
        string const elemRegionName = targetTokens[3];
        ElementRegionBase const & elemRegion = elemRegionManager.getRegion( elemRegionName );
        string const elemSubRegionName = targetTokens[4];
        ElementSubRegionBase const & elemSubRegion = elemRegion.getSubRegion( elemSubRegionName );
        targetGroup = &elemSubRegion;
      }
      else
      {
        for( localIndex pathLevel = 2; pathLevel < targetTokenLength; ++pathLevel )
        {
          dataRepository::Group const * const childGroup = targetGroup->getGroupPointer( targetTokens[pathLevel] );
          if( childGroup != nullptr )
          {
            targetGroup=childGroup;
          }
          else
          {
            string const targetTokensStr = stringutilities::join( targetTokens.begin(),
                                                                  targetTokens.begin()+pathLevel,
                                                                  '/' );
            GEOS_THROW( targetTokens[pathLevel] << " not found in path " <<
                        objectPath << std::endl << targetGroup->dumpSubGroupsNames(),
                        std::domain_error );
          }
        }
      }
      return targetGroup;
    }
  }
  catch( std::exception const & e )
  {
    throw InputError( e, getDataContext().toString() + " has a wrong objectPath: " + objectPath + "\n" );
  }
}

}
