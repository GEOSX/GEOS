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
 * @file StatOutputController.cpp
 */

#include "StatOutputController.hpp"
#include "events/EventManager.hpp"
#include "events/EventBase.hpp"
#include "fileIO/Outputs/TimeHistoryOutput.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

void StatOutputController::generatePackCollection( TasksManager & taskManager,
                                                   string_view regionName,
                                                   string_view path,
                                                   string_view fieldName,
                                                   string_array & packCollectionPaths )
{
  string const taskManagerKey = GEOS_FMT( "packCollection{}{}", regionName, fieldName );
  PackCollection * packCollection = &taskManager.registerGroup< PackCollection >( taskManagerKey );
  string & pcObjectPath = packCollection->getReference< string >( PackCollection::viewKeysStruct::objectPathString());
  pcObjectPath = path;
  string & pcName = packCollection->getReference< string >( PackCollection::viewKeysStruct::fieldNameString());
  pcName = fieldName;

  m_packCollections.push_back( packCollection );
  packCollectionPaths.emplace_back( GEOS_FMT( "{}", packCollection->getPath()) );
}

void StatOutputController::generateTimeHistory( OutputManager & outputManager,
                                                string_view regionName,
                                                string_array const & packCollectionPaths )
{
  string const outputManagerKey = GEOS_FMT( "compFlowHistory_{}", regionName );
  //we expect timeHistory will add the extension by itself
  makeDirsForPath( m_outputDir );
  string const filename = GEOS_FMT( "{}/{}", m_outputDir, outputManagerKey );

  TimeHistoryOutput * timeHistory = &outputManager.registerGroup< TimeHistoryOutput >( outputManagerKey );
  string_array & collectorPaths =  timeHistory->getReference< string_array >( TimeHistoryOutput::viewKeys::timeHistoryOutputTargetString() );
  collectorPaths = packCollectionPaths;
  string & outputFile =  timeHistory->getReference< string >( TimeHistoryOutput::viewKeys::timeHistoryOutputFilenameString() );
  outputFile = filename;

  m_timeHistories.push_back( timeHistory );
}

StatOutputController::StatOutputController( const string & name,
                                            Group * const parent ):
  TaskBase( name, parent ),
  m_statistics( nullptr )
{
  registerWrapper( viewKeyStruct::outputDirString(), &m_outputDir ).
    setApplyDefaultValue( "." ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Output directory for statistics HDF regions file" );
}

static std::set< string > const allowedChildTypes = { SinglePhaseStatistics::catalogName(),
                                                      CompositionalMultiphaseStatistics::catalogName() };

void StatOutputController::initializePreSubGroups()
{

  Group & problemManager = this->getGroupByPath( "/Problem" );
  DomainPartition & domain = problemManager.getGroup< DomainPartition >( "domain" );
  TasksManager & taskManager = this->getGroupByPath< TasksManager >( "/Tasks" );
  OutputManager & outputManager = this->getGroupByPath< OutputManager >( "/Outputs" );

  Group & meshBodies = domain.getMeshBodies();
  std::vector< string > const groupNames = this->getSubGroupsNames();

  GEOS_ERROR_IF( groupNames.empty(),
                 GEOS_FMT( "StatOutputController must have one of the following components : {}",
                           stringutilities::join( allowedChildTypes, "," ) ) );

  m_statistics = &this->getGroup< TaskBase >( groupNames[0] );

  // creating all needed components to read sub-groups statistics.
  forSubStats( [&]( auto & statistics ) {
    using STATSTYPE = typename TYPEOFREF( statistics );

    statistics.getSolver()->forDiscretizationOnMeshTargets( meshBodies,
                                                            [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      for( string const & regionName : regionNames )
      {
        ElementRegionBase & region = elemManager.getRegion( regionName );
        string const regionStatPath = GEOS_FMT( "{}/regionStatistics", region.getPath() );
        typename STATSTYPE::RegionStatistics & regionStats =
          this->getGroupByPath< typename STATSTYPE::RegionStatistics >( regionStatPath );

        string_array packCollectionPaths;
        regionStats.forWrappers( [&]( WrapperBase const & wrapper )
        {
          generatePackCollection( taskManager,
                                  regionName,
                                  regionStatPath,
                                  wrapper.getName(),
                                  packCollectionPaths );
        } );

        generateTimeHistory( outputManager, regionName, packCollectionPaths );
      }
    } );
  } );

}

Group * StatOutputController::createChild( string const & childKey, string const & childName )
{
  GEOS_ERROR_IF( allowedChildTypes.count( childKey )==0,
                 GEOS_FMT( "Wrong children \"{}\" in {}, allowed child types are {}",
                           childKey, getDataContext(), stringutilities::join( allowedChildTypes, "," )) );
  GEOS_LOG_RANK_0( "Adding Statistics: " << childKey << ", " << childName );

  std::unique_ptr< TaskBase > task = TaskBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< TaskBase >( childName, std::move( task ) );
}

template< typename LAMBDA >
void StatOutputController::forSubStats( LAMBDA lambda )
{
  forSubGroups< SinglePhaseStatistics,
                CompositionalMultiphaseStatistics >( lambda );
}

void StatOutputController::expandObjectCatalogs()
{
  for( auto const & key: allowedChildTypes )
  {
    createChild( key, key );
  }
}

bool StatOutputController::execute( real64 const time_n,
                                    real64 const dt,
                                    integer const cycleNumber,
                                    integer const eventCounter,
                                    real64 const eventProgress,
                                    DomainPartition & domain )
{
  m_statistics->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  for( PackCollection * packCollection : m_packCollections )
  {
    packCollection->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  for( TimeHistoryOutput * timeHistory : m_timeHistories )
  {
    timeHistory->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        StatOutputController,
                        string const &, dataRepository::Group * const )

}
