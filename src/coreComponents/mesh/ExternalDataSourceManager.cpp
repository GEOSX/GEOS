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


#include "ExternalDataSourceManager.hpp"
#include "ExternalDataSourceBase.hpp"


namespace geos
{

using namespace dataRepository;

ExternalDataSourceManager::ExternalDataSourceManager( string const & name,
                                                      Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );
}

ExternalDataSourceManager::~ExternalDataSourceManager()
{}

Group * ExternalDataSourceManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding External Data Source: " << childKey << ", " << childName );
  std::unique_ptr< ExternalDataSourceBase > externalDataSource = ExternalDataSourceBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< ExternalDataSourceBase >( childName, std::move( externalDataSource ) );
}


void ExternalDataSourceManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ExternalDataSourceBase here
  for( auto & catalogIter: ExternalDataSourceBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


void ExternalDataSourceManager::open( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  forSubGroups< ExternalDataSourceBase >( []( ExternalDataSourceBase & external )
  {
    external.open();
  } );
}


} /* namespace geos */
