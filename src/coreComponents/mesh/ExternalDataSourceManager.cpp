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


#include "ExternalDataRepositoryManager.hpp"
#include "ExternalDataRepositoryBase.hpp"


namespace geos
{

using namespace dataRepository;

ExternalDataRepositoryManager::ExternalDataRepositoryManager( string const & name,
                                                              Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );
}

ExternalDataRepositoryManager::~ExternalDataRepositoryManager()
{}

Group * ExternalDataRepositoryManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding External Data Repository: " << childKey << ", " << childName );
  std::unique_ptr< ExternalDataRepositoryBase > externalDataRepo = ExternalDataRepositoryBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< ExternalDataRepositoryBase >( childName, std::move( externalDataRepo ) );
}


void ExternalDataRepositoryManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ExternalDataRepositoryBase here
  for( auto & catalogIter: ExternalDataRepositoryBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


void ExternalDataRepositoryManager::open( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  forSubGroups< ExternalDataRepositoryBase >( []( ExternalDataRepositoryBase & external )
  {
    external.open();
  } );
}


} /* namespace geos */
