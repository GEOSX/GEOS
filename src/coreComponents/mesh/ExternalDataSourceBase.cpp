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

#include "ExternalDataRepositoryBase.hpp"

namespace geos
{
using namespace dataRepository;

ExternalDataRepositoryBase::ExternalDataRepositoryBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

Group * ExternalDataRepositoryBase::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding External Data Repository: " << childKey << ", " << childName );
  std::unique_ptr< ExternalDataRepositoryBase > event = ExternalDataRepositoryBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< ExternalDataRepositoryBase >( childName, std::move( event ) );
}

void ExternalDataRepositoryBase::expandObjectCatalogs()
{
  // Only add children if the parent is of type EventManager
  // otherwise, this would fall into a loop
  if( strcmp( this->getParent().getName().c_str(), "ExternalDataRepository" ) == 0 )
  {
    for( auto & catalogIter: ExternalDataRepositoryBase::getCatalog())
    {
      createChild( catalogIter.first, catalogIter.first );
    }
  }
}

ExternalDataRepositoryBase::CatalogInterface::CatalogType & ExternalDataRepositoryBase::getCatalog()
{
  static ExternalDataRepositoryBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


}
