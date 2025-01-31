/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ExternalDataSourceBase.hpp"

namespace geos
{
using namespace dataRepository;

ExternalDataSourceBase::ExternalDataSourceBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

Group * ExternalDataSourceBase::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< ExternalDataSourceBase > event =
    ExternalDataSourceBase::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &this->registerGroup( childName, std::move( event ) );
}

void ExternalDataSourceBase::expandObjectCatalogs()
{
  // Only add children if the parent is of type EventManager
  // otherwise, this would fall into a loop
  if( strcmp( this->getParent().getName().c_str(), "ExternalDataSource" ) == 0 )
  {
    for( auto & catalogIter: ExternalDataSourceBase::getCatalog())
    {
      createChild( catalogIter.first, catalogIter.first );
    }
  }
}

ExternalDataSourceBase::CatalogInterface::CatalogType & ExternalDataSourceBase::getCatalog()
{
  static ExternalDataSourceBase::CatalogInterface::CatalogType catalog;
  return catalog;
}


}
