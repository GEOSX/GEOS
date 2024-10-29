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
 * @file FiniteElementSpaceManager.cpp
 */

#include "FiniteElementDiscretizationManager.hpp"

#include "FiniteElementDiscretization.hpp"

namespace geos
{
using namespace dataRepository;

FiniteElementDiscretizationManager::FiniteElementDiscretizationManager( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FiniteElementDiscretizationManager::~FiniteElementDiscretizationManager()
{
  // TODO Auto-generated destructor stub
}


Group * FiniteElementDiscretizationManager::createChild( string const & childKey, string const & childName )
{
  // These objects should probably not be registered on managed group...
  std::unique_ptr< Group > fem = Group::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup( childName, std::move( fem ) );
}


void FiniteElementDiscretizationManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from Group here
  for( auto & catalogIter: Group::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geos */
