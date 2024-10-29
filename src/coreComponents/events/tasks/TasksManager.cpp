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
 * @file TasksManager.cpp
 */

#include "TasksManager.hpp"
#include "TaskBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace LvArray;

TasksManager::TasksManager( string const & name,
                            Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

TasksManager::~TasksManager()
{ }

Group * TasksManager::createChild( string const & childKey, string const & childName )
{
  std::unique_ptr< TaskBase > task = TaskBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< TaskBase >( childName, std::move( task ) );
}

void TasksManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from TaskBase here
  for( auto & catalogIter: TaskBase::getCatalog() )
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

} /* namespace geos */
