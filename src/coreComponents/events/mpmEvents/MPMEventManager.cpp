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

/**
 * @file MPMEventManager.cpp
 */

#include "MPMEventManager.hpp"

#include "common/TimingMacros.hpp"
#include "events/EventBase.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{

using namespace dataRepository;


MPMEventManager::MPMEventManager( string const & name,
                                  Group * const parent ):
                                  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );

  // This enables logLevel filtering
  enableLogLevelInput();
}


MPMEventManager::~MPMEventManager()
{}


Group * MPMEventManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding MPM Event: " << childKey << ", " << childName );
  std::unique_ptr< MPMEventBase > event = MPMEventBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< MPMEventBase >( childName, std::move( event ) );
}


void MPMEventManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from MPMEventBase here
  for( auto & catalogIter: MPMEventBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

} /* namespace geos */
