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
 * @file MPMEventBase.cpp
 */

#include "MPMEventBase.hpp"

namespace geos
{
  
  using namespace dataRepository;

  MPMEventBase::MPMEventBase( string const & name,
                              Group * const parent  ) :
                              Group( name, parent),
                              m_time( DBL_MAX ),
                              m_interval( DBL_MAX ),
                              m_isComplete( 0 )
  {
    registerWrapper( viewKeyStruct::timeString(), &m_time ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Time at which event starts" );

    registerWrapper( viewKeyStruct::intervalString(), &m_interval ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Time interval over which event is performed" );

    registerWrapper( viewKeyStruct::isCompleteString(), &m_isComplete ).
      setInputFlag( InputFlags::FALSE ).
      setDescription( "Flag for whether event has been completed" );
  }


  MPMEventBase::~MPMEventBase()
  {}


  MPMEventBase::CatalogInterface::CatalogType & MPMEventBase::getCatalog()
  {
    static MPMEventBase::CatalogInterface::CatalogType catalog;
    return catalog;
  }

  Group * MPMEventBase::createChild( string const & childKey, string const & childName )
  {
    GEOS_LOG_RANK_0( "Adding MPM Event: " << childKey << ", " << childName );
    std::unique_ptr< MPMEventBase > event = MPMEventBase::CatalogInterface::factory( childKey, childName, this );
    return &this->registerGroup< MPMEventBase >( childName, std::move( event ) );
  }


  void MPMEventBase::expandObjectCatalogs()
  {
    // Only add children if the parent is of type EventManager
    // otherwise, this would fall into a loop
    if( strcmp( this->getParent().getName().c_str(), "MPMEvents" ) == 0 )
    {
      for( auto & catalogIter: MPMEventBase::getCatalog() )
      {
        createChild( catalogIter.first, catalogIter.first );
      }
    }
  }


} /* namespace geos */
