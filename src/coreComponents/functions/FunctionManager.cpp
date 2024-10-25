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

/**
 * @file FunctionManager.cpp
 */

#include "FunctionManager.hpp"

namespace geos
{

FunctionManager * FunctionManager::m_instance = nullptr;

using namespace dataRepository;


FunctionManager::FunctionManager( const string & name,
                                  Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  GEOS_ERROR_IF( m_instance != nullptr, "Only one FunctionManager can exist at a time." );
  m_instance = this;
}

FunctionManager::~FunctionManager()
{
  GEOS_ERROR_IF( m_instance != this, "m_instance != this should not be possible." );
  m_instance = nullptr;
}

FunctionManager & FunctionManager::getInstance()
{
  GEOS_ERROR_IF( m_instance == nullptr,
                 "FunctionManager has not been constructed, or is already been destructed." );
  return *m_instance;
}

Group * FunctionManager::createChild( string const & functionCatalogKey,
                                      string const & functionName )
{
  GEOS_LOG_RANK_0( "   " << functionCatalogKey << ": " << functionName );
  std::unique_ptr< FunctionBase > function = FunctionBase::CatalogInterface::factory( functionCatalogKey, getDataContext(),
                                                                                      functionName, this );
  return &this->registerGroup< FunctionBase >( functionName, std::move( function ) );
}


void FunctionManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FunctionBase here
  for( auto & catalogIter: FunctionBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

} // end of namespace geos
