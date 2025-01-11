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

/**
 * @file ConstitutiveManager.cpp
 */

#include "ConstitutiveManager.hpp"

#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


ConstitutiveManager::ConstitutiveManager( string const & name,
                                          Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

ConstitutiveManager::~ConstitutiveManager()
{}


Group * ConstitutiveManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< ConstitutiveBase > material = ConstitutiveBase::CatalogInterface::factory( childKey, childName, this );
  return &registerGroup< ConstitutiveBase >( childName, std::move( material ) );
}


void ConstitutiveManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ConstitutiveBase here
  for( auto & catalogIter: ConstitutiveBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


void
ConstitutiveManager::hangConstitutiveRelation( string const & constitutiveRelationInstanceName,
                                               dataRepository::Group * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex ) const
{
  dataRepository::Group * constitutiveGroup = parent->getGroupPointer( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
  if( constitutiveGroup == nullptr )
  {
    constitutiveGroup = &parent->registerGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() ).setSizedFromParent( 1 );
    constitutiveGroup->resize( parent->size() );
  }

  // 1. Allocate constitutive relation
  // we only register the constitutive relation if it has not been registered yet.
  GEOS_ERROR_IF( constitutiveGroup->hasGroup( constitutiveRelationInstanceName ),
                 GEOS_FMT( "Error! The constitutive relation {} has already been registered on the subRegion {}. "
                           "Make sure that the same constitutive model is not listed as a material on a"
                           " region both as a stand-alone one and as part of a compound constitutive model.",
                           constitutiveRelationInstanceName, parent->getDataContext().toString() ) );

  ConstitutiveBase const & constitutiveRelation = getConstitutiveRelation( constitutiveRelationInstanceName );

  std::unique_ptr< ConstitutiveBase > material = constitutiveRelation.deliverClone( constitutiveRelationInstanceName, parent );

  material->allocateConstitutiveData( *parent,
                                      numConstitutivePointsPerParentIndex );
  GEOS_LOG_RANK_0( GEOS_FMT( "  {}/{} allocated",
                             constitutiveGroup->getName(),
                             constitutiveRelation.getName(),
                             numConstitutivePointsPerParentIndex ) );

  ConstitutiveBase &
  materialGroup = constitutiveGroup->registerGroup< ConstitutiveBase >( constitutiveRelationInstanceName, std::move( material ) );
  materialGroup.setSizedFromParent( 1 );
  materialGroup.resize( constitutiveGroup->size() );

  // 2. Allocate subrelations (for compound models)
  std::vector< string > const subRelationNames = constitutiveRelation.getSubRelationNames();
  for( string const & subRelationName : subRelationNames )
  {
    hangConstitutiveRelation( subRelationName, parent, numConstitutivePointsPerParentIndex );
  }
}

} /* namespace constitutive */

} /* namespace geos */
