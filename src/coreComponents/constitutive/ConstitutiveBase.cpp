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
 * @file ConstitutiveBase.cpp
 */



#include "ConstitutiveBase.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

ConstitutiveBase::ConstitutiveBase( string const & name,
                                    Group * const parent ):
  Group( name, parent ),
  m_numQuadraturePoints( 1 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

ConstitutiveBase::CatalogInterface::CatalogType & ConstitutiveBase::getCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void ConstitutiveBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  m_numQuadraturePoints = numConstitutivePointsPerParentIndex;

  for( auto & group : this->getSubGroups() )
  {
    for( auto & wrapper : group.second->wrappers() )
    {
      if( wrapper.second->sizedFromParent() )
      {
        string const wrapperName = makeFieldName( this->getName(), wrapper.first );
        parent.registerWrapper( wrapper.second->clone( wrapperName, parent ) ).
          setRestartFlags( RestartFlags::NO_WRITE );
      }
    }
  }

  for( auto & wrapper : this->wrappers() )
  {
    if( wrapper.second->sizedFromParent() )
    {
      string const wrapperName = makeFieldName( this->getName(), wrapper.first );
      parent.registerWrapper( wrapper.second->clone( wrapperName, parent ) ).
        setRestartFlags( RestartFlags::NO_WRITE );
    }
  }

  this->resize( parent.size() );
}

std::unique_ptr< ConstitutiveBase >
ConstitutiveBase::deliverClone( string const & name,
                                Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase >
  newModel = ConstitutiveBase::CatalogInterface::factory( this->getCatalogName(), name, parent );

  newModel->forWrappers( [&]( WrapperBase & wrapper )
  {
    wrapper.copyWrapper( this->getWrapperBase( wrapper.getName() ) );
  } );

  return newModel;
}


}
} /* namespace geos */
