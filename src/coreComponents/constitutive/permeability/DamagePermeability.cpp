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
 * @file DamagePermeability.cpp
 */

#include "DamagePermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


DamagePermeability::DamagePermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::bulkPermeabilityString(), &m_bulkPermeability ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Permeability of the intact bulk material" );

  registerWrapper( viewKeyStruct::damageDependenceConstantString(), &m_damageDependenceConstant ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Damage dependeny coefficient" );
}

std::unique_ptr< ConstitutiveBase >
DamagePermeability::deliverClone( string const & name,
                                  Group * const parent ) const
{
  return PermeabilityBase::deliverClone( name, parent );
}

void DamagePermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                   localIndex const numConstitutivePointsPerParentIndex )
{
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void DamagePermeability::postInputInitialization()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamagePermeability, string const &, Group * const )

}
} /* namespace geos */
