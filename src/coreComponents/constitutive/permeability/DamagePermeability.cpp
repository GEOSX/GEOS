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
    setDescription( "Permeability of an undamaged bulk rock." );
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
