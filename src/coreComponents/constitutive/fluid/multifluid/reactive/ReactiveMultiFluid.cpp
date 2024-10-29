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
 * @file ReactiveMultiFluid.cpp
 */
#include "ReactiveMultiFluid.hpp"
#include "ReactiveMultiFluidFields.hpp"


namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ReactiveMultiFluid::
  ReactiveMultiFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  // For now this is being hardcoded. We will see where this should come from.
  m_numPrimarySpecies = 7;
  m_numSecondarySpecies = 11;
  m_numKineticReactions = 2;

  registerField( fields::reactivefluid::primarySpeciesConcentration{}, &m_primarySpeciesConcentration );
  registerField( fields::reactivefluid::secondarySpeciesConcentration{}, &m_secondarySpeciesConcentration );
  registerField( fields::reactivefluid::primarySpeciesTotalConcentration{}, &m_primarySpeciesTotalConcentration );
  registerField( fields::reactivefluid::kineticReactionRates{}, &m_kineticReactionRates );
}

bool ReactiveMultiFluid::isThermal() const
{
  return true;
}

std::unique_ptr< ConstitutiveBase > ReactiveMultiFluid::
  deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  ReactiveMultiFluid & newConstitutiveRelation = dynamicCast< ReactiveMultiFluid & >( *clone );

  newConstitutiveRelation.createChemicalReactions();

  return clone;
}

void ReactiveMultiFluid::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  GEOS_THROW_IF_NE_MSG( numFluidPhases(), 1,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );

  createChemicalReactions();
}

void ReactiveMultiFluid::resizeFields( localIndex const size, localIndex const numPts )
{
  MultiFluidBase::resizeFields( size, numPts );

  integer const numPrimarySpecies = this->numPrimarySpecies();
  integer const numSecondarySpecies = this->numSecondarySpecies();
  integer const numKineticReactions = this->numKineticReactions();

  m_primarySpeciesConcentration.resize( size, numPrimarySpecies );
  m_secondarySpeciesConcentration.resize( size, numSecondarySpecies );
  m_primarySpeciesTotalConcentration.resize( size, numPrimarySpecies );
  m_kineticReactionRates.resize( size, numKineticReactions );
}

void ReactiveMultiFluid::createChemicalReactions()
{
  // instantiate reactions objects
  m_equilibriumReactions = std::make_unique< chemicalReactions::EquilibriumReactions >( getName() + "_equilibriumReactions", m_numPrimarySpecies, m_numSecondarySpecies );
  m_kineticReactions = std::make_unique< chemicalReactions::KineticReactions >( getName() + "_kineticReactions", m_numPrimarySpecies, m_numSecondarySpecies, m_numKineticReactions );
}

} //namespace constitutive

} //namespace geos
