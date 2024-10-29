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
 * @file ParticleFluidBase.cpp
 */

#include "ParticleFluidBase.hpp"

#include "ParticleFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ParticleFluidBase::ParticleFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).
    setApplyDefaultValue( 0.6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max proppant concentration" );

  registerWrapper( viewKeyStruct::isCollisionalSlipString(), &m_isCollisionalSlip ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether the collisional component of the slip velocity is considered" );

  registerField( fields::particlefluid::settlingFactor{}, &m_settlingFactor );
  registerField( fields::particlefluid::dSettlingFactor_dPressure{}, &m_dSettlingFactor_dPressure );
  registerField( fields::particlefluid::dSettlingFactor_dProppantConcentration{}, &m_dSettlingFactor_dProppantConcentration );
  registerField( fields::particlefluid::dSettlingFactor_dComponentConcentration{}, &m_dSettlingFactor_dComponentConcentration );

  registerField( fields::particlefluid::collisionFactor{}, &m_collisionFactor );
  registerField( fields::particlefluid::dCollisionFactor_dProppantConcentration{}, &m_dCollisionFactor_dProppantConcentration );

  registerField( fields::particlefluid::proppantPackPermeability{}, &m_proppantPackPermeability );
}

ParticleFluidBase::~ParticleFluidBase() = default;

void ParticleFluidBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();
}

void ParticleFluidBase::allocateConstitutiveData( Group & parent,
                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent.size() );
  m_dSettlingFactor_dComponentConcentration.resize( parent.size(), MAX_NUM_COMPONENTS );
}

} //namespace constitutive

} //namespace geos
