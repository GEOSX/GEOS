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
 * @file WillisRichardsPermeability.cpp
 */

#include "WillisRichardsPermeability.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


WillisRichardsPermeability::WillisRichardsPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::maxFracApertureString(), &m_maxFracAperture ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum fracture aperture at zero contact stress." );

  registerWrapper( viewKeyStruct::dilationCoefficientString(), &m_dilationCoefficient ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Dilation coefficient (tan of dilation angle)." );

  registerWrapper( viewKeyStruct::refClosureStressString(), &m_refClosureStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Effective normal stress causes 90% reduction in aperture." );

  registerField( fields::permeability::dPerm_dDispJump{}, &m_dPerm_dDispJump );
  registerField( fields::permeability::dPerm_dTraction{}, &m_dPerm_dTraction );
}

std::unique_ptr< ConstitutiveBase >
WillisRichardsPermeability::deliverClone( string const & name,
                                          Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void WillisRichardsPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                           localIndex const numConstitutivePointsPerParentIndex )
{
// NOTE: enforcing 1 quadrature point
  m_dPerm_dDispJump.resize( 0, 1, 3, 3 );
  m_dPerm_dTraction.resize( 0, 1, 3, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, WillisRichardsPermeability, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geos */
