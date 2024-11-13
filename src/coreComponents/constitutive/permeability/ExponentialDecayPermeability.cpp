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
 * @file ExponentialDecayPermeability.cpp
 */

#include "ExponentialDecayPermeability.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


ExponentialDecayPermeability::ExponentialDecayPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::empiricalConstantString(), &m_empiricalConstant ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "an empirical constant." );

  registerWrapper( viewKeyStruct::initialPermeabilityString(), &m_initialPermeability ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( " initial permeability of the fracture." );

  registerField( fields::permeability::dPerm_dTraction{}, &m_dPerm_dTraction );
  registerField( fields::permeability::dPerm_dDispJump{}, &m_dPerm_dDispJump );
}

std::unique_ptr< ConstitutiveBase >
ExponentialDecayPermeability::deliverClone( string const & name,
                                            Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void ExponentialDecayPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
// NOTE: enforcing 1 quadrature point
  m_dPerm_dTraction.resize( 0, 1, 3, 3 );
  m_dPerm_dDispJump.resize( 0, 1, 3, 3 );

  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ExponentialDecayPermeability, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geos */
