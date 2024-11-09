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
 * @file SinglePhaseThermalConductivityBase.cpp
 */

#include "SinglePhaseThermalConductivityBase.hpp"
#include "ThermalConductivityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SinglePhaseThermalConductivityBase::SinglePhaseThermalConductivityBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerField( fields::thermalconductivity::effectiveConductivity{}, &m_effectiveConductivity );
  registerField( fields::thermalconductivity::dEffectiveConductivity_dT{}, &m_dEffectiveConductivity_dT );
}

void SinglePhaseThermalConductivityBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  m_effectiveConductivity.resize( 0, 0, 3 );
  m_dEffectiveConductivity_dT.resize( 0, 0, 3 );
}

void SinglePhaseThermalConductivityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                                   localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_effectiveConductivity.resize( 0, 1, 3 );
  m_dEffectiveConductivity_dT.resize( 0, 1, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

} // namespace constitutive

} // namespace geos
