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
 * @file SinglePhaseConstantThermalConductivity.cpp
 */

#include "SinglePhaseConstantThermalConductivity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SinglePhaseConstantThermalConductivity::SinglePhaseConstantThermalConductivity( string const & name, Group * const parent ):
  SinglePhaseThermalConductivityBase( name, parent )
{
  registerWrapper( viewKeyStruct::thermalConductivityComponentsString(), &m_thermalConductivityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy, and zz components of a diagonal thermal conductivity tensor [J/(s.m.K)]" );
}

std::unique_ptr< ConstitutiveBase >
SinglePhaseConstantThermalConductivity::deliverClone( string const & name,
                                                      Group * const parent ) const
{
  return SinglePhaseThermalConductivityBase::deliverClone( name, parent );
}

void SinglePhaseConstantThermalConductivity::initializeRockFluidState( arrayView2d< real64 const > const & initialPorosity ) const
{
  for( localIndex ei = 0; ei < initialPorosity.size( 0 ); ++ei )
  {
    // NOTE: enforcing 1 quadrature point
    for( localIndex q = 0; q < 1; ++q )
    {
      m_effectiveConductivity[ei][q][0] = m_thermalConductivityComponents[0];
      m_effectiveConductivity[ei][q][1] = m_thermalConductivityComponents[1];
      m_effectiveConductivity[ei][q][2] = m_thermalConductivityComponents[2];
    }
  }
}

void SinglePhaseConstantThermalConductivity::update( arrayView2d< real64 const > const & initialPorosity ) const
{
  real64 thermalConductivityComponents[3];
  for( int i = 0; i<3; ++i )
  {
    thermalConductivityComponents[i] = m_thermalConductivityComponents[i];
  }
  arrayView3d< real64 > const effectiveConductivity = m_effectiveConductivity;

  forAll< parallelDevicePolicy<> >( initialPorosity.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    // NOTE: enforcing 1 quadrature point
    for( localIndex q = 0; q < 1; ++q )
    {
      effectiveConductivity[ei][q][0] = thermalConductivityComponents[0];
      effectiveConductivity[ei][q][1] = thermalConductivityComponents[1];
      effectiveConductivity[ei][q][2] = thermalConductivityComponents[2];
    }
  } );
}

void SinglePhaseConstantThermalConductivity::allocateConstitutiveData( dataRepository::Group & parent,
                                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  SinglePhaseThermalConductivityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void SinglePhaseConstantThermalConductivity::postInputInitialization()
{
  GEOS_THROW_IF( m_thermalConductivityComponents[0] <= 0 ||
                 m_thermalConductivityComponents[1] <= 0 ||
                 m_thermalConductivityComponents[2] <= 0,
                 GEOS_FMT( "{}: the components of the thermal conductivity tensor must be strictly positive",
                           getFullName() ),
                 InputError );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SinglePhaseConstantThermalConductivity, string const &, Group * const )

} // namespace constitutive

} // namespace geos
