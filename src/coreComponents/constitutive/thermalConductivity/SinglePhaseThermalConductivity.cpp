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
 * @file SinglePhaseThermalConductivity.cpp
 */

#include "SinglePhaseThermalConductivity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SinglePhaseThermalConductivity::SinglePhaseThermalConductivity( string const & name, Group * const parent ):
  SinglePhaseThermalConductivityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultThermalConductivityComponentsString(), &m_defaultThermalConductivityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy, and zz diagonal components of the default thermal conductivity tensor [J/(s.m.K)]" );

  registerWrapper( viewKeyStruct::thermalConductivityGradientComponentsString(), &m_thermalConductivityGradientComponents ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {0.0, 0.0, 0.0} ).
    setDescription( "xx, yy, and zz diagonal components of the thermal conductivity gradient tensor w.r.t. temperature [J/(s.m.K^2)]" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "The reference temperature at which the conductivity components are equal to the default values" );
}

std::unique_ptr< ConstitutiveBase >
SinglePhaseThermalConductivity::deliverClone( string const & name,
                                              Group * const parent ) const
{
  return SinglePhaseThermalConductivityBase::deliverClone( name, parent );
}

void SinglePhaseThermalConductivity::initializeRockFluidState( arrayView2d< real64 const > const & initialPorosity ) const
{
  arrayView3d< real64 > dEffectiveConductivity_dT = m_dEffectiveConductivity_dT.toView();
  arrayView3d< real64 > effectiveConductivity = m_effectiveConductivity.toView();
  R1Tensor const defaultThermalConductivityComponents = m_defaultThermalConductivityComponents;
  R1Tensor const thermalConductivityGradientComponents = m_thermalConductivityGradientComponents;

  forAll< parallelDevicePolicy<> >( initialPorosity.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    // NOTE: enforcing 1 quadrature point
    for( localIndex q = 0; q < 1; ++q )
    {
      effectiveConductivity[ei][q][0] = defaultThermalConductivityComponents[0];
      effectiveConductivity[ei][q][1] = defaultThermalConductivityComponents[1];
      effectiveConductivity[ei][q][2] = defaultThermalConductivityComponents[2];

      dEffectiveConductivity_dT[ei][q][0] = thermalConductivityGradientComponents[0];
      dEffectiveConductivity_dT[ei][q][1] = thermalConductivityGradientComponents[1];
      dEffectiveConductivity_dT[ei][q][2] = thermalConductivityGradientComponents[2];
    }
  } );
}

void SinglePhaseThermalConductivity::updateFromTemperature( arrayView1d< real64 const > const & temperature ) const
{
  arrayView3d< real64 > dEffectiveConductivity_dT = m_dEffectiveConductivity_dT.toView();
  arrayView3d< real64 > effectiveConductivity = m_effectiveConductivity.toView();
  R1Tensor const defaultThermalConductivityComponents = m_defaultThermalConductivityComponents;
  R1Tensor const thermalConductivityGradientComponents = m_thermalConductivityGradientComponents;
  real64 const referenceTemperature = m_referenceTemperature;

  forAll< parallelDevicePolicy<> >( temperature.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( localIndex q = 0; q < 1; ++q )
    {

      real64 const deltaTemperature = temperature[ei] - referenceTemperature;

      effectiveConductivity[ei][q][0] = defaultThermalConductivityComponents[0] + thermalConductivityGradientComponents[0] * deltaTemperature;
      effectiveConductivity[ei][q][1] = defaultThermalConductivityComponents[1] + thermalConductivityGradientComponents[1] * deltaTemperature;
      effectiveConductivity[ei][q][2] = defaultThermalConductivityComponents[2] + thermalConductivityGradientComponents[2] * deltaTemperature;

      for( localIndex i=0; i<=2; i++ )
      {
        if( effectiveConductivity[ei][q][i] <1e-2 )
        {
          effectiveConductivity[ei][q][i] = 1e-2; // W/m/K To avoid negative conductivity
        }
      }

      dEffectiveConductivity_dT[ei][q][0] = thermalConductivityGradientComponents[0];
      dEffectiveConductivity_dT[ei][q][1] = thermalConductivityGradientComponents[1];
      dEffectiveConductivity_dT[ei][q][2] = thermalConductivityGradientComponents[2];
    }
  } );
}

void SinglePhaseThermalConductivity::allocateConstitutiveData( dataRepository::Group & parent,
                                                               localIndex const numConstitutivePointsPerParentIndex )
{
  SinglePhaseThermalConductivityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void SinglePhaseThermalConductivity::postInputInitialization()
{
  GEOS_THROW_IF( m_defaultThermalConductivityComponents[0] <= 0 ||
                 m_defaultThermalConductivityComponents[1] <= 0 ||
                 m_defaultThermalConductivityComponents[2] <= 0,
                 GEOS_FMT( "{}: the components of the default thermal conductivity tensor must be strictly positive",
                           getFullName() ),
                 InputError );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SinglePhaseThermalConductivity, string const &, Group * const )

} // namespace constitutive

} // namespace geos
