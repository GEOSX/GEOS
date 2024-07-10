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
 * @file LohrenzBrayClarkViscosity.cpp
 */

#include "LohrenzBrayClarkViscosity.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

LohrenzBrayClarkViscosityUpdate::LohrenzBrayClarkViscosityUpdate( MixingType const mixing_type,
                                                                  arrayView1d< real64 const > const & componentCriticalVolume )
  : m_mixing_type( mixing_type ),
  m_componentCriticalVolume( componentCriticalVolume )
{}

LohrenzBrayClarkViscosity::LohrenzBrayClarkViscosity( string const & name,
                                                      ComponentProperties const & componentProperties,
                                                      integer const phaseIndex,
                                                      ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties ),
  m_parameters( modelParameters.get< Parameters >() )
{
  GEOS_UNUSED_VAR( phaseIndex );
}

LohrenzBrayClarkViscosity::KernelWrapper
LohrenzBrayClarkViscosity::createKernelWrapper() const
{
  auto const mixingType = EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::fromString( m_parameters->m_componentMixingType );
  return KernelWrapper( mixingType, m_parameters->m_componentCriticalVolume );
}

std::unique_ptr< ModelParameters >
LohrenzBrayClarkViscosity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< Parameters >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< Parameters >( std::move( parameters ) );
}

LohrenzBrayClarkViscosity::Parameters::Parameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{
  constexpr LohrenzBrayClarkViscosityUpdate::MixingType defaultMixing = LohrenzBrayClarkViscosityUpdate::MixingType::HERNING_ZIPPERER;
  m_componentMixingType = EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::toString( defaultMixing );
}

void LohrenzBrayClarkViscosity::Parameters::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::componentCriticalVolumeString(), &m_componentCriticalVolume ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Component critical volumes" );

  fluid->registerWrapper( viewKeyStruct::componentMixingTypeString(), &m_componentMixingType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_componentMixingType ).
    setDescription( "Viscosity mixing rule to be used for Lohrenz-Bray-Clark computation. Valid options:\n* " +
                    EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::concat( "\n* " ) );
}

void LohrenzBrayClarkViscosity::Parameters::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                                         ComponentProperties const & componentProperties )
{
  integer const numComponents = fluid->numFluidComponents();

  if( m_componentCriticalVolume.empty() )
  {
    m_componentCriticalVolume.resize( numComponents );

    arrayView1d< real64 > const & componentCriticalPressure = componentProperties.getComponentCriticalPressure();
    arrayView1d< real64 > const & componentCriticalTemperature = componentProperties.getComponentCriticalTemperature();

    calculateCriticalVolume( numComponents,
                             componentCriticalPressure,
                             componentCriticalTemperature,
                             m_componentCriticalVolume );
  }

  GEOS_THROW_IF_NE_MSG( m_componentCriticalVolume.size(), numComponents,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", fluid->getFullName(),
                                  viewKeyStruct::componentCriticalVolumeString() ),
                        InputError );

  // If the value is invalid, this will throw
  EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::fromString( m_componentMixingType );
}

void LohrenzBrayClarkViscosity::Parameters::calculateCriticalVolume(
  integer const numComponents,
  arrayView1d< const real64 > const criticalPressure,
  arrayView1d< const real64 > const criticalTemperature,
  arrayView1d< real64 > const criticalVolume )
{
  for( integer ic=0; ic<numComponents; ++ic )
  {
    criticalVolume[ic] = 2.215e-6 * criticalTemperature[ic] / (0.025 + 1e-6*criticalPressure[ic] );   // m^3/mol
  }
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
