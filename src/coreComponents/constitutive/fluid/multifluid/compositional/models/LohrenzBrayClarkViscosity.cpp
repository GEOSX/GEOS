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
 * @file LohrenzBrayClarkViscosity.cpp
 */

#include "LohrenzBrayClarkViscosity.hpp"
#include "CriticalVolume.hpp"
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
  m_parameters( modelParameters )
{
  GEOS_UNUSED_VAR( phaseIndex );
}

LohrenzBrayClarkViscosity::KernelWrapper
LohrenzBrayClarkViscosity::createKernelWrapper() const
{
  Parameters const * parameters = m_parameters.get< Parameters >();
  CriticalVolume const * criticalVolume = m_parameters.get< CriticalVolume >();
  auto const mixingType = EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::fromString( parameters->m_componentMixingType );
  return KernelWrapper( mixingType, criticalVolume->m_componentCriticalVolume );
}

std::unique_ptr< ModelParameters >
LohrenzBrayClarkViscosity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< Parameters >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< Parameters >( CriticalVolume::create( std::move( parameters )) );
}

LohrenzBrayClarkViscosity::Parameters::Parameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{
  constexpr LohrenzBrayClarkViscosityUpdate::MixingType defaultMixing = LohrenzBrayClarkViscosityUpdate::MixingType::HERNING_ZIPPERER;
  m_componentMixingType = EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::toString( defaultMixing );
}

void LohrenzBrayClarkViscosity::Parameters::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::componentMixingTypeString(), &m_componentMixingType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_componentMixingType ).
    setDescription( "Viscosity mixing rule to be used for Lohrenz-Bray-Clark computation. Valid options:\n* " +
                    EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::concat( "\n* " ) );
}

void LohrenzBrayClarkViscosity::Parameters::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                                         ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( fluid, componentProperties );
  // If the value is invalid, this will throw
  EnumStrings< LohrenzBrayClarkViscosityUpdate::MixingType >::fromString( m_componentMixingType );
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
