/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ImmiscibleWaterParameters.cpp
 */

#include "ImmiscibleWaterParameters.hpp"
#include "ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "dataRepository/InputFlags.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

ImmiscibleWaterParameters::ImmiscibleWaterParameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

std::unique_ptr< ModelParameters >
ImmiscibleWaterParameters::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< ImmiscibleWaterParameters >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< ImmiscibleWaterParameters >( std::move( parameters ) );
}

integer ImmiscibleWaterParameters::getWaterComponentIndex( ComponentProperties const & componentProperties )
{
  auto componentNames = componentProperties.getComponentName();
  integer const numComps = componentNames.size();
  for( integer ic = 0; ic < numComps; ++ic )
  {
    string const compName = stringutilities::toLower( componentNames[ic] );
    if( compName == waterComponentName )
    {
      return ic;
    }
  }
  return -1;
}

void ImmiscibleWaterParameters::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::waterReferencePressureString(), &m_waterReferencePressure ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "The reference pressure for water density and viscosity" );

  fluid->registerWrapper( viewKeyStruct::waterReferenceTemperatureString(), &m_waterReferenceTemperature ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_waterReferenceTemperature ).
    setDescription( "The reference temperature for water density and viscosity" );

  fluid->registerWrapper( viewKeyStruct::waterDensityString(), &m_waterDensity ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "The water density at the reference pressure and temperature" );

  fluid->registerWrapper( viewKeyStruct::waterViscosityString(), &m_waterViscosity ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "The water viscosity at the reference pressure and temperature" );

  fluid->registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterCompressibility ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "The compressibility of water" );

  fluid->registerWrapper( viewKeyStruct::waterViscosityCompressibilityString(), &m_waterViscosityCompressibility ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_waterViscosityCompressibility ).
    setDescription( "The compressibility (normalized derivative with respect to pressure) of the water viscosity" );

  fluid->registerWrapper( viewKeyStruct::waterExpansionCoefficientString(), &m_waterExpansionCoefficient ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_waterExpansionCoefficient ).
    setDescription( "The volumetric coefficient of thermal expansion of water" );

  fluid->registerWrapper( viewKeyStruct::waterViscosityExpansionCoefficientString(), &m_waterViscosityExpansionCoefficient ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDefaultValue( m_waterViscosityExpansionCoefficient ).
    setDescription( "The coefficient of thermal expansion (normalized derivative with respect to temperature) of water viscosity" );
}

void ImmiscibleWaterParameters::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                             ComponentProperties const & componentProperties )
{
  integer const waterIndex = fluid->getWaterPhaseIndex();
  GEOS_THROW_IF_LT_MSG( waterIndex, 0,
                        GEOS_FMT( "{}: water phase not found '{}'", fluid->getFullName(),
                                  MultiFluidBase::viewKeyStruct::phaseNamesString() ),
                        InputError );

  integer const h2oIndex = getWaterComponentIndex( componentProperties );
  GEOS_THROW_IF_LT_MSG( h2oIndex, 0,
                        GEOS_FMT( "{}: water component not found '{}'", fluid->getFullName(),
                                  MultiFluidBase::viewKeyStruct::componentNamesString() ),
                        InputError );

  // Pretty much everything should be positive
  auto const checkLowerBound = [&]( real64 const & value, real64 const & bound, string const & attribute )
  {
    GEOS_THROW_IF_LT_MSG( value, bound,
                          GEOS_FMT( "{}: invalid number of value in attribute '{}'. Should be greater than {}",
                                    fluid->getFullName(), bound, attribute ),
                          InputError );
  };

  real64 constexpr epsilon = MultiFluidConstants::epsilon;

  checkLowerBound( m_waterDensity, epsilon, viewKeyStruct::waterDensityString());
  checkLowerBound( m_waterViscosity, epsilon, viewKeyStruct::waterViscosityString());
  checkLowerBound( m_waterCompressibility, 0.0, viewKeyStruct::waterCompressibilityString());
  checkLowerBound( m_waterViscosityCompressibility, 0.0, viewKeyStruct::waterViscosityCompressibilityString());
  checkLowerBound( m_waterExpansionCoefficient, 0.0, viewKeyStruct::waterExpansionCoefficientString());
  checkLowerBound( m_waterViscosityExpansionCoefficient, 0.0, viewKeyStruct::waterViscosityExpansionCoefficientString());
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
