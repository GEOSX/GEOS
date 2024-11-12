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
 * @file CriticalVolume.cpp
 */

#include "CriticalVolume.hpp"
#include "ComponentProperties.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

CriticalVolume::CriticalVolume( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

std::unique_ptr< ModelParameters > CriticalVolume::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< CriticalVolume >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< CriticalVolume >( std::move( parameters ) );
}

void CriticalVolume::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::componentCriticalVolumeString(), &m_componentCriticalVolume ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Component critical volumes" );
}

void CriticalVolume::postInputInitializationImpl( MultiFluidBase const * fluid,
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
}

void CriticalVolume::calculateCriticalVolume( integer const numComponents,
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

} // end namespace constitutive

} // end namespace geos
