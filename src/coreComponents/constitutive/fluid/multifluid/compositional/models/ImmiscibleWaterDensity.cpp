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
 * @file ImmiscibleWaterDensity.cpp
 */

#include "ImmiscibleWaterDensity.hpp"
#include "ImmiscibleWaterParameters.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

ImmiscibleWaterDensityUpdate::ImmiscibleWaterDensityUpdate( real64 const waterMolecularWeight,
                                                            real64 const referencePressure,
                                                            real64 const referenceTemperature,
                                                            real64 const density,
                                                            real64 const compressibility,
                                                            real64 const expansionCoefficient ):
  m_waterMolecularWeight( waterMolecularWeight ),
  m_referencePressure( referencePressure ),
  m_referenceTemperature( referenceTemperature ),
  m_density( density ),
  m_compressibility( compressibility ),
  m_expansionCoefficient( expansionCoefficient )
{}

ImmiscibleWaterDensity::ImmiscibleWaterDensity( string const & name,
                                                ComponentProperties const & componentProperties,
                                                integer const phaseIndex,
                                                ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties ),
  m_parameters( modelParameters )
{
  GEOS_UNUSED_VAR( phaseIndex );
  integer const h2oIndex = ImmiscibleWaterParameters::getWaterComponentIndex( componentProperties );
  GEOS_THROW_IF_LT_MSG( h2oIndex, 0, "Water component not found", InputError );
  m_waterMolecularWeight = componentProperties.getComponentMolarWeight()[h2oIndex];
}

ImmiscibleWaterDensity::KernelWrapper
ImmiscibleWaterDensity::createKernelWrapper() const
{
  ImmiscibleWaterParameters const * waterParameters = m_parameters.get< ImmiscibleWaterParameters >();
  return KernelWrapper( m_waterMolecularWeight,
                        waterParameters->m_waterReferencePressure,
                        waterParameters->m_waterReferenceTemperature,
                        waterParameters->m_waterDensity,
                        waterParameters->m_waterCompressibility,
                        waterParameters->m_waterExpansionCoefficient );
}

std::unique_ptr< ModelParameters >
ImmiscibleWaterDensity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  return ImmiscibleWaterParameters::create( std::move( parameters ) );
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
