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
 * @file ImmiscibleWaterViscosity.cpp
 */

#include "ImmiscibleWaterViscosity.hpp"
#include "ImmiscibleWaterParameters.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

ImmiscibleWaterViscosityUpdate::ImmiscibleWaterViscosityUpdate( real64 const referencePressure,
                                                                real64 const referenceTemperature,
                                                                real64 const viscosity,
                                                                real64 const compressibility,
                                                                real64 const expansionCoefficient ):
  m_referencePressure( referencePressure ),
  m_referenceTemperature( referenceTemperature ),
  m_viscosity( viscosity ),
  m_compressibility( compressibility ),
  m_expansionCoefficient( expansionCoefficient )
{}

ImmiscibleWaterViscosity::ImmiscibleWaterViscosity( string const & name,
                                                    ComponentProperties const & componentProperties,
                                                    integer const phaseIndex,
                                                    ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties ),
  m_parameters( modelParameters )
{
  GEOS_UNUSED_VAR( phaseIndex );
}

ImmiscibleWaterViscosity::KernelWrapper
ImmiscibleWaterViscosity::createKernelWrapper() const
{
  ImmiscibleWaterParameters const * waterParameters = m_parameters.get< ImmiscibleWaterParameters >();
  return KernelWrapper( waterParameters->m_waterReferencePressure,
                        waterParameters->m_waterReferenceTemperature,
                        waterParameters->m_waterViscosity,
                        waterParameters->m_waterViscosityCompressibility,
                        waterParameters->m_waterViscosityExpansionCoefficient );
}

std::unique_ptr< ModelParameters >
ImmiscibleWaterViscosity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  return ImmiscibleWaterParameters::create( std::move( parameters ) );
}

} // namespace compositional

} // namespace constitutive

} // end namespace geos
