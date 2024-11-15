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
  GEOS_UNUSED_VAR( fluid );
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
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
