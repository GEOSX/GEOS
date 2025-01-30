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
 * @file ComponentProperties.cpp
 */

#include "ComponentProperties.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

void ComponentProperties::classifyComponents( string_array const & componentNames, array1d< integer > & componentType )
{
  integer const numComps = componentNames.size();
  componentType.resize( numComps );
  std::unordered_map< string, ComponentType > const nameDict =
  {
    { "n2", ComponentType::Nitrogen },
    { "h2o", ComponentType::Water },
    { "h2s", ComponentType::HydrogenSulphide },
    { "co2", ComponentType::CarbonDioxide }
  };

  for( integer ic = 0; ic < numComps; ++ic )
  {
    string const name = stringutilities::toLower( componentNames[ic] );
    auto const it = nameDict.find( name );
    if( it == nameDict.end())
    {
      // Default is hydrocarbon
      componentType[ic] = static_cast< integer >(ComponentType::HydroCarbon);
    }
    else
    {
      componentType[ic] = static_cast< integer >(it->second);
    }
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
