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
 * @file WaterDensity.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/WaterDensity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{


WaterDensity::WaterDensity( string const & name,
                            string_array const & inputParams,
                            string_array const & componentNames,
                            array1d< real64 > const & componentMolarWeight,
                            bool const printTable ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  GEOS_UNUSED_VAR( inputParams );
  m_waterDensityTable = PureWaterProperties::makeSaturationDensityTable( m_functionName, FunctionManager::getInstance() );
  if( printTable )
    m_waterDensityTable->print( m_waterDensityTable->getName() );
}

void WaterDensity::checkTablesParameters( real64 const pressure,
                                          real64 const temperature ) const
{
  m_waterDensityTable->checkCoord( pressure, 0 );
  m_waterDensityTable->checkCoord( temperature, 1 );
}

WaterDensity::KernelWrapper
WaterDensity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_waterDensityTable );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, WaterDensity, string const &, string_array const &, string_array const &, array1d< real64 > const &, bool const )

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
