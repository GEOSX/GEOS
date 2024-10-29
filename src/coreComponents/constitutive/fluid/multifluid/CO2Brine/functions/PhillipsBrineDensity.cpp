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
 * @file PhillipsBrineDensity.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

namespace
{

void calculateBrineDensity( PTTableCoordinates const & tableCoords,
                            real64 const & salinity,
                            array1d< real64 > const & densities )
{
  // these coefficients come from Phillips et al. (1981), equations (4) and (5), pages 14 and 15
  constexpr real64 c1 = -9.9595;
  constexpr real64 c2 = 7.0845;
  constexpr real64 c3 = 3.9093;

  constexpr real64 a1 = -0.004539;
  constexpr real64 a2 = -0.0001638;
  constexpr real64 a3 = 0.00002551;

  constexpr real64 AA = -3.033405;
  constexpr real64 BB = 10.128163;
  constexpr real64 CC = -8.750567;
  constexpr real64 DD = 2.663107;

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const P = tableCoords.getPressure( i ) / 1e5;

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const T = tableCoords.getTemperature( j );
      // see Phillips et al. (1981), equations (4) and (5), pages 14 and 15
      real64 const x = c1 * exp( a1 * salinity )
                       + c2 * exp( a2 * T )
                       + c3 * exp( a3 * P );
      densities[j*nPressures+i] = (AA + BB * x + CC * x * x + DD * x * x * x) * 1000.0;
    }
  }
}

void calculatePureWaterDensity( PTTableCoordinates const & tableCoords,
                                string const & functionName,
                                array1d< real64 > const & densities )
{
  // if no salinity, we fall back to the standard approach in three steps
  // 1- Get the saturation density as a function of temperature
  // 2- Get the saturation pressure as a function of temperature
  // 3- Get the pure water density

  TableFunction const * waterSatDensityTable =
    PureWaterProperties::makeSaturationDensityTable( functionName, FunctionManager::getInstance() );
  TableFunction const * waterSatPressureTable =
    PureWaterProperties::makeSaturationPressureTable( functionName, FunctionManager::getInstance() );

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const P = tableCoords.getPressure( i );

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const T = tableCoords.getTemperature( j );

      // Step 1: get the saturation density
      real64 const waterSatDensity = waterSatDensityTable->evaluate( &T );
      // Step 2: get the saturation pressure
      real64 const waterSatPressure = waterSatPressureTable->evaluate( &T );
      // Step 3: get the pure water density
      // Note: for now, we keep a constant water compressibility for consistency with the Ezrokhi model
      // In the future, we should query the water compressibility as a function of pressure and temperature in a table
      real64 const waterCompressibility = 4.5e-10; // Pa-1 // TODO: consolidate to a unique file as Dick started doing
      densities[j*nPressures+i] = waterSatDensity * exp( waterCompressibility * ( P - waterSatPressure ) );
    }
  }
}

TableFunction const * makeDensityTable( string_array const & inputParams,
                                        string const & functionName,
                                        FunctionManager & functionManager )
{
  string const tableName = functionName + "_table";

  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    TableFunction * const densityTable = functionManager.getGroupPointer< TableFunction >( tableName );
    densityTable->initializeFunction();
    densityTable->setDimUnits( { units::Pressure, units::TemperatureInC } );
    densityTable->setValueUnits( units::Density );
    return densityTable;
  }
  else
  {
    GEOS_THROW_IF_LT_MSG( inputParams.size(), 9,
                          GEOS_FMT( "{}: insufficient number of model parameters", functionName ),
                          InputError );

    // initialize the (p,T) coordinates
    PTTableCoordinates tableCoords;
    PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

    // initialize salinity
    real64 salinity;
    try
    {
      salinity = stod( inputParams[8] );
    }
    catch( std::invalid_argument const & e )
    {
      GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
    }

    array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
    if( !isZero( salinity ) )
    {
      // if we are in the range of validity of the Phillips method, everything is good
      // if we are not, we issue a warning message
      calculateBrineDensity( tableCoords, salinity, densities );
      GEOS_LOG_RANK_0_IF( salinity < 0.25,
                          GEOS_FMT( "{}: Warning! The salinity value of {} is below the range of validity of the Phillips model, results may be inaccurate",
                                    functionName, salinity ) );
    }
    else
    {
      // the Phillips correlation is inaccurate in the absence of salinity.
      // since this is a very frequent case, we implement an alternate (more accurate) method below
      calculatePureWaterDensity( tableCoords, functionName, densities );
    }

    TableFunction * const densityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );
    densityTable->setTableCoordinates( tableCoords.getCoords(),
                                       { units::Pressure, units::TemperatureInC } );
    densityTable->setTableValues( densities, units::Density );
    densityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return densityTable;
  }
}

} // namespace

PhillipsBrineDensity::PhillipsBrineDensity( string const & name,
                                            string_array const & inputParams,
                                            string_array const & componentNames,
                                            array1d< real64 > const & componentMolarWeight,
                                            bool const printTable ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames, "componentNames" );

  string const expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames, "componentNames" );

  m_brineDensityTable = makeDensityTable( inputParams, m_functionName, FunctionManager::getInstance() );
  if( printTable )
    m_brineDensityTable->print( m_brineDensityTable->getName() );
}

PhillipsBrineDensity::KernelWrapper
PhillipsBrineDensity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_brineDensityTable,
                        m_CO2Index,
                        m_waterIndex );
}

void PhillipsBrineDensity::checkTablesParameters( real64 const pressure,
                                                  real64 const temperature ) const
{
  m_brineDensityTable->checkCoord( pressure, 0 );
  m_brineDensityTable->checkCoord( temperature, 1 );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, PhillipsBrineDensity, string const &, string_array const &, string_array const &, array1d< real64 > const &, bool const )

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
