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
 * @file BrineEnthalpy.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/BrineEnthalpy.hpp"

#include "functions/FunctionManager.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Enthalpy.hpp"


namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

namespace
{

real64 michaelidesBrineEnthalpy( real64 const & T,
                                 real64 const & m )
{

  static real64 a[4][3] = {
    {-9633.6, -4080.0, 286.49},
    {166.58, 68.577, -4.6856},
    {-0.90963, -0.36524, 0.0249667},
    {1.7965e-3, 7.1924e-4, -4.9e-5}
  };

  real64 x1, x2, h1, h2, dh;

  x1 = 1000.0 / (1000.0 + 58.44 * m);
  x2 = 58.44 * m  / (1000.0 + 58.44 * m);

  dh = 0.0;

  for( localIndex i = 0; i < 4; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      dh += a[i][j] * pow( T, real64( i )) * pow( m, real64( j ));
    }
  }

  dh *= 4.184 / (1000.0 + 58.44 * m);

  h1 = 0.12453e-4 * pow( T, 3.0 ) - 0.45137e-2 * pow( T, 2.0 ) + 4.81155 * T - 29.578;
  h2 = (-0.83624e-3 * pow( T, 3.0 ) + 0.16792 * pow( T, 2.0 ) - 25.9293 * T) * 4.184 / 58.44;

  return ( (x1 * h1 + x2 * h2 + m * dh) * 1000.0 );
}

void calculateBrineEnthalpy( PTTableCoordinates const & tableCoords,
                             real64 const & m,
                             array1d< real64 > const & enthalpies )
{

  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nTemperatures; ++i )
  {
    real64 const TC = tableCoords.getTemperature( i );
    enthalpies[i] = michaelidesBrineEnthalpy( TC, m );
  }
}

TableFunction const * makeCO2EnthalpyTable( string_array const & inputParams,
                                            string const & functionName,
                                            FunctionManager & functionManager )
{
  string const tableName = functionName + "_CO2_enthalpy_table";

  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    TableFunction * const enthalpyTable = functionManager.getGroupPointer< TableFunction >( tableName );
    enthalpyTable->initializeFunction();
    enthalpyTable->setDimUnits( PTTableCoordinates::coordsUnits );
    enthalpyTable->setValueUnits( units::Enthalpy );
    return enthalpyTable;
  }
  else
  {
    // initialize the (p,T) coordinates
    PTTableCoordinates tableCoords;
    PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

    real64 tolerance = 1e-10;

    try
    {
      if( inputParams.size() >= 10 )
      {
        tolerance = stod( inputParams[9] );
      }
    }
    catch( const std::invalid_argument & e )
    {
      GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
    }

    array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
    array1d< real64 > enthalpies( tableCoords.nPressures() * tableCoords.nTemperatures() );

    SpanWagnerCO2Density::calculateCO2Density( functionName, tolerance, tableCoords, densities );

    CO2Enthalpy::calculateCO2Enthalpy( tableCoords, densities, enthalpies );

    TableFunction * const enthalpyTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    enthalpyTable->setTableCoordinates( tableCoords.getCoords(), tableCoords.coordsUnits );
    enthalpyTable->setTableValues( enthalpies, units::Enthalpy );
    enthalpyTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return enthalpyTable;
  }
}

TableFunction const * makeBrineEnthalpyTable( string_array const & inputParams,
                                              string const & functionName,
                                              FunctionManager & functionManager )
{
  string const tableName = functionName + "_brine_enthalpy_table";

  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    TableFunction * const enthalpyTable = functionManager.getGroupPointer< TableFunction >( tableName );
    enthalpyTable->initializeFunction();
    enthalpyTable->setDimUnits( { PTTableCoordinates::coordsUnits[1] } );
    enthalpyTable->setValueUnits( units::Enthalpy );
    return enthalpyTable;
  }
  else
  {
    // initialize the (p,T) coordinates
    PTTableCoordinates tableCoords;
    PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

    // initialize salinity
    GEOS_THROW_IF_LT_MSG( inputParams.size(), 9,
                          GEOS_FMT( "{}: insufficient number of model parameters", functionName ),
                          InputError );
    real64 salinity;

    try
    {
      salinity = stod( inputParams[8] );
    }
    catch( std::invalid_argument const & e )
    {
      GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
    }

    array1d< real64 > enthalpies( tableCoords.nTemperatures() );
    array1d< array1d< real64 > > temperatures;
    temperatures.resize( 1 );
    temperatures[0] = tableCoords.getTemperatures();

    calculateBrineEnthalpy( tableCoords, salinity, enthalpies );


    TableFunction * const enthalpyTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    enthalpyTable->setTableCoordinates( temperatures, { tableCoords.coordsUnits[1] } );
    enthalpyTable->setTableValues( enthalpies, units::Enthalpy );
    enthalpyTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return enthalpyTable;
  }
}

} // namespace

BrineEnthalpy::BrineEnthalpy( string const & name,
                              string_array const & inputParams,
                              string_array const & componentNames,
                              array1d< real64 > const & componentMolarWeight,
                              TableFunction::OutputOptions const pvtOutputOpts ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames, "componentNames" );

  string const expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames, "componentNames" );

  m_CO2EnthalpyTable = makeCO2EnthalpyTable( inputParams, m_functionName, FunctionManager::getInstance() );
  m_brineEnthalpyTable = makeBrineEnthalpyTable( inputParams, m_functionName, FunctionManager::getInstance() );

  m_CO2EnthalpyTable->outputPVTTableData( pvtOutputOpts );
  m_brineEnthalpyTable->outputPVTTableData( pvtOutputOpts );
}


void BrineEnthalpy::checkTablesParameters( real64 const pressure,
                                           real64 const temperature ) const
{
  m_brineEnthalpyTable->checkCoord( temperature, 0 );

  m_CO2EnthalpyTable->checkCoord( pressure, 0 );
  m_CO2EnthalpyTable->checkCoord( temperature, 1 );
}



BrineEnthalpy::KernelWrapper
BrineEnthalpy::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_CO2EnthalpyTable,
                        *m_brineEnthalpyTable,
                        m_CO2Index,
                        m_waterIndex );
}

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
