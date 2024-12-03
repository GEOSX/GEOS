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
 * @file CO2Solubility.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Solubility.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2SolubilitySpycherPruess.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2SolubilityDuanSun.hpp"

#include "functions/FunctionManager.hpp"
#include "common/Units.hpp"

namespace geos
{

using namespace stringutilities;

namespace
{

TableFunction const * makeTable( string const & tableName,
                                 constitutive::PVTProps::PTTableCoordinates const & tableCoords,
                                 array1d< real64 > && values,
                                 FunctionManager & functionManager )
{
  TableFunction * tableFunction = nullptr;
  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    tableFunction = functionManager.getGroupPointer< TableFunction >( tableName );
  }
  else
  {
    tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );
    tableFunction->setTableCoordinates( tableCoords.getCoords(), { units::Pressure, units::TemperatureInC } );
    tableFunction->setTableValues( values, units::Solubility );
    tableFunction->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  }
  tableFunction->initializeFunction();
  return tableFunction;
}

std::pair< TableFunction const *, TableFunction const * >
makeSolubilityTables( string const & functionName,
                      string_array const & inputParams,
                      constitutive::PVTProps::CO2Solubility::SolubilityModel const & solubilityModel )
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  constitutive::PVTProps::PTTableCoordinates tableCoords;

  // Check solubility model for explicit table input
  if( solubilityModel == constitutive::PVTProps::CO2Solubility::SolubilityModel::Tables )
  {
    // The default table is a table with all zeros unless the name is explicitly provided
    // The pressure and temperature values below will be used only to create the zero table so they
    // simply give a range large enough to cover most values.
    tableCoords.appendPressure( 1.0e5 ).appendPressure( 1.0e8 )
      .appendTemperature( 0.0 ).appendTemperature( 800.0 );

    TableFunction const * tables[2] = { nullptr, nullptr };
    for( integer tableIndex : { 0, 1 } )
    {
      array1d< real64 > values( 4 );
      values.zero();

      string inputTableName = inputParams[2 + tableIndex];
      if( inputTableName.empty() )
      {
        inputTableName = GEOS_FMT( "{}_zeroDissolution_table", constitutive::PVTProps::CO2Solubility::catalogName() );
      }
      else
      {
        // If a name is explicitly given, then check that it exists
        GEOS_THROW_IF( !functionManager.hasGroup< TableFunction >( inputTableName ),
                       GEOS_FMT( "{}: Could not find TableFunction with name {}", functionName, inputTableName ),
                       InputError );
      }
      tables[tableIndex] = makeTable( inputTableName, tableCoords, std::move ( values ), functionManager );
    }
    return { tables[0], tables[1] };
  }

  // If the tables have already been created, then simply retrieve them
  string const co2TableName = functionName + "_co2Dissolution_table";
  string const h2oTableName = functionName + "_waterVaporization_table";

  if( functionManager.hasGroup< TableFunction >( co2TableName ) && functionManager.hasGroup< TableFunction >( h2oTableName ))
  {
    TableFunction const * co2SolubilityTable = functionManager.getGroupPointer< TableFunction >( co2TableName );
    TableFunction const * h2oSolubilityTable = functionManager.getGroupPointer< TableFunction >( h2oTableName );
    return {co2SolubilityTable, h2oSolubilityTable};
  }

  // Initialize the (p,T) coordinates
  constitutive::PVTProps::PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

  // Initialize salinity and tolerance
  GEOS_THROW_IF_LT_MSG( inputParams.size(), 9,
                        GEOS_FMT( "{}: insufficient number of model parameters", functionName ),
                        InputError );

  real64 tolerance = 1e-9;
  real64 salinity = 0.0;
  try
  {
    salinity = stod( inputParams[8] );
    if( inputParams.size() >= 10 )
    {
      tolerance = stod( inputParams[9] );
    }
  }
  catch( const std::invalid_argument & e )
  {
    GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
  }

  integer const nPressures = tableCoords.nPressures();
  integer const nTemperatures = tableCoords.nTemperatures();

  array1d< real64 > co2Solubility( nPressures * nTemperatures );
  array1d< real64 > h2oSolubility( nPressures * nTemperatures );

  if( solubilityModel == constitutive::PVTProps::CO2Solubility::SolubilityModel::DuanSun )
  {
    constitutive::PVTProps::CO2SolubilityDuanSun::populateSolubilityTables(
      functionName,
      tableCoords,
      salinity,
      tolerance,
      co2Solubility,
      h2oSolubility );
  }
  else if( solubilityModel == constitutive::PVTProps::CO2Solubility::SolubilityModel::SpycherPruess )
  {
    constitutive::PVTProps::CO2SolubilitySpycherPruess::populateSolubilityTables(
      functionName,
      tableCoords,
      salinity,
      tolerance,
      co2Solubility,
      h2oSolubility );
  }

  // Truncate negative solubility and warn
  integer constexpr maxBad = 5;     // Maximum number of bad values to report
  stackArray2d< real64, maxBad *4 > badValues( maxBad, 4 );
  integer badCount = 0;
  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const P = tableCoords.getPressure( i );
    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const T = tableCoords.getTemperature( j );
      if( co2Solubility[j*nPressures+i] < 0.0 || h2oSolubility[j*nPressures+i] < 0.0 )
      {
        badValues( badCount % maxBad, 0 ) = P;
        badValues( badCount % maxBad, 1 ) = T;
        badValues( badCount % maxBad, 2 ) = co2Solubility[j*nPressures+i];
        badValues( badCount % maxBad, 3 ) = h2oSolubility[j*nPressures+i];
        ++badCount;
      }
      if( co2Solubility[j*nPressures+i] < 0.0 )
      {
        co2Solubility[j*nPressures+i] = 0.0;
      }
      if( h2oSolubility[j*nPressures+i] < 0.0 )
      {
        h2oSolubility[j*nPressures+i] = 0.0;
      }
    }
  }

  if( 0 < badCount )
  {
    std::ostringstream badValueTable;
    badValueTable
      << std::setw( 15 ) << "Pressure (Pa)" << " "
      << std::setw( 15 ) << "Temperature (C)" << " "
      << std::setw( 23 ) << "CO2 solubility (mol/kg)" << " "
      << std::setw( 23 ) << "H2O solubility (mol/kg)" << " ";
    for( integer row = 0; row < LvArray::math::min( maxBad, badCount ); ++row )
    {
      badValueTable
        << "\n"
        << std::setw( 15 ) << badValues( row, 0 ) << " "
        << std::setw( 15 ) << badValues( row, 1 ) << " "
        << std::setw( 23 ) << badValues( row, 2 ) << " "
        << std::setw( 23 ) << badValues( row, 3 ) << " ";
    }
    GEOS_LOG_RANK_0( GEOS_FMT( "CO2Solubility: {} negative solubility values encountered. These will be truncated to zero.\nCheck out report table with max {} values.\n{}",
                               badCount, maxBad, badValueTable.str() ) );
  }

  TableFunction const * co2SolubilityTable = makeTable( co2TableName, tableCoords, std::move( co2Solubility ), functionManager );
  TableFunction const * h2oSolubilityTable = makeTable( h2oTableName, tableCoords, std::move( h2oSolubility ), functionManager );
  return {co2SolubilityTable, h2oSolubilityTable};
}

} // namespace

namespace constitutive
{
namespace PVTProps
{

CO2Solubility::CO2Solubility( string const & name,
                              string_array const & inputParams,
                              string_array const & phaseNames,
                              string_array const & componentNames,
                              array1d< real64 > const & componentMolarWeight,
                              TableFunction::OutputOptions const pvtOutputOpts ):
  FlashModelBase( name,
                  componentNames,
                  componentMolarWeight )
{
  GEOS_THROW_IF_NE_MSG( phaseNames.size(), 2,
                        "The CO2Solubility model is a two-phase model",
                        InputError );
  GEOS_THROW_IF_NE_MSG( componentNames.size(), 2,
                        "The CO2Solubility model is a two-component model",
                        InputError );

  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames, "componentNames" );

  string const expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames, "componentNames" );

  string const expectedGasPhaseNames[] = { "CO2", "co2", "gas", "Gas" };
  m_phaseGasIndex = PVTFunctionHelpers::findName( phaseNames, expectedGasPhaseNames, "phaseNames" );

  string const expectedWaterPhaseNames[] = { "Water", "water", "Liquid", "liquid" };
  m_phaseLiquidIndex = PVTFunctionHelpers::findName( phaseNames, expectedWaterPhaseNames, "phaseNames" );

  SolubilityModel solubilityModel = SolubilityModel::DuanSun;   // Default solubility model
  if( inputParams[1] == EnumStrings< SolubilityModel >::toString( SolubilityModel::Tables ) )
  {
    solubilityModel = SolubilityModel::Tables;
  }
  else if( 11 <= inputParams.size() )
  {
    solubilityModel = EnumStrings< SolubilityModel >::fromString( inputParams[10] );
  }

  std::tie( m_CO2SolubilityTable, m_WaterVapourisationTable ) = makeSolubilityTables( m_modelName, inputParams, solubilityModel );

  m_CO2SolubilityTable->outputPVTTableData( pvtOutputOpts );
  m_WaterVapourisationTable->outputPVTTableData( pvtOutputOpts );

}

void CO2Solubility::checkTablesParameters( real64 const pressure,
                                           real64 const temperature ) const
{
  m_CO2SolubilityTable->checkCoord( pressure, 0 );
  m_CO2SolubilityTable->checkCoord( temperature, 1 );
  m_WaterVapourisationTable->checkCoord( pressure, 0 );
  m_WaterVapourisationTable->checkCoord( temperature, 1 );
}

CO2Solubility::KernelWrapper CO2Solubility::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_CO2SolubilityTable,
                        *m_WaterVapourisationTable,
                        m_CO2Index,
                        m_waterIndex,
                        m_phaseGasIndex,
                        m_phaseLiquidIndex );
}

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geos
