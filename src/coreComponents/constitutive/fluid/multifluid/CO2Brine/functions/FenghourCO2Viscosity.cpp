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
 * @file FenghourCO2Viscosity.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/FenghourCO2Viscosity.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp"
#include "functions/FunctionManager.hpp"
#include "common/Units.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

namespace
{

real64 fenghourCO2ViscosityFunction( real64 const & temperatureCent,
                                     real64 const & density )
{
  constexpr real64 espar = 251.196;
  constexpr real64 esparInv = 1.0 / espar;
  // coefficients from Table (1) of Fenghour and Wakeham (1998)
  constexpr real64 aa[5] = { 0.235156, -0.491266, 5.211155e-2, 5.347906e-2, -1.537102e-2 };
  // coefficients from Table (3) of Fenghour and Wakeham (1998)
  constexpr real64 d11 =  0.4071119e-2;
  constexpr real64 d21 =  0.7198037e-4;
  constexpr real64 d64 =  0.2411697e-16;
  constexpr real64 d81 =  0.2971072e-22;
  constexpr real64 d82 = -0.1627888e-22;
  // we neglect critical viscosity
  constexpr real64 vcrit = 0.0;

  // temperature in Kelvin
  real64 const temperatureKelvin = units::convertCToK( temperatureCent );
  // equation (5) of Fenghour and Wakeham (1998)
  real64 const Tred = temperatureKelvin * esparInv;
  real64 const x = log( Tred );
  // equation (4) of Fenghour and Wakeham (1998)
  real64 const lnGfun = aa[0] + x * (aa[1] + x * (aa[2] + x *(aa[3] + x * aa[4])));
  real64 const GfunInv = exp( -lnGfun );
  // equation (3) of Fenghour and Wakeham (1998)
  real64 const vlimit = 1.00697 * sqrt( temperatureKelvin ) * GfunInv;

  // equation (8) of Fenghour and Wakeham (1998)
  real64 const d2 = density * density;
  real64 const vxcess = density * (d11 + density * (d21 + d2*d2*(d64 / (Tred*Tred*Tred) + d2*(d81 + d82/Tred))));

  // equation (1) of Fenghour and Wakeham (1998)
  return 1e-6 * (vlimit + vxcess + vcrit);
}

void calculateCO2Viscosity( PTTableCoordinates const & tableCoords,
                            array1d< real64 > const & densities,
                            array1d< real64 > const & viscosities )
{

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const T = tableCoords.getTemperature( j );
      viscosities[j*nPressures+i] = fenghourCO2ViscosityFunction( T, densities[j*nPressures+i] );
    }
  }
}

TableFunction const * makeViscosityTable( string_array const & inputParams,
                                          string const & functionName,
                                          FunctionManager & functionManager )
{
  string const tableName = functionName + "_table";

  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    TableFunction * const viscosityTable = functionManager.getGroupPointer< TableFunction >( tableName );
    viscosityTable->initializeFunction();
    viscosityTable->setDimUnits( { units::Pressure, units::TemperatureInC } );
    viscosityTable->setValueUnits( units::Viscosity );
    return viscosityTable;
  }
  else
  {
    PTTableCoordinates tableCoords;
    PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

    real64 tolerance = 1e-10;
    try
    {
      if( inputParams.size() >= 9 )
      {
        tolerance = stod( inputParams[8] );
      }
    }
    catch( const std::invalid_argument & e )
    {
      GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
    }

    localIndex const nP = tableCoords.nPressures();
    localIndex const nT = tableCoords.nTemperatures();
    array1d< real64 > density( nP * nT );
    array1d< real64 > viscosity( nP * nT );
    SpanWagnerCO2Density::calculateCO2Density( functionName, tolerance, tableCoords, density );
    calculateCO2Viscosity( tableCoords, density, viscosity );

    TableFunction * const viscosityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );
    viscosityTable->setTableCoordinates( tableCoords.getCoords(),
                                         { units::Pressure, units::TemperatureInC } );
    viscosityTable->setTableValues( viscosity, units::Viscosity );
    viscosityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return viscosityTable;
  }
}

} // namespace

FenghourCO2Viscosity::FenghourCO2Viscosity( string const & name,
                                            string_array const & inputParams,
                                            string_array const & componentNames,
                                            array1d< real64 > const & componentMolarWeight,
                                            TableFunction::OutputOptions const pvtOutputOpts )
  : PVTFunctionBase( name,
                     componentNames,
                     componentMolarWeight )
{
  m_CO2ViscosityTable = makeViscosityTable( inputParams, m_functionName, FunctionManager::getInstance() );

  m_CO2ViscosityTable->outputPVTTableData( pvtOutputOpts );
}

void FenghourCO2Viscosity::checkTablesParameters( real64 const pressure,
                                                  real64 const temperature ) const
{
  m_CO2ViscosityTable->checkCoord( pressure, 0 );
  m_CO2ViscosityTable->checkCoord( temperature, 1 );
}

FenghourCO2Viscosity::KernelWrapper
FenghourCO2Viscosity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_CO2ViscosityTable );
}

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geos
