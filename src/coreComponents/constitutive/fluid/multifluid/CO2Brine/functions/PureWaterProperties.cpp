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
 * @file PureWaterProperties.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

TableFunction const *
PureWaterProperties::makeSaturationViscosityTable( string const & functionName,
                                                   FunctionManager & functionManager,
                                                   bool const printTable )
{
  array1d< array1d< real64 > > temperatures;
  array1d< real64 > viscosities;

  integer const nValues = 26;
  temperatures.resize( 1 );
  temperatures[0].resize( nValues );
  viscosities.resize( nValues );

  temperatures[0][0] = 0.01;
  temperatures[0][1] = 10;
  temperatures[0][2] = 20;
  temperatures[0][3] = 25;
  temperatures[0][4] = 30;
  temperatures[0][5] = 40;
  temperatures[0][6] = 50;
  temperatures[0][7] = 60;
  temperatures[0][8] = 70;
  temperatures[0][9] = 80;
  temperatures[0][10] = 90;
  temperatures[0][11] = 100;
  temperatures[0][12] = 110;
  temperatures[0][13] = 120;
  temperatures[0][14] = 140;
  temperatures[0][15] = 160;
  temperatures[0][16] = 180;
  temperatures[0][17] = 200;
  temperatures[0][18] = 220;
  temperatures[0][19] = 240;
  temperatures[0][20] = 260;
  temperatures[0][21] = 280;
  temperatures[0][22] = 300;
  temperatures[0][23] = 320;
  temperatures[0][24] = 340;
  temperatures[0][25] = 360;

  viscosities[0] = 0.0017914;
  viscosities[1] = 0.0013060;
  viscosities[2] = 0.0010016;
  viscosities[3] = 0.0008900;
  viscosities[4] = 0.0007972;
  viscosities[5] = 0.0006527;
  viscosities[6] = 0.0005465;
  viscosities[7] = 0.0004660;
  viscosities[8] = 0.0004035;
  viscosities[9] = 0.0003540;
  viscosities[10] = 0.0003142;
  viscosities[11] = 0.0002816;
  viscosities[12] = 0.0002546;
  viscosities[13] = 0.0002320;
  viscosities[14] = 0.0001966;
  viscosities[15] = 0.0001704;
  viscosities[16] = 0.0001504;
  viscosities[17] = 0.0001346;
  viscosities[18] = 0.0001218;
  viscosities[19] = 0.0001111;
  viscosities[20] = 0.0001018;
  viscosities[21] = 0.0000936;
  viscosities[22] = 0.0000859;
  viscosities[23] = 0.0000783;
  viscosities[24] = 0.0000703;
  viscosities[25] = 0.0000603;

  if( printTable )
  {
    std::ofstream table_file( "PureWaterSaturationViscosity.csv" );
    table_file << "T[C]" << std::endl;
    for( localIndex i = 0; i < nValues; ++i )
    {
      table_file << temperatures[0][i] << "," <<viscosities[i] << std::endl;
    }
    table_file.close();
  }

  string const tableName = functionName +  "_table";
  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< TableFunction >( tableName );
  }
  else
  {
    TableFunction * const viscosityTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    viscosityTable->setTableCoordinates( temperatures, { units::TemperatureInC } );
    viscosityTable->setTableValues( viscosities, units::Viscosity );
    viscosityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return viscosityTable;
  }
}

TableFunction const *
PureWaterProperties::makeSaturationDensityTable( string const & functionName,
                                                 FunctionManager & functionManager,
                                                 bool const printTable )
{
  array1d< array1d< real64 > > temperatures;
  array1d< real64 > densities;

  integer const nValues = 26;
  temperatures.resize( 1 );
  temperatures[0].resize( nValues );
  densities.resize( nValues );

  temperatures[0][0] = 0.01;
  temperatures[0][1] = 10;
  temperatures[0][2] = 20;
  temperatures[0][3] = 25;
  temperatures[0][4] = 30;
  temperatures[0][5] = 40;
  temperatures[0][6] = 50;
  temperatures[0][7] = 60;
  temperatures[0][8] = 70;
  temperatures[0][9] = 80;
  temperatures[0][10] = 90;
  temperatures[0][11] = 100;
  temperatures[0][12] = 110;
  temperatures[0][13] = 120;
  temperatures[0][14] = 140;
  temperatures[0][15] = 160;
  temperatures[0][16] = 180;
  temperatures[0][17] = 200;
  temperatures[0][18] = 220;
  temperatures[0][19] = 240;
  temperatures[0][20] = 260;
  temperatures[0][21] = 280;
  temperatures[0][22] = 300;
  temperatures[0][23] = 320;
  temperatures[0][24] = 340;
  temperatures[0][25] = 360;

  densities[0] = 999.85;
  densities[1] = 999.7;
  densities[2] = 998.21;
  densities[3] = 997.05;
  densities[4] = 995.65;
  densities[5] = 992.25;
  densities[6] = 988.04;
  densities[7] = 983.2;
  densities[8] = 977.76;
  densities[9] = 971.79;
  densities[10] = 965.31;
  densities[11] = 958.35;
  densities[12] = 950.95;
  densities[13] = 943.11;
  densities[14] = 926.13;
  densities[15] = 907.45;
  densities[16] = 887;
  densities[17] = 864.66;
  densities[18] = 840.22;
  densities[19] = 813.37;
  densities[20] = 783.63;
  densities[21] = 750.28;
  densities[22] = 712.14;
  densities[23] = 667.09;
  densities[24] = 610.67;
  densities[25] = 527.59;

  if( printTable )
  {
    std::ofstream table_file( "PureWaterSaturationDensity.csv" );
    table_file << "T[C]" << std::endl;
    for( localIndex i = 0; i < nValues; ++i )
    {
      table_file << temperatures[0][i] << "," << densities[i] << std::endl;
    }
    table_file.close();
  }

  string const tableName = functionName +  "_sat_density_table";
  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< TableFunction >( tableName );
  }
  else
  {
    TableFunction * const densityTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    densityTable->setTableCoordinates( temperatures, { units::TemperatureInC } );
    densityTable->setTableValues( densities, units::Density );
    densityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return densityTable;
  }
}

TableFunction const *
PureWaterProperties::makeSaturationPressureTable( string const & functionName,
                                                  FunctionManager & functionManager,
                                                  bool const printTable )
{
  array1d< array1d< real64 > > temperatures;
  array1d< real64 > pressures;

  integer const nValues = 26;
  temperatures.resize( 1 );
  temperatures[0].resize( nValues );
  pressures.resize( nValues );

  temperatures[0][0] = 0.01;
  temperatures[0][1] = 10;
  temperatures[0][2] = 20;
  temperatures[0][3] = 25;
  temperatures[0][4] = 30;
  temperatures[0][5] = 40;
  temperatures[0][6] = 50;
  temperatures[0][7] = 60;
  temperatures[0][8] = 70;
  temperatures[0][9] = 80;
  temperatures[0][10] = 90;
  temperatures[0][11] = 100;
  temperatures[0][12] = 110;
  temperatures[0][13] = 120;
  temperatures[0][14] = 140;
  temperatures[0][15] = 160;
  temperatures[0][16] = 180;
  temperatures[0][17] = 200;
  temperatures[0][18] = 220;
  temperatures[0][19] = 240;
  temperatures[0][20] = 260;
  temperatures[0][21] = 280;
  temperatures[0][22] = 300;
  temperatures[0][23] = 320;
  temperatures[0][24] = 340;
  temperatures[0][25] = 360;

  pressures[0] =  612;
  pressures[1] =  1200;
  pressures[2] =  2300;
  pressures[3] =  3200;
  pressures[4] =  4200;
  pressures[5] =  7400;
  pressures[6] =  12400;
  pressures[7] =  19900;
  pressures[8] =  31200;
  pressures[9] =  47400;
  pressures[10] = 70200;
  pressures[11] = 101000;
  pressures[12] = 143000;
  pressures[13] = 199000;
  pressures[14] = 362000;
  pressures[15] = 618000;
  pressures[16] = 1000000;
  pressures[17] = 1550000;
  pressures[18] = 2320000;
  pressures[19] = 3350000;
  pressures[20] = 4690000;
  pressures[21] = 6420000;
  pressures[22] = 8590000;
  pressures[23] = 11300000;
  pressures[24] = 14600000;
  pressures[25] = 18700000;

  if( printTable )
  {
    std::ofstream table_file( "PureWaterSaturationPressure.csv" );
    table_file << "T[C]" << std::endl;
    for( localIndex i = 0; i < nValues; ++i )
    {
      table_file << temperatures[0][i] << "," << pressures[i] << std::endl;
    }
    table_file.close();
  }

  string const tableName = functionName +  "_sat_pressure_table";
  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< TableFunction >( tableName );
  }
  else
  {
    TableFunction * const pressureTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    pressureTable->setTableCoordinates( temperatures, { units::TemperatureInC } );
    pressureTable->setTableValues( pressures, units::Pressure );
    pressureTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return pressureTable;
  }
}

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geos
