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

// Source includes
#include "constitutiveTestHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::constitutive::cappres;
using namespace geos::dataRepository;

CapillaryPressureBase & makeBrooksCoreyCapPressureTwoPhase( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "oil";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1; phaseMinSat[1] = 0.05;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = 4;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 2 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = 1;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postInputInitializationRecursive();
  return capPressure;
}


CapillaryPressureBase & makeBrooksCoreyCapPressureThreePhase( string const & name, Group & parent )
{
  BrooksCoreyCapillaryPressure & capPressure = parent.registerGroup< BrooksCoreyCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "water"; phaseNames[1] = "oil"; phaseNames[2] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.02; phaseMinSat[2] = 0.1;

  array1d< real64 > & phaseCapPressureExpInv =
    capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 2; phaseCapPressureExpInv[1] = -3; phaseCapPressureExpInv[2] = 2.5;

  array1d< real64 > & phaseEntryPressure = capPressure.getReference< array1d< real64 > >( BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString() );
  phaseEntryPressure.resize( 3 );
  phaseEntryPressure[0] = 1; phaseEntryPressure[1] = -1; phaseEntryPressure[2] = 2;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-7;

  capPressure.postInputInitializationRecursive();
  return capPressure;
}


CapillaryPressureBase & makeVanGenuchtenCapPressureTwoPhase( string const & name, Group & parent )
{
  VanGenuchtenCapillaryPressure & capPressure = parent.registerGroup< VanGenuchtenCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1;

  array1d< real64 > & phaseCapPressureExpInv = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 0.4; phaseCapPressureExpInv[1] = 0.5;

  array1d< real64 > & phaseCapPressureMultiplier = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString() );
  phaseCapPressureMultiplier.resize( 2 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postInputInitializationRecursive();
  return capPressure;
}

CapillaryPressureBase & makeVanGenuchtenCapPressureThreePhase( string const & name, Group & parent )
{
  VanGenuchtenCapillaryPressure & capPressure = parent.registerGroup< VanGenuchtenCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & phaseMinSat = capPressure.getReference< array1d< real64 > >( VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString() );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04; phaseMinSat[1] = 0.1; phaseMinSat[2] = 0.;

  array1d< real64 > & phaseCapPressureExpInv = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString() );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 0.33; phaseCapPressureExpInv[1] = 0.4; phaseCapPressureExpInv[ 2 ] = 0.5;

  array1d< real64 > & phaseCapPressureMultiplier = capPressure.getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString() );
  phaseCapPressureMultiplier.resize( 3 );
  phaseCapPressureMultiplier[0] = 0.5; phaseCapPressureMultiplier[1] = 1; phaseCapPressureMultiplier[2] = 0.2;

  real64 & capPressureEpsilon = capPressure.getReference< real64 >( VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString() );
  capPressureEpsilon = 1e-4;

  capPressure.postInputInitializationRecursive();
  return capPressure;
}

CapillaryPressureBase & makeTableCapPressureTwoPhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex Naxis = 6;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = 0.22;
  coordinates[0][1] = 0.3;
  coordinates[0][2] = 0.5;
  coordinates[0][3] = 0.6;
  coordinates[0][4] = 0.8;
  coordinates[0][5] = 1.0;

  real64_array values( Naxis );
  values[0] = 48263.3;
  values[1] = 27579;
  values[2] = 20684.3;
  values[3] = 13789.5;
  values[4] = 6894.76;
  values[5] = 3204.28;

  TableFunction & table_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_pc" ) );
  table_w.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_w.setTableValues( values, units::Pressure );
  table_w.reInitializeFunction();

  table_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  TableCapillaryPressure & capPressure = parent.registerGroup< TableCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  auto & waterTableName = capPressure.getReference< string >( TableCapillaryPressure::viewKeyStruct::wettingNonWettingCapPresTableNameString() );
  waterTableName = "water_pc";

  capPressure.postInputInitializationRecursive();
  capPressure.initialize(); // to test all the checks
  return capPressure;
}

CapillaryPressureBase & makeTableCapPressureThreePhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex const Naxis = 6;

  // 1.a) First pair of phases (ow)

  // Setup table
  array1d< real64_array > coordinates_w;
  coordinates_w.resize( 1 );
  coordinates_w[0].resize( Naxis );
  coordinates_w[0][0] = 0.22;
  coordinates_w[0][1] = 0.3;
  coordinates_w[0][2] = 0.5;
  coordinates_w[0][3] = 0.6;
  coordinates_w[0][4] = 0.8;
  coordinates_w[0][5] = 1.0;

  real64_array values_w( Naxis );
  values_w[0] = 48263.3;
  values_w[1] = 27579;
  values_w[2] = 20684.3;
  values_w[3] = 13789.5;
  values_w[4] = 6894.76;
  values_w[5] = 3294.76;

  TableFunction & table_ow_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_ow_pc" ) );
  table_ow_w.setTableCoordinates( coordinates_w, { units::Dimensionless } );
  table_ow_w.setTableValues( values_w, units::Pressure );
  table_ow_w.reInitializeFunction();

  table_ow_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 1.a) Second pair of phases (og)

  array1d< real64_array > coordinates_g;
  coordinates_g.resize( 1 );
  coordinates_g[0].resize( Naxis );
  coordinates_g[0][0] = 0.0;
  coordinates_g[0][1] = 0.04;
  coordinates_g[0][2] = 0.1;
  coordinates_g[0][3] = 0.2;
  coordinates_g[0][4] = 0.7;
  coordinates_g[0][5] = 0.78;

  real64_array values_g( Naxis );
  values_g[0] = 0.0;
  values_g[1] = 1723.689;
  values_g[2] = 3447.38;
  values_g[3] = 6894.76;
  values_g[4] = 24131.7;
  values_g[5] = 26889.6;

  TableFunction & table_og_g = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "gas_og_pc" ) );
  table_og_g.setTableCoordinates( coordinates_g, { units::Dimensionless } );
  table_og_g.setTableValues( values_g, units::Pressure );
  table_og_g.reInitializeFunction();

  table_og_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  TableCapillaryPressure & capPressure = parent.registerGroup< TableCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & waterTableName = capPressure.getReference< string >( TableCapillaryPressure::viewKeyStruct::wettingIntermediateCapPresTableNameString() );
  waterTableName = "water_ow_pc";

  auto & gasTableName = capPressure.getReference< string >( TableCapillaryPressure::viewKeyStruct::nonWettingIntermediateCapPresTableNameString() );
  gasTableName = "gas_og_pc";

  capPressure.postInputInitializationRecursive();
  capPressure.initialize(); // to test all the checks
  return capPressure;
}

CapillaryPressureBase & makeJFunctionCapPressureTwoPhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex const Naxis = 10;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );

  coordinates[0][0] = 0.2924;
  coordinates[0][1] = 0.3639;
  coordinates[0][2] = 0.4354;
  coordinates[0][3] = 0.5068;
  coordinates[0][4] = 0.5783;
  coordinates[0][5] = 0.6498;
  coordinates[0][6] = 0.7213;
  coordinates[0][7] = 0.7927;
  coordinates[0][8] = 0.8642;
  coordinates[0][9] = 0.9357;

  real64_array values( Naxis );
  values[0] = 0.3772;
  values[1] = 0.0908;
  values[2] = 0.0607;
  values[3] = 0.0468;
  values[4] = 0.0381;
  values[5] = 0.0318;
  values[6] = 0.0267;
  values[7] = 0.0222;
  values[8] = 0.0178;
  values[9] = 0.0127;


  TableFunction & table_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_jFunction" ) );
  table_w.setTableCoordinates( coordinates, { units::Dimensionless } );
  table_w.setTableValues( values, units::Dimensionless );
  table_w.reInitializeFunction();

  table_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  JFunctionCapillaryPressure & capPressure = parent.registerGroup< JFunctionCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "water"; phaseNames[1] = "gas";

  auto & waterTableName = capPressure.getReference< string >( JFunctionCapillaryPressure::viewKeyStruct::wettingNonWettingJFuncTableNameString() );
  waterTableName = "water_jFunction";

  auto & surfaceTension = capPressure.getReference< real64 >( JFunctionCapillaryPressure::viewKeyStruct::wettingNonWettingSurfaceTensionString() );
  surfaceTension = 23.86955676433857e-3;

  auto & permeabilityDirection =
    capPressure.getReference< JFunctionCapillaryPressure::PermeabilityDirection >( JFunctionCapillaryPressure::viewKeyStruct::permeabilityDirectionString() );
  permeabilityDirection = JFunctionCapillaryPressure::PermeabilityDirection::XY;

  capPressure.postInputInitializationRecursive();
  capPressure.initialize(); // to test all the checks

  return capPressure;
}

CapillaryPressureBase & makeJFunctionCapPressureThreePhase( string const & name, Group & parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables

  // 1D table, various interpolation methods
  localIndex const Naxis = 6;

  // 1.a) First pair of phases (ow)

  // Setup table
  array1d< real64_array > coordinates_w;
  coordinates_w.resize( 1 );
  coordinates_w[0].resize( Naxis );
  coordinates_w[0][0] = 0.22;
  coordinates_w[0][1] = 0.3;
  coordinates_w[0][2] = 0.5;
  coordinates_w[0][3] = 0.6;
  coordinates_w[0][4] = 0.8;
  coordinates_w[0][5] = 1.0;

  real64_array values_w( Naxis );
  values_w[0] = 0.40203;
  values_w[1] = 0.31311;
  values_w[2] = 0.22423;
  values_w[3] = 0.15011;
  values_w[4] = 0.04224;
  values_w[5] = 0.00603;

  TableFunction & table_ow_w = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "water_ow_jFunction" ) );
  table_ow_w.setTableCoordinates( coordinates_w, { units::Dimensionless } );
  table_ow_w.setTableValues( values_w, units::Dimensionless );
  table_ow_w.reInitializeFunction();

  table_ow_w.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 1.a) Second pair of phases (og)

  array1d< real64_array > coordinates_g;
  coordinates_g.resize( 1 );
  coordinates_g[0].resize( Naxis );
  coordinates_g[0][0] = 0.0;
  coordinates_g[0][1] = 0.04;
  coordinates_g[0][2] = 0.1;
  coordinates_g[0][3] = 0.2;
  coordinates_g[0][4] = 0.7;
  coordinates_g[0][5] = 0.78;

  real64_array values_g( Naxis );
  values_g[0] = 0.0;
  values_g[1] = 0.02101;
  values_g[2] = 0.11607;
  values_g[3] = 0.19103;
  values_g[4] = 0.21033;
  values_g[5] = 0.27203;

  TableFunction & table_og_g = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "gas_og_jFunction" ) );
  table_og_g.setTableCoordinates( coordinates_g, { units::Dimensionless } );
  table_og_g.setTableValues( values_g, units::Dimensionless );
  table_og_g.reInitializeFunction();

  table_og_g.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then set up the constitutive model

  JFunctionCapillaryPressure & capPressure = parent.registerGroup< JFunctionCapillaryPressure >( name );

  string_array & phaseNames = capPressure.getReference< string_array >( CapillaryPressureBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  auto & waterTableName = capPressure.getReference< string >( JFunctionCapillaryPressure::viewKeyStruct::wettingIntermediateJFuncTableNameString() );
  waterTableName = "water_ow_jFunction";

  auto & gasTableName = capPressure.getReference< string >( JFunctionCapillaryPressure::viewKeyStruct::nonWettingIntermediateJFuncTableNameString() );
  gasTableName = "gas_og_jFunction";

  auto & waterSurfaceTension = capPressure.getReference< real64 >( JFunctionCapillaryPressure::viewKeyStruct::wettingIntermediateSurfaceTensionString() );
  waterSurfaceTension = 23.86955676433857e-3;

  auto & gasSurfaceTension = capPressure.getReference< real64 >( JFunctionCapillaryPressure::viewKeyStruct::nonWettingIntermediateSurfaceTensionString() );
  gasSurfaceTension = 14.24643678933543e-3;

  auto & permeabilityDirection =
    capPressure.getReference< JFunctionCapillaryPressure::PermeabilityDirection >( JFunctionCapillaryPressure::viewKeyStruct::permeabilityDirectionString() );
  permeabilityDirection = JFunctionCapillaryPressure::PermeabilityDirection::Z;

  capPressure.postInputInitializationRecursive();
  capPressure.initialize(); // to test all the checks
  return capPressure;
}


class CapillaryPressureTest : public ConstitutiveTestBase< CapillaryPressureBase >
{
public:
  void test( arraySlice1d< real64 const > const sat,
             real64 const eps,
             real64 const tol )
  {
    arrayView3d< real64 const, USD_CAPPRES > phaseCapPressure;
    arrayView4d< real64 const, USD_CAPPRES_DS > dPhaseCapPressure_dPhaseVolFraction;
    testNumericalDerivatives( m_parent,
                              *m_model,
                              sat,
                              eps,
                              tol,
                              "phaseCapPressure",
                              [&phaseCapPressure] ( CapillaryPressureBase & relPerm )
    {
      phaseCapPressure = relPerm.phaseCapPressure();
      return phaseCapPressure[ 0 ][ 0 ];
    },
                              [&dPhaseCapPressure_dPhaseVolFraction] ( CapillaryPressureBase & relPerm )
    {
      dPhaseCapPressure_dPhaseVolFraction = relPerm.dPhaseCapPressure_dPhaseVolFraction();
      return dPhaseCapPressure_dPhaseVolFraction[ 0 ][ 0 ];
    }
                              );
  }
};

TEST_F( CapillaryPressureTest, numericalDerivatives_brooksCoreyCapPressureTwoPhase )
{
  initialize( makeBrooksCoreyCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const startSat = 0.4;
  real64 const endSat = 0.6;
  real64 const dS = 1e-1;

  array1d< real64 > sat( 2 );

  sat[0] = startSat;
  sat[1] = 1.0-sat[0];

  while( sat[0] <= endSat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_brooksCoreyCapPressureThreePhase )
{
  initialize( makeBrooksCoreyCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const startSat = 0.4;
  real64 const endSat = 0.6;
  real64 const dS = 1e-1;

  array1d< real64 > sat( 3 );

  sat[0] = startSat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];

  while( sat[0] <= endSat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5 * ( 1-sat[0] );
    sat[2] = 1.0 - sat[0] - sat[1];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_vanGenuchtenCapPressureTwoPhase )
{
  initialize( makeVanGenuchtenCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const startSat = 0.4;
  real64 const endSat = 0.6;
  real64 const dS = 1e-1;

  array1d< real64 > sat( 2 );

  sat[0] = startSat;
  sat[1] = 1-sat[1];

  while( sat[0] <= endSat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }

}


TEST_F( CapillaryPressureTest, numericalDerivatives_vanGenuchtenCapPressureThreePhase )
{
  initialize( makeVanGenuchtenCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const startSat = 0.4;
  real64 const endSat = 0.6;
  real64 const dS = 1e-1;

  array1d< real64 > sat( 3 );

  sat[0] = startSat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];

  while( sat[0] <= endSat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5*(1-sat[0]);
    sat[2] = 1 - sat[0] - sat[1];
  }
}


TEST_F( CapillaryPressureTest, numericalDerivatives_tableCapPressureTwoPhase )
{
  initialize( makeTableCapPressureTwoPhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const startSat = 0.25;
  real64 const endSat   = 0.75;
  real64 const dS        = 1e-1;

  array1d< real64 > sat( 2 );
  sat[0] = startSat; sat[1] = 1-sat[0];
  while( sat[0] <= endSat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }

}


TEST_F( CapillaryPressureTest, numericalDerivatives_tableCapPressureThreePhase )
{
  initialize( makeTableCapPressureThreePhase( "capPressure", m_parent ) );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const startSat = 0.25;
  real64 const endSat = 0.75;
  real64 const dS = 1e-1;

  array1d< real64 > sat( 3 );

  sat[0] = startSat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];

  while( sat[0] <= endSat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5*(1-sat[0]);
    sat[2] = 1 - sat[0] - sat[1];
  }
}

TEST_F( CapillaryPressureTest, numericalDerivatives_jFunctionCapPressureTwoPhase )
{
  initialize( makeJFunctionCapPressureTwoPhase( "capPressure", m_parent ) );

  // here, we have to apply a special treatment to this test
  // to make sure that the J-function multiplier is initialized using initializeRockState
  // this requires calling allocateConstitutiveData in advance (it will be called again later, in the "test" function)

  // setup some values for porosity and permeability
  array2d< real64 > porosity;
  porosity.resize( 1, 1 );
  porosity[0][0] = 0.13496794266569806;
  array3d< real64 > permeability;
  permeability.resize( 1, 1, 3 );
  permeability[0][0][0] = 0.1722194e-15;
  permeability[0][0][1] = 0.3423156e-15;
  permeability[0][0][2] = 0.2324191e-15;

  // initialize the J-function multiplier (done on GPU if GPU is available)
  m_model->allocateConstitutiveData( m_parent, 1 );
  m_model->initializeRockState( porosity.toViewConst(), permeability.toViewConst() );

  // move the multiplier back to the CPU since the test is performed on the CPU
  auto & jFuncMultiplier =
    m_model->getReference< array2d< real64 > >( fields::cappres::jFuncMultiplier::key() );
  jFuncMultiplier.move( hostMemorySpace, false );

  // we are ready to proceed to the test

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.3;
  real64 const end_sat   = 0.9;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat; sat[1] = 1-sat[0];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }

}

TEST_F( CapillaryPressureTest, numericalDerivatives_jFunctionCapPressureThreePhase )
{
  initialize( makeJFunctionCapPressureThreePhase( "capPressure", m_parent ) );

  // here, we have to apply a special treatment to this test
  // to make sure that the J-function multiplier is initialized using initializeRockState
  // this requires calling allocateConstitutiveData in advance (it will be called again later, in the "test" function)

  // setup some values for porosity and permeability
  array2d< real64 > porosity;
  porosity.resize( 1, 1 );
  porosity[0][0] = 0.03496794266569806;
  array3d< real64 > permeability;
  permeability.resize( 1, 1, 3 );
  permeability[0][0][0] = 0.3722194e-15;
  permeability[0][0][1] = 0.4423156e-15;
  permeability[0][0][2] = 0.2324191e-15;

  // initialize the J-function multiplier (done on the GPU if GPU is available)
  m_model->allocateConstitutiveData( m_parent, 1 );
  m_model->initializeRockState( porosity.toViewConst(), permeability.toViewConst() );

  // move the multiplier back to the CPU since the test is performed on the CPU
  auto & jFuncMultiplier =
    m_model->getReference< array2d< real64 > >( fields::cappres::jFuncMultiplier::key() );
  jFuncMultiplier.move( hostMemorySpace, false );

  // we are ready to proceed to the test

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.25;
  real64 const end_sat   = 0.75;
  real64 const dS        = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5*(1-sat[0]);
  sat[2] = 1.0-sat[0]-sat[1];
  while( sat[0] <= end_sat )
  {
    test( sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5*(1-sat[0]);
    sat[2] = 1 - sat[0] - sat[1];
  }
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
