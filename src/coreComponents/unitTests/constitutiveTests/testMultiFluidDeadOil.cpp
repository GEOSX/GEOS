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
 * @file testMultiFluidDeadOil.cpp
 */

#include "MultiFluidTest.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

static constexpr char const * pvdgTableContent = "# Pg(Pa) Bg(m3/sm3) Visc(Pa.s)\n"
                                                 "3000000  0.04234  0.00001344\n"
                                                 "6000000  0.02046  0.0000142\n"
                                                 "9000000  0.01328  0.00001526\n"
                                                 "12000000 0.00977  0.0000166\n"
                                                 "15000000 0.00773  0.00001818\n"
                                                 "18000000 0.006426 0.00001994\n"
                                                 "21000000 0.005541 0.00002181\n"
                                                 "24000000 0.004919 0.0000237\n"
                                                 "27000000 0.004471 0.00002559 -- this is a comment\n"
                                                 "29500000 0.004194 0.00002714\n"
                                                 "31000000 0.004031 0.00002806\n"
                                                 "33000000 0.00391  0.00002832\n"
                                                 "53000000 0.003868 0.00002935";

static constexpr char const * pvdoTableContent = "#P[Pa] Bo[m3/sm3] Visc(Pa.s)\n"
                                                 "10000000.0 1.23331 0.00015674\n"
                                                 "12500000.0 1.21987 0.00016570\n"
                                                 "15000000.0 1.20802 0.00017445\n"
                                                 "20000000.0 1.18791 0.00019143\n"
                                                 "25000000.0 1.17137 0.00020779\n"
                                                 "30000000.0 1.15742 0.00022361\n"
                                                 "33200000.3 1.14946 0.00023359\n"
                                                 "35000000.0 1.14543 0.00023894\n"
                                                 "40000000.0 1.13498 0.00025383 -- this is a comment\n"
                                                 "50000000.0 1.11753 0.00028237\n"
                                                 "60000000.0 1.10346 0.00030941\n"
                                                 "70000000.0 1.09180 0.00033506\n"
                                                 "80000000.0 1.08194 0.00035945\n"
                                                 "90000000.0 1.07347 0.00038266\n"
                                                 "95000000.0 1.06966 0.00039384\n"
                                                 "100000000.0 1.06610 0.00040476\n"
                                                 "110000000.0 1.05961 0.00042584\n"
                                                 "112500000.0 1.05811 0.00043096\n"
                                                 "115000000.0 1.05665 0.00043602\n"
                                                 "117500000.0 1.05523 0.00044102\n"
                                                 "120000000.0 1.05385 0.00044596\n";

static constexpr char const * pvdwTableContent = "# Pref[Pa] Bw[m3/sm3] Cp[1/Pa]     Visc[Pa.s]\n"
                                                 " 30600000.1 1.03  0.00000000041 0.0003";

template< bool FROM_TABLE >
class MultiFluidDeadOilTest : public MultiFluidTest< DeadOilFluid, 3, 3 >
{
public:
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-4;
public:
  MultiFluidDeadOilTest()
  {
    if constexpr (!FROM_TABLE)
    {
      writeTableToFile( "pvdo.txt", pvdoTableContent );
      writeTableToFile( "pvdg.txt", pvdgTableContent );
      writeTableToFile( "pvdw.txt", pvdwTableContent );
    }

    m_parent.resize( 1 );
    string const fluidName = GEOS_FMT( "fluid{}", (FROM_TABLE ? "Tables" : "Files"));
    m_model = makeDeadOilFluid( fluidName, &m_parent );

    m_parent.initialize();
    m_parent.initializePostInitialConditions();
  }

  ~MultiFluidDeadOilTest() override
  {
    if constexpr (!FROM_TABLE)
    {
      removeFile( "pvdo.txt" );
      removeFile( "pvdg.txt" );
      removeFile( "pvdw.txt" );
    }
  }

private:
  static DeadOilFluid * makeDeadOilFluid( string const & name, Group * parent );
  static void fillPhysicalProperties( DeadOilFluid & fluid );
};

template< bool FROM_TABLE >
void MultiFluidDeadOilTest< FROM_TABLE >::fillPhysicalProperties( DeadOilFluid & fluid )
{
  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  fill< 3 >( phaseNames, {"oil", "water", "gas"} );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  fill< 3 >( compNames, {"oil", "water", "gas"} );

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  fill< 3 >( molarWgt, {114e-3, 18e-3, 16e-3} );

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluidBase::viewKeyStruct::surfacePhaseMassDensitiesString() );
  fill< 3 >( surfaceDens, {800.0, 1022.0, 0.9907} );
}

template<>
DeadOilFluid * MultiFluidDeadOilTest< false >::makeDeadOilFluid( string const & name, Group * parent )
{
  DeadOilFluid & fluid = parent->registerGroup< DeadOilFluid >( name );

  fillPhysicalProperties( fluid );

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluidBase::viewKeyStruct::tableFilesString() );
  fill< 3 >( tableNames, {"pvdo.txt", "pvdw.txt", "pvdg.txt"} );

  fluid.postInputInitializationRecursive();
  return &fluid;
}

template<>
DeadOilFluid * MultiFluidDeadOilTest< true >::makeDeadOilFluid( string const & name, Group * parent )
{
  // 1) First, define the tables (PVDO, PVDG)

  // 1D table with linear interpolation
  integer constexpr NaxisPVDO = 21;
  integer constexpr NaxisPVDG = 13;

  array1d< real64_array > coordinatesPVDO( 1 );
  real64_array valuesPVDO_Bo;
  real64_array valuesPVDO_visc;
  fill< NaxisPVDO >( coordinatesPVDO[0], {
    1.000e+07, 1.250e+07, 1.500e+07, 2.000e+07, 2.500e+07, 3.000e+07, 3.320e+07, 3.500e+07,
    4.000e+07, 5.000e+07, 6.000e+07, 7.000e+07, 8.000e+07, 9.000e+07, 9.500e+07, 1.000e+08,
    1.100e+08, 1.125e+08, 1.150e+08, 1.175e+08, 1.200e+08 } );
  fill< NaxisPVDO >( valuesPVDO_Bo, {
    1.23331, 1.21987, 1.20802, 1.18791, 1.17137, 1.15742, 1.14946, 1.14543,
    1.13498, 1.11753, 1.10346, 1.09180, 1.08194, 1.07347, 1.06966, 1.06610,
    1.05961, 1.05811, 1.05665, 1.05523, 1.05385 } );
  fill< NaxisPVDO >( valuesPVDO_visc, {
    1.56740e-04, 1.65700e-04, 1.74450e-04, 1.91430e-04, 2.07790e-04, 2.23610e-04, 2.33590e-04, 2.38940e-04,
    2.53830e-04, 2.82370e-04, 3.09410e-04, 3.35060e-04, 3.59450e-04, 3.82660e-04, 3.93840e-04, 4.04760e-04,
    4.25840e-04, 4.30960e-04, 4.36020e-04, 4.41020e-04, 4.45960e-04 } );

  array1d< real64_array > coordinatesPVDG( 1 );
  real64_array valuesPVDG_Bg;
  real64_array valuesPVDG_visc;
  fill< NaxisPVDG >( coordinatesPVDG[0], {
    3.000e+06, 6.000e+06, 9.000e+06, 1.200e+07, 1.500e+07, 1.800e+07, 2.100e+07, 2.400e+07,
    2.700e+07, 2.950e+07, 3.100e+07, 3.300e+07, 5.300e+07 } );
  fill< NaxisPVDG >( valuesPVDG_Bg, {
    4.23400e-02, 2.04600e-02, 1.32800e-02, 9.77000e-03, 7.73000e-03, 6.42600e-03, 5.54100e-03, 4.91900e-03,
    4.47100e-03, 4.19400e-03, 4.03100e-03, 3.91000e-03, 3.86800e-03 } );
  fill< NaxisPVDG >( valuesPVDG_visc, {
    1.34400e-05, 1.42000e-05, 1.52600e-05, 1.66000e-05, 1.81800e-05, 1.99400e-05, 2.18100e-05, 2.37000e-05,
    2.55900e-05, 2.71400e-05, 2.80600e-05, 2.83200e-05, 2.93500e-05 } );

  initializeTable( "PVDO_Bo", coordinatesPVDO, valuesPVDO_Bo );
  initializeTable( "PVDO_visc", coordinatesPVDO, valuesPVDO_visc );
  initializeTable( "PVDG_Bg", coordinatesPVDG, valuesPVDG_Bg );
  initializeTable( "PVDG_visc", coordinatesPVDG, valuesPVDG_visc );

  // 2) Then, define the Dead-Oil constitutive model

  DeadOilFluid & fluid = parent->registerGroup< DeadOilFluid >( name );

  fillPhysicalProperties( fluid );

  string_array & FVFTableNames = fluid.getReference< string_array >( DeadOilFluid::viewKeyStruct::formationVolumeFactorTableNamesString() );
  fill< 2 >( FVFTableNames, {"PVDG_Bg", "PVDO_Bo"} );

  string_array & viscosityTableNames = fluid.getReference< string_array >( DeadOilFluid::viewKeyStruct::viscosityTableNamesString() );
  fill< 2 >( viscosityTableNames, {"PVDG_visc", "PVDO_visc"} );

  real64 & waterRefPressure = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterRefPressureString() );
  waterRefPressure = 30600000.1;
  real64 & waterFormationVolumeFactor = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterFormationVolumeFactorString() );
  waterFormationVolumeFactor = 1.03;
  real64 & waterCompressibility = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterCompressibilityString() );
  waterCompressibility = 0.00000000041;
  real64 & waterViscosity = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterViscosityString() );
  waterViscosity = 0.0003;

  fluid.postInputInitializationRecursive();
  return &fluid;
}

using MultiFluidDeadOilTestFromFiles = MultiFluidDeadOilTest< false >;
using MultiFluidDeadOilTestFromTables = MultiFluidDeadOilTest< true >;

TEST_F( MultiFluidDeadOilTestFromFiles, numericalDerivativesMolar )
{
  auto & fluid = getFluid();
  fluid.setMassFlag( false );

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());

  for( real64 const pressure : { 1.24e7, 3.21e7, 5.01e7 } )
  {
    TestData data ( pressure, 297.15, { 0.1, 0.3, 0.6 } );
    testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
  }
}

TEST_F( MultiFluidDeadOilTestFromFiles, numericalDerivativesMass )
{
  auto & fluid = getFluid();
  fluid.setMassFlag( true );

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());

  for( real64 const pressure : { 1.24e7, 3.21e7, 5.01e7 } )
  {
    TestData data ( pressure, 297.15, { 0.1, 0.3, 0.6 } );
    testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
  }
}

TEST_F( MultiFluidDeadOilTestFromTables, numericalDerivativesMolar )
{
  auto & fluid = getFluid();
  fluid.setMassFlag( false );

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());

  for( real64 const pressure : { 1.24e7, 3.21e7, 5.01e7 } )
  {
    TestData data ( pressure, 297.15, { 0.1, 0.3, 0.6 } );
    testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
  }
}

TEST_F( MultiFluidDeadOilTestFromTables, numericalDerivativesMass )
{
  auto & fluid = getFluid();
  fluid.setMassFlag( true );

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());

  for( real64 const pressure : { 1.24e7, 3.21e7, 5.01e7 } )
  {
    TestData data ( pressure, 297.15, { 0.1, 0.3, 0.6 } );
    testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
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
