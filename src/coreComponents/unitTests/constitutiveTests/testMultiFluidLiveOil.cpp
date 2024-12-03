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
 * @file testMultiFluidLiveOil.cpp
 */

#include "MultiFluidTest.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::dataRepository;
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

static const char * pvtoTableContent = "# Rs[sm3/sm3]\tPbub[Pa]\tBo[m3/sm3]\tVisc(Pa.s)\n"
                                       "\n"
                                       "  2\t            2000000\t    1.02\t    0.000975\n"
                                       "  5\t            5000000\t    1.03\t    0.00091\n"
                                       " 10\t            10000000\t1.04\t    0.00083\n"
                                       " 15\t            20000000\t1.05\t    0.000695\n"
                                       "                90000000\t1.03\t    0.000985  -- some line comment\n"
                                       " 30\t            30000000\t1.07\t    0.000594\n"
                                       " 40\t            40000000\t1.08\t    0.00051\n"
                                       "                50000000\t1.07\t    0.000549  -- another one\n"
                                       "                90000000\t1.06\t    0.00074\n"
                                       " 50\t            50000000.7\t1.09\t    0.000449\n"
                                       "                90000000.7\t1.08\t    0.000605";

static const char * pvtwTableContent = "#\tPref[Pa]\tBw[m3/sm3]\tCp[1/Pa]\t    Visc[Pa.s]\n"
                                       "\t30600000.1\t1.03\t\t0.00000000041\t0.0003";

class MultiFluidLiveOilTest : public MultiFluidTest< BlackOilFluid, 3, 3 >
{
public:
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-4;
public:
  MultiFluidLiveOilTest()
  {
    writeTableToFile( pvtoFileName, pvtoTableContent );
    writeTableToFile( pvdgFileName, pvdgTableContent );
    writeTableToFile( pvtwFileName, pvtwTableContent );

    m_parent.resize( 1 );
    m_model = makeLiveOilFluid( "fluid", &m_parent );

    m_parent.initialize();
    m_parent.initializePostInitialConditions();
  }

  ~MultiFluidLiveOilTest() override
  {
    removeFile( pvtoFileName );
    removeFile( pvdgFileName );
    removeFile( pvtwFileName );
  }

private:
  static BlackOilFluid * makeLiveOilFluid( string const & name, Group * parent );
  static constexpr const char * pvtoFileName = "pvto.txt";
  static constexpr const char * pvdgFileName = "pvdg.txt";
  static constexpr const char * pvtwFileName = "pvtw.txt";
};

BlackOilFluid * MultiFluidLiveOilTest::makeLiveOilFluid( string const & name, Group * parent )
{
  BlackOilFluid & fluid = parent->registerGroup< BlackOilFluid >( name );

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  fill< 3 >( phaseNames, {"oil", "gas", "water"} );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  fill< 3 >( compNames, {"oil", "gas", "water"} );

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  fill< 3 >( molarWgt, {114e-3, 16e-3, 18e-3} );

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluidBase::viewKeyStruct::surfacePhaseMassDensitiesString() );
  fill< 3 >( surfaceDens, {800.0, 0.9907, 1022.0} );

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluidBase::viewKeyStruct::tableFilesString() );
  fill< 3 >( tableNames, {pvtoFileName, pvdgFileName, pvtwFileName} );

  fluid.postInputInitializationRecursive();
  return &fluid;
}

TEST_F( MultiFluidLiveOilTest, numericalDerivativesMolar )
{
  auto & fluid = getFluid();
  fluid.setMassFlag( false );

  real64 constexpr eps = 1.0e-6;

  array2d< real64 > samples( 2, 3 );
  fill< 3 >( samples[0], { 0.10000, 0.3, 0.60000 } );
  fill< 3 >( samples[1], { 0.79999, 0.2, 0.00001 } );

  for( real64 const pressure : { 1.24e7, 3.21e7, 5.01e7 } )
  {
    for( integer sampleIndex = 0; sampleIndex < samples.size( 0 ); sampleIndex++ )
    {
      TestData data ( pressure, 297.15, samples[sampleIndex] );
      testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
    }
  }
}

TEST_F( MultiFluidLiveOilTest, numericalDerivativesMass )
{
  auto & fluid = getFluid();
  fluid.setMassFlag( true );

  real64 constexpr eps = 1.0e-6;

  array2d< real64 > samples( 2, 3 );
  fill< 3 >( samples[0], { 0.10000, 0.3, 0.60000 } );
  fill< 3 >( samples[1], { 0.79999, 0.2, 0.00001 } );

  for( real64 const pressure : { 1.24e7, 3.21e7, 5.01e7 } )
  {
    for( integer sampleIndex = 0; sampleIndex < samples.size( 0 ); sampleIndex++ )
    {
      TestData data ( pressure, 297.15, samples[sampleIndex] );
      testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
    }
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
