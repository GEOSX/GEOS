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
 * @file testMultiFluidCompositionalMultiphasePVTPackage.cpp
 */

#include "MultiFluidTest.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

namespace geos
{
namespace testing
{

template< >
struct UsePVTPackage< CompositionalMultiphaseFluidPVTPackage >
{
  static constexpr bool value = true;
};

} // namespace testing

} // namespace geos

enum class EOS_TYPE : int { PR, SRK };
ENUM_STRINGS( EOS_TYPE, "PR", "SRK" );

template< integer NUM_COMP >
struct Fluid
{};

template< EOS_TYPE EOS, integer NUM_COMP >
class MultiFluidCompositionalMultiphasePVTPackageTest : public MultiFluidTest< CompositionalMultiphaseFluidPVTPackage, 2, NUM_COMP >
{
public:
  using Base = MultiFluidTest< CompositionalMultiphaseFluidPVTPackage, 2, NUM_COMP >;
  static constexpr real64 relTol = 1.0e-4;
public:
  MultiFluidCompositionalMultiphasePVTPackageTest()
  {
    auto & parent = this->m_parent;
    parent.resize( 1 );

    string fluidName = GEOS_FMT( "fluid{}{}", EnumStrings< EOS_TYPE >::toString( EOS ), NUM_COMP );
    this->m_model = makeFluid( fluidName, &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~MultiFluidCompositionalMultiphasePVTPackageTest() override = default;

  void testNumericatDerivatives( const bool useMass )
  {
    auto & fluid = this->getFluid();
    fluid.setMassFlag( useMass );

    auto * parent = &(this->getParent());

    array2d< real64 > samples;
    Fluid< Base::numComp >::getSamples( samples );
    integer const sampleCount = samples.size( 0 );

    real64 constexpr eps = 1.0e-7;

    constexpr real64 pressures[] = { 1.0e5, 50.0e5, 100.0e5, 600.0e5 };
    constexpr real64 temperatures[] = { 15.5, 24.0, 40.0, 80.0 };

    for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
    {
      for( real64 const pressure : pressures )
      {
        for( real64 const temperature : temperatures )
        {
          typename Base::TestData data ( pressure, units::convertCToK( temperature ), samples[sampleIndex].toSliceConst() );
          Base::testNumericalDerivatives( fluid, parent, data, eps, relTol );
        }
      }
    }
  }

private:
  static CompositionalMultiphaseFluidPVTPackage * makeFluid( string const & name, Group * parent );
};

template< integer NUM_COMP >
static void fillBinaryCoeffs( array2d< real64 > & binaryCoeff, std::array< real64 const, NUM_COMP *(NUM_COMP-1)/2 > const data )
{
  auto bic = data.begin();
  binaryCoeff.resize( NUM_COMP, NUM_COMP );
  for( integer i = 0; i < NUM_COMP; ++i )
  {
    binaryCoeff( i, i ) = 0.0;
    for( integer j = i+1; j < NUM_COMP; ++j )
    {
      binaryCoeff( i, j ) = *bic++;
      binaryCoeff( j, i ) = binaryCoeff( i, j );
    }
  }
}

template< >
struct Fluid< 4 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 4 >( componentNames, {"N2", "C10", "C20", "H20"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 4 >( molarWeight, {28e-3, 134e-3, 275e-3, 18e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentCriticalPressureString() );
    fill< 4 >( criticalPressure, {34e5, 25.3e5, 14.6e5, 220.5e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 4 >( criticalTemperature, {126.2, 622.0, 782.0, 647.0} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentAcentricFactorString() );
    fill< 4 >( acentricFactor, {0.04, 0.443, 0.816, 0.344} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 2, 4 );
    fill< 4 >( samples[0], {0.099, 0.300, 0.600, 0.001} );
    fill< 4 >( samples[1], {0.000, 0.000, 0.000, 1.000} );
  }
};

template< EOS_TYPE EOS, integer NUM_COMP >
CompositionalMultiphaseFluidPVTPackage *
MultiFluidCompositionalMultiphasePVTPackageTest< EOS, NUM_COMP >::
makeFluid( string const & name, Group * parent )
{
  CompositionalMultiphaseFluidPVTPackage & fluid = parent->registerGroup< CompositionalMultiphaseFluidPVTPackage >( name );

  string_array & phaseNames = fluid.getReference< string_array >( string( MultiFluidBase::viewKeyStruct::phaseNamesString()) );
  fill< 2 >( phaseNames, {"oil", "gas"} );

  string const eosName = EnumStrings< EOS_TYPE >::toString( EOS );
  string_array & equationOfState = fluid.getReference< string_array >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::equationsOfStateString() );
  fill< 2 >( equationOfState, {eosName, eosName} );

  Fluid< NUM_COMP >::fillProperties( fluid );

  fluid.postInputInitializationRecursive();
  return &fluid;
}

using PengRobinson4Test = MultiFluidCompositionalMultiphasePVTPackageTest< EOS_TYPE::PR, 4 >;
using SoaveRedlichKwong4Test = MultiFluidCompositionalMultiphasePVTPackageTest< EOS_TYPE::SRK, 4 >;

TEST_F( PengRobinson4Test, numericalDerivativesMolar )
{
  testNumericatDerivatives( false );
}
TEST_F( PengRobinson4Test, numericalDerivativesMass )
{
  testNumericatDerivatives( true );
}

TEST_F( SoaveRedlichKwong4Test, numericalDerivativesMolar )
{
  testNumericatDerivatives( false );
}
TEST_F( SoaveRedlichKwong4Test, numericalDerivativesMass )
{
  testNumericatDerivatives( true );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
