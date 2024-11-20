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
 * @file testMultiFluidCompositionalMultiphase.cpp
 */

#include "MultiFluidTest.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

enum class VISCOSITY_TYPE : int { CONSTANT, LBC };
ENUM_STRINGS( VISCOSITY_TYPE, "Constant", "LohrenzBrayClark" );

template< VISCOSITY_TYPE VISCOSITY >
struct Viscosity {};

template<>
struct Viscosity< VISCOSITY_TYPE::CONSTANT >
{
  using FluidType = CompositionalTwoPhaseConstantViscosity;
};
template<>
struct Viscosity< VISCOSITY_TYPE::LBC >
{
  using FluidType = CompositionalTwoPhaseLohrenzBrayClarkViscosity;
};

template< typename FluidModel, integer NUM_COMP >
struct Fluid
{};

template< EquationOfStateType EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
class MultiFluidCompositionalMultiphaseTest : public MultiFluidTest< typename Viscosity< VISCOSITY >::FluidType, 2, NUM_COMP >
{
public:
  using FluidModel = typename Viscosity< VISCOSITY >::FluidType;
  using Base = MultiFluidTest< FluidModel, 2, NUM_COMP >;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-4;
public:
  MultiFluidCompositionalMultiphaseTest()
  {
    auto & parent = this->m_parent;
    parent.resize( 1 );

    string fluidName = GEOS_FMT( "fluid{}{}{}",
                                 NUM_COMP,
                                 EnumStrings< EquationOfStateType >::toString( EOS ),
                                 EnumStrings< VISCOSITY_TYPE >::toString( VISCOSITY ) );
    this->m_model = makeFluid( fluidName, &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~MultiFluidCompositionalMultiphaseTest() override = default;

  void testNumericalDerivatives( const bool useMass )
  {
    auto & fluid = this->getFluid();
    fluid.setMassFlag( useMass );

    auto * parent = &(this->getParent());

    array2d< real64 > samples;
    Fluid< FluidModel, Base::numComp >::getSamples( samples );
    integer const sampleCount = samples.size( 0 );

    real64 constexpr eps = 1.0e-6;

    constexpr real64 pressures[] = { 1.0e5, 50.0e5, 100.0e5, 600.0e5 };
    constexpr real64 temperatures[] = { 15.5, 24.0, 40.0, 80.0 };

    for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
    {
      for( real64 const pressure : pressures )
      {
        for( real64 const temperature : temperatures )
        {
          typename Base::TestData data ( pressure, units::convertCToK( temperature ), samples[sampleIndex].toSliceConst() );
          Base::testNumericalDerivatives( fluid, parent, data, eps, relTol, absTol );
        }
      }
    }
  }

private:
  static FluidModel * makeFluid( string const & name, Group * parent );
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

template< typename FluidModel >
struct Fluid< FluidModel, 4 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 4 >( componentNames, {"N2", "C10", "C20", "H20"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 4 >( molarWeight, {28e-3, 134e-3, 275e-3, 18e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalPressureString() );
    fill< 4 >( criticalPressure, {34e5, 25.3e5, 14.6e5, 220.5e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 4 >( criticalTemperature, {126.2, 622.0, 782.0, 647.0} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentAcentricFactorString() );
    fill< 4 >( acentricFactor, {0.04, 0.443, 0.816, 0.344} );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( FluidModel::viewKeyStruct::componentBinaryCoeffString() );
    fillBinaryCoeffs< 4 >( binaryCoeff, {0.0, 0.1, 0.0, 0.0, 0.0, 0.0} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 3, 4 );
    fill< 4 >( samples[0], {0.099, 0.300, 0.600, 0.001} );
    fill< 4 >( samples[1], {0.350, 0.350, 0.200, 0.100} );
    fill< 4 >( samples[2], {0.000, 0.000, 0.000, 1.000} );
  }
};

template< typename FluidModel >
struct Fluid< FluidModel, 5 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 5 >( componentNames, {"CO2", "N2", "C1", "C2", "C4"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 5 >( molarWeight, {44.0098e-3, 28.0135e-3, 16.0428e-3, 30.0700e-3, 82.4191e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalPressureString() );
    fill< 5 >( criticalPressure, {73.77300e5, 33.95800e5, 45.99200e5, 48.71800e5, 33.20710e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 5 >( criticalTemperature, {304.1280, 126.1920, 190.5640, 305.3300, 504.2160} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentAcentricFactorString() );
    fill< 5 >( acentricFactor, {0.223000, 0.037200, 0.010400, 0.099100, 0.250274} );
    array1d< real64 > & volumeShift = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentVolumeShiftString() );
    fill< 5 >( volumeShift, {1.845465e-01, -1.283880e-01, 9.225800e-02, 6.458060e-02, 0.000000e+00} );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( FluidModel::viewKeyStruct::componentBinaryCoeffString() );
    fillBinaryCoeffs< 5 >( binaryCoeff, {0.0, 0.1, 0.03, 0.139, 0.032, 0.0, 0.12, 0.03, 0.0, 0.0} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 1, 5 );
    fill< 5 >( samples[0], {0.050, 0.150, 0.550, 0.150, 0.100} );
  }
};
template< EquationOfStateType EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
typename MultiFluidCompositionalMultiphaseTest< EOS, VISCOSITY, NUM_COMP >::FluidModel *
MultiFluidCompositionalMultiphaseTest< EOS, VISCOSITY, NUM_COMP >::
makeFluid( string const & name, Group * parent )
{
  FluidModel & fluid = parent->registerGroup< FluidModel >( name );

  string_array & phaseNames = fluid.template getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  fill< 2 >( phaseNames, {"oil", "gas"} );

  string const eosName = EnumStrings< EquationOfStateType >::toString( EOS );
  string_array & equationOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  fill< 2 >( equationOfState, {eosName, eosName} );

  Fluid< FluidModel, NUM_COMP >::fillProperties( fluid );

  fluid.postInputInitializationRecursive();
  return &fluid;
}

using PengRobinson4Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::CONSTANT, 4 >;
using PengRobinsonLBC4Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::LBC, 4 >;
using SoaveRedlichKwong4Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::CONSTANT, 4 >;
using SoaveRedlichKwongLBC4Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::LBC, 4 >;
TEST_F( PengRobinson4Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( PengRobinsonLBC4Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( SoaveRedlichKwong4Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( SoaveRedlichKwongLBC4Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}

using PengRobinsonLBC5Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::LBC, 5 >;
using SoaveRedlichKwongLBC5Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::LBC, 5 >;
TEST_F( PengRobinsonLBC5Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( SoaveRedlichKwongLBC5Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
