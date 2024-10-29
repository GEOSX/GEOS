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
 * @file testMultiFluidCO2Brine.cpp
 */

#include "MultiFluidTest.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "common/GeosxConfig.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

enum class BrineModelType : int {Phillips, Ezrokhi};
enum class FlashType : int {DuanSun, SpycherPruess};

ENUM_STRINGS( BrineModelType, "Phillips", "Ezrokhi" );
ENUM_STRINGS( FlashType, "DuanSun", "SpycherPruess" );

template< BrineModelType BRINE, bool THERMAL >
struct FluidType {};

template<>
struct FluidType< BrineModelType::Phillips, false >
{
  using type = CO2BrinePhillipsFluid;
  static constexpr const char * brineContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                               "ViscosityFun PhillipsBrineViscosity 0.1";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1";
};
template<>
struct FluidType< BrineModelType::Phillips, true >
{
  using type = CO2BrinePhillipsThermalFluid;
  static constexpr const char * brineContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                               "ViscosityFun PhillipsBrineViscosity 0.1\n"
                                               "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 299.15 369.15 10 0";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "EnthalpyFun CO2Enthalpy 1e6 1.5e7 5e4 367.15 369.15 1";
};
template<>
struct FluidType< BrineModelType::Ezrokhi, false >
{
  using type = CO2BrineEzrokhiFluid;
  static constexpr const char * brineContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                               "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1";
};
template<>
struct FluidType< BrineModelType::Ezrokhi, true >
{
  using type = CO2BrineEzrokhiThermalFluid;
  static constexpr const char * brineContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                               "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4\n"
                                               "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 299.15 369.15 10 0";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "EnthalpyFun CO2Enthalpy 1e6 1.5e7 5e4 367.15 369.15 1";
};

template< FlashType FLASH >
struct FlashModel {};

template<>
struct FlashModel< FlashType::DuanSun >
{
  static constexpr const char * flashContent = "FlashModel CO2Solubility 1e6 1.5e7 5e4 367.15 369.15 1 0.15";
};
template<>
struct FlashModel< FlashType::SpycherPruess >
{
  static constexpr const char * flashContent = "FlashModel CO2Solubility 1e6 1.5e7 5e4 367.15 369.15 1 0.15 1.0e-10 SpycherPruess";
};

template< BrineModelType BRINE, FlashType FLASH, bool THERMAL >
class MultiFluidCO2BrineTest : public MultiFluidTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >,
  public ::testing::WithParamInterface< typename MultiFluidTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >::TestData >
{
public:
  using CO2BrineFluid = typename FluidType< BRINE, THERMAL >::type;
  using Base = MultiFluidTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-4;

public:
  MultiFluidCO2BrineTest()
  {
    Base::writeTableToFile( pvtGasFileName, FluidType< BRINE, THERMAL >::gasContent );
    Base::writeTableToFile( pvtLiquidFileName, FluidType< BRINE, THERMAL >::brineContent );
    Base::writeTableToFile( pvtFlashFileName, FlashModel< FLASH >::flashContent );

    auto & parent = this->m_parent;
    parent.resize( 1 );
    string const fluidName = GEOS_FMT( "fluid{}{}{}",
                                       EnumStrings< BrineModelType >::toString( BRINE ),
                                       EnumStrings< FlashType >::toString( FLASH ),
                                       (THERMAL ? "Thermal" : ""));
    this->m_model = makeCO2BrineFluid( fluidName, &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~MultiFluidCO2BrineTest() override
  {
    Base::removeFile( pvtGasFileName );
    Base::removeFile( pvtLiquidFileName );
    Base::removeFile( pvtFlashFileName );
  }

  // Test numerical derivatives at selected data points
  void testNumericalDerivatives( bool const useMass )
  {
    auto * parent = &(this->getParent());
    auto & fluid = this->getFluid();
    fluid.setMassFlag( useMass );

    real64 constexpr eps = 1.0e-5;

    // Some of the functions are simply table lookups. We need to keep the test points away from
    // the table nodes because the kink in the linear interpolation might cause numerical derivative
    // mismatches. Some of these values have been manually inspected and the differences, although
    // not meeting the tolerance here, are small as expected.
    constexpr real64 temperatures[] = { 367.65, 368.00, 368.75 };
    constexpr real64 pressures[] = { 5.001e6, 7.501e6, 1.201e7 };

    for( real64 const pressure : pressures )
    {
      for( real64 const temperature : temperatures )
      {
        typename Base::TestData data ( pressure, temperature, { 0.3, 0.7 } );
        Base::testNumericalDerivatives( fluid, parent, data, eps, relTol, absTol );
      }
    }
  }

private:
  static CO2BrineFluid * makeCO2BrineFluid( string const & name, Group * parent );
  static constexpr const char * pvtGasFileName = "pvtgas.txt";
  static constexpr const char * pvtLiquidFileName = "pvtliquid.txt";
  static constexpr const char * pvtFlashFileName = "co2flash.txt";
};

template< BrineModelType BRINE, FlashType FLASH, bool THERMAL >
typename MultiFluidCO2BrineTest< BRINE, FLASH, THERMAL >::CO2BrineFluid *
MultiFluidCO2BrineTest< BRINE, FLASH, THERMAL >::makeCO2BrineFluid( string const & name, Group * parent )
{
  CO2BrineFluid & co2BrineFluid = parent->registerGroup< CO2BrineFluid >( name );

  Group & fluid = co2BrineFluid;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  fill< 2 >( phaseNames, {"gas", "liquid"} );

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  fill< 2 >( compNames, {"co2", "water"} );

  auto & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  fill< 2 >( molarWeight, {44e-3, 18e-3} );

  auto & phasePVTParaFileNames = fluid.getReference< path_array >( CO2BrineFluid::viewKeyStruct::phasePVTParaFilesString() );
  fill< 2 >( phasePVTParaFileNames, {pvtGasFileName, pvtLiquidFileName} );

  auto & flashModelParaFileName = fluid.getReference< Path >( CO2BrineFluid::viewKeyStruct::flashModelParaFileString() );
  flashModelParaFileName = pvtFlashFileName;

  co2BrineFluid.postInputInitializationRecursive();

  return &co2BrineFluid;
}

using MultiFluidCO2BrineBrinePhillipsTest = MultiFluidCO2BrineTest< BrineModelType::Phillips,
                                                                    FlashType::DuanSun,
                                                                    false >;
using MultiFluidCO2BrineBrineEzrokhiTest = MultiFluidCO2BrineTest< BrineModelType::Ezrokhi,
                                                                   FlashType::DuanSun,
                                                                   false >;
using MultiFluidCO2BrineBrinePhillipsThermalTest = MultiFluidCO2BrineTest< BrineModelType::Phillips,
                                                                           FlashType::DuanSun,
                                                                           true >;
using MultiFluidCO2BrineBrinePhillipsSpycherPruessTest = MultiFluidCO2BrineTest< BrineModelType::Phillips,
                                                                                 FlashType::SpycherPruess,
                                                                                 false >;
#if !defined(GEOS_DEVICE_COMPILE)
using MultiFluidCO2BrineBrineEzrokhiThermalTest = MultiFluidCO2BrineTest< BrineModelType::Ezrokhi,
                                                                          FlashType::DuanSun,
                                                                          true >;
#endif

TEST_F( MultiFluidCO2BrineBrinePhillipsTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrinePhillipsTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

TEST_F( MultiFluidCO2BrineBrineEzrokhiTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrineEzrokhiTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

TEST_F( MultiFluidCO2BrineBrinePhillipsThermalTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrinePhillipsThermalTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

#if !defined(GEOS_DEVICE_COMPILE)
TEST_F( MultiFluidCO2BrineBrineEzrokhiThermalTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrineEzrokhiThermalTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}
#endif

TEST_F( MultiFluidCO2BrineBrinePhillipsSpycherPruessTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrinePhillipsSpycherPruessTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

TEST_P( MultiFluidCO2BrineBrinePhillipsTest, checkAgainstPreviousImplementationMass )
{
  auto & fluid = dynamicCast< CO2BrineFluid & >( this->getFluid());
  auto fluidWrapper = fluid.createKernelWrapper();
  testValuesAgainstPreviousImplementation( fluidWrapper, GetParam(), relTol );
}

TEST_P( MultiFluidCO2BrineBrineEzrokhiTest, checkAgainstPreviousImplementationMass )
{
  auto & fluid = dynamicCast< CO2BrineFluid & >( this->getFluid());
  auto fluidWrapper = fluid.createKernelWrapper();
  testValuesAgainstPreviousImplementation( fluidWrapper, GetParam(), relTol );
}

TEST_P( MultiFluidCO2BrineBrinePhillipsThermalTest, checkAgainstPreviousImplementationMass )
{
  auto & fluid = dynamicCast< CO2BrineFluid & >( this->getFluid());
  auto fluidWrapper = fluid.createKernelWrapper();
  testValuesAgainstPreviousImplementation( fluidWrapper, GetParam(), relTol );
}

TEST_P( MultiFluidCO2BrineBrinePhillipsSpycherPruessTest, checkAgainstPreviousImplementationMass )
{
  auto & fluid = dynamicCast< CO2BrineFluid & >( this->getFluid());
  auto fluidWrapper = fluid.createKernelWrapper();
  testValuesAgainstPreviousImplementation( fluidWrapper, GetParam(), relTol );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

#define D( ... ) MultiFluidCO2BrineBrinePhillipsTest::TestData{ __VA_ARGS__ }
INSTANTIATE_TEST_SUITE_P(
  MultiFluidCO2Brine, MultiFluidCO2BrineBrinePhillipsTest,
  ::testing::Values(
     //| pressure  | temp       | composition               | density    | phase fraction            |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    |
     D( 5.00100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 5.88317e+03, {2.94136e-01, 7.05864e-01}, {1.87668e+03, 5.32967e+04}, {8.25738e+01, 9.70853e+02}, {1.90427e-05, 3.03214e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 5.00100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 5.87456e+03, {2.94150e-01, 7.05850e-01}, {1.87383e+03, 5.32812e+04}, {8.24487e+01, 9.70542e+02}, {1.90569e-05, 3.02064e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 5.00100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 5.85630e+03, {2.94181e-01, 7.05819e-01}, {1.86780e+03, 5.32478e+04}, {8.21833e+01, 9.69875e+02}, {1.90872e-05, 2.99598e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 9.18273e+03, {2.91939e-01, 7.08061e-01}, {3.05367e+03, 5.32486e+04}, {1.34361e+02, 9.74235e+02}, {2.00622e-05, 3.03214e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 9.16610e+03, {2.91956e-01, 7.08044e-01}, {3.04766e+03, 5.32333e+04}, {1.34097e+02, 9.73923e+02}, {2.00722e-05, 3.02064e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 9.13094e+03, {2.91992e-01, 7.08008e-01}, {3.03496e+03, 5.32005e+04}, {1.33538e+02, 9.73254e+02}, {2.00938e-05, 2.99598e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 1.57730e+04, {2.89059e-01, 7.10941e-01}, {5.77607e+03, 5.32330e+04}, {2.54147e+02, 9.79494e+02}, {2.39022e-05, 3.03214e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 1.57318e+04, {2.89075e-01, 7.10925e-01}, {5.75768e+03, 5.32179e+04}, {2.53338e+02, 9.79185e+02}, {2.38854e-05, 3.02064e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 1.56452e+04, {2.89109e-01, 7.10891e-01}, {5.71917e+03, 5.31853e+04}, {2.51643e+02, 9.78520e+02}, {2.38514e-05, 2.99598e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} )
  )
);
#undef D

#define D( ... ) MultiFluidCO2BrineBrineEzrokhiTest::TestData{ __VA_ARGS__ }
INSTANTIATE_TEST_SUITE_P(
  MultiFluidCO2Brine, MultiFluidCO2BrineBrineEzrokhiTest,
  ::testing::Values(
     //| pressure  | temp       | composition               | density    | phase fraction            |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    |
     D( 5.00100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 5.89876e+03, {2.94136e-01, 7.05864e-01}, {1.87668e+03, 5.51674e+04}, {8.25738e+01, 1.00493e+03}, {1.90427e-05, 3.12148e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 5.00100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 5.89023e+03, {2.94150e-01, 7.05850e-01}, {1.87383e+03, 5.51663e+04}, {8.24487e+01, 1.00488e+03}, {1.90569e-05, 3.11023e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 5.00100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 5.87213e+03, {2.94181e-01, 7.05819e-01}, {1.86780e+03, 5.51642e+04}, {8.21833e+01, 1.00478e+03}, {1.90872e-05, 3.08611e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 9.23469e+03, {2.91939e-01, 7.08061e-01}, {3.05367e+03, 5.58209e+04}, {1.34361e+02, 1.02130e+03}, {2.00622e-05, 3.16876e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 9.21829e+03, {2.91956e-01, 7.08044e-01}, {3.04766e+03, 5.58254e+04}, {1.34097e+02, 1.02135e+03}, {2.00722e-05, 3.15764e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 9.18360e+03, {2.91992e-01, 7.08008e-01}, {3.03496e+03, 5.58353e+04}, {1.33538e+02, 1.02146e+03}, {2.00938e-05, 3.13381e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 1.59791e+04, {2.89059e-01, 7.10941e-01}, {5.77607e+03, 5.67057e+04}, {2.54147e+02, 1.04339e+03}, {2.39022e-05, 3.23075e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 1.59385e+04, {2.89075e-01, 7.10925e-01}, {5.75768e+03, 5.67188e+04}, {2.53338e+02, 1.04360e+03}, {2.38854e-05, 3.21991e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 1.58533e+04, {2.89109e-01, 7.10891e-01}, {5.71917e+03, 5.67472e+04}, {2.51643e+02, 1.04405e+03}, {2.38514e-05, 3.19665e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} )
  )
);
#undef D

#define D( ... ) MultiFluidCO2BrineBrinePhillipsThermalTest::TestData{ __VA_ARGS__ }
INSTANTIATE_TEST_SUITE_P(
  MultiFluidCO2Brine, MultiFluidCO2BrineBrinePhillipsThermalTest,
  ::testing::Values(
     //| pressure  | temp       | composition               | density    | phase fraction            |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    |
     D( 5.00100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 5.88317e+03, {2.94136e-01, 7.05864e-01}, {1.87668e+03, 5.32967e+04}, {8.25738e+01, 9.70853e+02}, {1.90427e-05, 3.03214e-04}, {1.21447e+07, 2.18796e+07}, {1.20841e+07, 2.18744e+07} ),
     D( 5.00100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 5.87456e+03, {2.94150e-01, 7.05850e-01}, {1.87383e+03, 5.32812e+04}, {8.24487e+01, 9.70542e+02}, {1.90569e-05, 3.02064e-04}, {1.21537e+07, 2.19628e+07}, {1.20931e+07, 2.19576e+07} ),
     D( 5.00100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 5.85630e+03, {2.94181e-01, 7.05819e-01}, {1.86780e+03, 5.32478e+04}, {8.21833e+01, 9.69875e+02}, {1.90872e-05, 2.99598e-04}, {1.21731e+07, 2.21410e+07}, {1.21123e+07, 2.21359e+07} ),
     D( 7.50100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 9.18273e+03, {2.91939e-01, 7.08061e-01}, {3.05367e+03, 5.32486e+04}, {1.34361e+02, 9.74235e+02}, {2.00622e-05, 3.03214e-04}, {1.17163e+07, 2.18445e+07}, {1.16605e+07, 2.18368e+07} ),
     D( 7.50100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 9.16610e+03, {2.91956e-01, 7.08044e-01}, {3.04766e+03, 5.32333e+04}, {1.34097e+02, 9.73923e+02}, {2.00722e-05, 3.02064e-04}, {1.17269e+07, 2.19275e+07}, {1.16709e+07, 2.19198e+07} ),
     D( 7.50100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 9.13094e+03, {2.91992e-01, 7.08008e-01}, {3.03496e+03, 5.32005e+04}, {1.33538e+02, 9.73254e+02}, {2.00938e-05, 2.99598e-04}, {1.17494e+07, 2.21054e+07}, {1.16932e+07, 2.20977e+07} ),
     D( 1.20100e+07, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 1.57730e+04, {2.89059e-01, 7.10941e-01}, {5.77607e+03, 5.32330e+04}, {2.54147e+02, 9.79494e+02}, {2.39022e-05, 3.03214e-04}, {1.08306e+07, 2.17898e+07}, {1.07833e+07, 2.17775e+07} ),
     D( 1.20100e+07, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 1.57318e+04, {2.89075e-01, 7.10925e-01}, {5.75768e+03, 5.32179e+04}, {2.53338e+02, 9.79185e+02}, {2.38854e-05, 3.02064e-04}, {1.08453e+07, 2.18726e+07}, {1.07979e+07, 2.18603e+07} ),
     D( 1.20100e+07, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 1.56452e+04, {2.89109e-01, 7.10891e-01}, {5.71917e+03, 5.31853e+04}, {2.51643e+02, 9.78520e+02}, {2.38514e-05, 2.99598e-04}, {1.08767e+07, 2.20500e+07}, {1.08290e+07, 2.20378e+07} )
  )
);
#undef D

#define D( ... ) MultiFluidCO2BrineBrinePhillipsSpycherPruessTest::TestData{ __VA_ARGS__ }
INSTANTIATE_TEST_SUITE_P(
  MultiFluidCO2Brine, MultiFluidCO2BrineBrinePhillipsSpycherPruessTest,
  ::testing::Values(
     //| pressure  | temp       | composition               | density    | phase fraction            |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    |
     D( 5.00100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 5.77217e+03, {3.00488e-01, 6.99512e-01}, {1.87668e+03, 5.32842e+04}, {8.25738e+01, 9.70985e+02}, {1.90427e-05, 3.03214e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 5.00100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 5.76239e+03, {3.00580e-01, 6.99420e-01}, {1.87383e+03, 5.32686e+04}, {8.24487e+01, 9.70677e+02}, {1.90569e-05, 3.02064e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 5.00100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 5.74154e+03, {3.00781e-01, 6.99219e-01}, {1.86780e+03, 5.32348e+04}, {8.21833e+01, 9.70013e+02}, {1.90872e-05, 2.99598e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 9.05928e+03, {2.96731e-01, 7.03269e-01}, {3.05367e+03, 5.32295e+04}, {1.34361e+02, 9.74440e+02}, {2.00622e-05, 3.03214e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 9.04143e+03, {2.96804e-01, 7.03196e-01}, {3.04766e+03, 5.32141e+04}, {1.34097e+02, 9.74130e+02}, {2.00722e-05, 3.02064e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 7.50100e+06, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 9.00358e+03, {2.96962e-01, 7.03038e-01}, {3.03496e+03, 5.31808e+04}, {1.33538e+02, 9.73465e+02}, {2.00938e-05, 2.99598e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.67650e+02, {3.00000e-01, 7.00000e-01}, 1.56195e+04, {2.93050e-01, 7.06950e-01}, {5.77607e+03, 5.32048e+04}, {2.54147e+02, 9.79798e+02}, {2.39022e-05, 3.03214e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.68000e+02, {3.00000e-01, 7.00000e-01}, 1.55771e+04, {2.93105e-01, 7.06895e-01}, {5.75768e+03, 5.31894e+04}, {2.53338e+02, 9.79492e+02}, {2.38854e-05, 3.02064e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} ),
     D( 1.20100e+07, 3.68750e+02, {3.00000e-01, 7.00000e-01}, 1.54878e+04, {2.93225e-01, 7.06775e-01}, {5.71917e+03, 5.31561e+04}, {2.51643e+02, 9.78834e+02}, {2.38514e-05, 2.99598e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00} )
  )
);
#undef D

/* UNCRUSTIFY-ON */

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
