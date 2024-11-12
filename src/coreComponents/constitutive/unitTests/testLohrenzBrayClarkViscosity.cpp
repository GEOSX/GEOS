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
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CompositionalDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/models/LohrenzBrayClarkViscosity.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
using ViscosityData = std::tuple<
  integer const,            // Mixing type
  real64 const,             // pressure
  real64 const,             // temperature
  Feed< NC > const,         // phase composition
  real64 const              // expected viscosity
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 9 >
{
  static std::unique_ptr< TestFluid< 9 > > createFluid()
  {
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::C5, Fluid::C2, Fluid::C3, Fluid::C4, Fluid::C5, Fluid::C10} );
    const std::array< real64 const, 36 > bics = {
      0.01, 0, 0.003732, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.01, 0, 0.028, 0.01, 0.01, 0, 0, 0.01, 0, 0.04532, 0.01, 0.01, 0, 0, 0
    };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< int NC >
class LohrenzBrayClarkViscosityTestFixture :  public ::testing::TestWithParam< ViscosityData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  LohrenzBrayClarkViscosityTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = CompositionalDensity::createParameters( std::make_unique< ModelParameters >() );
    m_parameters = LohrenzBrayClarkViscosity::createParameters( std::move( m_parameters ) );

    auto * parameters = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    parameters->m_componentCriticalVolume.resize( NC );
    TestFluid< 9 >::populateArray( parameters->m_componentCriticalVolume, this->m_fluid->criticalVolume );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    m_density = std::make_unique< CompositionalDensity >( "PhaseDensity", componentProperties, 0, *m_parameters );
    m_viscosity = std::make_unique< LohrenzBrayClarkViscosity >( "PhaseViscosity", componentProperties, 0, *m_parameters );
  }

  ~LohrenzBrayClarkViscosityTestFixture() = default;

  void testViscosity( ViscosityData< NC > const & data )
  {
    auto const mixing_type = static_cast< LohrenzBrayClarkViscosityUpdate::MixingType >(std::get< 0 >( data ));
    real64 const pressure = std::get< 1 >( data );
    real64 const temperature = std::get< 2 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 3 >( data ));
    real64 const expectedViscosity = std::get< 4 >( data );

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    viscosityKernelWrapper.setMixingType( mixing_type );

    densityKernelWrapper.compute( componentProperties,
                                  pressure,
                                  temperature,
                                  phaseComposition.toSliceConst(),
                                  molarDensity,
                                  tempDerivs.toSlice(),
                                  massDensity,
                                  tempDerivs.toSlice(),
                                  false );

    viscosityKernelWrapper.compute( componentProperties,
                                    pressure,
                                    temperature,
                                    phaseComposition.toSliceConst(),
                                    massDensity,
                                    tempDerivs.toSliceConst(),
                                    viscosity,
                                    tempDerivs.toSlice(),
                                    false );

    checkRelativeError( viscosity, expectedViscosity, relTol, absTol );
  }

  void testViscosityDerivatives( ViscosityData< NC > const & data )
  {
    auto const mixing_type = static_cast< LohrenzBrayClarkViscosityUpdate::MixingType >(std::get< 0 >( data ));
    real64 const pressure = std::get< 1 >( data );
    real64 const temperature = std::get< 2 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 3 >( data ));

    auto componentProperties = m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    viscosityKernelWrapper.setMixingType( mixing_type );

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > molarDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > viscosityDerivs( numDofs );

    densityKernelWrapper.compute( componentProperties,
                                  pressure,
                                  temperature,
                                  phaseComposition.toSliceConst(),
                                  molarDensity,
                                  molarDensityDerivs.toSlice(),
                                  massDensity,
                                  massDensityDerivs.toSlice(),
                                  false );

    viscosityKernelWrapper.compute( componentProperties,
                                    pressure,
                                    temperature,
                                    phaseComposition.toSliceConst(),
                                    massDensity,
                                    massDensityDerivs.toSliceConst(),
                                    viscosity,
                                    viscosityDerivs.toSlice(),
                                    false );

    auto calculateViscosity = [&]( real64 const p, real64 const t, auto const & zmf ) -> real64 {
      real64 densityMolar = 0.0;
      real64 densityMass = 0.0;
      real64 phaseViscosity = 0.0;
      stackArray1d< real64, numDofs > tempDerivs( numDofs );
      densityKernelWrapper.compute( componentProperties, p, t, zmf.toSliceConst(),
                                    densityMolar, tempDerivs.toSlice(), densityMass, tempDerivs.toSlice(), false );
      viscosityKernelWrapper.compute( componentProperties, p, t, zmf.toSliceConst(),
                                      densityMass, tempDerivs.toSliceConst(), phaseViscosity, tempDerivs.toSlice(), false );
      return phaseViscosity;
    };

    // Viscosity values are very small so we will inflate the values to avoid false positives due
    // to the absolute value check
    real64 constexpr scale = 1.0e6;

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, scale*viscosityDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return scale*calculateViscosity( p, temperature, phaseComposition );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, scale*viscosityDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return scale*calculateViscosity( pressure, t, phaseComposition );
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < 1; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, scale*viscosityDerivs[Deriv::dC + ic],
        [&]( real64 const z ) -> real64 {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = phaseComposition[jc];
        }
        zmf[ic] += z;
        return scale*calculateViscosity( pressure, temperature, zmf );
      } );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< CompositionalDensity > m_density{};
  std::unique_ptr< LohrenzBrayClarkViscosity > m_viscosity{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using LohrenzBrayClarkViscosity9 = LohrenzBrayClarkViscosityTestFixture< 9 >;

TEST_P( LohrenzBrayClarkViscosity9, testViscosity )
{
  testViscosity( GetParam() );
}

TEST_P( LohrenzBrayClarkViscosity9, testViscosityDerivatives )
{
  testViscosityDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(
  LohrenzBrayClarkViscosity,
  LohrenzBrayClarkViscosity9,
  ::testing::ValuesIn( {
      ViscosityData< 9 >{ 0, 1.839590e+06, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.041140e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+06, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.001152e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+06, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.064158e-04 },
      ViscosityData< 9 >{ 0, 1.839590e+06, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.761220e-05 },
      ViscosityData< 9 >{ 1, 1.839590e+06, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.659626e-05 },
      ViscosityData< 9 >{ 2, 1.839590e+06, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.844604e-05 },
      ViscosityData< 9 >{ 0, 1.839590e+06, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.062011e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+06, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.021978e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+06, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.085108e-04 },
      ViscosityData< 9 >{ 0, 1.839590e+06, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.681288e-05 },
      ViscosityData< 9 >{ 1, 1.839590e+06, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.230528e-05 },
      ViscosityData< 9 >{ 2, 1.839590e+06, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.959149e-05 },
      ViscosityData< 9 >{ 0, 1.839590e+06, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 2.033595e-05 },
      ViscosityData< 9 >{ 1, 1.839590e+06, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.923634e-05 },
      ViscosityData< 9 >{ 2, 1.839590e+06, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 2.133312e-05 },
      ViscosityData< 9 >{ 0, 1.839590e+06, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.681118e-05 },
      ViscosityData< 9 >{ 1, 1.839590e+06, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.229779e-05 },
      ViscosityData< 9 >{ 2, 1.839590e+06, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.959939e-05 },
      ViscosityData< 9 >{ 0, 1.839590e+08, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 6.123583e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+08, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 6.083595e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+08, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 6.146601e-04 },
      ViscosityData< 9 >{ 0, 1.839590e+08, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.451048e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+08, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.440889e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+08, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.459387e-04 },
      ViscosityData< 9 >{ 0, 1.839590e+08, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 6.170310e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+08, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 6.130277e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+08, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 6.193407e-04 },
      ViscosityData< 9 >{ 0, 1.839590e+08, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 4.720985e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+08, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 4.675909e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+08, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 4.748771e-04 },
      ViscosityData< 9 >{ 0, 1.839590e+08, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.078437e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+08, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.067441e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+08, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 1.088409e-04 },
      ViscosityData< 9 >{ 0, 1.839590e+08, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 4.763030e-04 },
      ViscosityData< 9 >{ 1, 1.839590e+08, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 4.717896e-04 },
      ViscosityData< 9 >{ 2, 1.839590e+08, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 4.790912e-04 }
    } )
  );

} // testing

} // geos
