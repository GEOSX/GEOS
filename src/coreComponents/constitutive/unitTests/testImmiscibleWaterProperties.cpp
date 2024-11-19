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
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterViscosity.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
using TestData = std::tuple<
  real64 const,             // pressure
  real64 const,             // temperature
  Feed< NC > const,         // phase composition
  real64 const,             // expected molar density
  real64 const,             // expected mass density
  real64 const              // expected viscosity
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 3 >
{
  static std::unique_ptr< TestFluid< 3 > > createFluid()
  {
    auto fluid = TestFluid< 3 >::create( {Fluid::C1, Fluid::C10, Fluid::H2O} );
    std::array< real64, 3 > const bics = {0.25, 0.0, 0.0};
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< int NC >
class ImmiscibleWaterPropertiesTestFixture :  public ::testing::TestWithParam< TestData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  ImmiscibleWaterPropertiesTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = ImmiscibleWaterDensity::createParameters( std::make_unique< ModelParameters >() );
    m_parameters = ImmiscibleWaterViscosity::createParameters( std::move( m_parameters ) );

    auto * waterParameters = const_cast< ImmiscibleWaterParameters * >(m_parameters->get< ImmiscibleWaterParameters >());
    waterParameters->m_waterReferencePressure = 215.0e5;
    waterParameters->m_waterDensity = 1020.0;
    waterParameters->m_waterCompressibility = 4.1483E-10;
    waterParameters->m_waterViscosityCompressibility = 2.0E-11;
    waterParameters->m_waterViscosity = 0.32929e-3;

    m_density = std::make_unique< ImmiscibleWaterDensity >( "PhaseDensity", componentProperties, 0, *m_parameters );
    m_viscosity = std::make_unique< ImmiscibleWaterViscosity >( "PhaseViscosity", componentProperties, 0, *m_parameters );
  }

  ~ImmiscibleWaterPropertiesTestFixture() = default;

  void testProperties( TestData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));
    real64 const expectedMolarDensity = std::get< 3 >( data );
    real64 const expectedMassDensity = std::get< 4 >( data );
    real64 const expectedViscosity = std::get< 5 >( data );

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

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

    checkRelativeError( molarDensity, expectedMolarDensity, relTol, absTol );
    checkRelativeError( massDensity, expectedMassDensity, relTol, absTol );
    checkRelativeError( viscosity, expectedViscosity, relTol, absTol );
  }

  void testPropertyDerivatives( TestData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));

    auto componentProperties = m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > molarDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > viscosityDerivs( numDofs );
    stackArray1d< real64, 3 > derivatives( 3 );

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

    // Viscosity values are very small so we will inflate the values to avoid false positives due
    // to the absolute value check
    real64 constexpr viscosityScale = 1.0e6;

    auto calculateProperties = [&]( real64 const p, real64 const t, auto const & zmf, auto & values ) {
      stackArray1d< real64, numDofs > tempDerivs( numDofs );
      densityKernelWrapper.compute( componentProperties, p, t, zmf.toSliceConst(),
                                    values[0], tempDerivs.toSlice(), values[1], tempDerivs.toSlice(), false );
      viscosityKernelWrapper.compute( componentProperties, p, t, zmf.toSliceConst(),
                                      values[1], tempDerivs.toSliceConst(), values[2], tempDerivs.toSlice(), false );
      values[2] *= viscosityScale;
    };

    auto concatDerivatives = [&]( int idof ){
      derivatives[0] = molarDensityDerivs[idof];
      derivatives[1] = massDensityDerivs[idof];
      derivatives[2] = viscosityScale * viscosityDerivs[idof];
    };

    // Compare against numerical derivatives
    // -- Pressure derivative
    concatDerivatives( Deriv::dP );
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative< 3 >( pressure, dp, derivatives.toSliceConst(),
                                            [&]( real64 const p, auto & values ) {
      calculateProperties( p, temperature, phaseComposition, values );
    }, absTol, relTol );

    // -- Temperature derivative
    concatDerivatives( Deriv::dT );
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative< 3 >( temperature, dT, derivatives.toSliceConst(),
                                            [&]( real64 const t, auto & values ) {
      calculateProperties( pressure, t, phaseComposition, values );
    }, absTol, relTol );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < NC; ++ic )
    {
      concatDerivatives( Deriv::dC+ic );
      internal::testNumericalDerivative< 3 >( 0.0, dz, derivatives.toSliceConst(),
                                              [&]( real64 const z, auto & values ) {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = phaseComposition[jc];
        }
        zmf[ic] += z;
        calculateProperties( pressure, temperature, zmf, values );
      }, absTol, relTol );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< ImmiscibleWaterDensity > m_density{};
  std::unique_ptr< ImmiscibleWaterViscosity > m_viscosity{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using ImmiscibleWaterProperties3 = ImmiscibleWaterPropertiesTestFixture< 3 >;

TEST_P( ImmiscibleWaterProperties3, testProperties )
{
  testProperties( GetParam() );
}

TEST_P( ImmiscibleWaterProperties3, testPropertyDerivatives )
{
  testPropertyDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterProperties, ImmiscibleWaterProperties3,
  ::testing::Values( 
    TestData< 3 >( 2.0e+06, 293.15, { 0.30, 0.30, 0.40 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 293.15, { 0.30, 0.30, 0.40 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 293.15, { 0.30, 0.30, 0.40 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 353.15, { 0.30, 0.30, 0.40 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 353.15, { 0.30, 0.30, 0.40 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 353.15, { 0.30, 0.30, 0.40 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 573.15, { 0.30, 0.30, 0.40 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 573.15, { 0.30, 0.30, 0.40 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 573.15, { 0.30, 0.30, 0.40 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 293.15, { 0.00, 0.00, 1.00 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 293.15, { 0.00, 0.00, 1.00 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 293.15, { 0.00, 0.00, 1.00 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 353.15, { 0.00, 0.00, 1.00 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 353.15, { 0.00, 0.00, 1.00 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 353.15, { 0.00, 0.00, 1.00 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 573.15, { 0.00, 0.00, 1.00 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 573.15, { 0.00, 0.00, 1.00 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 573.15, { 0.00, 0.00, 1.00 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 293.15, { 0.20, 0.80, 0.00 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 293.15, { 0.20, 0.80, 0.00 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 293.15, { 0.20, 0.80, 0.00 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 353.15, { 0.20, 0.80, 0.00 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 353.15, { 0.20, 0.80, 0.00 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 353.15, { 0.20, 0.80, 0.00 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 ),
    TestData< 3 >( 2.0e+06, 573.15, { 0.20, 0.80, 0.00 }, 5.616333e+04, 1.011782e+03, 3.291616e-04 ),
    TestData< 3 >( 1.5e+07, 573.15, { 0.20, 0.80, 0.00 }, 5.646702e+04, 1.017253e+03, 3.292472e-04 ),
    TestData< 3 >( 6.0e+07, 573.15, { 0.20, 0.80, 0.00 }, 5.753101e+04, 1.036421e+03, 3.295437e-04 )
  )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
