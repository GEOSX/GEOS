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
#include "constitutive/fluid/multifluid/compositional/functions/CompositionalProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
using CompositionalPropertiesTestData = std::tuple<
  real64 const,     // 0 - pressure
  real64 const,     // 1 - temperature
  Feed< NC > const, // 2 - composition
  real64 const,     // 3 - expected molar density
  real64 const,     // 4 - expected molecular weight
  real64 const      // 5 - expected mass density
  >;

template< typename EOS_TYPE, integer NC >
std::vector< CompositionalPropertiesTestData< NC > > generateTestData();

template< integer NC >
std::unique_ptr< TestFluid< NC > > createFluid();
template<>
std::unique_ptr< TestFluid< 4 > > createFluid< 4 >()
{
  return TestFluid< 4 >::create( {Fluid::N2, Fluid::C8, Fluid::C10, Fluid::H2O} );
}

template< typename EOS_TYPE, integer NC >
class CompositionalPropertiesTestDataTestFixture : public ::testing::TestWithParam< CompositionalPropertiesTestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  static constexpr integer numDof = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  CompositionalPropertiesTestDataTestFixture()
    : m_fluid( createFluid< NC >() )
  {}
  ~CompositionalPropertiesTestDataTestFixture() = default;

  // Compares the calculated molar density against the expected value from PVT package
  void testMolarDensity( CompositionalPropertiesTestData< NC > const & data )
  {
    const auto [pressure, temperature, composition] = getInputData( data );
    real64 const expectedMolarDensity = std::get< 3 >( data );

    real64 const molarDensity = computeMolarDensity( pressure,
                                                     temperature,
                                                     composition.toSliceConst() );
    checkRelativeError( molarDensity, expectedMolarDensity, internal::relTol, internal::absTol );
  }

  // Compares the calculated molar density derivatives against numerical calculated
  // finite difference values
  void testMolarDensityDerivative( CompositionalPropertiesTestData< NC > const & data )
  {
    const auto [pressure, temperature, composition] = getInputData( data );

    real64 molarDensity = 0.0;
    stackArray1d< real64, numDof > molarDensityDerivs( numDof );

    computeMolarDensity( pressure,
                         temperature,
                         composition.toSliceConst(),
                         molarDensity,
                         molarDensityDerivs.toSlice() );

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, molarDensityDerivs[Deriv::dP],
      [this, & t=temperature, & zmf=composition]( real64 const p ) -> real64 {
      return computeMolarDensity( p, t, zmf );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, molarDensityDerivs[Deriv::dT],
      [this, & p=pressure, & zmf=composition]( real64 const t ) -> real64 {
      return computeMolarDensity( p, t, zmf );
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, molarDensityDerivs[Deriv::dC + ic],
        [this, & p=pressure, & t=temperature, zmf=composition, ic]( real64 const z ) -> real64 {
        zmf[ic] += z;
        real64 const density = computeMolarDensity( p, t, zmf );
        zmf[ic] -= z;
        return density;
      } );
    }
  }

  // Compares the calculated mass density against the expected value from PVT package
  void testMassDensity( CompositionalPropertiesTestData< NC > const & data )
  {
    const auto [pressure, temperature, composition] = getInputData( data );
    real64 const expectedMassDensity = std::get< 5 >( data );

    real64 const massDensity = computeMassDensity( pressure,
                                                   temperature,
                                                   composition );
    checkRelativeError( massDensity, expectedMassDensity, internal::relTol, internal::absTol );
  }

  // Compares the calculated mass density derivatives against numerical calculated
  // finite difference values
  void testMassDensityDerivative( CompositionalPropertiesTestData< NC > const & data )
  {
    const auto [pressure, temperature, composition] = getInputData( data );

    real64 massDensity = 0.0;
    array1d< real64 > massDensityDerivs( numDof );

    computeMassDensity( pressure,
                        temperature,
                        composition,
                        massDensity,
                        massDensityDerivs );

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, massDensityDerivs[Deriv::dP],
      [this, & t=temperature, & zmf=composition]( real64 const p ) -> real64 {
      return computeMassDensity( p, t, zmf );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, massDensityDerivs[Deriv::dT],
      [this, & p=pressure, & zmf=composition]( real64 const t ) -> real64 {
      return computeMassDensity( p, t, zmf );
    } );

    // -- Composition derivatives
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, massDensityDerivs[Deriv::dC + ic],
        [this, & p=pressure, & t=temperature, & zmf=composition, ic]( real64 const z ) -> real64 {
        real64 const z0 = zmf[ic];
        zmf[ic] += z;
        real64 const density = computeMassDensity( p, t, zmf );
        zmf[ic] = z0;
        return density;
      } );
    }
  }

private:
  std::tuple< real64 const, real64 const, array1d< real64 > >
  getInputData( CompositionalPropertiesTestData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    array1d< real64 > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    return {pressure, temperature, composition};
  }

  real64 computeMolarDensity( real64 const pressure, real64 const temperature,
                              arraySlice1d< real64 const > const & composition ) const
  {
    real64 molarDensity = 0.0;
    stackArray1d< real64, numDof > molarDensityDerivs( numDof );
    computeMolarDensity( pressure, temperature, composition, molarDensity, molarDensityDerivs.toSlice() );
    return molarDensity;
  }

  void computeMolarDensity( real64 const pressure, real64 const temperature,
                            arraySlice1d< real64 const > const & composition,
                            real64 & molarDensity,
                            arraySlice1d< real64 > const molarDensityDerivs ) const
  {
    auto const componentProperties = this->m_fluid->createKernelWrapper();
    auto const binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;
    auto const volumeShift = componentProperties.m_componentVolumeShift;

    real64 compressibilityFactor = 0.0;
    stackArray1d< real64, numComps > aPureCoefficient( numComps );
    stackArray1d< real64, numComps > bPureCoefficient( numComps );
    real64 aMixtureCoefficient = 0.0;
    real64 bMixtureCoefficient = 0.0;

    CubicEOSPhaseModel< EOS_TYPE >::
    computeMixtureCoefficients( numComps,
                                pressure,
                                temperature,
                                composition,
                                componentProperties,
                                aPureCoefficient.toSlice(),
                                bPureCoefficient.toSlice(),
                                aMixtureCoefficient,
                                bMixtureCoefficient );

    stackArray1d< real64, numDof > aMixtureCoefficientDerivs( numDof );
    stackArray1d< real64, numDof > bMixtureCoefficientDerivs( numDof );

    CubicEOSPhaseModel< EOS_TYPE >::
    computeMixtureCoefficients( numComps,
                                pressure,
                                temperature,
                                composition,
                                componentProperties,
                                aPureCoefficient.toSliceConst(),
                                bPureCoefficient.toSliceConst(),
                                aMixtureCoefficient,
                                bMixtureCoefficient,
                                aMixtureCoefficientDerivs.toSlice(),
                                bMixtureCoefficientDerivs.toSlice() );

    CubicEOSPhaseModel< EOS_TYPE >::
    computeCompressibilityFactor( numComps,
                                  composition,
                                  binaryInteractionCoefficients,
                                  aPureCoefficient.toSliceConst(),
                                  bPureCoefficient.toSliceConst(),
                                  aMixtureCoefficient,
                                  bMixtureCoefficient,
                                  compressibilityFactor );

    stackArray1d< real64, numDof > compressibilityFactorDerivs( numDof );

    CubicEOSPhaseModel< EOS_TYPE >::
    computeCompressibilityFactor( numComps,
                                  aMixtureCoefficient,
                                  bMixtureCoefficient,
                                  compressibilityFactor,
                                  aMixtureCoefficientDerivs.toSliceConst(),
                                  bMixtureCoefficientDerivs.toSliceConst(),
                                  compressibilityFactorDerivs.toSlice() );

    CompositionalProperties::
      computeMolarDensity( numComps,
                           pressure,
                           temperature,
                           composition,
                           volumeShift,
                           compressibilityFactor,
                           compressibilityFactorDerivs.toSliceConst(),
                           molarDensity,
                           molarDensityDerivs );
  }

  real64 computeMassDensity( real64 const pressure, real64 const temperature,
                             arraySlice1d< real64 const > const & composition ) const
  {
    real64 massDensity = 0.0;
    stackArray1d< real64, numDof > massDensityDerivs( numDof );
    computeMassDensity( pressure, temperature, composition, massDensity, massDensityDerivs.toSlice() );
    return massDensity;
  }

  void computeMassDensity( real64 const pressure, real64 const temperature,
                           arraySlice1d< real64 const > const & composition,
                           real64 & massDensity,
                           arraySlice1d< real64 > const massDensityDerivs ) const
  {
    auto const componentProperties = this->m_fluid->createKernelWrapper();
    auto const molecularWeight = componentProperties.m_componentMolarWeight;

    real64 molarDensity = 0.0;
    stackArray1d< real64, numDof > molarDensityDerivs( numDof );

    computeMolarDensity( pressure, temperature,
                         composition,
                         molarDensity,
                         molarDensityDerivs.toSlice() );

    CompositionalProperties::
      computeMassDensity( numComps,
                          composition,
                          molecularWeight,
                          molarDensity,
                          molarDensityDerivs.toSliceConst(),
                          massDensity,
                          massDensityDerivs );
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid;
};

using PR4Comp = CompositionalPropertiesTestDataTestFixture< PengRobinsonEOS, 4 >;

TEST_P( PR4Comp, testMolarDensity )
{
  testMolarDensity( GetParam() );
}

TEST_P( PR4Comp, testMolarDensityDerivative )
{
  testMolarDensityDerivative( GetParam() );
}

TEST_P( PR4Comp, testMassDensity )
{
  testMassDensity( GetParam() );
}

TEST_P( PR4Comp, testMassDensityDerivative )
{
  testMassDensityDerivative( GetParam() );
}

// Test data generated from PVT package
// All compositions are single phase
template<>
std::vector< CompositionalPropertiesTestData< 4 > > generateTestData< PengRobinsonEOS, 4 >()
{
  return {
    { 1.000000e+05, 2.771500e+02, { 0.000000, 0.495099, 0.495118, 0.009783 }, 3.733061e+03, 1.271768e+02, 4.747588e+02 },
    { 1.000000e+05, 2.771500e+02, { 0.000652, 0.128231, 0.128281, 0.742836 }, 1.119134e+04, 4.630013e+01, 5.181608e+02 },
    { 1.000000e+05, 2.771500e+02, { 0.855328, 0.000205, 0.000000, 0.144467 }, 4.348717e+01, 2.658628e+01, 1.156162e+00 },
    { 1.000000e+05, 2.886500e+02, { 0.000507, 0.112984, 0.113029, 0.773480 }, 1.214893e+04, 4.293636e+01, 5.216311e+02 },
    { 1.000000e+05, 2.886500e+02, { 0.777870, 0.000520, 0.000000, 0.221610 }, 4.177763e+01, 2.584219e+01, 1.079625e+00 },
    { 1.000000e+05, 2.886500e+02, { 0.985235, 0.000000, 0.000000, 0.014765 }, 4.169466e+01, 2.786538e+01, 1.161838e+00 },
    { 1.000000e+05, 2.981500e+02, { 0.653033, 0.000901, 0.000000, 0.346066 }, 4.049414e+01, 2.463074e+01, 9.974005e-01 },
    { 1.000000e+05, 2.981500e+02, { 0.000506, 0.143046, 0.143326, 0.713122 }, 1.013248e+04, 4.959370e+01, 5.025071e+02 },
    { 1.000000e+05, 2.981500e+02, { 0.582848, 0.000748, 0.000000, 0.416404 }, 4.053125e+01, 2.391430e+01, 9.692764e-01 },
    { 1.000000e+05, 2.981500e+02, { 0.972848, 0.000000, 0.000000, 0.027152 }, 4.036515e+01, 2.774153e+01, 1.119791e+00 },
    { 1.000000e+05, 3.331500e+02, { 0.000000, 0.477146, 0.477164, 0.045691 }, 3.754190e+03, 1.232183e+02, 4.625850e+02 },
    { 1.000000e+05, 3.331500e+02, { 0.210877, 0.008984, 0.000001, 0.780137 }, 3.640845e+01, 2.098792e+01, 7.641376e-01 },
    { 1.000000e+05, 3.331500e+02, { 0.818043, 0.000000, 0.000000, 0.181957 }, 3.614852e+01, 2.619379e+01, 9.468669e-01 },
    { 1.000000e+05, 3.721500e+02, { 0.000000, 0.104818, 0.104822, 0.790360 }, 3.312721e+01, 4.112577e+01, 1.362382e+00 },
    { 1.000000e+05, 3.721500e+02, { 0.000117, 0.356347, 0.549688, 0.093848 }, 3.581375e+03, 1.206090e+02, 4.319462e+02 },
    { 1.013250e+05, 2.771500e+02, { 0.000000, 0.495099, 0.495118, 0.009783 }, 3.733064e+03, 1.271768e+02, 4.747592e+02 },
    { 1.013250e+05, 2.886500e+02, { 0.000516, 0.112950, 0.112994, 0.773540 }, 1.215157e+04, 4.292888e+01, 5.216533e+02 },
    { 1.013250e+05, 2.886500e+02, { 0.780815, 0.000514, 0.000000, 0.218672 }, 4.233123e+01, 2.587101e+01, 1.095152e+00 },
    { 1.013250e+05, 2.886500e+02, { 0.000599, 0.132633, 0.132752, 0.734016 }, 1.081051e+04, 4.727865e+01, 5.111063e+02 },
    { 1.013250e+05, 2.886500e+02, { 0.743752, 0.000426, 0.000000, 0.255823 }, 4.234963e+01, 2.549200e+01, 1.079576e+00 },
    { 1.013250e+05, 2.981500e+02, { 0.000448, 0.114592, 0.114680, 0.770280 }, 1.192251e+04, 4.329570e+01, 5.161936e+02 },
    { 1.013250e+05, 2.981500e+02, { 0.657746, 0.000890, 0.000000, 0.341364 }, 4.103031e+01, 2.467683e+01, 1.012498e+00 },
    { 1.013250e+05, 2.981500e+02, { 0.000516, 0.142463, 0.142736, 0.714285 }, 1.016359e+04, 4.946432e+01, 5.027349e+02 },
    { 1.013250e+05, 3.331500e+02, { 0.000000, 0.477146, 0.477164, 0.045690 }, 3.754193e+03, 1.232184e+02, 4.625856e+02 },
    { 1.013250e+05, 3.331500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.593026e+04, 1.801500e+01, 8.274336e+02 },
    { 1.013250e+05, 3.331500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.593025e+04, 1.801500e+01, 8.274336e+02 },
    { 1.013250e+05, 3.331500e+02, { 0.820407, 0.000000, 0.000000, 0.179593 }, 3.662747e+01, 2.621743e+01, 9.602781e-01 },
    { 1.013250e+05, 3.721500e+02, { 0.080408, 0.000000, 0.000000, 0.919592 }, 3.299996e+01, 1.881892e+01, 6.210236e-01 },
    { 5.000000e+06, 2.771500e+02, { 0.000000, 0.495216, 0.495235, 0.009549 }, 3.742134e+03, 1.272026e+02, 4.760091e+02 },
    { 5.000000e+06, 2.771500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.772799e+04, 1.801500e+01, 8.598197e+02 },
    { 5.000000e+06, 2.771500e+02, { 0.029212, 0.107914, 0.107918, 0.754956 }, 1.264269e+04, 4.210047e+01, 5.322630e+02 },
    { 5.000000e+06, 2.886500e+02, { 0.000000, 0.493093, 0.493111, 0.013796 }, 3.739164e+03, 1.267344e+02, 4.738808e+02 },
    { 5.000000e+06, 2.886500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.739475e+04, 1.801500e+01, 8.538165e+02 },
    { 5.000000e+06, 2.886500e+02, { 0.999549, 0.000000, 0.000000, 0.000451 }, 2.122222e+03, 2.800849e+01, 5.944023e+01 },
    { 5.000000e+06, 2.981500e+02, { 0.000000, 0.490855, 0.490873, 0.018272 }, 3.739300e+03, 1.262410e+02, 4.720529e+02 },
    { 5.000000e+06, 2.981500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.710908e+04, 1.801500e+01, 8.486701e+02 },
    { 5.000000e+06, 2.981500e+02, { 0.030096, 0.107834, 0.107839, 0.754231 }, 1.241638e+04, 4.209177e+01, 5.226274e+02 },
    { 5.000000e+06, 3.331500e+02, { 0.031858, 0.107709, 0.107720, 0.752713 }, 1.198645e+04, 4.208259e+01, 5.044209e+02 },
    { 5.000000e+06, 3.331500e+02, { 0.964491, 0.000264, 0.000000, 0.035245 }, 1.824920e+03, 2.768343e+01, 5.052004e+01 },
    { 5.000000e+06, 3.331500e+02, { 0.035989, 0.120533, 0.120574, 0.722904 }, 1.106526e+04, 4.495504e+01, 4.974394e+02 },
    { 5.000000e+06, 3.331500e+02, { 0.961874, 0.000247, 0.000000, 0.037879 }, 1.826197e+03, 2.765560e+01, 5.050457e+01 },
    { 5.000000e+06, 3.721500e+02, { 0.037408, 0.121753, 0.121923, 0.718916 }, 1.046451e+04, 4.525417e+01, 4.735627e+02 },
    { 5.000000e+06, 3.721500e+02, { 0.889694, 0.001014, 0.000003, 0.109288 }, 1.643713e+03, 2.700817e+01, 4.439368e+01 },
    { 1.000000e+07, 2.771500e+02, { 0.000004, 0.000000, 0.000000, 0.999996 }, 4.775121e+04, 1.801504e+01, 8.602398e+02 },
    { 1.000000e+07, 2.771500e+02, { 0.999835, 0.000000, 0.000000, 0.000165 }, 4.468580e+03, 2.801135e+01, 1.251710e+02 },
    { 1.000000e+07, 2.886500e+02, { 0.000000, 0.493260, 0.493279, 0.013460 }, 3.748279e+03, 1.267714e+02, 4.751745e+02 },
    { 1.000000e+07, 2.886500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.742059e+04, 1.801500e+01, 8.542819e+02 },
    { 1.000000e+07, 2.886500e+02, { 0.055934, 0.104932, 0.104936, 0.734198 }, 1.273695e+04, 4.171005e+01, 5.312588e+02 },
    { 1.000000e+07, 2.981500e+02, { 0.000000, 0.491084, 0.491102, 0.017814 }, 3.748930e+03, 1.262914e+02, 4.734578e+02 },
    { 1.000000e+07, 2.981500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.713720e+04, 1.801500e+01, 8.491766e+02 },
    { 1.000000e+07, 2.981500e+02, { 0.056951, 0.104818, 0.104822, 0.733409 }, 1.263544e+04, 4.169517e+01, 5.268366e+02 },
    { 1.000000e+07, 2.981500e+02, { 0.065451, 0.116299, 0.116311, 0.701940 }, 1.172878e+04, 4.431245e+01, 5.197309e+02 },
    { 1.000000e+07, 2.981500e+02, { 0.991760, 0.000074, 0.000000, 0.008166 }, 4.103168e+03, 2.793775e+01, 1.146333e+02 },
    { 1.000000e+07, 3.331500e+02, { 0.000000, 0.478424, 0.478442, 0.043135 }, 3.777896e+03, 1.235001e+02, 4.665706e+02 },
    { 1.000000e+07, 3.331500e+02, { 0.000000, 0.000000, 0.000000, 1.000000 }, 4.600669e+04, 1.801500e+01, 8.288105e+02 },
    { 1.000000e+07, 3.331500e+02, { 0.996818, 0.000000, 0.000000, 0.003182 }, 3.589847e+03, 2.798118e+01, 1.004482e+02 },
    { 1.000000e+07, 3.721500e+02, { 0.000000, 0.453197, 0.453215, 0.093588 }, 3.878777e+03, 1.179381e+02, 4.574555e+02 },
    { 1.000000e+08, 2.771500e+02, { 0.000000, 0.496608, 0.496627, 0.006766 }, 3.834889e+03, 1.275094e+02, 4.889844e+02 },
    { 1.000000e+08, 2.886500e+02, { 0.156951, 0.104818, 0.104822, 0.633409 }, 1.310108e+04, 4.269497e+01, 5.593502e+02 },
    { 1.000000e+08, 2.886500e+02, { 0.000031, 0.000000, 0.000000, 0.999969 }, 4.783830e+04, 1.801531e+01, 8.618217e+02 },
    { 1.000000e+08, 2.886500e+02, { 0.999447, 0.000000, 0.000000, 0.000553 }, 2.330195e+04, 2.800748e+01, 6.526287e+02 },
    { 1.000000e+08, 2.981500e+02, { 0.000000, 0.493623, 0.493642, 0.012735 }, 3.840586e+03, 1.268513e+02, 4.871835e+02 },
    { 1.000000e+08, 3.331500e+02, { 0.010000, 0.000000, 0.000000, 0.990000 }, 4.638822e+04, 1.811498e+01, 8.403216e+02 },
    { 1.000000e+08, 3.721500e+02, { 0.156951, 0.104818, 0.104822, 0.633409 }, 1.247160e+04, 4.269497e+01, 5.324747e+02 },
    { 1.000000e+08, 3.721500e+02, { 0.000671, 0.000000, 0.000000, 0.999329 }, 4.538427e+04, 1.802171e+01, 8.179019e+02 },
    { 1.000000e+08, 3.721500e+02, { 0.991955, 0.000000, 0.000000, 0.008045 }, 1.979639e+04, 2.793257e+01, 5.529641e+02 }
  };
}

INSTANTIATE_TEST_SUITE_P(
  CompositionalPropertiesTest,
  PR4Comp,
  ::testing::ValuesIn( generateTestData< PengRobinsonEOS, 4 >())
  );

} // testing

} // geos
