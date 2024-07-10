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

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
struct FluidData
{};

template<>
struct FluidData< 2 >
{
  static std::unique_ptr< TestFluid< 2 > > createFluid()
  {
    return TestFluid< 2 >::create( {Fluid::CO2, Fluid::C5} );
  }
};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    return TestFluid< 4 >::create( {Fluid::N2, Fluid::C8, Fluid::C10, Fluid::H2O} );
  }
};

template< int NC >
using FlashData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const,   // total composition
  bool,               // expected flash status (success/failure)
  real64 const,       // expected vapour fraction
  Feed< NC > const,   // expected liquid composition
  Feed< NC > const    // expected vapour composition
  >;

template< int NC, EquationOfStateType EOS_TYPE >
class NegativeTwoPhaseFlashTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  NegativeTwoPhaseFlashTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {}

  ~NegativeTwoPhaseFlashTestFixture() = default;

  void testFlash( FlashData< NC > const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    bool const expectedStatus = std::get< 3 >( data );
    real64 const expectedVapourFraction = std::get< 4 >( data );

    stackArray1d< real64, numComps > expectedLiquidComposition;
    TestFluid< NC >::createArray( expectedLiquidComposition, std::get< 5 >( data ));
    stackArray1d< real64, numComps > expectedVapourComposition;
    TestFluid< NC >::createArray( expectedVapourComposition, std::get< 6 >( data ));

    real64 vapourFraction = -1.0;
    stackArray1d< real64, numComps > liquidComposition( numComps );
    stackArray1d< real64, numComps > vapourComposition( numComps );
    stackArray2d< real64, numComps > kValues( 1, numComps );
    kValues.zero();

    bool status = NegativeTwoPhaseFlash::compute(
      numComps,
      pressure,
      temperature,
      composition.toSliceConst(),
      componentProperties,
      EOS_TYPE,
      EOS_TYPE,
      kValues.toSlice(),
      vapourFraction,
      liquidComposition.toSlice(),
      vapourComposition.toSlice() );

    // Check the flash success result
    ASSERT_EQ( expectedStatus, status );

    if( !expectedStatus )
    {
      return;
    }

    // Check the vaopur fraction
    checkRelativeError( expectedVapourFraction, vapourFraction, relTol, absTol );

    // Check liquid composition
    if( expectedVapourFraction < 1.0 - absTol )
    {
      for( integer ic=0; ic<numComps; ++ic )
      {
        checkRelativeError( expectedLiquidComposition[ic], liquidComposition[ic], relTol, absTol );
      }
    }

    // Check vapour composition
    if( absTol < expectedVapourFraction )
    {
      for( integer ic=0; ic<numComps; ++ic )
      {
        checkRelativeError( expectedVapourComposition[ic], vapourComposition[ic], relTol, absTol );
      }
    }
  }

  void testFlashDerivatives( FlashData< NC > const & data )
  {
    // Number of output values from each flash calculation
    constexpr integer numValues = 1 + 2*numComps;

    auto componentProperties = this->m_fluid->createKernelWrapper();

    bool const expectedStatus = std::get< 3 >( data );
    if( !expectedStatus ) return;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    real64 vapourFraction = -1.0;
    stackArray1d< real64, numComps > liquidComposition( numComps );
    stackArray1d< real64, numComps > vapourComposition( numComps );
    stackArray2d< real64, numComps > kValues( 1, numComps );
    kValues.zero();

    stackArray1d< real64, numDofs > vapourFractionDerivs( numDofs );
    stackArray2d< real64, numComps * numDofs > liquidCompositionDerivs( numComps, numDofs );
    stackArray2d< real64, numComps * numDofs > vapourCompositionDerivs( numComps, numDofs );
    stackArray1d< real64, numValues > derivatives( numValues );

    // Combine values and derivatives into a single output
    auto const concatDerivatives = []( integer const kc, auto & derivs, auto const & v, auto const & xmf, auto const & ymf ){
      derivs[0] = v[kc];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        derivs[1+ic] = xmf( ic, kc );
        derivs[1+ic+numComps] = ymf( ic, kc );
      }
    };
    std::cout << std::scientific << std::setprecision( 8 );

    auto const evaluateFlash = [&]( real64 const p, real64 const t, auto const & zmf, auto & values ){
      stackArray1d< real64, numComps > displacedLiquidComposition( numComps );
      stackArray1d< real64, numComps > displacedVapourComposition( numComps );
      kValues.zero();

      NegativeTwoPhaseFlash::compute(
        numComps,
        p,
        t,
        zmf.toSliceConst(),
        componentProperties,
        EOS_TYPE,
        EOS_TYPE,
        kValues.toSlice(),
        values[0],
        displacedLiquidComposition.toSlice(),
        displacedVapourComposition.toSlice() );
      for( integer ic = 0; ic < numComps; ++ic )
      {
        values[1+ic] = displacedLiquidComposition[ic];
        values[1+ic+numComps] = displacedVapourComposition[ic];
      }
    };

    NegativeTwoPhaseFlash::compute(
      numComps,
      pressure,
      temperature,
      composition.toSliceConst(),
      componentProperties,
      EOS_TYPE,
      EOS_TYPE,
      kValues.toSlice(),
      vapourFraction,
      liquidComposition.toSlice(),
      vapourComposition.toSlice() );

    NegativeTwoPhaseFlash::computeDerivatives(
      numComps,
      pressure,
      temperature,
      composition.toSliceConst(),
      componentProperties,
      EOS_TYPE,
      EOS_TYPE,
      vapourFraction,
      liquidComposition.toSliceConst(),
      vapourComposition.toSliceConst(),
      vapourFractionDerivs.toSlice(),
      liquidCompositionDerivs.toSlice(),
      vapourCompositionDerivs.toSlice() );

    // Test against numerically calculated values
    // --- Pressure derivatives ---
    concatDerivatives( Deriv::dP, derivatives, vapourFractionDerivs, liquidCompositionDerivs, vapourCompositionDerivs );
    real64 const dp = 1.0e-4 * pressure;
    geos::testing::internal::testNumericalDerivative< numValues >(
      pressure, dp, derivatives,
      [&]( real64 const p, auto & values ) {
      evaluateFlash( p, temperature, composition, values );
    } );

    // --- Temperature derivatives ---
    concatDerivatives( Deriv::dT, derivatives, vapourFractionDerivs, liquidCompositionDerivs, vapourCompositionDerivs );
    real64 const dT = 1.0e-6 * temperature;
    geos::testing::internal::testNumericalDerivative< numValues >(
      temperature, dT, derivatives,
      [&]( real64 const t, auto & values ) {
      evaluateFlash( pressure, t, composition, values );
    } );

    // --- Composition derivatives ---
    real64 constexpr dz = 1.0e-7;
    for( integer jc = 0; jc < numComps; ++jc )
    {
      if( composition[jc] < 1.0e-6 ) continue;
      integer const kc = Deriv::dC + jc;
      concatDerivatives( kc, derivatives, vapourFractionDerivs, liquidCompositionDerivs, vapourCompositionDerivs );
      geos::testing::internal::testNumericalDerivative< numValues >(
        0.0, dz, derivatives,
        [&]( real64 const z, auto & values ) {
        real64 const originalFraction = composition[jc];
        composition[jc] += z;
        evaluateFlash( pressure, temperature, composition, values );
        composition[jc] = originalFraction;
      }, relTol, absTol );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

using NegativeTwoPhaseFlash2CompPR = NegativeTwoPhaseFlashTestFixture< 2, EquationOfStateType::PengRobinson >;
using NegativeTwoPhaseFlash2CompSRK = NegativeTwoPhaseFlashTestFixture< 2, EquationOfStateType::SoaveRedlichKwong >;
using NegativeTwoPhaseFlash4CompPR = NegativeTwoPhaseFlashTestFixture< 4, EquationOfStateType::PengRobinson >;
using NegativeTwoPhaseFlash4CompSRK = NegativeTwoPhaseFlashTestFixture< 4, EquationOfStateType::SoaveRedlichKwong >;

TEST_P( NegativeTwoPhaseFlash2CompPR, testNegativeFlash )
{
  testFlash( GetParam() );
}

TEST_P( NegativeTwoPhaseFlash2CompSRK, testNegativeFlash )
{
  testFlash( GetParam() );
}

TEST_P( NegativeTwoPhaseFlash4CompPR, testNegativeFlash )
{
  testFlash( GetParam() );
}

TEST_P( NegativeTwoPhaseFlash4CompSRK, testNegativeFlash )
{
  testFlash( GetParam() );
}

TEST_P( NegativeTwoPhaseFlash2CompPR, testNegativeFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

TEST_P( NegativeTwoPhaseFlash2CompSRK, testNegativeFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

TEST_P( NegativeTwoPhaseFlash4CompPR, testNegativeFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

TEST_P( NegativeTwoPhaseFlash4CompSRK, testNegativeFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data generated by PVTPackage
//-------------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlash,
  NegativeTwoPhaseFlash2CompPR,
  ::testing::Values(
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.10000000, 0.90000000 }, true, 0.89038113, { 0.89566514, 0.10433486 }, { 0.00204205, 0.99795795 } ),
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.75000000, 0.25000000 }, true, 0.16300512, { 0.89566514, 0.10433486 }, { 0.00204205, 0.99795795 } ),
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.50000000, 0.50000000 }, true, 0.44276513, { 0.89566514, 0.10433486 }, { 0.00204205, 0.99795795 } ),
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.00204205, 0.99795795 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.10000000, 0.90000000 }, true, 0.89010742, { 0.89366817, 0.10633183 }, { 0.00201380, 0.99798620 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.75000000, 0.25000000 }, true, 0.16112541, { 0.89366817, 0.10633183 }, { 0.00201380, 0.99798620 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.50000000, 0.50000000 }, true, 0.44150310, { 0.89366817, 0.10633183 }, { 0.00201380, 0.99798620 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.00201380, 0.99798620 } ),
    FlashData< 2 >( 1.000000e+05, 1.931500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+05, 1.931500e+02, { 0.75000000, 0.25000000 }, true, 1.00000000, { 0.75000000, 0.25000000 }, { 0.75000000, 0.25000000 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.88618228, 0.11381772 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.75000000, 0.25000000 }, true, 0.18928601, { 0.88618228, 0.11381772 }, { 0.16672986, 0.83327014 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.50000000, 0.50000000 }, true, 0.53677251, { 0.88618228, 0.11381772 }, { 0.16672986, 0.83327014 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.16672986, 0.83327014 } ),
    FlashData< 2 >( 1.000000e+08, 1.931500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+05, 2.771500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 5.000000e+06, 2.771500e+02, { 0.50000000, 0.50000000 }, true, 1.00000000, { 0.93119724, 0.06880276 }, { 0.50000000, 0.50000000 } ),
    FlashData< 2 >( 5.000000e+06, 2.771500e+02, { 0.90000000, 0.10000000 }, true, 0.27263042, { 0.93119724, 0.06880276 }, { 0.81676673, 0.18323327 } ),
    FlashData< 2 >( 1.000000e+07, 2.771500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.55705424, 0.44294576 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+07, 2.771500e+02, { 0.75000000, 0.25000000 }, true, 0.00000000, { 0.75000000, 0.25000000 }, { 0.55705421, 0.44294579 } ),
    FlashData< 2 >( 1.000000e+08, 3.331500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+05, 3.731500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 5.000000e+06, 3.731500e+02, { 0.90000000, 0.10000000 }, true, 1.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+07, 3.731500e+02, { 0.90000000, 0.10000000 }, true, 1.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+08, 3.731500e+02, { 0.10000000, 0.90000000 }, true, 0.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+08, 3.731500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } )
    )
  );

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlash,
  NegativeTwoPhaseFlash2CompSRK,
  ::testing::Values(
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.10000000, 0.90000000 }, true, 0.89111708, { 0.90429170, 0.09570830 }, { 0.00172601, 0.99827399 } ),
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.75000000, 0.25000000 }, true, 0.17094789, { 0.90429170, 0.09570830 }, { 0.00172601, 0.99827399 } ),
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.50000000, 0.50000000 }, true, 0.44793604, { 0.90429170, 0.09570830 }, { 0.00172601, 0.99827399 } ),
    FlashData< 2 >( 1.000000e+05, 1.231500e+02, { 0.90000000, 0.10000000 }, true, 0.00475500, { 0.90429170, 0.09570830 }, { 0.00172601, 0.99827399 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.10000000, 0.90000000 }, true, 0.89087481, { 0.90248098, 0.09751902 }, { 0.00170237, 0.99829763 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.75000000, 0.25000000 }, true, 0.16927687, { 0.90248098, 0.09751902 }, { 0.00170237, 0.99829763 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.50000000, 0.50000000 }, true, 0.44681454, { 0.90248098, 0.09751902 }, { 0.00170237, 0.99829763 } ),
    FlashData< 2 >( 1.013250e+05, 1.231500e+02, { 0.90000000, 0.10000000 }, true, 0.00275426, { 0.90248098, 0.09751902 }, { 0.00170237, 0.99829763 } ),
    FlashData< 2 >( 1.000000e+06, 1.231500e+02, { 0.10000000, 0.90000000 }, true, 0.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+06, 1.231500e+02, { 0.50000000, 0.50000000 }, true, 0.00000000, { 0.50000000, 0.50000000 }, { 0.49950202, 0.50049798 } ),
    FlashData< 2 >( 1.000000e+06, 1.231500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 5.000000e+06, 1.231500e+02, { 0.75000000, 0.25000000 }, true, 0.00000000, { 0.75000000, 0.25000000 }, { 0.74999999, 0.25000001 } ),
    FlashData< 2 >( 5.000000e+06, 1.231500e+02, { 0.50000000, 0.50000000 }, true, 0.00000000, { 0.50000000, 0.50000000 }, { 0.49979984, 0.50020016 } ),
    FlashData< 2 >( 5.000000e+06, 1.231500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+08, 1.231500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+05, 1.931500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+05, 1.931500e+02, { 0.75000000, 0.25000000 }, true, 1.00000000, { 0.75000000, 0.25000000 }, { 0.75000000, 0.25000000 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.89288212, 0.10711788 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.75000000, 0.25000000 }, true, 0.19641336, { 0.89288212, 0.10711788 }, { 0.16542588, 0.83457412 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.50000000, 0.50000000 }, true, 0.54007664, { 0.89288212, 0.10711788 }, { 0.16542588, 0.83457412 } ),
    FlashData< 2 >( 1.000000e+06, 1.931500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.16542588, 0.83457412 } ),
    FlashData< 2 >( 1.000000e+08, 1.931500e+02, { 0.50000000, 0.50000000 }, true, 0.00000000, { 0.50000000, 0.50000000 }, { 0.49999999, 0.50000001 } ),
    FlashData< 2 >( 1.000000e+08, 1.931500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+05, 2.771500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+05, 2.771500e+02, { 0.90000000, 0.10000000 }, true, 1.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 5.000000e+06, 2.771500e+02, { 0.50000000, 0.50000000 }, true, 1.00000000, { 0.93366481, 0.06633519 }, { 0.50000000, 0.50000000 } ),
    FlashData< 2 >( 5.000000e+06, 2.771500e+02, { 0.90000000, 0.10000000 }, true, 0.30050187, { 0.93366481, 0.06633519 }, { 0.82163618, 0.17836382 } ),
    FlashData< 2 >( 1.000000e+07, 2.771500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.55571828, 0.44428172 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+07, 2.771500e+02, { 0.50000000, 0.50000000 }, true, 1.00000000, { 0.55571830, 0.44428170 }, { 0.50000000, 0.50000000 } ),
    FlashData< 2 >( 1.000000e+07, 2.771500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.55571831, 0.44428169 } ),
    FlashData< 2 >( 1.000000e+07, 2.886500e+02, { 0.10000000, 0.90000000 }, true, 1.00000000, { 0.66024784, 0.33975216 }, { 0.10000000, 0.90000000 } ),
    FlashData< 2 >( 1.000000e+07, 2.886500e+02, { 0.75000000, 0.25000000 }, true, 0.00000000, { 0.75000000, 0.25000000 }, { 0.66024775, 0.33975225 } ),
    FlashData< 2 >( 1.000000e+07, 2.886500e+02, { 0.50000000, 0.50000000 }, true, 1.00000000, { 0.66024781, 0.33975219 }, { 0.50000000, 0.50000000 } ),
    FlashData< 2 >( 1.000000e+07, 2.886500e+02, { 0.90000000, 0.10000000 }, true, 0.00000000, { 0.90000000, 0.10000000 }, { 0.66024783, 0.33975217 } ),
    FlashData< 2 >( 1.000000e+07, 3.731500e+02, { 0.90000000, 0.10000000 }, true, 1.00000000, { 0.90000000, 0.10000000 }, { 0.90000000, 0.10000000 } ),
    FlashData< 2 >( 1.000000e+08, 3.731500e+02, { 0.10000000, 0.90000000 }, true, 0.00000000, { 0.10000000, 0.90000000 }, { 0.10000000, 0.90000000 } )
    )
  );

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlash,
  NegativeTwoPhaseFlash4CompPR,
  ::testing::Values(
    FlashData< 4 >( 1.000000e+05, 1.931500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.05629583, { 0.00070397, 0.11107082, 0.11107506, 0.77715015 },
                    { 0.99983707, 0.00000005, 0.00000000, 0.00016289 } ),
    FlashData< 4 >( 1.000000e+05, 1.931500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.15624142, { 0.00088099, 0.12422747, 0.12423222, 0.75065932 },
                    { 0.99978392, 0.00000003, 0.00000000, 0.00021605 } ),
    FlashData< 4 >( 5.000000e+06, 1.931500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.02878730, { 0.02899884, 0.10792486, 0.10792898, 0.75514731 },
                    { 0.99998806, 0.00000004, 0.00000000, 0.00001190 } ),
    FlashData< 4 >( 5.000000e+06, 1.931500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.12514257, { 0.03636069, 0.11981152, 0.11981609, 0.72401170 },
                    { 0.99998417, 0.00000003, 0.00000000, 0.00001580 } ),
    FlashData< 4 >( 5.000000e+06, 1.931500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.79032837, { 0.00000000, 0.49991503, 0.49993411, 0.00015087 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 5.000000e+06, 1.931500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.09334204, { 0.01265794, 0.00000000, 0.11561361, 0.87172845 },
                    { 0.99999506, 0.00000000, 0.00000000, 0.00000494 } ),
    FlashData< 4 >( 1.000000e+07, 1.931500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.10262297, { 0.06054432, 0.11680480, 0.11680932, 0.70584155 },
                    { 0.99997026, 0.00000052, 0.00000000, 0.00002922 } ),
    FlashData< 4 >( 1.000000e+08, 1.931500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.26661844, { 0.00000001, 0.00000000, 0.00000000, 0.99999999 },
                    { 0.21360485, 0.39313860, 0.39315361, 0.00010294 } ),
    FlashData< 4 >( 1.000000e+06, 2.771500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.10252763, { 0.00359465, 0.00000000, 0.11679690, 0.87960845 },
                    { 0.99087343, 0.00000000, 0.00000008, 0.00912649 } ),
    FlashData< 4 >( 1.000000e+06, 2.771500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 0.99071680, { 0.00000042, 0.00000000, 0.00000000, 0.99999958 },
                    { 0.99927648, 0.00000000, 0.00000000, 0.00072352 } ),
    FlashData< 4 >( 5.000000e+06, 2.771500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.02869272, { 0.02921217, 0.10791374, 0.10791847, 0.75495561 },
                    { 0.99596684, 0.00002082, 0.00000001, 0.00401233 } ),
    FlashData< 4 >( 5.000000e+06, 2.771500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.12793659, { 0.03394197, 0.12019282, 0.12019997, 0.72566524 },
                    { 0.99542635, 0.00001751, 0.00000000, 0.00455613 } ),
    FlashData< 4 >( 5.000000e+06, 2.771500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.78833888, { 0.00000000, 0.49521612, 0.49523502, 0.00954887 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 5.000000e+06, 2.771500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.01000018, { 0.00000196, 0.00000000, 0.00000000, 0.99999804 },
                    { 0.99978790, 0.00000000, 0.00000000, 0.00021210 } ),
    FlashData< 4 >( 1.000000e+07, 2.771500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.00244308, { 0.05464929, 0.10507462, 0.10507872, 0.73519737 },
                    { 0.99678327, 0.00003656, 0.00000004, 0.00318014 } ),
    FlashData< 4 >( 1.000000e+07, 2.771500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.10021652, { 0.06345880, 0.11648904, 0.11649691, 0.70355525 },
                    { 0.99636092, 0.00003083, 0.00000002, 0.00360823 } ),
    FlashData< 4 >( 1.000000e+07, 2.771500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.78838676, { 0.00000000, 0.49532818, 0.49534708, 0.00932474 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 1.000000e+05, 2.886500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 1.00000000, { 0.00000007, 0.00000000, 0.00000000, 0.99999993 },
                    { 0.99000000, 0.00000000, 0.00000000, 0.01000000 } ),
    FlashData< 4 >( 1.000000e+05, 2.886500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.01014979, { 0.00000007, 0.00000000, 0.00000000, 0.99999993 },
                    { 0.98523527, 0.00000000, 0.00000000, 0.01476473 } ),
    FlashData< 4 >( 1.000000e+07, 2.886500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 0.99033139, { 0.00000617, 0.00000000, 0.00000000, 0.99999383 },
                    { 0.99966531, 0.00000000, 0.00000000, 0.00033469 } ),
    FlashData< 4 >( 1.000000e+07, 2.886500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.00999724, { 0.00000617, 0.00000000, 0.00000000, 0.99999383 },
                    { 0.99966531, 0.00000000, 0.00000000, 0.00033469 } ),
    FlashData< 4 >( 1.000000e+08, 2.886500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.26901755, { 0.00001146, 0.00000000, 0.00000000, 0.99998854 },
                    { 0.21166879, 0.38963257, 0.38964744, 0.00905121 } ),
    FlashData< 4 >( 1.000000e+05, 3.331500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.86540955, { 0.00020370, 0.00000000, 0.77880508, 0.22099122 },
                    { 0.12108785, 0.00000000, 0.00000263, 0.87890952 } ),
    FlashData< 4 >( 1.000000e+05, 3.331500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 1.00000000, { 0.00000033, 0.00000000, 0.00000000, 0.99999967 },
                    { 0.99000000, 0.00000000, 0.00000000, 0.01000000 } ),
    FlashData< 4 >( 1.000000e+05, 3.331500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.01222390, { 0.00000033, 0.00000000, 0.00000000, 0.99999967 },
                    { 0.81804281, 0.00000000, 0.00000000, 0.18195719 } ),
    FlashData< 4 >( 1.000000e+07, 3.731500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 1.00000000, { 0.00012025, 0.00000000, 0.00000000, 0.99987975 },
                    { 0.99000000, 0.00000000, 0.00000000, 0.01000000 } ),
    FlashData< 4 >( 1.000000e+08, 3.731500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.00000000, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 },
                    { 0.93121567, 0.00696526, 0.00174080, 0.06007828 } ),
    FlashData< 4 >( 1.000000e+08, 3.731500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.77592410, { 0.00000000, 0.46777900, 0.46779685, 0.06442416 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 1.013250e+05, 4.731500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.91046790, { 0.00021364, 0.00000000, 0.97656019, 0.02322617 },
                    { 0.11510441, 0.00000000, 0.01909844, 0.86579715 } ),
    FlashData< 4 >( 1.013250e+05, 4.731500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 1.00000000, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 },
                    { 0.01000000, 0.00000000, 0.00000000, 0.99000000 } ),
    FlashData< 4 >( 5.000000e+06, 4.731500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.01474382, { 0.00043465, 0.00000000, 0.00000000, 0.99956535 },
                    { 0.64920515, 0.00000000, 0.00000000, 0.35079485 } ),
    FlashData< 4 >( 1.000000e+07, 4.731500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.00000000, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 },
                    { 0.60848558, 0.01033586, 0.00049071, 0.38068785 } ),
    FlashData< 4 >( 1.000000e+07, 4.731500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.13867033, { 0.08254496, 0.11986286, 0.12161793, 0.67597425 },
                    { 0.61911282, 0.01136912, 0.00049663, 0.36902144 } ),
    FlashData< 4 >( 1.000000e+07, 4.731500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.66978458, { 0.00000000, 0.31742202, 0.31743521, 0.36514277 },
                    { 0.00000000, 0.00000053, 0.00000000, 0.99999947 } ),
    FlashData< 4 >( 1.000000e+07, 4.731500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.08507006, { 0.05879909, 0.00000000, 0.11449503, 0.82670588 },
                    { 0.59975210, 0.00000000, 0.00078842, 0.39945949 } ),
    FlashData< 4 >( 1.000000e+07, 4.731500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.01121076, { 0.00102378, 0.00000000, 0.00000000, 0.99897622 },
                    { 0.80170284, 0.00000000, 0.00000000, 0.19829716 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.00000000, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 },
                    { 0.72438623, 0.02564505, 0.01665696, 0.23331176 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.00000000, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 },
                    { 0.73612750, 0.02738195, 0.01777184, 0.21871871 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.72801768, { 0.00000000, 0.38538472, 0.38540005, 0.22921523 },
                    { 0.00000000, 0.00000023, 0.00000000, 0.99999977 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.00000000, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 },
                    { 0.74504275, 0.00000000, 0.01613702, 0.23882023 } )
    )
  );

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlash,
  NegativeTwoPhaseFlash4CompSRK,
  ::testing::Values(
    FlashData< 4 >( 1.000000e+05, 1.931500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.05645957, { 0.00052723, 0.11109010, 0.11109434, 0.77728833 },
                    { 0.99989306, 0.00000003, 0.00000000, 0.00010691 } ),
    FlashData< 4 >( 1.000000e+05, 1.931500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.15640660, { 0.00067175, 0.12425180, 0.12425654, 0.75081991 },
                    { 0.99985756, 0.00000002, 0.00000000, 0.00014242 } ),
    FlashData< 4 >( 1.000000e+05, 1.931500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.79033121, { 0.00000000, 0.49992180, 0.49994088, 0.00013733 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 1.000000e+05, 1.931500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.10462698, { 0.00021850, 0.00000000, 0.11707076, 0.88271074 },
                    { 0.99995588, 0.00000000, 0.00000000, 0.00004412 } ),
    FlashData< 4 >( 5.000000e+06, 1.931500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.03575702, { 0.02198019, 0.10870496, 0.10870911, 0.76060573 },
                    { 0.99999256, 0.00000002, 0.00000000, 0.00000742 } ),
    FlashData< 4 >( 5.000000e+06, 1.931500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.13260833, { 0.02806574, 0.12084275, 0.12084737, 0.73024414 },
                    { 0.99999008, 0.00000002, 0.00000000, 0.00000990 } ),
    FlashData< 4 >( 1.000000e+05, 2.771500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.17767058, { 0.00055814, 0.12742716, 0.12746960, 0.74454511 },
                    { 0.88079876, 0.00017392, 0.00000000, 0.11902732 } ),
    FlashData< 4 >( 1.000000e+05, 2.771500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.78822635, { 0.00000000, 0.49495299, 0.49497188, 0.01007514 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 1.000000e+05, 2.771500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.11241759, { 0.00027478, 0.00000000, 0.11809830, 0.88162692 },
                    { 0.93022908, 0.00000000, 0.00000025, 0.06977066 } ),
    FlashData< 4 >( 1.000000e+08, 2.771500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.00000000, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 },
                    { 0.13237045, 0.10086691, 0.43406051, 0.33270214 } ),
    FlashData< 4 >( 1.000000e+08, 2.771500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.00000000, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 },
                    { 0.49094516, 0.08091778, 0.17876080, 0.24937626 } ),
    FlashData< 4 >( 1.000000e+08, 2.771500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.78896344, { 0.00000000, 0.49668171, 0.49670067, 0.00661762 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 1.013250e+05, 2.886500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 1.00000000, { 0.00000004, 0.00000000, 0.00000000, 0.99999996 },
                    { 0.99000000, 0.00000000, 0.00000000, 0.01000000 } ),
    FlashData< 4 >( 1.000000e+06, 2.886500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 0.99135986, { 0.00000042, 0.00000000, 0.00000000, 0.99999958 },
                    { 0.99862829, 0.00000000, 0.00000000, 0.00137171 } ),
    FlashData< 4 >( 1.000000e+06, 2.886500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.01001332, { 0.00000042, 0.00000000, 0.00000000, 0.99999958 },
                    { 0.99862829, 0.00000000, 0.00000000, 0.00137171 } ),
    FlashData< 4 >( 5.000000e+06, 2.886500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.03305889, { 0.02489908, 0.10840062, 0.10840577, 0.75829452 },
                    { 0.99443907, 0.00002967, 0.00000001, 0.00553126 } ),
    FlashData< 4 >( 1.000000e+07, 2.981500e+02, { 0.05695100, 0.10481800, 0.10482200, 0.73340900 }, true, 0.00905374, { 0.04838984, 0.10577509, 0.10577970, 0.74005537 },
                    { 0.99398426, 0.00006297, 0.00000006, 0.00595271 } ),
    FlashData< 4 >( 1.000000e+07, 2.981500e+02, { 0.15695100, 0.10481800, 0.10482200, 0.63340900 }, true, 0.10758197, { 0.05612207, 0.11744726, 0.11745840, 0.70897227 },
                    { 0.99335094, 0.00005531, 0.00000004, 0.00659371 } ),
    FlashData< 4 >( 1.000000e+07, 2.981500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.78640424, { 0.00000000, 0.49073071, 0.49074943, 0.01851986 },
                    { 0.00000000, 0.00000000, 0.00000000, 1.00000000 } ),
    FlashData< 4 >( 1.000000e+05, 3.331500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 1.00000000, { 0.00000022, 0.00000000, 0.00000000, 0.99999978 },
                    { 0.99000000, 0.00000000, 0.00000000, 0.01000000 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.00000000, 0.10481800, 0.10482200, 0.79036000 }, true, 0.72419599, { 0.00000000, 0.38004425, 0.38005974, 0.23989601 },
                    { 0.00000000, 0.00000038, 0.00000000, 0.99999962 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 }, true, 0.00000000, { 0.10481800, 0.00000000, 0.10482200, 0.79036000 },
                    { 0.84227242, 0.00000000, 0.00524021, 0.15248737 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.99000000, 0.00000000, 0.00000000, 0.01000000 }, true, 1.00000000, { 0.00566116, 0.00000000, 0.00000000, 0.99433884 },
                    { 0.99000000, 0.00000000, 0.00000000, 0.01000000 } ),
    FlashData< 4 >( 1.000000e+08, 4.731500e+02, { 0.01000000, 0.00000000, 0.00000000, 0.99000000 }, true, 0.00460663, { 0.00566116, 0.00000000, 0.00000000, 0.99433884 },
                    { 0.94752980, 0.00000000, 0.00000000, 0.05247020 } )
    )
  );

} // testing

} // geos
