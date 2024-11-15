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
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterFlashModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

constexpr integer numTestComps = 3;

template< integer NC >
using FlashData = std::tuple<
  real64 const,                 // pressure
  real64 const,                 // temperature
  Feed< NC > const,             // phase composition
  real64 const,                 // expected liquid fraction
  real64 const,                 // expected vapour fraction
  real64 const,                 // expected aqueous fraction
  Feed< numTestComps > const,   // expected liquid mole fractions
  Feed< numTestComps > const,   // expected vapour mole fractions
  Feed< numTestComps > const    // expected aqueous mole fractions
  >;

template< integer NC >
struct FluidData {};

template<>
struct FluidData< 3 >
{
  static constexpr integer testComponents[numTestComps] = {0, 1, 2};
  static std::unique_ptr< TestFluid< 3 > > createFluid()
  {
    auto fluid = TestFluid< 3 >::create( {Fluid::C1, Fluid::C10, Fluid::H2O} );
    const std::array< real64 const, 3 > bics = { 0.25, 0.0, 0.0 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template<>
struct FluidData< 9 >
{
  static constexpr integer testComponents[numTestComps] = {0, 2, 8};
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

template< integer NC >
class ImmiscibleWaterFlashModelTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr integer numPhases = 3;
  static constexpr integer numComps = NC;
  static constexpr integer numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using PhasePropSlice = ImmiscibleWaterFlashModelUpdate::PhaseProp::SliceType;
  using PhaseCompSlice = ImmiscibleWaterFlashModelUpdate::PhaseComp::SliceType;

public:
  ImmiscibleWaterFlashModelTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = ImmiscibleWaterFlashModel::createParameters( std::move( m_parameters ) );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    m_flash = std::make_unique< ImmiscibleWaterFlashModel >( "FlashModel", componentProperties, *m_parameters );
  }

  ~ImmiscibleWaterFlashModelTestFixture() = default;

  void testFlash( FlashData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    real64 const expectedPhaseFraction[numPhases] = {std::get< 3 >( data ), std::get< 4 >( data ), std::get< 5 >( data )};
    Feed< numTestComps > const expectedPhaseComponentFraction[numPhases] = {std::get< 6 >( data ), std::get< 7 >( data ), std::get< 8 >( data )    };

    stackArray2d< real64, (numPhases-1)*numComps > kValues( numPhases-1, numComps );
    LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

    StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
    StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDofs );
    StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > phaseComponentFractionData( 1, 1, numPhases, numComps );
    StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

    auto phaseFraction = phaseFractionData[0][0];
    auto dPhaseFraction = dPhaseFractionData[0][0];
    auto phaseComponentFraction = phaseComponentFractionData[0][0];
    auto dPhaseComponentFraction = dPhaseComponentFractionData[0][0];

    auto componentProperties = m_fluid->createKernelWrapper();
    auto flashKernelWrapper = m_flash->createKernelWrapper();

    flashKernelWrapper.compute( componentProperties,
                                pressure,
                                temperature,
                                composition.toSliceConst(),
                                kValues.toSlice(),
                                PhasePropSlice( phaseFraction, dPhaseFraction ),
                                PhaseCompSlice( phaseComponentFraction, dPhaseComponentFraction ) );

    for( integer ip = 0; ip < numPhases; ip++ )
    {
      checkRelativeError( phaseFraction[ip], expectedPhaseFraction[ip], relTol, absTol );
      for( integer i = 0; i < numTestComps; ++i )
      {
        integer const ic = FluidData< numComps >::testComponents[i];
        checkRelativeError( phaseComponentFraction[ip][ic], expectedPhaseComponentFraction[ip][i], relTol, absTol );
      }
    }
  }

  void testFlashDerivatives( FlashData< NC > const & data )
  {
    // Number of output values from each flash calculation
    constexpr integer numValues = numPhases * (1 + numComps);

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    stackArray2d< real64, (numPhases-1)*numComps > kValues( numPhases-1, numComps );
    LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

    StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
    StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDofs );
    StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > phaseComponentFractionData( 1, 1, numPhases, numComps );
    StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

    auto phaseFraction = phaseFractionData[0][0];
    auto dPhaseFraction = dPhaseFractionData[0][0];
    auto phaseComponentFraction = phaseComponentFractionData[0][0];
    auto dPhaseComponentFraction = dPhaseComponentFractionData[0][0];

    stackArray1d< real64, numValues > derivatives( numValues );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto flashKernelWrapper = m_flash->createKernelWrapper();

    flashKernelWrapper.compute( componentProperties,
                                pressure,
                                temperature,
                                composition.toSliceConst(),
                                kValues.toSlice(),
                                PhasePropSlice( phaseFraction, dPhaseFraction ),
                                PhaseCompSlice( phaseComponentFraction, dPhaseComponentFraction ) );

    // Combine derivatives into a single output
    auto const concatDerivatives = []( integer const kc, auto & derivs, auto const & phaseFractionDerivs, auto const & phaseComponentFractionDerivs ){
      integer j = 0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        derivs[j++] = phaseFractionDerivs( ip, kc );
        for( integer ic = 0; ic < numComps; ++ic )
        {
          derivs[j++] = phaseComponentFractionDerivs( ip, ic, kc );
        }
      }
    };

    auto const evaluateFlash = [&]( real64 const p, real64 const t, auto const & zmf, auto & values ){
      StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > displacedPhaseFractionData( 1, 1, numPhases );
      StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > displacedPhaseFractionDerivsData( 1, 1, numPhases, numDofs );
      StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > displacedPhaseComponentFractionData( 1, 1, numPhases, numComps );
      StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > displacedPhaseComponentFractionDerivsData( 1, 1, numPhases, numComps, numDofs );

      auto displacedPhaseFraction = displacedPhaseFractionData[0][0];
      auto displacedPhaseFractionDerivs = displacedPhaseFractionDerivsData[0][0];
      auto displacedPhaseComponentFraction = displacedPhaseComponentFractionData[0][0];
      auto displacedPhaseComponentFractionDerivs = displacedPhaseComponentFractionDerivsData[0][0];

      flashKernelWrapper.compute( componentProperties,
                                  p,
                                  t,
                                  zmf,
                                  kValues.toSlice(),
                                  PhasePropSlice( displacedPhaseFraction, displacedPhaseFractionDerivs ),
                                  PhaseCompSlice( displacedPhaseComponentFraction, displacedPhaseComponentFractionDerivs ) );
      integer j = 0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        values[j++] = displacedPhaseFraction[ip];
        for( integer ic = 0; ic < numComps; ++ic )
        {
          values[j++] = displacedPhaseComponentFraction( ip, ic );
        }
      }
    };

    // Test against numerically calculated values
    // --- Pressure derivatives ---
    concatDerivatives( Deriv::dP, derivatives, dPhaseFraction, dPhaseComponentFraction );
    real64 const dp = 1.0e-4 * pressure;
    geos::testing::internal::testNumericalDerivative< numValues >(
      pressure, dp, derivatives,
      [&]( real64 const p, auto & values ) {
      evaluateFlash( p, temperature, composition.toSliceConst(), values );
    } );

    // -- Temperature derivative
    concatDerivatives( Deriv::dT, derivatives, dPhaseFraction, dPhaseComponentFraction );
    real64 const dT = 1.0e-6 * temperature;
    geos::testing::internal::testNumericalDerivative< numValues >(
      temperature, dT, derivatives,
      [&]( real64 const t, auto & values ) {
      evaluateFlash( pressure, t, composition.toSliceConst(), values );
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer const ic : FluidData< numComps >::testComponents )
    {
      real64 sumZ = 0.0;
      for( integer jc = 0; jc < numComps; ++jc )
      {
        sumZ += composition[jc];
      }
      sumZ -= composition[ic];
      if( sumZ < absTol )
      {
        continue;
      }
      concatDerivatives( Deriv::dC+ic, derivatives, dPhaseFraction, dPhaseComponentFraction );
      geos::testing::internal::testNumericalDerivative< numValues >(
        0.0, dz, derivatives,
        [&]( real64 const z, auto & values ) {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = composition[jc];
        }
        zmf[ic] += z;
        evaluateFlash( pressure, temperature, zmf.toSliceConst(), values );
      } );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< ImmiscibleWaterFlashModel > m_flash{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using ImmiscibleWaterFlashModel3 = ImmiscibleWaterFlashModelTestFixture< 3 >;
using ImmiscibleWaterFlashModel9 = ImmiscibleWaterFlashModelTestFixture< 9 >;

TEST_P( ImmiscibleWaterFlashModel3, testFlash )
{
  testFlash( GetParam() );
}
TEST_P( ImmiscibleWaterFlashModel3, testFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

TEST_P( ImmiscibleWaterFlashModel9, testFlash )
{
  testFlash( GetParam() );
}
TEST_P( ImmiscibleWaterFlashModel9, testFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterFlashModel, ImmiscibleWaterFlashModel3,
  ::testing::Values(
    FlashData<3>{1.0e+05, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 293.15, {0.300000, 0.300000, 0.400000}, 0.300217, 0.299783, 0.400000, {0.000723, 0.999277, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 293.15, {0.300000, 0.300000, 0.400000}, 0.302162, 0.297838, 0.400000, {0.007157, 0.992843, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 293.15, {0.300000, 0.300000, 0.400000}, 0.320811, 0.279189, 0.400000, {0.064871, 0.935129, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 313.15, {0.300000, 0.300000, 0.400000}, 0.300236, 0.299764, 0.400000, {0.000787, 0.999213, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 313.15, {0.300000, 0.300000, 0.400000}, 0.302354, 0.297646, 0.400000, {0.007785, 0.992215, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 313.15, {0.300000, 0.300000, 0.400000}, 0.322797, 0.277203, 0.400000, {0.070623, 0.929377, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 353.15, {0.300000, 0.300000, 0.400000}, 0.300271, 0.299729, 0.400000, {0.000919, 0.999081, 0.000000}, {0.999982, 0.000018, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 353.15, {0.300000, 0.300000, 0.400000}, 0.302753, 0.297247, 0.400000, {0.009094, 0.990906, 0.000000}, {0.999998, 0.000002, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 353.15, {0.300000, 0.300000, 0.400000}, 0.326941, 0.273059, 0.400000, {0.082404, 0.917596, 0.000000}, {0.999999, 0.000001, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 293.15, {0.200000, 0.800000, 0.000000}, 0.800579, 0.199421, 0.000000, {0.000723, 0.999277, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 293.15, {0.200000, 0.800000, 0.000000}, 0.805767, 0.194233, 0.000000, {0.007157, 0.992843, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 293.15, {0.200000, 0.800000, 0.000000}, 0.855497, 0.144503, 0.000000, {0.064871, 0.935129, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 313.15, {0.200000, 0.800000, 0.000000}, 0.800630, 0.199370, 0.000000, {0.000787, 0.999213, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 313.15, {0.200000, 0.800000, 0.000000}, 0.806277, 0.193723, 0.000000, {0.007785, 0.992215, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 313.15, {0.200000, 0.800000, 0.000000}, 0.860791, 0.139209, 0.000000, {0.070623, 0.929377, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 353.15, {0.200000, 0.800000, 0.000000}, 0.800732, 0.199268, 0.000000, {0.000919, 0.999081, 0.000000}, {0.999982, 0.000018, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 353.15, {0.200000, 0.800000, 0.000000}, 0.807341, 0.192659, 0.000000, {0.009094, 0.990906, 0.000000}, {0.999998, 0.000002, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 353.15, {0.200000, 0.800000, 0.000000}, 0.871843, 0.128157, 0.000000, {0.082404, 0.917596, 0.000000}, {0.999999, 0.000001, 0.000000}, {0.000000, 0.000000, 1.000000}}
  )
);

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterFlashModel, ImmiscibleWaterFlashModel9,
  ::testing::Values(
    FlashData<9>{1.0e+06, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.198828, 0.792172, 0.009000, {0.000000, 0.011551, 0.850986}, {0.000000, 0.672081, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.307007, 0.683993, 0.009000, {0.000000, 0.111647, 0.551128}, {0.000000, 0.731621, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.539556, 0.170737}, {0.000000, 0.539556, 0.170737}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.190493, 0.800507, 0.009000, {0.000000, 0.010985, 0.888223}, {0.000000, 0.665337, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.288620, 0.702380, 0.009000, {0.000000, 0.106713, 0.586237}, {0.000000, 0.717418, 0.000001}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.539556, 0.170737}, {0.000000, 0.539556, 0.170737}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.181787, 0.809213, 0.009000, {0.000000, 0.010453, 0.930750}, {0.000000, 0.658417, 0.000003}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.261660, 0.729340, 0.009000, {0.000000, 0.101277, 0.646622}, {0.000000, 0.696794, 0.000007}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.947424, 0.043576, 0.009000, {0.000000, 0.529993, 0.178339}, {0.000000, 0.747468, 0.005438}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.120388, 0.479612, 0.400000, {0.000000, 0.011551, 0.850996}, {0.000000, 0.672084, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.185888, 0.414112, 0.400000, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.539550, 0.170750}, {0.000000, 0.539550, 0.170750}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.115342, 0.484658, 0.400000, {0.000000, 0.010985, 0.888229}, {0.000000, 0.665341, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.174755, 0.425245, 0.400000, {0.000000, 0.106714, 0.586247}, {0.000000, 0.717425, 0.000001}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.539550, 0.170750}, {0.000000, 0.539550, 0.170750}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.110071, 0.489929, 0.400000, {0.000000, 0.010453, 0.930753}, {0.000000, 0.658421, 0.000003}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.158432, 0.441568, 0.400000, {0.000000, 0.101278, 0.646630}, {0.000000, 0.696800, 0.000007}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.573644, 0.026356, 0.400000, {0.000000, 0.529997, 0.178345}, {0.000000, 0.747479, 0.005437}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}}
  )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
