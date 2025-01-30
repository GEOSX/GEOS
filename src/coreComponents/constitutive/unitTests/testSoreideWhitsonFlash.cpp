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

// Source includes
#include "constitutive/fluid/multifluid/compositional/functions/SoreideWhitsonEOSModel.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

namespace geos
{

using namespace constitutive;
using namespace constitutive::compositional;

namespace testing
{

template< int NC >
using TestData = std::tuple<
  real64 const,         // Pressure
  real64 const,         // Temperature
  real64 const,         // Salinity
  Feed< NC > const      // Input composition
  >;

template< integer NC >
class SoreideWhitsonFlashTestFixture : public ::testing::TestWithParam< TestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  static constexpr integer numDof = NC + 2;
  static constexpr real64 absTol = 1.0e-4;
  static constexpr real64 relTol = 1.0e-5;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using ParamType = TestData< NC >;
public:
  SoreideWhitsonFlashTestFixture();
  ~SoreideWhitsonFlashTestFixture() = default;

  void testFlash( ParamType const & testData )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();
    real64 const pressure = std::get< 0 >( testData );
    real64 const temperature = std::get< 1 >( testData );
    real64 const salinity = std::get< 2 >( testData );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 3 >( testData ));
    
    GEOS_UNUSED_VAR(salinity);

    stackArray2d< real64, 2*numComps > phaseComposition(2, numComps);
    real64 vapourFraction = -1.0;
    arraySlice1d< real64 > liquidComposition = phaseComposition[0];
    arraySlice1d< real64 > vapourComposition = phaseComposition[1];

    stackArray2d< real64, numComps > kValues(1, numComps);

    auto printComposition = [](auto const & mf){
      std::ostringstream os;
      os << std::fixed << std::setprecision(5);
      os << "{" << mf[0];
      for (int ic = 1; ic < numComps; ic++)
      {
        os << ", " << mf[ic];
      }
      os << "}";
      return os.str();
    };

    bool flashStatus = NegativeTwoPhaseFlash::compute(numComps,
                       pressure,
                       temperature,
                       composition.toSliceConst(),
                       componentProperties,
                       EquationOfStateType::PengRobinson,
                       EquationOfStateType::SoreideWhitson,
                       kValues.toSlice(),
                       vapourFraction,
                       liquidComposition,
                       vapourComposition );

    std::cout << "TestData<" << numComps << ">{"
        << std::scientific << std::setprecision(2) << pressure << ", "
        << std::fixed << std::setprecision(2) << temperature << ", "
        << printComposition(composition) << ", "
        << flashStatus << ", "
        << std::scientific << std::setprecision(5) << vapourFraction << ", "
        << printComposition(liquidComposition) << ", "
        << printComposition(vapourComposition)
        << "},\n";
    flashStatus = NegativeTwoPhaseFlash::compute(numComps,
                       pressure,
                       temperature,
                       composition.toSliceConst(),
                       componentProperties,
                       EquationOfStateType::PengRobinson,
                       EquationOfStateType::PengRobinson,
                       kValues.toSlice(),
                       vapourFraction,
                       liquidComposition,
                       vapourComposition );
    std::cout << "TestData<" << numComps << ">{"
        << std::scientific << std::setprecision(2) << pressure << ", "
        << std::fixed << std::setprecision(2) << temperature << ", "
        << printComposition(composition) << ", "
        << flashStatus << ", "
        << std::scientific << std::setprecision(5) << vapourFraction << ", "
        << printComposition(liquidComposition) << ", "
        << printComposition(vapourComposition)
        << "},\n";
    std::cout << "TestData -----------------------------------------------------\n";
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

template< int NC >
struct FluidData {};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > create()
  {
    return TestFluid< 4 >::create( {Fluid::N2, Fluid::C1, Fluid::CO2, Fluid::H2O} );
  }
};

template< integer NC >
SoreideWhitsonFlashTestFixture< NC >::SoreideWhitsonFlashTestFixture():
  m_fluid( FluidData< NC >::create() )
{}

using SoreideWhitsonFlash4 = SoreideWhitsonFlashTestFixture< 4 >;

TEST_P( SoreideWhitsonFlash4, testFlash )
{
  testFlash( GetParam() );
}

template< int NC >
struct TestFeed {};

template<>
struct TestFeed< 4 >
{
  static std::array< Feed< 4 >, 1 > constexpr feeds = {
    Feed< 4 >{0., 0., 0.5, 0.5}
  };
};

template< int NC >
std::vector< TestData< NC > > generateTestData()
{
  std::array< real64 const, 2 > pressures( {1.83959e+06, 1.83959e+08} );
  std::array< real64 const, 2 > temperatures( {2.97150e+02, 3.63000e+02} );
  std::array< real64 const, 2 > salinities( {0.0} );
  std::vector< TestData< NC > > testData;
  for( const auto & composition : TestFeed< NC >::feeds )
  {
    for( const real64 pressure : pressures )
    {
      for( const real64 temperature : temperatures )
      {
        for( const real64 salinity : salinities )
        {
          testData.emplace_back( pressure, temperature, salinity, composition );
        }
      }
    }
  }
  return testData;
}

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonFlashTest, SoreideWhitsonFlash4, ::testing::ValuesIn( generateTestData< 4 >()) );

} // namespace testing

} // namespace geos
