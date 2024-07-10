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
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/StabilityTest.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

static constexpr integer numComps = 2;

using StabilityData = std::tuple<
  real64 const,   // pressure
  real64 const,   // temperature
  real64 const,   // CH4 mole fraction
  real64 const    // expected tangent plane distance
  >;

template< EquationOfStateType EOS_TYPE >
class StabilityTestTest2CompFixture :  public ::testing::TestWithParam< StabilityData >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;

public:
  StabilityTestTest2CompFixture()
    : m_fluid( createFluid() )
  {}

  ~StabilityTestTest2CompFixture() = default;

  void testStability( StabilityData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    real64 const zCH4 = std::get< 2 >( data );
    real64 const expectedTangentPlaneDistance = std::get< 3 >( data );

    stackArray1d< real64, numComps > composition( numComps );
    composition[0] = zCH4;
    composition[1] = 1.0 - zCH4;

    real64 tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;
    stackArray1d< real64, numComps > kValues( numComps );

    bool const stabilityStatus = StabilityTest::compute( numComps,
                                                         pressure,
                                                         temperature,
                                                         composition.toSliceConst(),
                                                         componentProperties,
                                                         EOS_TYPE,
                                                         tangentPlaneDistance,
                                                         kValues.toSlice() );

    // Expect this to succeed
    ASSERT_EQ( stabilityStatus, true );

    // Check the tanget plane distance
    checkRelativeError( expectedTangentPlaneDistance, tangentPlaneDistance, relTol, absTol );
  }

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};

private:
  static std::unique_ptr< TestFluid< numComps > > createFluid()
  {
    std::unique_ptr< TestFluid< numComps > > fluid = TestFluid< numComps >::create( {Fluid::C1, Fluid::C3} );
    fluid->setBinaryCoefficients( Feed< 1 >{ 0.1 } );
    return fluid;
  }
};

using PengRobinson = StabilityTestTest2CompFixture< EquationOfStateType::PengRobinson >;
using SoaveRedlichKwong = StabilityTestTest2CompFixture< EquationOfStateType::SoaveRedlichKwong >;
TEST_P( PengRobinson, testStabilityTest )
{
  testStability( GetParam() );
}
TEST_P( SoaveRedlichKwong, testStabilityTest )
{
  testStability( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  StabilityTest, PengRobinson,
  ::testing::Values(
    StabilityData{ 1.00000e+06, 297.15, 0.2,   1.1102230e-16 },
    StabilityData{ 1.00000e+06, 353.15, 0.2,  -2.2204460e-16 },
    StabilityData{ 5.00000e+06, 297.15, 0.2,  -1.0438517e+00 },
    StabilityData{ 5.00000e+06, 353.15, 0.2,  -2.0228439e-03 },
    StabilityData{ 2.00000e+07, 297.15, 0.2,  -3.3306691e-16 },
    StabilityData{ 2.00000e+07, 353.15, 0.2,  -6.6613381e-16 }    
  )
);

INSTANTIATE_TEST_SUITE_P(
  StabilityTest, SoaveRedlichKwong,
  ::testing::Values(
    StabilityData{ 1.00000e+06, 297.15, 0.2,  -2.2204460e-16 },
    StabilityData{ 1.00000e+06, 353.15, 0.2,  -3.3306691e-16 },
    StabilityData{ 5.00000e+06, 297.15, 0.2,  -1.0966021e+00 },
    StabilityData{ 5.00000e+06, 353.15, 0.2,  -3.6652014e-03 },
    StabilityData{ 2.00000e+07, 297.15, 0.2,  -2.4424907e-15 },
    StabilityData{ 2.00000e+07, 353.15, 0.2,  -7.7715612e-16 }
  )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
