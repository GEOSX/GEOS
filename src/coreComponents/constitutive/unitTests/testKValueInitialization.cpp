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
#include "common/DataTypes.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/KValueInitialization.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{
static constexpr real64 relTol = 1.0e-4;

// Parameters are
// 0 - pressure
// 1 - temperature
// 2 - index of k-value
// 3 - expected k-value
template< int NC >
class WilsonKValueInitializationTestFixture :
  public ::testing::TestWithParam< std::tuple< real64 const, real64 const, integer const, real64 const > >
{
public:
  using Feed = std::initializer_list< real64 const >;
public:
  WilsonKValueInitializationTestFixture( Feed const Pc, Feed const Tc, Feed const omega )
    : numComps( NC )
  {
    assign( criticalPressure, Pc );
    assign( criticalTemperature, Tc );
    assign( acentricFactor, omega );
  }
  virtual ~WilsonKValueInitializationTestFixture() = default;

  void testKValues()
  {
    real64 const pressure = std::get< 0 >( GetParam() );
    real64 const temperature = std::get< 1 >( GetParam() );
    integer const compIndex = std::get< 2 >( GetParam() );
    real64 const expectedKValue = std::get< 3 >( GetParam() );

    stackArray1d< real64, NC > kValues( numComps );

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        createKernelWrapper(),
                                                        kValues.toSlice() );

    ASSERT_EQ( kValues.size(), NC );
    checkRelativeError( expectedKValue, kValues[compIndex], relTol );
  }

private:
  void assign( array1d< real64 > & values, Feed const data )
  {
    for( Feed::const_iterator it = data.begin(); it != data.end(); ++it )
    {
      values.emplace_back( *it );
    }
  }

  ComponentProperties::KernelWrapper createKernelWrapper() const
  {
    return ComponentProperties::KernelWrapper(
      discarded,
      criticalPressure,
      criticalTemperature,
      acentricFactor,
      discarded,
      discarded2d );
  }

protected:
  const integer numComps;
  array1d< real64 > criticalPressure;
  array1d< real64 > criticalTemperature;
  array1d< real64 > acentricFactor;
  array1d< real64 > discarded;
  array2d< real64 > discarded2d;
};

class WilsonKValues2CompFixture : public WilsonKValueInitializationTestFixture< 2 >
{
public:
  WilsonKValues2CompFixture():
    WilsonKValueInitializationTestFixture< 2 >( {12.96e5, 45.99e5}, {33.15, 190.6}, {-0.219, 0.0114} )
  {}
};

class WilsonKValues9CompFixture : public WilsonKValueInitializationTestFixture< 9 >
{
public:
  WilsonKValues9CompFixture():
    WilsonKValueInitializationTestFixture< 9 >(
      {46.32700e5, 73.70000e5, 48.83900e5, 40.48100e5, 32.29200e5, 24.84900e5, 17.01500e5, 10.45300e5, 13.78000e5},
      {190.1580, 304.2020, 305.4000, 388.7960, 531.0700, 668.5570, 815.4070, 968.7040, 1105.1500},
      {0.0117, 0.2389, 0.0986, 0.1655, 0.2642, 0.3974, 0.6231, 1.0804, 1.15} )
  {}
};

TEST_P( WilsonKValues2CompFixture, testKValues )
{
  testKValues();
}

TEST_P( WilsonKValues9CompFixture, testKValues )
{
  testKValues();
}

// 2-component fluid test
INSTANTIATE_TEST_SUITE_P(
  WilsonKValues2CompTest,
  WilsonKValues2CompFixture,
  ::testing::Values(
    std::make_tuple( 1.01325e+5, 2.88650e+2, 0, 5.2375203367823417e+2 ),
    std::make_tuple( 1.01325e+5, 2.88650e+2, 1, 2.8719541295945191e+2 ),
    std::make_tuple( 1.00000e+6, 2.88650e+2, 0, 5.3069174812447081e+1 ),
    std::make_tuple( 1.00000e+6, 2.88650e+2, 1, 2.9100075218116466e+1 ),
    std::make_tuple( 5.00000e+6, 2.88650e+2, 0, 1.0613834962489415e+1 ),
    std::make_tuple( 5.00000e+6, 2.88650e+2, 1, 5.8200150436232931e+0 ),
    std::make_tuple( 1.00000e+8, 2.88650e+2, 0, 5.3069174812447081e-1 ),
    std::make_tuple( 1.00000e+8, 2.88650e+2, 1, 2.9100075218116466e-1 ),
    std::make_tuple( 1.01325e+5, 3.50000e+2, 0, 5.6989140865114678e+2 ),
    std::make_tuple( 1.01325e+5, 3.50000e+2, 1, 5.3850288275236221e+2 ),
    std::make_tuple( 1.00000e+6, 3.50000e+2, 0, 5.7744246981577454e+1 ),
    std::make_tuple( 1.00000e+6, 3.50000e+2, 1, 5.4563804594883109e+1 ),
    std::make_tuple( 5.00000e+6, 3.50000e+2, 0, 1.1548849396315489e+1 ),
    std::make_tuple( 5.00000e+6, 3.50000e+2, 1, 1.0912760918976621e+1 ),
    std::make_tuple( 1.00000e+8, 3.50000e+2, 0, 5.7744246981577454e-1 ),
    std::make_tuple( 1.00000e+8, 3.50000e+2, 1, 5.4563804594883109e-1 )
    ));

// 9-component fluid test
INSTANTIATE_TEST_SUITE_P(
  WilsonKValues9CompTest,
  WilsonKValues9CompFixture,
  ::testing::Values(
    std::make_tuple( 1.01325e+05, 2.88650e+02, 0, 2.91876325e+02 ),
    std::make_tuple( 1.01325e+05, 2.88650e+02, 1, 5.08252148e+01 ),
    std::make_tuple( 1.01325e+05, 2.88650e+02, 8, 8.91433739e-14 ),
    std::make_tuple( 1.00000e+06, 2.88650e+02, 0, 2.95743686e+01 ),
    std::make_tuple( 1.00000e+06, 2.88650e+02, 1, 5.14986489e+00 ),
    std::make_tuple( 1.00000e+06, 2.88650e+02, 8, 9.03245236e-15 ),
    std::make_tuple( 5.00000e+06, 2.88650e+02, 0, 5.91487373e+00 ),
    std::make_tuple( 5.00000e+06, 2.88650e+02, 1, 1.02997298e+00 ),
    std::make_tuple( 5.00000e+06, 2.88650e+02, 8, 1.80649047e-15 ),
    std::make_tuple( 1.00000e+08, 2.88650e+02, 0, 2.95743686e-01 ),
    std::make_tuple( 1.00000e+08, 2.88650e+02, 1, 5.14986489e-02 ),
    std::make_tuple( 1.00000e+08, 2.88650e+02, 8, 9.03245236e-17 ),
    std::make_tuple( 1.01325e+05, 3.50000e+02, 0, 5.46584215e+02 ),
    std::make_tuple( 1.01325e+05, 3.50000e+02, 1, 1.73708806e+02 ),
    std::make_tuple( 1.01325e+05, 3.50000e+02, 8, 2.06610518e-10 ),
    std::make_tuple( 1.00000e+06, 3.50000e+02, 0, 5.53826456e+01 ),
    std::make_tuple( 1.00000e+06, 3.50000e+02, 1, 1.76010447e+01 ),
    std::make_tuple( 1.00000e+06, 3.50000e+02, 8, 2.09348108e-11 ),
    std::make_tuple( 5.00000e+06, 3.50000e+02, 0, 1.10765291e+01 ),
    std::make_tuple( 5.00000e+06, 3.50000e+02, 1, 3.52020894e+00 ),
    std::make_tuple( 5.00000e+06, 3.50000e+02, 8, 4.18696215e-12 ),
    std::make_tuple( 1.00000e+08, 3.50000e+02, 0, 5.53826456e-01 ),
    std::make_tuple( 1.00000e+08, 3.50000e+02, 1, 1.76010447e-01 ),
    std::make_tuple( 1.00000e+08, 3.50000e+02, 8, 2.09348108e-13 )
    ));


// Parameters are
// 0 - pressure
// 1 - temperature
// 2 - expected k-value
class GasWaterKValueInitializationTestFixture :
  public ::testing::TestWithParam< std::tuple< real64 const, real64 const, real64 const > >
{
public:
  GasWaterKValueInitializationTestFixture() = default;
  virtual ~GasWaterKValueInitializationTestFixture() = default;

  void testKValues()
  {
    real64 const pressure = std::get< 0 >( GetParam() );
    real64 const temperature = std::get< 1 >( GetParam() );
    real64 const expectedKValue = std::get< 2 >( GetParam() );

    real64 calculatedKValue = KValueInitialization::computeWaterGasKvalue( pressure, temperature );
    checkRelativeError( expectedKValue, calculatedKValue, relTol );
  }
};

// gas-water test
TEST_P( GasWaterKValueInitializationTestFixture, testKValues )
{
  testKValues();
}

// 9-component fluid test
INSTANTIATE_TEST_SUITE_P(
  GasWaterKValueInitializationTest,
  GasWaterKValueInitializationTestFixture,
  ::testing::Values(
    std::make_tuple( 9.50000e+04, 2.88650e+02, 2.23565965e-02 ),
    std::make_tuple( 1.01325e+05, 2.88650e+02, 2.09610330e-02 ),
    std::make_tuple( 1.00000e+06, 2.88650e+02, 2.12387667e-03 ),
    std::make_tuple( 5.00000e+06, 2.88650e+02, 4.24775333e-04 ),
    std::make_tuple( 1.00000e+08, 2.88650e+02, 2.12387667e-05 ),
    std::make_tuple( 9.50000e+04, 3.50000e+02, 4.23601370e-01 ),
    std::make_tuple( 1.01325e+05, 3.50000e+02, 3.97158945e-01 ),
    std::make_tuple( 1.00000e+06, 3.50000e+02, 4.02421301e-02 ),
    std::make_tuple( 5.00000e+06, 3.50000e+02, 8.04842603e-03 ),
    std::make_tuple( 1.00000e+08, 3.50000e+02, 4.02421301e-04 ),
    std::make_tuple( 9.50000e+04, 3.73150e+02, 9.99692622e-01 ),
    std::make_tuple( 1.01325e+05, 3.73150e+02, 9.37288912e-01 ),
    std::make_tuple( 1.00000e+06, 3.73150e+02, 9.49707991e-02 ),
    std::make_tuple( 5.00000e+06, 3.73150e+02, 1.89941598e-02 ),
    std::make_tuple( 1.00000e+08, 3.73150e+02, 9.49707991e-04 )
    ));

} // testing

} // geos
