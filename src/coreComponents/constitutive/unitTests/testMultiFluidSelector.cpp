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
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

template< typename FluidType >
class MultiFluidSelectorTest : public ::testing::Test
{
public:
  MultiFluidSelectorTest():
    m_node(),
    m_parent( "parent", m_node )
  {
    m_model = &m_parent.registerGroup< FluidType >( "fluid" );
  }
  ~MultiFluidSelectorTest() override = default;

  MultiFluidBase & getFluid() const { return *m_model; }

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  FluidType * m_model{};
};

using MultiFluidSelectorTestDeadOilFluid = MultiFluidSelectorTest< DeadOilFluid >;
using MultiFluidSelectorTestCO2BrinePhillipsThermalFluid = MultiFluidSelectorTest< CO2BrinePhillipsThermalFluid >;
using MultiFluidSelectorTestCompositionalTwoPhaseConstantViscosity = MultiFluidSelectorTest< CompositionalTwoPhaseConstantViscosity >;

TEST_F( MultiFluidSelectorTestDeadOilFluid, testValidComponents )
{
  bool isExecuted = false;
  constitutiveComponentUpdatePassThru( getFluid(), 2, [&]( auto &, auto NC )
  {
    integer constexpr numComps = NC();
    EXPECT_EQ( numComps, 2 );
    isExecuted = true;
  } );
  EXPECT_TRUE( isExecuted );

  isExecuted = false;
  constitutiveComponentUpdatePassThru( getFluid(), 3, [&]( auto &, auto NC )
  {
    integer constexpr numComps = NC();
    EXPECT_EQ( numComps, 3 );
    isExecuted = true;
  } );
  EXPECT_TRUE( isExecuted );
}

TEST_F( MultiFluidSelectorTestDeadOilFluid, testInvalidComponents )
{
  EXPECT_THROW( constitutiveComponentUpdatePassThru( getFluid(), 1, []( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );

  EXPECT_THROW( constitutiveComponentUpdatePassThru( getFluid(), 4, []( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );
}

TEST_F( MultiFluidSelectorTestDeadOilFluid, testThermal )
{
  EXPECT_THROW( constitutiveComponentUpdatePassThru< true >( getFluid(), 2, []( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );
}

TEST_F( MultiFluidSelectorTestCO2BrinePhillipsThermalFluid, testValidComponents )
{
  bool isExecuted = false;
  constitutiveComponentUpdatePassThru( getFluid(), 2, [&]( auto &, auto NC )
  {
    integer constexpr numComps = NC();
    EXPECT_EQ( numComps, 2 );
    isExecuted = true;
  } );
  EXPECT_TRUE( isExecuted );
}

TEST_F( MultiFluidSelectorTestCO2BrinePhillipsThermalFluid, testInvalidComponents )
{
  EXPECT_THROW( constitutiveComponentUpdatePassThru( getFluid(), 1, []( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );

  EXPECT_THROW( constitutiveComponentUpdatePassThru( getFluid(), 3, []( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );
}

TEST_F( MultiFluidSelectorTestCO2BrinePhillipsThermalFluid, testThermal )
{
  bool isExecuted = false;
  constitutiveComponentUpdatePassThru< true >( getFluid(), 2, [&]( auto &, auto )
  {
    isExecuted = true;
  } );
  EXPECT_TRUE( isExecuted );
}

TEST_F( MultiFluidSelectorTestCompositionalTwoPhaseConstantViscosity, testValidComponents )
{
  for( integer nc = 2; nc <= 5; nc++ )
  {
    bool isExecuted = false;
    constitutiveComponentUpdatePassThru( getFluid(), nc, [&]( auto &, auto NC )
    {
      integer constexpr numComps = NC();
      EXPECT_EQ( numComps, nc );
      isExecuted = true;
    } );
    EXPECT_TRUE( isExecuted );
  }
}

TEST_F( MultiFluidSelectorTestCompositionalTwoPhaseConstantViscosity, testInvalidComponents )
{
  EXPECT_THROW( constitutiveComponentUpdatePassThru( getFluid(), 1, []( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );

  EXPECT_THROW( constitutiveComponentUpdatePassThru( getFluid(), 6, []( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );
}

TEST_F( MultiFluidSelectorTestCompositionalTwoPhaseConstantViscosity, testThermal )
{
  EXPECT_THROW( constitutiveComponentUpdatePassThru< true >( getFluid(), 2, [&]( auto &, auto )
  {
    FAIL(); // Shouldn't be called
  } ), InputError );
}
