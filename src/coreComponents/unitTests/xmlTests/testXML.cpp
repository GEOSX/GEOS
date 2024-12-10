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
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

TEST( testXML, testXMLString )
{
  char const * xmlInput =
    R"xml(
    <?xml version="1.0" ?>
    <Problem>
      <Solvers>
        <SolidMechanics_LagrangianFEM
          name="lagsolve"
          cflFactor="0.25"
          discretization="FE1"
          targetRegions="{ Region2 }"/>
      </Solvers>
      <Mesh>
        <InternalMesh
          name="mesh1"
          elementTypes="{ C3D8 }"
          xCoords="{ 0, 3 }"
          yCoords="{ 0, 1 }"
          zCoords="{ 0, 1 }"
          nx="{ 4 }"
          ny="{ 1 }"
          nz="{ 1 }"
          cellBlockNames="{ cb1 }"/>
      </Mesh>
      <Events
        maxTime="1.0e-3">
        <PeriodicEvent
          name="solverApplications"
          forceDt="1.0e-3"
          target="/Solvers/lagsolve"/>
      </Events>
      <NumericalMethods>
        <FiniteElements>
          <FiniteElementSpace
            name="FE1"
            order="1"/>
        </FiniteElements>
      </NumericalMethods>
      <ElementRegions>
        <CellElementRegion
          name="Region2"
          cellBlocks="{ * }"
          materialList="{ shale }"/>
      </ElementRegions>
      <Constitutive>
        <ElasticIsotropic
          name="shale"
          defaultDensity="2700"
          defaultBulkModulus="5.5556e9"
          defaultShearModulus="4.16667e9"/>
      </Constitutive>
    </Problem>
    )xml";

  geos::getGlobalState().getProblemManager().parseInputString( xmlInput );
}

TEST( testXML, testXMLStringExpectedFail )
{
  char const * xmlInput =
    R"xml(
    <?xml version="1.0" ?>
    <Problem>
      <Solvers>
    </Problem>
    )xml";

  EXPECT_THROW( geos::getGlobalState().getProblemManager().parseInputString( xmlInput ), geos::InputError );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv, false ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
