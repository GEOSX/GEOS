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
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "unitTests/testingUtilities/TestingTasks.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseStatistics.hpp"
#include <hdf5.h>
#include "dataRepository/Group.hpp"
#include <filesystem>
#include "fileIO/Outputs/TimeHistoryOutput.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

namespace ComponentsGeneration
{


TEST( testStatOutputController, checkSinglePhaseFluxStatistics )
{
  string const testSet =
    R"xml(
<Problem>

  <Solvers>
    <SinglePhaseFVM name="testSolver"
                    discretization="singlePhaseTPFA"
                    targetRegions="{ reservoir }" >

      <NonlinearSolverParameters newtonMaxIter="40"
                                 allowNonConverged="1" />
      <LinearSolverParameters solverType="gmres"
                              preconditionerType="iluk"
                              krylovTol="1.0e-6" />

    </SinglePhaseFVM>
  </Solvers>

   <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="singlePhaseTPFA" />
    </FiniteVolume>
  </NumericalMethods>

  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="{ C3D8 }"
                  xCoords="{   0, 10 }"
                  yCoords="{   0, 10 }"
                  zCoords="{ -10,  0 }"
                  nx="{ 10 }"
                  ny="{ 10 }"
                  nz="{ 10 }"
                  cellBlockNames="{ cellBlock }" />
  </Mesh>

  <ElementRegions>
    <CellElementRegion name="reservoir"
                       cellBlocks="{ cellBlock }"
                       materialList="{ water, rock }" />
  </ElementRegions>

   <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  defaultDensity="1000"
                                  defaultViscosity="0.001"
                                  referencePressure="0.0"
                                  compressibility="5e-10"
                                  viscosibility="0.0" />

    <CompressibleSolidConstantPermeability name="rock"
                                           solidModelName="nullSolid"
                                           porosityModelName="rockPorosity"
                                           permeabilityModelName="rockPerm" />
    <NullModel name="nullSolid" />
    <PressurePorosity name="rockPorosity"
                      defaultReferencePorosity="0.05"
                      referencePressure="0.0"
                      compressibility="1.0e-9" />
    <ConstantPermeability name="rockPerm"
                          permeabilityComponents="{ 1.0e-12, 1.0e-12, 1.0e-15 }" />
  </Constitutive>

  <Geometry>
    <!-- source selects 2 elements -->
    <Box name="sourceBox"
         xMin="{ -0.01, -0.01, -10.01 }"
         xMax="{  2.01,  1.01,  -8.99 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{  4.99, 8.99, -1.01 }"
         xMax="{ 10.01, 10.01, 0.01 }" />
  </Geometry>

  <Events
    maxTime="1000.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />
    <PeriodicEvent
      name="statistics"
      timeFrequency="500"
      target="/Tasks/statController"/>
  </Events>

  <Functions>
    <!-- Unscaled injection / production rate in mol/s -->
    <TableFunction
      name="FluxRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0 }"
      values="{       0.000,  1.000,  0.767 }"
    />
  </Functions>

  <Tasks>

    <StatOutputController 
      name="statController"
      outputDirectory="./hdf/statistics" >
      <SinglePhaseStatistics
        name="singFlowStatistics"
        flowSolverName="testSolver"
        logLevel="1"
      />
    </StatOutputController>

  </Tasks>
  
  <!-- SPHINX_MESH -->  
  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="{ C3D8 }"
                  xCoords="{   0, 10 }"
                  yCoords="{   0, 10 }"
                  zCoords="{ -10,  0 }"
                  nx="{ 10 }"
                  ny="{ 10 }"
                  nz="{ 10 }"
                  cellBlockNames="{ cellBlock }" />
  </Mesh>
  <!-- SPHINX_MESH_END -->  
  
  <Geometry>
    <!-- source selects 2 elements -->
    <Box name="sourceBox"
         xMin="{ -0.01, -0.01, -10.01 }"
         xMax="{  2.01,  1.01,  -8.99 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{  4.99, 8.99, -1.01 }"
         xMax="{ 10.01, 10.01, 0.01 }" />
  </Geometry>
  
</Problem>

)xml";

  OutputBase::setOutputDirectory( "." );
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  setupProblemFromXML( problem, testSet.data() );

  std::vector< string > const refCollectorPaths =
  {
    "/Tasks/packCollectionreservoiraveragePressure",
    "/Tasks/packCollectionreservoirminPressure",
    "/Tasks/packCollectionreservoirmaxPressure",
    "/Tasks/packCollectionreservoirminDeltaPressure",
    "/Tasks/packCollectionreservoirmaxDeltaPressure",
    "/Tasks/packCollectionreservoirtotalMass",
    "/Tasks/packCollectionreservoiraverageTemperature",
    "/Tasks/packCollectionreservoirminTemperature",
    "/Tasks/packCollectionreservoirmaxTemperature",
    "/Tasks/packCollectionreservoirtotalPoreVolume",
    "/Tasks/packCollectionreservoirtotalUncompactedPoreVolume",
  };
  string const outputPath = "/Outputs/compFlowHistory_reservoir";

  { // verify component creation

    for( string const & path : refCollectorPaths )
    {
      EXPECT_NO_THROW( {
          Group const & group = problem.getGroupByPath( path );
          ASSERT_STREQ( path.c_str(), group.getPath().c_str() );
        } );
    }

    EXPECT_NO_THROW( {
        Group const & group = problem.getGroupByPath( outputPath );
        ASSERT_STREQ( outputPath.c_str(), group.getPath().c_str() );
      } );
  }

  { // check all timeHistory paths
    TimeHistoryOutput & timeHistory = problem.getGroupByPath< TimeHistoryOutput >( "/Outputs/compFlowHistory_reservoir" );
    string_array & collectorPaths =  timeHistory.getReference< string_array >( TimeHistoryOutput::viewKeys::timeHistoryOutputTargetString() );

    for( size_t idxPath = 0; idxPath < refCollectorPaths.size(); idxPath++ )
    {
      ASSERT_STREQ( refCollectorPaths[idxPath].c_str(), collectorPaths[idxPath].c_str() );
    }
  }


  // run simulation
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  string fileFromPath = "./hdf/statistics/compFlowHistory_reservoir.hdf5";
  int64_t fileId = H5Fopen( fileFromPath.c_str(), H5F_ACC_RDWR, 0 );
  EXPECT_TRUE( fileId != -1 );
  H5Fclose( fileId );

  integer status = std::remove( fileFromPath.c_str() );
  EXPECT_TRUE( status == 0 );

}
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
