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

#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "unitTests/testingUtilities/TestingTasks.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/fluidFlow/SourceFluxStatistics.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseStatistics.hpp"

#include <gtest/gtest.h>


using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;


//////////////////////////////// Test base utilities ////////////////////////////////


/**
 * @brief this struct is used to provide the input data of each flow tests
 */
struct TestInputs
{
  string xmlInput;
  std::map< string, string > tableFiles;

  string sourceFluxName;
  string sinkFluxName;
  string timeStepCheckerPath;
  string timeStepFluxStatsPath;
  string wholeSimFluxStatsPath;
  string flowSolverPath;

  // rates for each timesteps, for each phases
  array2d< real64 > sourceRates;
  array2d< real64 > sinkRates;

  // parameters for precomputing results
  real64 dt;
  real64 sourceRateFactor;
  real64 sinkRateFactor;
  integer sourceElementsCount;
  integer sinkElementsCount;

  /// In order to be sure that sub-timestepping is supported, requires the test to have at least one sub-timestep.
  /// At least one simulation should test the timestep cuts !
  integer requiredSubTimeStep = 0;
};

/**
 * @brief this struct computes from the test inputs the values to expect from the simulation.
 */
struct TestSet
{
  TestInputs const inputs;

  integer timestepCount;
  integer totalElementsCount;
  integer phaseCount;

  // stats for each timesteps, for each phases
  array2d< real64 > sourceRates;
  array2d< real64 > sinkRates;
  array2d< real64 > sourceMassProd;
  array2d< real64 > sinkMassProd;
  array2d< real64 > massDeltas;

  // overall simulation stats for each phases
  array1d< real64 > sourceMeanRate;
  array1d< real64 > sinkMeanRate;
  array1d< real64 > totalSourceMassProd;
  array1d< real64 > totalSinkMassProd;
  array1d< real64 > totalMassProd;
  array1d< real64 > totalMeanRate;

  /**
   * @brief Compute the expected statistics set for the tested simulation.
   * @param inputParams the test simulation input parameters
   */
  TestSet( TestInputs const & inputParams ):
    inputs( inputParams )
  {
    // tables must provide the same timestep & phase rates
    EXPECT_EQ( inputs.sourceRates.size( 0 ), inputs.sinkRates.size( 0 ));
    EXPECT_EQ( inputs.sourceRates.size( 1 ), inputs.sinkRates.size( 1 ));

    timestepCount = inputs.sourceRates.size( 0 );
    phaseCount = inputs.sourceRates.size( 1 );
    totalElementsCount = inputs.sourceElementsCount + inputs.sinkElementsCount;

    sourceRates.resize( timestepCount, phaseCount );
    sinkRates.resize( timestepCount, phaseCount );
    sourceMassProd.resize( timestepCount, phaseCount );
    sinkMassProd.resize( timestepCount, phaseCount );
    massDeltas.resize( timestepCount, phaseCount );
    sourceMeanRate.resize( phaseCount );
    sinkMeanRate.resize( phaseCount );
    totalSourceMassProd.resize( phaseCount );
    totalSinkMassProd.resize( phaseCount );
    totalMassProd.resize( phaseCount );
    totalMeanRate.resize( phaseCount );

    for( integer ip = 0; ip < phaseCount; ++ip )
    {
      for( integer timestepId = 0; timestepId < timestepCount; ++timestepId )
      {
        // mass production / injection calculation
        sourceRates[timestepId][ip] = inputs.sourceRates[timestepId][ip] * inputs.sourceRateFactor;
        sourceMassProd[timestepId][ip] = inputs.sourceRates[timestepId][ip] * inputs.dt * inputs.sourceRateFactor;
        totalSourceMassProd[ip] += sourceMassProd[timestepId][ip];
        sinkRates[timestepId][ip] = inputs.sinkRates[timestepId][ip] * inputs.sinkRateFactor;
        sinkMassProd[timestepId][ip] = inputs.sinkRates[timestepId][ip] * inputs.dt * inputs.sinkRateFactor;
        massDeltas[timestepId][ip] = -( sourceMassProd[timestepId][ip] + sinkMassProd[timestepId][ip] );
        totalSinkMassProd[ip] += sinkMassProd[timestepId][ip];
        // rates accumulations
        sourceMeanRate[ip] += inputs.sourceRates[timestepId][ip] * inputs.sourceRateFactor;
        sinkMeanRate[ip] += inputs.sinkRates[timestepId][ip] * inputs.sinkRateFactor;
      }
      // mean rates calculation
      real64 const ratesMeanDivisor = 1.0 / double( timestepCount - 1 );
      sourceMeanRate[ip] *= ratesMeanDivisor;
      sinkMeanRate[ip] *= ratesMeanDivisor;
      // totals
      totalMassProd[ip] = totalSinkMassProd[ip] + totalSourceMassProd[ip];
      totalMeanRate[ip] = sinkMeanRate[ip] + sourceMeanRate[ip];
    }
  }
};


class FlowStatisticsTest : public ::testing::Test
{
public:

  void writeTableFiles( std::map< string, string > const & files )
  {
    for( auto const & [fileName, content] : files )
    {
      std::ofstream os( fileName );
      ASSERT_TRUE( os.is_open() );
      os << content;
      os.close();

      m_tableFileNames.push_back( fileName );
    }
  }

  void TearDown() override
  {
    // removing temp table files
    for( string const & fileName : m_tableFileNames )
    {
      ASSERT_TRUE( std::remove( fileName.c_str() ) == 0 );
    }
    m_tableFileNames.clear();
  }

private:
  std::vector< string > m_tableFileNames;
};


void setRateTable( array2d< real64 > & rateTable, std::initializer_list< std::initializer_list< real64 > > timestepPhaseValues )
{
  rateTable.resize( timestepPhaseValues.size(), timestepPhaseValues.begin()->size() );
  integer timestepId = 0;
  for( auto const & phaseValues : timestepPhaseValues )
  {
    integer ip = 0;
    for( auto const & phaseValue : phaseValues )
    {
      rateTable[timestepId][ip++] = phaseValue;
    }
    ++timestepId;
  }
}

real64 getTotalFluidMass( ProblemManager & problem, string_view flowSolverPath )
{
  real64 totalMass = 0.0;
  PhysicsSolverBase const & solver = problem.getGroupByPath< PhysicsSolverBase >( string( flowSolverPath ) );
  solver.forDiscretizationOnMeshTargets( problem.getDomainPartition().getMeshBodies(),
                                         [&] ( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & )
  {
    mesh.getElemManager().forElementRegions( [&]( ElementRegionBase & region )
    {
      SinglePhaseStatistics::RegionStatistics & regionStats = region.getReference< SinglePhaseStatistics::RegionStatistics >(
        SinglePhaseStatistics::viewKeyStruct::regionStatisticsString() );

      totalMass += regionStats.totalMass;
    } );
  } );
  return totalMass;
}


/**
 * @brief Verification that the source flux statistics are correct for the current timestep
 * @param expectedMasses the expected mass values per phase
 * @param expectedRates the expected rate values per phase
 * @param expectedElementCount the number of expected targeted elements
 * @param stats the timestep stats
 * @param context a context string to provide in any error message.
 */
void checkFluxStats( arraySlice1d< real64 > const & expectedMasses,
                     arraySlice1d< real64 > const & expectedRates,
                     integer const expectedElementCount,
                     SourceFluxStatsAggregator::WrappedStats const & stats,
                     string_view context )
{
  for( int ip = 0; ip < stats.stats().getPhaseCount(); ++ip )
  {
    EXPECT_DOUBLE_EQ( stats.stats().m_producedMass[ip], expectedMasses[ip] ) << GEOS_FMT( "The flux named '{}' did not produce the expected mass ({}, phase = {}).",
                                                                                          stats.getFluxName(), context, ip );
    EXPECT_DOUBLE_EQ( stats.stats().m_productionRate[ip], expectedRates[ip] ) << GEOS_FMT( "The flux named '{}' did not produce at the expected rate ({}, phase = {}).",
                                                                                           stats.getFluxName(), context, ip );
  }
  EXPECT_DOUBLE_EQ( stats.stats().m_elementCount, expectedElementCount ) << GEOS_FMT( "The flux named '{}' did not produce in the expected elements ({}).",
                                                                                      stats.getFluxName(), context );
}

/**
 * @brief Verification that the source flux statistics are correct over the whole simulation
 * @param problem the simulation ProblemManager
 * @param testSet the simulation TestSet
 */
void checkWholeSimFluxStats( ProblemManager & problem, TestSet const & testSet )
{
  DomainPartition & domain = problem.getDomainPartition();
  SourceFluxStatsAggregator & wholeSimStats =
    problem.getGroupByPath< SourceFluxStatsAggregator >( testSet.inputs.wholeSimFluxStatsPath );

  wholeSimStats.forMeshLevelStatsWrapper( domain,
                                          [&] ( MeshLevel & meshLevel,
                                                SourceFluxStatsAggregator::WrappedStats & meshLevelStats )
  {
    wholeSimStats.forAllFluxStatsWrappers( meshLevel,
                                           [&] ( MeshLevel &,
                                                 SourceFluxStatsAggregator::WrappedStats & fluxStats )
    {
      if( fluxStats.getFluxName() == testSet.inputs.sourceFluxName )
      {
        checkFluxStats( testSet.totalSourceMassProd,
                        testSet.sourceMeanRate,
                        testSet.inputs.sourceElementsCount,
                        fluxStats, "over whole simulation" );
      }
      else if( fluxStats.getFluxName() == testSet.inputs.sinkFluxName )
      {
        checkFluxStats( testSet.totalSinkMassProd,
                        testSet.sinkMeanRate,
                        testSet.inputs.sinkElementsCount,
                        fluxStats, "over whole simulation" );
      }
      else
      {
        FAIL() << "Unexpected SourceFlux found!";
      }
    } );

    for( int ip = 0; ip < meshLevelStats.stats().getPhaseCount(); ++ip )
    {
      EXPECT_DOUBLE_EQ( meshLevelStats.stats().m_producedMass[ip], testSet.totalMassProd[ip] ) << "The fluxes did not produce the expected total mass (over whole simulation, phase = " << ip << ").";
      EXPECT_DOUBLE_EQ( meshLevelStats.stats().m_productionRate[ip], testSet.totalMeanRate[ip] ) << "The fluxes did not produce at the expected rate (over whole simulation, phase = " << ip << ").";
    }
  } );
}

/**
 * @brief Verification that the source flux statistics are correct for a given timestep
 * @param problem the simulation ProblemManager
 * @param testSet the simulation TestSet
 * @param time_n the current timestep start
 * @param timestepId the current timestep id (= cycle)
 */
void checkTimeStepFluxStats( ProblemManager & problem, TestSet const & testSet,
                             real64 const time_n, integer const timestepId )
{
  DomainPartition & domain = problem.getDomainPartition();
  SourceFluxStatsAggregator & timestepStats =
    problem.getGroupByPath< SourceFluxStatsAggregator >( testSet.inputs.timeStepFluxStatsPath );

  timestepStats.forMeshLevelStatsWrapper( domain,
                                          [&] ( MeshLevel & meshLevel,
                                                SourceFluxStatsAggregator::WrappedStats & )
  {
    timestepStats.forAllFluxStatsWrappers( meshLevel,
                                           [&] ( MeshLevel &,
                                                 SourceFluxStatsAggregator::WrappedStats & fluxStats )
    {
      if( fluxStats.getFluxName() == testSet.inputs.sourceFluxName )
      {
        checkFluxStats( testSet.sourceMassProd[timestepId],
                        testSet.sourceRates[timestepId],
                        testSet.inputs.sourceElementsCount,
                        fluxStats, GEOS_FMT( "for timestep at t = {} s", time_n ) );
      }
      else if( fluxStats.getFluxName() == testSet.inputs.sinkFluxName )
      {
        checkFluxStats( testSet.sinkMassProd[timestepId],
                        testSet.sinkRates[timestepId],
                        testSet.inputs.sinkElementsCount,
                        fluxStats, GEOS_FMT( "for timestep at t = {} s", time_n ) );
      }
      else
      {
        FAIL() << "Unexpected SourceFlux found!";
      }
    } );
  } );
}

void checkTimeStepStats( TestSet const & testSet,
                         real64 const time_n,
                         integer const timestepId )
{
  EXPECT_LT( timestepId, testSet.timestepCount ) << GEOS_FMT( "The tested time-step count were higher than expected (t = {} s).",
                                                              time_n );
}

void checkWholeSimTimeStepStats( ProblemManager & problem,
                                 TestSet const & testSet,
                                 TimeStepChecker const & timeStepChecker )
{
  EXPECT_EQ( timeStepChecker.getTestedTimeStepCount(), testSet.timestepCount ) << "The tested time-step were different than expected.";

  PhysicsSolverBase const & solver = problem.getGroupByPath< PhysicsSolverBase >( testSet.inputs.flowSolverPath );
  SolverStatistics const & solverStats = solver.getSolverStatistics();
  EXPECT_GE( solverStats.getNumTimeStepCuts(), testSet.inputs.requiredSubTimeStep ) << "The test did not encountered any timestep cut, but were expected to. "
                                                                                       "Consider adapting the simulation so a timestep cut occurs to check they work as expected.";
}


//////////////////////////////// SinglePhase Flux Statistics Test ////////////////////////////////
namespace SinglePhaseFluxStatisticsTest
{


TestSet getTestSet()
{
  TestInputs testInputs;

  testInputs.xmlInput =
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

  <FieldSpecifications>
    <!-- The rates are positive, but we need negative values for injection (scale = -1) -->
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="-1"
                functionName="FluxRate"
                setNames="{ sourceBox }" />
    <!-- sink is producing 3x source rate -->
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="3"
                functionName="FluxRate"
                setNames="{ sinkBox }"/>

    <HydrostaticEquilibrium name="equil"
                            objectPath="ElementRegions"
                            maxNumberOfEquilibrationIterations="100"
                            datumElevation="-5"
                            datumPressure="1.895e7" />
  </FieldSpecifications>

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

  <!-- We are adding 500s to the whole sim time to force the wholeSimStatsEvent to be executed -->
  <Events maxTime="5500.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />

    <PeriodicEvent name="timestepStatsEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepFluxStats" />
    <PeriodicEvent name="timestepReservoirStatsEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepReservoirStats" />

    <PeriodicEvent name="timestepsCheckEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepChecker" />

    <PeriodicEvent name="wholeSimStatsEvent"
                   timeFrequency="5000.0"
                   target="/Tasks/wholeSimFluxStats" />
  </Events>

  <Tasks>
    <SourceFluxStatistics name="timeStepFluxStats"
                          flowSolverName="testSolver"
                          fluxNames="{*}"
                          logLevel="3" />
    <SourceFluxStatistics name="wholeSimFluxStats"
                          fluxNames="{*}"
                          flowSolverName="testSolver"
                          logLevel="3" />

    <SinglePhaseStatistics name="timeStepReservoirStats"
                           flowSolverName="testSolver"
                           logLevel="1" />

    <TimeStepChecker name="timeStepChecker" />
  </Tasks>

  <Functions>
    <!-- Unscaled injection / production rate in mol/s -->
    <TableFunction
      name="FluxRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.767,  0.894,  0.561,  0.234,  0.194,  0.178,  0.162,  0.059,  0.000,  0.000 }"
    />
  </Functions>

</Problem>
)xml";


  testInputs.sourceFluxName = "sourceFlux";
  testInputs.sinkFluxName = "sinkFlux";
  testInputs.timeStepCheckerPath = "/Tasks/timeStepChecker";
  testInputs.timeStepFluxStatsPath = "/Tasks/timeStepFluxStats";
  testInputs.wholeSimFluxStatsPath = "/Tasks/wholeSimFluxStats";
  testInputs.flowSolverPath = "/Solvers/testSolver";

  testInputs.dt = 500.0;
  testInputs.sourceElementsCount = 2;
  testInputs.sinkElementsCount = 5;

  // FluxRate table from 0.0s to 5000.0s
  setRateTable( testInputs.sourceRates,
                { { 0.000 },
                  { 0.000 },
                  { 0.767 },
                  { 0.894 },
                  { 0.561 },
                  { 0.234 },
                  { 0.194 },
                  { 0.178 },
                  { 0.162 },
                  { 0.059 },
                  { 0.000 } } );
  testInputs.sinkRates=testInputs.sourceRates;

  // sink is 3x source production
  testInputs.sourceRateFactor = -1.0;
  testInputs.sinkRateFactor = 3.0;

  return TestSet( testInputs );
}

TEST_F( FlowStatisticsTest, checkSinglePhaseFluxStatistics )
{
  TestSet const testSet = getTestSet();

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();

  setupProblemFromXML( problem, testSet.inputs.xmlInput.data() );

  real64 firstMass;

  TimeStepChecker & timeStepChecker = problem.getGroupByPath< TimeStepChecker >( testSet.inputs.timeStepCheckerPath );
  timeStepChecker.setTimeStepCheckingFunction( [&]( real64 const time_n )
  {
    integer const timestepId = timeStepChecker.getTestedTimeStepCount();
    checkTimeStepStats( testSet, time_n, timestepId );
    checkTimeStepFluxStats( problem, testSet, time_n, timestepId );

    static bool passedFirstTimeStep = false;
    if( !passedFirstTimeStep )
    {
      passedFirstTimeStep = true;
      firstMass = getTotalFluidMass( problem, testSet.inputs.flowSolverPath );
    }
  } );

  // run simulation
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  checkWholeSimFluxStats( problem, testSet );
  checkWholeSimTimeStepStats( problem, testSet, timeStepChecker );

  // check singlephasestatistics results
  real64 const lastMass = getTotalFluidMass( problem, testSet.inputs.flowSolverPath );
  real64 const massDiffTol = 1e-7;
  EXPECT_NEAR( lastMass - firstMass,
               -testSet.totalMassProd[0],
               massDiffTol * std::abs( testSet.totalMassProd[0] ) ) << GEOS_FMT( "{} total mass difference from start to end is not consistent with fluxes production.",
                                                                                 SinglePhaseStatistics::catalogName() );
}


} /* namespace SinglePhaseFluxStatisticsTest */


//////////////////////////////// Multiphase Flux Statistics Test ////////////////////////////////
namespace MultiPhaseFluxStatisticsTestMass
{


TestSet getTestSet()
{
  TestInputs testInputs;

  testInputs.xmlInput =
    R"xml(
<Problem>

  <Solvers>
    <CompositionalMultiphaseFVM name="testSolver"
                                discretization="fluidTPFA"
                                targetRegions="{ reservoir }"
                                temperature="366.483"
                                useMass="1"
                                logLevel="1" >
      <NonlinearSolverParameters newtonMaxIter="8"
                                 maxTimeStepCuts="8"
                                 allowNonConverged="1" />
      <LinearSolverParameters solverType="gmres"
                              preconditionerType="iluk"
                              krylovTol="1.0e-6" />
    </CompositionalMultiphaseFVM>
  </Solvers>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="fluidTPFA" />
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
                       materialList="{ fluid, rock, relperm }" />
  </ElementRegions>

  <Constitutive>
    <CO2BrineEzrokhiFluid name="fluid"
                          phaseNames="{ gas, water }"
                          componentNames="{ co2, water }"
                          componentMolarWeight="{ 44e-3, 18e-3 }"
                          phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
                          flashModelParaFile="co2flash.txt" />

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
    <TableRelativePermeability name="relperm"
                               phaseNames="{ gas, water }"
                               wettingNonWettingRelPermTableNames="{ gasRelativePermeabilityTable, waterRelativePermeabilityTable }" />
  </Constitutive>

  <FieldSpecifications>
    <!-- We are injecting CO2 (negative production values), scaling to convert from mol/s to kg/s -->
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="0"
                scale="-44e-3"
                functionName="FluxInjectionRate"
                setNames="{ sourceBox }" />
    <!-- We are depleting water, scaling to convert from mol/s to kg/s -->
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="18e-3"
                functionName="FluxProductionRate"
                setNames="{ sinkBox }" />

    <HydrostaticEquilibrium name="equil"
                            objectPath="ElementRegions"
                            maxNumberOfEquilibrationIterations="100"
                            datumElevation="-5"
                            datumPressure="1.895e7"
                            initialPhaseName="water"
                            componentNames="{ co2, water }"
                            componentFractionVsElevationTableNames="{ initGasCompFracTable, initWaterCompFracTable }"
                            temperatureVsElevationTableName="initTempTable" />
  </FieldSpecifications>

  <Geometry>
    <!-- source selects 2 elements -->
    <Box name="sourceBox"
         xMin="{ -0.01, -0.01, -10.01 }"
         xMax="{  1.01,  1.01,  -8.99 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{  8.99,  8.99, -1.01 }"
         xMax="{ 10.01, 10.01,  0.01 }" />
  </Geometry>

  <!-- We are adding 500s to the whole sim time to force the wholeSimStatsEvent to be executed -->
  <Events maxTime="5500.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />

    <PeriodicEvent name="timestepFluxStatsEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepFluxStats" />
    <PeriodicEvent name="timestepReservoirStatsEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepReservoirStats" />

    <PeriodicEvent name="timestepsCheckEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepChecker" />

    <PeriodicEvent name="wholeSimStatsEvent"
                   timeFrequency="5000.0"
                   target="/Tasks/wholeSimFluxStats" />
  </Events>

  <Tasks>
    <SourceFluxStatistics name="timeStepFluxStats"
                          flowSolverName="testSolver"
                          fluxNames="{*}"
                          logLevel="2" />
    <SourceFluxStatistics name="wholeSimFluxStats"
                          flowSolverName="testSolver"
                          fluxNames="{*}"
                          logLevel="2" />

    <CompositionalMultiphaseStatistics name="timeStepReservoirStats"
                                       flowSolverName="testSolver"
                                       logLevel="1" />

    <TimeStepChecker name="timeStepChecker" />
  </Tasks>

  <Functions>
    <!-- CO2 injection rates in mol/s -->
    <TableFunction
      name="FluxInjectionRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.267,  0.561,  0.194,  0.102,  0.059,  0.000,  0.000,  0.000,  0.000,  0.000 }"
    />
    <!-- water depletion rates in mol/s -->
    <TableFunction
      name="FluxProductionRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.003,  0.062,  0.121,  0.427,  0.502,  0.199,  0.083,  0.027,  0.000,  0.000 }"
    />

    <TableFunction name="initGasCompFracTable"
                   coordinates="{ -10.0, -7.0, -3.0, -1.0 }"
                   values="{      0.001,  0.001,  0.001,  0.001 }" />
    <TableFunction name="initWaterCompFracTable"
                   coordinates="{ -10.0, -7.0, -3.0, -1.0 }"
                   values="{      0.999,  0.999,  0.999,  0.999 }" />
    <TableFunction name="initTempTable"
                   coordinates="{ -10.0,   -7.0,   -3.0,   -1.0 }"
                   values="{     395.15, 389.15, 382.15, 378.15 }" />
    <TableFunction name="waterRelativePermeabilityTable"
                   coordinates="{ 0.3000,          0.3175,         0.3350,         0.3525,        0.3700,        0.3875,       0.4050,       0.4225,       0.4400,       0.4575,       0.4750,      0.4925,      0.5100,      0.5275,      0.5450,      0.5625,      0.5800,      0.5975,     0.6150,     0.6325,     0.6500,     0.6675,     0.6850,     0.7025,     0.7200,     0.7375,     0.7550,     0.7725,    0.7900,    0.8054,    0.8209,    0.8404,    0.8600,    0.8775,    0.8950,    0.9125,    0.9300,    0.9475,    0.9650,    0.9825, 1.0000   }"
                   values="{         0.0, 0.0000001069690, 0.000001523818, 0.000007304599, 0.00002242961, 0.00005398050, 0.0001113999, 0.0002068239, 0.0003554932, 0.0005762517, 0.0008921512, 0.001331180, 0.001927144, 0.002720726, 0.003760776, 0.005105868, 0.006826186, 0.009005830, 0.01174561, 0.01516648, 0.01941368, 0.02466185, 0.03112128, 0.03904542, 0.04874017, 0.06057494, 0.07499593, 0.09254174, 0.1138611, 0.1364565, 0.1632363, 0.2042135, 0.2547712, 0.3097943, 0.3755964, 0.4536528, 0.5451093, 0.6502388, 0.7674166, 0.8909226, 1.000000 }" />
    <TableFunction name="gasRelativePermeabilityTable"
                   coordinates="{     0.0000,       0.0175,      0.0350,      0.0525,      0.0700,     0.0875,     0.1050,     0.1225,     0.1400,     0.1595,     0.1790,     0.1945,     0.2100,     0.2275,     0.2450,     0.2625,     0.2800,     0.2975,     0.3150,    0.3325,    0.3500,    0.3675,    0.3850,    0.4025,    0.4200,    0.4375,    0.4550,    0.4725,    0.4900,    0.5075,    0.5250,    0.5425,    0.5600,    0.5775,    0.5950,    0.6125,    0.6300,    0.6475,    0.6650,    0.6825,    0.7000 }"
                   values="{     0.000000000, 0.0008885248, 0.002483741, 0.004583224, 0.007135315, 0.01012132, 0.01353719, 0.01738728, 0.02168159, 0.02701850, 0.03295183, 0.03808925, 0.04363513, 0.05042783, 0.05779578, 0.06577020, 0.07438478, 0.08367565, 0.09368138, 0.1044429, 0.1160032, 0.1284076, 0.1417029, 0.1559376, 0.1711607, 0.1874214, 0.2047679, 0.2232459, 0.2428968, 0.2637550, 0.2858446, 0.3091747, 0.3337331, 0.3594782, 0.3863263, 0.4141347, 0.4426735, 0.4715782, 0.5002513, 0.5275887, 0.5500000 }" />
  </Functions>

</Problem>
)xml";

  testInputs.tableFiles["pvtgas.txt"] = "DensityFun SpanWagnerCO2Density 1.5e7 2.5e7 1e5 370.15 400.15 2\n"
                                        "ViscosityFun FenghourCO2Viscosity 1.5e7 2.5e7 1e5 370.15 400.15 2\n";

  testInputs.tableFiles["pvtliquid.txt"] = "DensityFun EzrokhiBrineDensity 0.1033 -2.2991e-5 -2.3658e-6\n"
                                           "ViscosityFun EzrokhiBrineViscosity 0 0 0\n";

  testInputs.tableFiles["co2flash.txt"] = "FlashModel CO2Solubility 1.5e7 2.5e7 1e5 370.15 400.15 2 0\n";


  testInputs.sourceFluxName = "sourceFlux";
  testInputs.sinkFluxName = "sinkFlux";
  testInputs.timeStepCheckerPath = "/Tasks/timeStepChecker";
  testInputs.timeStepFluxStatsPath = "/Tasks/timeStepFluxStats";
  testInputs.wholeSimFluxStatsPath = "/Tasks/wholeSimFluxStats";
  testInputs.flowSolverPath = "/Solvers/testSolver";

  testInputs.dt = 500.0;
  testInputs.sourceElementsCount = 1;
  testInputs.sinkElementsCount = 1;

  // FluxInjectionRate & FluxProductionRate table from 0.0s to 5000.0s
  setRateTable( testInputs.sourceRates,
                { { 0.000, 0.0 },
                  { 0.000, 0.0 },
                  { 0.267, 0.0 },
                  { 0.561, 0.0 },
                  { 0.194, 0.0 },
                  { 0.102, 0.0 },
                  { 0.059, 0.0 },
                  { 0.000, 0.0 },
                  { 0.000, 0.0 },
                  { 0.000, 0.0 },
                  { 0.000, 0.0 } } );
  setRateTable( testInputs.sinkRates,
                { { 0.0, 0.000 },
                  { 0.0, 0.000 },
                  { 0.0, 0.003 },
                  { 0.0, 0.062 },
                  { 0.0, 0.121 },
                  { 0.0, 0.427 },
                  { 0.0, 0.502 },
                  { 0.0, 0.199 },
                  { 0.0, 0.083 },
                  { 0.0, 0.027 },
                  { 0.0, 0.000 } } );

  testInputs.sourceRateFactor = -44e-3;
  testInputs.sinkRateFactor = 18e-3;

  return TestSet( testInputs );
}


TEST_F( FlowStatisticsTest, checkMultiPhaseFluxStatisticsMass )
{
  TestSet const testSet = getTestSet();
  writeTableFiles( testSet.inputs.tableFiles );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();

  setupProblemFromXML( problem, testSet.inputs.xmlInput.data() );

  TimeStepChecker & timeStepChecker = problem.getGroupByPath< TimeStepChecker >( testSet.inputs.timeStepCheckerPath );
  timeStepChecker.setTimeStepCheckingFunction( [&]( real64 const time_n )
  {
    integer const timestepId = timeStepChecker.getTestedTimeStepCount();
    checkTimeStepStats( testSet, time_n, timestepId );
    checkTimeStepFluxStats( problem, testSet, time_n, timestepId );
  } );

  // run simulation
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  checkWholeSimFluxStats( problem, testSet );
  checkWholeSimTimeStepStats( problem, testSet, timeStepChecker );
}


}   /* namespace MultiPhaseFluxStatisticsTest */


//////////////////////////////// Multiphase Flux Statistics Test ////////////////////////////////
namespace MultiPhaseFluxStatisticsTestMol
{


TestSet getTestSet()
{
  TestInputs testInputs;

  testInputs.xmlInput =
    R"xml(
<Problem>

  <Solvers>
    <CompositionalMultiphaseFVM name="testSolver"
                                discretization="fluidTPFA"
                                targetRegions="{ reservoir }"
                                temperature="366.483"
                                useMass="0"
                                logLevel="1" >
      <NonlinearSolverParameters newtonMaxIter="8"
                                 maxTimeStepCuts="8"
                                 allowNonConverged="1" />
      <LinearSolverParameters solverType="gmres"
                              preconditionerType="iluk"
                              krylovTol="1.0e-6" />
    </CompositionalMultiphaseFVM>
  </Solvers>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="fluidTPFA" />
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
                       materialList="{ fluid, rock, relperm }" />
  </ElementRegions>

  <Constitutive>
    <CO2BrineEzrokhiFluid name="fluid"
                          phaseNames="{ gas, water }"
                          componentNames="{ co2, water }"
                          componentMolarWeight="{ 44e-3, 18e-3 }"
                          phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
                          flashModelParaFile="co2flash.txt" />

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
    <TableRelativePermeability name="relperm"
                               phaseNames="{ gas, water }"
                               wettingNonWettingRelPermTableNames="{ gasRelativePermeabilityTable, waterRelativePermeabilityTable }" />
  </Constitutive>

  <FieldSpecifications>
    <!-- We are injecting CO2 (negative production values) -->
    <SourceFlux name="sourceFlux"
                objectPath="ElementRegions/reservoir"
                component="0"
                scale="-8"
                functionName="FluxInjectionRate"
                setNames="{ sourceBox }" />
    <!-- We are depleting water -->
    <SourceFlux name="sinkFlux"
                objectPath="ElementRegions/reservoir"
                component="1"
                scale="8"
                functionName="FluxProductionRate"
                setNames="{ sinkBox }" />

    <HydrostaticEquilibrium name="equil"
                            objectPath="ElementRegions"
                            maxNumberOfEquilibrationIterations="100"
                            datumElevation="-5"
                            datumPressure="1.895e7"
                            initialPhaseName="water"
                            componentNames="{ co2, water }"
                            componentFractionVsElevationTableNames="{ initGasCompFracTable, initWaterCompFracTable }"
                            temperatureVsElevationTableName="initTempTable" />
  </FieldSpecifications>

  <Geometry>
    <!-- source selects 2 elements -->
    <Box name="sourceBox"
         xMin="{ -0.01, -0.01, -10.01 }"
         xMax="{  1.01,  1.01,  -8.99 }" />
    <!-- sink selects 2 elements -->
    <Box name="sinkBox"
         xMin="{  8.99,  8.99, -1.01 }"
         xMax="{ 10.01, 10.01,  0.01 }" />
  </Geometry>

  <!-- We are adding 500s to the whole sim time to force the wholeSimStatsEvent to be executed -->
  <Events maxTime="5500.0">
    <PeriodicEvent name="solverApplications"
                   forceDt="500.0"
                   target="/Solvers/testSolver" />

    <PeriodicEvent name="timestepFluxStatsEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepFluxStats" />
    <PeriodicEvent name="timestepReservoirStatsEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepReservoirStats" />

    <PeriodicEvent name="timestepsCheckEvent"
                   timeFrequency="500.0"
                   target="/Tasks/timeStepChecker" />

    <PeriodicEvent name="wholeSimStatsEvent"
                   timeFrequency="5000.0"
                   target="/Tasks/wholeSimFluxStats" />
  </Events>

  <Tasks>
    <SourceFluxStatistics name="timeStepFluxStats"
                          flowSolverName="testSolver"
                          fluxNames="{*}"
                          logLevel="1" />
    <SourceFluxStatistics name="wholeSimFluxStats"
                          flowSolverName="testSolver"
                          fluxNames="{*}"
                          logLevel="1" />

    <CompositionalMultiphaseStatistics name="timeStepReservoirStats"
                                       flowSolverName="testSolver"
                                       logLevel="1" />

    <TimeStepChecker name="timeStepChecker" />
  </Tasks>

  <Functions>
    <!-- CO2 injection rates in mol/s -->
    <TableFunction
      name="FluxInjectionRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.267,  0.561,  0.194,  0.102,  0.059,  0.000,  0.000,  0.000,  0.000,  0.000 }"
    />
    <!-- water depletion rates in mol/s -->
    <TableFunction
      name="FluxProductionRate"
      inputVarNames="{ time }"
      interpolation="lower"
      coordinates="{    0.0,  500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0 }"
      values="{       0.000,  0.000,  0.003,  0.062,  0.121,  0.427,  0.502,  0.199,  0.083,  0.027,  0.000,  0.000 }"
    />

    <TableFunction name="initGasCompFracTable"
                   coordinates="{ -10.0, -7.0, -3.0, -1.0 }"
                   values="{      0.001,  0.001,  0.001,  0.001 }" />
    <TableFunction name="initWaterCompFracTable"
                   coordinates="{ -10.0, -7.0, -3.0, -1.0 }"
                   values="{      0.999,  0.999,  0.999,  0.999 }" />
    <TableFunction name="initTempTable"
                   coordinates="{ -10.0,   -7.0,   -3.0,   -1.0 }"
                   values="{     395.15, 389.15, 382.15, 378.15 }" />
    <TableFunction name="waterRelativePermeabilityTable"
                   coordinates="{ 0.3000,          0.3175,         0.3350,         0.3525,        0.3700,        0.3875,       0.4050,       0.4225,       0.4400,       0.4575,       0.4750,      0.4925,      0.5100,      0.5275,      0.5450,      0.5625,      0.5800,      0.5975,     0.6150,     0.6325,     0.6500,     0.6675,     0.6850,     0.7025,     0.7200,     0.7375,     0.7550,     0.7725,    0.7900,    0.8054,    0.8209,    0.8404,    0.8600,    0.8775,    0.8950,    0.9125,    0.9300,    0.9475,    0.9650,    0.9825, 1.0000   }"
                   values="{         0.0, 0.0000001069690, 0.000001523818, 0.000007304599, 0.00002242961, 0.00005398050, 0.0001113999, 0.0002068239, 0.0003554932, 0.0005762517, 0.0008921512, 0.001331180, 0.001927144, 0.002720726, 0.003760776, 0.005105868, 0.006826186, 0.009005830, 0.01174561, 0.01516648, 0.01941368, 0.02466185, 0.03112128, 0.03904542, 0.04874017, 0.06057494, 0.07499593, 0.09254174, 0.1138611, 0.1364565, 0.1632363, 0.2042135, 0.2547712, 0.3097943, 0.3755964, 0.4536528, 0.5451093, 0.6502388, 0.7674166, 0.8909226, 1.000000 }" />
    <TableFunction name="gasRelativePermeabilityTable"
                   coordinates="{     0.0000,       0.0175,      0.0350,      0.0525,      0.0700,     0.0875,     0.1050,     0.1225,     0.1400,     0.1595,     0.1790,     0.1945,     0.2100,     0.2275,     0.2450,     0.2625,     0.2800,     0.2975,     0.3150,    0.3325,    0.3500,    0.3675,    0.3850,    0.4025,    0.4200,    0.4375,    0.4550,    0.4725,    0.4900,    0.5075,    0.5250,    0.5425,    0.5600,    0.5775,    0.5950,    0.6125,    0.6300,    0.6475,    0.6650,    0.6825,    0.7000 }"
                   values="{     0.000000000, 0.0008885248, 0.002483741, 0.004583224, 0.007135315, 0.01012132, 0.01353719, 0.01738728, 0.02168159, 0.02701850, 0.03295183, 0.03808925, 0.04363513, 0.05042783, 0.05779578, 0.06577020, 0.07438478, 0.08367565, 0.09368138, 0.1044429, 0.1160032, 0.1284076, 0.1417029, 0.1559376, 0.1711607, 0.1874214, 0.2047679, 0.2232459, 0.2428968, 0.2637550, 0.2858446, 0.3091747, 0.3337331, 0.3594782, 0.3863263, 0.4141347, 0.4426735, 0.4715782, 0.5002513, 0.5275887, 0.5500000 }" />
  </Functions>

</Problem>
)xml";

  testInputs.tableFiles["pvtgas.txt"] = "DensityFun SpanWagnerCO2Density 1.5e7 2.5e7 1e5 370.15 400.15 2\n"
                                        "ViscosityFun FenghourCO2Viscosity 1.5e7 2.5e7 1e5 370.15 400.15 2\n";

  testInputs.tableFiles["pvtliquid.txt"] = "DensityFun EzrokhiBrineDensity 0.1033 -2.2991e-5 -2.3658e-6\n"
                                           "ViscosityFun EzrokhiBrineViscosity 0 0 0\n";

  testInputs.tableFiles["co2flash.txt"] = "FlashModel CO2Solubility 1.5e7 2.5e7 1e5 370.15 400.15 2 0\n";


  testInputs.sourceFluxName = "sourceFlux";
  testInputs.sinkFluxName = "sinkFlux";
  testInputs.timeStepCheckerPath = "/Tasks/timeStepChecker";
  testInputs.timeStepFluxStatsPath = "/Tasks/timeStepFluxStats";
  testInputs.wholeSimFluxStatsPath = "/Tasks/wholeSimFluxStats";
  testInputs.flowSolverPath = "/Solvers/testSolver";

  testInputs.dt = 500.0;
  testInputs.sourceElementsCount = 1;
  testInputs.sinkElementsCount = 1;

  // FluxInjectionRate & FluxProductionRate table from 0.0s to 5000.0s
  setRateTable( testInputs.sourceRates,
                { { 0.000, 0.0 },
                  { 0.000, 0.0 },
                  { 0.267, 0.0 },
                  { 0.561, 0.0 },
                  { 0.194, 0.0 },
                  { 0.102, 0.0 },
                  { 0.059, 0.0 },
                  { 0.000, 0.0 },
                  { 0.000, 0.0 },
                  { 0.000, 0.0 },
                  { 0.000, 0.0 } } );
  setRateTable( testInputs.sinkRates,
                { { 0.0, 0.000 },
                  { 0.0, 0.000 },
                  { 0.0, 0.003 },
                  { 0.0, 0.062 },
                  { 0.0, 0.121 },
                  { 0.0, 0.427 },
                  { 0.0, 0.502 },
                  { 0.0, 0.199 },
                  { 0.0, 0.083 },
                  { 0.0, 0.027 },
                  { 0.0, 0.000 } } );

  // scale is set to high values to make the solver generate timestep cuts
  testInputs.sourceRateFactor = -8.0;
  testInputs.sinkRateFactor = 8.0;

  // this simulation is set-up to have at least one timestep cut.
  testInputs.requiredSubTimeStep = 2;

  return TestSet( testInputs );
}


TEST_F( FlowStatisticsTest, checkMultiPhaseFluxStatisticsMol )
{
  TestSet const testSet = getTestSet();
  writeTableFiles( testSet.inputs.tableFiles );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();

  setupProblemFromXML( problem, testSet.inputs.xmlInput.data() );

  TimeStepChecker & timeStepChecker = problem.getGroupByPath< TimeStepChecker >( testSet.inputs.timeStepCheckerPath );
  timeStepChecker.setTimeStepCheckingFunction( [&]( real64 const time_n )
  {
    integer const timestepId = timeStepChecker.getTestedTimeStepCount();
    checkTimeStepStats( testSet, time_n, timestepId );
    checkTimeStepFluxStats( problem, testSet, time_n, timestepId );
  } );

  // run simulation
  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  checkWholeSimFluxStats( problem, testSet );
  checkWholeSimTimeStepStats( problem, testSet, timeStepChecker );
}


}   /* namespace MultiPhaseFluxStatisticsTest */


//////////////////////////////// Main ////////////////////////////////


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
