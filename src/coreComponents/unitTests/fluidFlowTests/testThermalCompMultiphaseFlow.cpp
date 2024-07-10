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

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::constitutive::multifluid;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * pvtLiquid = "DensityFun PhillipsBrineDensity 1e6 7.5e7 5e5 295.15 370.15 25 0\n"
                         "ViscosityFun PhillipsBrineViscosity 0\n"
                         "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 295.15 370.15 25 0\n";

char const * pvtGas = "DensityFun SpanWagnerCO2Density 1e6 7.5e7 5e5 295.15 370.15 25\n"
                      "ViscosityFun FenghourCO2Viscosity 1e6 7.5e7 5e5 295.15 370.15 25\n"
                      "EnthalpyFun CO2Enthalpy 1e6 7.5e7 5e5 295.15 370.15 25\n";

char const * co2flash = "FlashModel CO2Solubility  1e6 7.5e7 5e5 295.15 370.15 25 0";

char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers>
      <CompositionalMultiphaseFVM name="compflow"
                                  logLevel="1"
                                  discretization="fluidTPFA"
                                  temperature="368.15"
                                  useMass="1"
                                  isThermal="1"
                                  initialDt="1000"
                                  maxCompFractionChange="0.5"
                                  targetRegions="{ region }">
        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="100"
                                   lineSearchAction="None"
                                   maxTimeStepCuts="5" />
        <LinearSolverParameters directParallel="0" />
      </CompositionalMultiphaseFVM>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{ C3D8 }"
                    xCoords="{ 0, 20 }"
                    yCoords="{ 0, 1 }"
                    zCoords="{ 0, 1 }"
                    nx="{ 5 }"
                    ny="{ 1 }"
                    nz="{ 1 }"
                    cellBlockNames="{ cb }" />
    </Mesh>
    <Geometry>
      <Box name="sink"
           xMin="{ -0.01, -0.01, -0.01 }"
           xMax="{ 4.01, 1.01, 1.01 }" />
      <Box name="source"
           xMin="{ -0.01, -0.01, -0.01 }"
           xMax="{ 4.01, 1.01, 1.01 }" />
    </Geometry>
    <Events maxTime="1000">
      <PeriodicEvent name="solverApplications"
                     maxEventDt="1000"
                     target="/Solvers/compflow" />
    </Events>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA" />
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region"
                         cellBlocks="{ cb }"
                         materialList="{ fluid, rock, relperm, thermalCond }" />
    </ElementRegions>
    <Constitutive>
      <CompressibleSolidConstantPermeability name="rock"
                                             solidModelName="nullSolid"
                                             porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm"
                                             solidInternalEnergyModelName="rockInternalEnergy" />
      <NullModel name="nullSolid" />
      <PressurePorosity name="rockPorosity"
                        defaultReferencePorosity="0.2"
                        referencePressure="0.0"
                        compressibility="1.0e-9" />
      <SolidInternalEnergy name="rockInternalEnergy"
                           referenceVolumetricHeatCapacity="1.95e6"
                           referenceTemperature="368.15"
                           referenceInternalEnergy="0" />
      <ConstantPermeability name="rockPerm"
                            permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }" />
      <CO2BrinePhillipsThermalFluid name="fluid"
                                    phaseNames="{ gas, water }"
                                    componentNames="{ co2, water }"
                                    componentMolarWeight="{ 44e-3, 18e-3 }"
                                    phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
                                    flashModelParaFile="co2flash.txt" />
      <BrooksCoreyRelativePermeability name="relperm"
                                       phaseNames="{ gas, water }"
                                       phaseMinVolumeFraction="{ 0.0, 0.0 }"
                                       phaseRelPermExponent="{ 1.5, 1.5 }"
                                       phaseRelPermMaxValue="{ 0.9, 0.9 }" />
      <MultiPhaseConstantThermalConductivity name="thermalCond"
                                             phaseNames="{ gas, water }"
                                             thermalConductivityComponents="{ 0.6, 0.6, 0.6 }" />
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="9e6" />
      <FieldSpecification name="initialTemperature"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="temperature"
                          scale="368.15" />
      <FieldSpecification name="initialComposition_co2"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="0"
                          scale="0.1" />
      <FieldSpecification name="initialComposition_water"
                          initialCondition="1"
                          setNames="{ all }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="1"
                          scale="0.9" />
      <FieldSpecification name="sinkPressure"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="7e6" />
      <FieldSpecification name="sinkTemperature"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="temperature"
                          scale="368.15" />
      <FieldSpecification name="sinkTermComposition_co2"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="0"
                          scale="0.1" />
      <FieldSpecification name="sinkTermComposition_water"
                          setNames="{ sink }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="1"
                          scale="0.9" />
      <FieldSpecification name="sourcePressure"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="pressure"
                          scale="1.45e7" />
      <FieldSpecification name="sourceTemperature"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="temperature"
                          scale="300.15" />
      <FieldSpecification name="sourceTermComposition_co2"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="0"
                          scale="0.9" />
      <FieldSpecification name="sourceTermComposition_water"
                          setNames="{ source }"
                          objectPath="ElementRegions/region/cb"
                          fieldName="globalCompFraction"
                          component="1"
                          scale="0.1" />
    </FieldSpecifications>
  </Problem>
  )xml";

// Sphinx end before input XML

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseFVM & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );

  // assemble the analytical residual
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();

  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
  residual.move( hostMemorySpace, false );

  // copy the analytical residual
  array1d< real64 > residualOrig( residual );

  // create the numerical jacobian
  jacobian.move( hostMemorySpace );
  CRSMatrix< real64, globalIndex > jacobianFD( jacobian );
  jacobianFD.zero();

  // fill jacobian FD
  fillCellCenteredNumericalJacobian( solver,
                                     domain,
                                     true,
                                     perturbParameter,
                                     residual.toView(),
                                     residualOrig.toView(),
                                     jacobian.toView(),
                                     jacobianFD.toView(),
                                     assembleFunction );

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol, 1e-6 );
}

void writeTableToFile( string const & filename, char const * str )
{
  std::ofstream os( filename );
  ASSERT_TRUE( os.is_open() );
  os << str;
  os.close();
}

void removeFile( string const & filename )
{
  int const ret = std::remove( filename.c_str() );
  ASSERT_TRUE( ret == 0 );
}

class ThermalCompositionalMultiphaseFlowTest : public ::testing::Test
{
public:

  ThermalCompositionalMultiphaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {
    writeTableToFile( pvtLiquidFilename, pvtLiquid );
    writeTableToFile( pvtGasFilename, pvtGas );
    writeTableToFile( co2flashFilename, co2flash );
  }

  ~ThermalCompositionalMultiphaseFlowTest() override
  {
    removeFile( pvtLiquidFilename );
    removeFile( pvtGasFilename );
    removeFile( co2flashFilename );
  }

protected:

  void SetUp() override
  {

    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseFVM >( "compflow" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution() );

    solver->implicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  CompositionalMultiphaseFVM * solver;

  string const pvtLiquidFilename = "pvtliquid.txt";
  string const pvtGasFilename = "pvtgas.txt";
  string const co2flashFilename = "co2flash.txt";
};

real64 constexpr ThermalCompositionalMultiphaseFlowTest::time;
real64 constexpr ThermalCompositionalMultiphaseFlowTest::dt;
real64 constexpr ThermalCompositionalMultiphaseFlowTest::eps;

TEST_F( ThermalCompositionalMultiphaseFlowTest, derivativeNumericalCheck_composition )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-4;

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testCompositionNumericalDerivatives( *solver, domain, perturb, tol );
}

TEST_F( ThermalCompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseVolumeFraction )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  testPhaseVolumeFractionNumericalDerivatives( *solver, domain, true, perturb, tol );
}

TEST_F( ThermalCompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseMobility )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testPhaseMobilityNumericalDerivatives( *solver, domain, true, perturb, tol );
}

TEST_F( ThermalCompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#if 0
TEST_F( ThermalCompositionalMultiphaseFlowTest, jacobianNumericalCheck_accumulationVolumeBalance )
{
  real64 const perturb = sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleAccumulationAndVolumeBalanceTerms( domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}
#endif

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
