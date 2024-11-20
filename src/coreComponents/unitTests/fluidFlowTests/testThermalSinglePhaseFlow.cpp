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

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "unitTests/fluidFlowTests/testSingleFlowUtils.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML

char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers>
      <SinglePhaseFVM name="singleflow"
                      logLevel="1"
                      discretization="fluidTPFA"
                      temperature="368.15"
                      isThermal="1"
                      targetRegions="{ region }">
        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="100" />
        <LinearSolverParameters solverType="gmres"
                                krylovTol="1.0e-10" />
      </SinglePhaseFVM>
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
                     target="/Solvers/singleflow" />
    </Events>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA" />
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region"
                         cellBlocks="{ * }"
                         materialList="{ water, rock, thermalCond }" />
    </ElementRegions>
    <Constitutive>
      <CompressibleSolidConstantPermeability name="rock"
                                             solidModelName="nullSolid"
                                             porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm"
                                             solidInternalEnergyModelName="rockInternalEnergy" />
      <NullModel name="nullSolid" />
      <PressurePorosity name="rockPorosity"
                        defaultReferencePorosity="0.05"
                        referencePressure="0.0"
                        compressibility="1.0e-9" />
      <SolidInternalEnergy name="rockInternalEnergy"
                           referenceVolumetricHeatCapacity="1.95e6"
                           referenceTemperature="368.15"
                           referenceInternalEnergy="0" />
      <ConstantPermeability name="rockPerm"
                            permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }" />
      <ThermalCompressibleSinglePhaseFluid name="water"
                                           defaultDensity="1000"
                                           defaultViscosity="0.001"
                                           referencePressure="0.0"
                                           referenceTemperature="0.0"
                                           compressibility="5e-10"
                                           thermalExpansionCoeff="7e-4"
                                           viscosibility="0.0"
                                           specificHeatCapacity="4.5e3" />
      <SinglePhaseThermalConductivity name="thermalCond"
                                      defaultThermalConductivityComponents="{ 0.6, 0.6, 0.6 }" />
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
    </FieldSpecifications>
  </Problem>
  )xml";
// Sphinx end before input XML

template< typename LAMBDA >
void testNumericalJacobian( SinglePhaseFVM<> & solver,
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

class ThermalSinglePhaseFlowTest : public ::testing::Test
{
public:

  ThermalSinglePhaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {

    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseFVM<> >( "singleflow" );

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
  SinglePhaseFVM<> * solver;
};

real64 constexpr ThermalSinglePhaseFlowTest::time;
real64 constexpr ThermalSinglePhaseFlowTest::dt;
real64 constexpr ThermalSinglePhaseFlowTest::eps;

TEST_F( ThermalSinglePhaseFlowTest, derivativeNumericalCheck_mobility )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testMobilityNumericalDerivatives( *solver, domain, true, perturb, tol );
}

TEST_F( ThermalSinglePhaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 2e-2; // 2% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    // The first input parameter denotes t_n, which is unused. Just input something here.
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#if 1
TEST_F( ThermalSinglePhaseFlowTest, jacobianNumericalCheck_accumulationBalance )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 2e-2; // 2% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleAccumulationTerms( domain, solver->getDofManager(), localMatrix, localRhs );
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
