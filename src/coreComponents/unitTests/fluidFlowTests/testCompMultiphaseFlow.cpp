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

char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <CompositionalMultiphaseFVM name="compflow"
                                   logLevel="0"
                                   discretization="fluidTPFA"
                                   targetRegions="{region}"
                                   temperature="297.15"
                                   useMass="1">

        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="2"/>
        <LinearSolverParameters solverType="gmres"
                                krylovTol="1.0e-10"/>
      </CompositionalMultiphaseFVM>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{C3D8}"
                    xCoords="{0, 3}"
                    yCoords="{0, 1}"
                    zCoords="{0, 1}"
                    nx="{3}"
                    ny="{1}"
                    nz="{1}"
                    cellBlockNames="{cb1}"/>
    </Mesh>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA"/>
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="region" cellBlocks="{cb1}" materialList="{fluid, rock, relperm, cappressure}" />
    </ElementRegions>
    <Constitutive>
      <CompositionalMultiphaseFluid name="fluid"
                                    phaseNames="{oil, gas}"
                                    equationsOfState="{PR, PR}"
                                    componentNames="{N2, C10, C20, H2O}"
                                    componentCriticalPressure="{34e5, 25.3e5, 14.6e5, 220.5e5}"
                                    componentCriticalTemperature="{126.2, 622.0, 782.0, 647.0}"
                                    componentAcentricFactor="{0.04, 0.443, 0.816, 0.344}"
                                    componentMolarWeight="{28e-3, 134e-3, 275e-3, 18e-3}"
                                    componentVolumeShift="{0, 0, 0, 0}"
                                    componentBinaryCoeff="{ {0, 0, 0, 0},
                                                            {0, 0, 0, 0},
                                                            {0, 0, 0, 0},
                                                            {0, 0, 0, 0} }"/>
      <CompressibleSolidConstantPermeability name="rock"
          solidModelName="nullSolid"
          porosityModelName="rockPorosity"
          permeabilityModelName="rockPerm"/>
     <NullModel name="nullSolid"/>
     <PressurePorosity name="rockPorosity"
                       defaultReferencePorosity="0.05"
                       referencePressure = "0.0"
                       compressibility="1.0e-9"/>
      <BrooksCoreyRelativePermeability name="relperm"
                                       phaseNames="{oil, gas}"
                                       phaseMinVolumeFraction="{0.1, 0.15}"
                                       phaseRelPermExponent="{2.0, 2.0}"
                                       phaseRelPermMaxValue="{0.8, 0.9}"/>
      <BrooksCoreyCapillaryPressure name="cappressure"
                                    phaseNames="{oil, gas}"
                                    phaseMinVolumeFraction="{0.2, 0.05}"
                                    phaseCapPressureExponentInv="{4.25, 3.5}"
                                    phaseEntryPressure="{0., 1e8}"
                                    capPressureEpsilon="0.0"/>
    <ConstantPermeability name="rockPerm"
                          permeabilityComponents="{2.0e-16, 2.0e-16, 2.0e-16}"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="pressure"
                 functionName="initialPressureFunc"
                 scale="5e6"/>
      <FieldSpecification name="initialComposition_N2"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="0"
                 scale="0.099"/>
      <FieldSpecification name="initialComposition_C10"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="1"
                 scale="0.3"/>
      <FieldSpecification name="initialComposition_C20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="2"
                 scale="0.6"/>
      <FieldSpecification name="initialComposition_H20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/region/cb1"
                 fieldName="globalCompFraction"
                 component="3"
                 scale="0.001"/>
    </FieldSpecifications>
    <Functions>
      <TableFunction name="initialPressureFunc"
                     inputVarNames="{elementCenter}"
                     coordinates="{0.0, 3.0}"
                     values="{1.0, 0.5}"/>
    </Functions>
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
                                     false,
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

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class CompositionalMultiphaseFlowTest : public ::testing::Test
{
public:

  CompositionalMultiphaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

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
};

real64 constexpr CompositionalMultiphaseFlowTest::time;
real64 constexpr CompositionalMultiphaseFlowTest::dt;
real64 constexpr CompositionalMultiphaseFlowTest::eps;

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_composition )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-4;

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testCompositionNumericalDerivatives( *solver, domain, perturb, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseVolumeFraction )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  testPhaseVolumeFractionNumericalDerivatives( *solver, domain, false, perturb, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseMobility )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testPhaseMobilityNumericalDerivatives( *solver, domain, false, perturb, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux )
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

/*
 * Accumulation numerical test not passing due to some numerical catastrophic cancellation
 * happenning in the kernel for the particular set of initial conditions we're running.
 * The test should be re-enabled and fixed at some point.
 */
#if 0
TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_accumulationVolumeBalance )
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
