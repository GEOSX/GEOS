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

#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <CompositionalMultiphaseReservoir name="reservoirSystem"
                 flowSolverName="compositionalMultiphaseFlow"
                 wellSolverName="compositionalMultiphaseWell"
                 logLevel="1"
                 targetRegions="{Region1,wellRegion1,wellRegion2}">
        <NonlinearSolverParameters newtonMaxIter="40"/>
        <LinearSolverParameters solverType="direct"
                                logLevel="2"/>
      </CompositionalMultiphaseReservoir>
      <CompositionalMultiphaseFVM name="compositionalMultiphaseFlow"
                                  logLevel="1"
                                  discretization="fluidTPFA"
                                  targetRegions="{Region1}"
                                  temperature="297.15"
                                  useMass="0">
      </CompositionalMultiphaseFVM>
      <CompositionalMultiphaseWell name="compositionalMultiphaseWell"
                                   logLevel="1"
                                   targetRegions="{wellRegion1,wellRegion2}"
                                   useMass="0">
          <WellControls name="wellControls1"
                        type="producer"
                        referenceElevation="1.25"
                        control="BHP"
                        targetBHP="2e6"
                        targetPhaseRate="1"
                        targetPhaseName="oil"/>
          <WellControls name="wellControls2"
                        type="injector"
                        referenceElevation="1.25"
                        control="totalVolRate"
                        targetBHP="6e7"
                        targetTotalRate="1e-5"
                        injectionTemperature="297.15"
                        injectionStream="{0.1, 0.1, 0.1, 0.7}"/>
      </CompositionalMultiphaseWell>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh1"
                    elementTypes="{C3D8}"
                    xCoords="{0, 5}"
                    yCoords="{0, 1}"
                    zCoords="{0, 1}"
                    nx="{3}"
                    ny="{1}"
                    nz="{1}"
                    cellBlockNames="{cb1}">
        <InternalWell name="well_producer1"
                      wellRegionName="wellRegion1"
                      wellControlsName="wellControls1"
                      polylineNodeCoords="{ {4.5, 0,  2  },
                                             {4.5, 0,  0.5} }"
                      polylineSegmentConn="{ {0, 1} }"
                      radius="0.1"
                      numElementsPerSegment="1">
            <Perforation name="producer1_perf1"
                         distanceFromHead="1.45"/>
        </InternalWell>
        <InternalWell name="well_injector1"
                      wellRegionName="wellRegion2"
                      wellControlsName="wellControls2"
                      polylineNodeCoords="{ {0.5, 0, 2  },
                                             {0.5, 0, 0.5} }"
                      polylineSegmentConn="{ {0, 1} }"
                      radius="0.1"
                      numElementsPerSegment="1">
            <Perforation name="injector1_perf1"
                         distanceFromHead="1.45"/>
        </InternalWell>
      </InternalMesh>
    </Mesh>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="fluidTPFA"/>
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="Region1"
                         cellBlocks="{cb1}"
                         materialList="{fluid1, rock, relperm}"/>
      <WellElementRegion name="wellRegion1"
                         materialList="{fluid1, relperm}"/>
      <WellElementRegion name="wellRegion2"
                         materialList="{fluid1, relperm}"/>
    </ElementRegions>
    <Constitutive>
      <CompositionalMultiphaseFluid name="fluid1"
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
    <ConstantPermeability name="rockPerm"
                          permeabilityComponents="{2.0e-16, 2.0e-16, 2.0e-16}"/>
      <BrooksCoreyRelativePermeability name="relperm"
                                       phaseNames="{oil, gas}"
                                       phaseMinVolumeFraction="{0.1, 0.15}"
                                       phaseRelPermExponent="{2.0, 2.0}"
                                       phaseRelPermMaxValue="{0.8, 0.9}"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="pressure"
                 scale="5e6"/>
      <FieldSpecification name="initialComposition_N2"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="0"
                 scale="0.099"/>
      <FieldSpecification name="initialComposition_C10"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="1"
                 scale="0.3"/>
      <FieldSpecification name="initialComposition_C20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="2"
                 scale="0.6"/>
      <FieldSpecification name="initialComposition_H20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region1/cb1"
                 fieldName="globalCompFraction"
                 component="3"
                 scale="0.001"/>
    </FieldSpecifications>
  </Problem>
  )xml";

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseReservoirAndWells<> & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA && assembleFunction )
{
  CompositionalMultiphaseWell & wellSolver = *solver.wellSolver();
  CompositionalMultiphaseFVM & flowSolver = dynamicCast< CompositionalMultiphaseFVM & >( *solver.reservoirSolver() );

  localIndex const NC = flowSolver.numFluidComponents();

  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows() );
  DofManager const & dofManager = solver.getDofManager();

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

  string const resDofKey  = dofManager.getKey( wellSolver.resElementDofName() );
  string const wellDofKey = dofManager.getKey( wellSolver.wellElementDofName() );

  // at this point we start assembling the finite-difference block by block

  ////////////////////////////////////////////////
  // Step 1) Compute the terms in J_RR and J_WR //
  ////////////////////////////////////////////////
  domain.forMeshBodies( [&] ( MeshBody & meshBody )
  {
    meshBody.forMeshLevels( [&] ( MeshLevel & mesh )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      for( localIndex er = 0; er < elemManager.numRegions(); ++er )
      {
        ElementRegionBase & elemRegion = elemManager.getRegion( er );
        elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const, CellElementSubRegion & subRegion )
        {
          // get the degrees of freedom and ghosting information
          arrayView1d< globalIndex const > const & dofNumber =
            subRegion.getReference< array1d< globalIndex > >( resDofKey );

          // get the primary variables on the reservoir elements
          arrayView1d< real64 > const & pres =
            subRegion.getField< fields::well::pressure >();
          pres.move( hostMemorySpace, false );

          arrayView2d< real64, compflow::USD_COMP > const & compDens =
            subRegion.getField< fields::well::globalCompDensity >();
          compDens.move( hostMemorySpace, false );

          // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
          for( localIndex ei = 0; ei < subRegion.size(); ++ei )
          {
            real64 totalDensity = 0.0;
            for( localIndex ic = 0; ic < NC; ++ic )
            {
              totalDensity += compDens[ei][ic];
            }

            {
              solver.resetStateToBeginningOfStep( domain );

              // here is the perturbation in the pressure of the element
              real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
              pres.move( hostMemorySpace, true );
              pres[ei] += dP;

              // after perturbing, update the pressure-dependent quantities in the reservoir
              flowSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                       MeshLevel & mesh2,
                                                                                       arrayView1d< string const > const & regionNames2 )
              {
                mesh2.getElemManager().forElementSubRegions( regionNames2,
                                                             [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion2 )
                {
                  flowSolver.updateFluidState( subRegion2 );
                } );
              } );

              wellSolver.updateState( domain );

              residual.zero();
              jacobian.zero();
              assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

              fillNumericalJacobian( residual.toViewConst(),
                                     residualOrig.toViewConst(),
                                     dofNumber[ei],
                                     dP,
                                     jacobianFD.toViewConstSizes() );
            }

            for( localIndex jc = 0; jc < NC; ++jc )
            {
              solver.resetStateToBeginningOfStep( domain );

              real64 const dRho = perturbParameter * totalDensity;
              compDens.move( hostMemorySpace, true );
              compDens[ei][jc] += dRho;

              flowSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                       MeshLevel & mesh2,
                                                                                       arrayView1d< string const > const & regionNames2 )
              {
                mesh2.getElemManager().forElementSubRegions( regionNames2,
                                                             [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion2 )
                {
                  flowSolver.updateFluidState( subRegion2 );
                } );
              } );

              residual.zero();
              jacobian.zero();
              assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

              fillNumericalJacobian( residual.toViewConst(),
                                     residualOrig.toViewConst(),
                                     dofNumber[ei] + jc + 1,
                                     dRho,
                                     jacobianFD.toViewConstSizes() );
            }
          }
        } );
      }
    } );
  } );

  /////////////////////////////////////////////////
  // Step 2) Compute the terms in J_RW and J_WW //
  /////////////////////////////////////////////////

  // loop over the wells
  wellSolver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                           MeshLevel & mesh,
                                                                           arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {
      // get the degrees of freedom, ghosting info and next well elem index
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the primary variables on the well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getField< fields::well::pressure >();
      wellElemPressure.move( hostMemorySpace, false );

      arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getField< fields::well::globalCompDensity >();
      wellElemCompDens.move( hostMemorySpace, false );

      arrayView1d< real64 > const & connRate =
        subRegion.getField< fields::well::mixtureConnectionRate >();
      connRate.move( hostMemorySpace, false );

      // a) compute all the derivatives wrt to the pressure in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
        real64 wellElemTotalDensity = 0.0;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          wellElemTotalDensity += wellElemCompDens[iwelem][ic];
        }

        {
          solver.resetStateToBeginningOfStep( domain );

          // here is the perturbation in the pressure of the well element
          real64 const dP = perturbParameter * ( wellElemPressure[iwelem] + perturbParameter );
          wellElemPressure.move( hostMemorySpace, true );
          wellElemPressure[iwelem] += dP;

          // after perturbing, update the pressure-dependent quantities in the well
          wellSolver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DPRES,
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          solver.resetStateToBeginningOfStep( domain );

          real64 const dRho = perturbParameter * wellElemTotalDensity;
          wellElemCompDens.move( hostMemorySpace, true );
          wellElemCompDens[iwelem][jc] += dRho;

          wellSolver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + jc,
                                 dRho,
                                 jacobianFD.toViewConstSizes() );
        }
      }

      // b) compute all the derivatives wrt to the connection in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
        {
          solver.resetStateToBeginningOfStep( domain );

          // here is the perturbation in the rate of the well element
          real64 const dRate = perturbParameter * ( connRate[iwelem] + perturbParameter );
          connRate.move( hostMemorySpace, true );
          connRate[iwelem] += dRate;

          wellSolver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC,
                                 dRate,
                                 jacobianFD.toViewConstSizes() );
        }
      }
    } );
  } );

  // assemble the analytical jacobian
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();
  assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class CompositionalMultiphaseReservoirSolverTest : public ::testing::Test
{
public:

  CompositionalMultiphaseReservoirSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseReservoirAndWells<> >( "reservoirSystem" );

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
  CompositionalMultiphaseReservoirAndWells<> * solver;
};

real64 constexpr CompositionalMultiphaseReservoirSolverTest::time;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::dt;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::eps;

#if 0

TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Perforation )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = *state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleCouplingTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#endif


TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}


TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_VolumeBalance )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assembleVolumeBalanceTerms( domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_PressureRel )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assemblePressureRelations( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
