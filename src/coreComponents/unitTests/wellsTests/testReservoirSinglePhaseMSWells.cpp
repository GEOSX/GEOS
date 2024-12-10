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

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/SinglePhaseReservoirAndWells.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"

#include "tests/meshDirName.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

char const * PreXmlInput =
  R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <SinglePhaseReservoir name="reservoirSystem"
                 flowSolverName="singlePhaseFlow"
                 wellSolverName="singlePhaseWell"
                 logLevel="1"
                 targetRegions="{Region1,wellRegion1,wellRegion2}">
        <NonlinearSolverParameters newtonMaxIter="40"/>
        <LinearSolverParameters solverType="direct"
                                logLevel="2"/>
      </SinglePhaseReservoir>
      <SinglePhaseFVM name="singlePhaseFlow"
                               logLevel="1"
                               discretization="singlePhaseTPFA"
                               targetRegions="{Region1}">
      </SinglePhaseFVM>
      <SinglePhaseWell name="singlePhaseWell"
                       logLevel="1"
                       targetRegions="{wellRegion1,wellRegion2}">
          <WellControls name="wellControls1"
                        type="producer"
                        referenceElevation="2"
                        control="BHP"
                        targetBHP="5e5"
                        targetTotalRate="1e-3"/>
          <WellControls name="wellControls2"
                        type="injector"
                        referenceElevation="2"
                        control="totalVolRate"
                        targetBHP="2e7"
                        targetTotalRate="1e-4"/>
      </SinglePhaseWell>
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
                    cellBlockNames="{cb1}">)xml";

char const * PostXmlInput =
  R"xml(
      </InternalMesh>
    </Mesh>
    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="singlePhaseTPFA"/>
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="Region1"
                         cellBlocks="{ * }"
                         materialList="{water, rock}"/>
      <WellElementRegion name="wellRegion1"
                         materialList="{water}"/>
      <WellElementRegion name="wellRegion2"
                         materialList="{water}"/>
    </ElementRegions>
    <Constitutive>
      <CompressibleSinglePhaseFluid name="water"
                                    defaultDensity="1000"
                                    defaultViscosity="0.001"
                                    referencePressure="0.0"
                                    referenceDensity="1000"
                                    compressibility="5e-10"
                                    referenceViscosity="0.001"
                                    viscosibility="0.0"/>
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
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification name="initialPressure"
                          initialCondition="1"
                          setNames="{all}"
                          objectPath="ElementRegions/Region1/cb1"
                          fieldName="pressure"
                          scale="5e6"/>
    </FieldSpecifications>
  </Problem>
  )xml";

template< typename LAMBDA >
void testNumericalJacobian( SinglePhaseReservoirAndWells<> & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA && assembleFunction )
{
  SinglePhaseWell & wellSolver = *solver.wellSolver();
  SinglePhaseFVM< SinglePhaseBase > & flowSolver = dynamicCast< SinglePhaseFVM< SinglePhaseBase > & >( *solver.reservoirSolver() );

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

  string const & resDofKey  = dofManager.getKey( wellSolver.resElementDofName() );
  string const & wellDofKey = dofManager.getKey( wellSolver.wellElementDofName() );

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
        elemRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
        {
          // get the dof numbers and ghosting information
          arrayView1d< globalIndex const > const & dofNumber =
            subRegion.getReference< array1d< globalIndex > >( resDofKey );

          // get the primary variables on reservoir elements
          arrayView1d< real64 > const & pres =
            subRegion.getField< fields::well::pressure >();
          pres.move( hostMemorySpace, false );

          // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
          for( localIndex ei = 0; ei < subRegion.size(); ++ei )
          {
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

      // get the degrees of freedom and ghosting information
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the primary variables on well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getField< fields::well::pressure >();
      wellElemPressure.move( hostMemorySpace, false );

      arrayView1d< real64 > const & connRate =
        subRegion.getField< fields::well::connectionRate >();
      connRate.move( hostMemorySpace, false );

      // a) compute all the derivatives wrt to the pressure in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
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

          // consider mass balance eq lid in RESERVOIR elems and WELL elems
          //      this is computing J_RW and J_WW
          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] + singlePhaseWellKernels::ColOffset::DPRES,
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }
      }

      // b) compute all the derivatives wrt to the connection in WELL elem iwelem
      for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
      {
        {
          solver.resetStateToBeginningOfStep( domain );

          // here is the perturbation in the pressure of the well element
          real64 const dRate = perturbParameter * ( connRate[iwelem] + perturbParameter );
          connRate.move( hostMemorySpace, true );
          connRate[iwelem] += dRate;

          // after perturbing, update the rate-dependent quantities in the well (well controls)
          wellSolver.updateState( domain );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          // consider mass balance eq lid in RESERVOIR elems and WELL elems
          //      this is computing J_RW and J_WW
          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 wellElemDofNumber[iwelem] + singlePhaseWellKernels::ColOffset::DRATE,
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


class SinglePhaseReservoirSolverTest : public ::testing::Test
{
public:

  SinglePhaseReservoirSolverTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< SinglePhaseReservoirAndWells<> >( "reservoirSystem" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution() );

    solver->implicitStepSetup( TIME, DT, domain );
  }

  void TestAssembleCouplingTerms()
  {
    real64 const perturb = std::sqrt( EPS );
    real64 const tol = 1e-1; // 10% error margin

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    testNumericalJacobian( *solver, domain, perturb, tol,
                           [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
    {
      solver->assembleCouplingTerms( TIME, DT, domain, solver->getDofManager(), localMatrix, localRhs );
    } );
  }

  void TestAssembleFluxTerms()
  {
    real64 const perturb = std::sqrt( EPS );
    real64 const tol = 1e-1; // 10% error margin

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    testNumericalJacobian( *solver, domain, perturb, tol,
                           [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
    {
      solver->wellSolver()->assembleFluxTerms( TIME, DT, domain, solver->getDofManager(), localMatrix, localRhs );
    } );
  }

  void TestAssemblePressureRelations()
  {
    real64 const perturb = std::sqrt( EPS );
    real64 const tol = 1e-1; // 10% error margin

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    testNumericalJacobian( *solver, domain, perturb, tol,
                           [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
    {
      solver->wellSolver()->assemblePressureRelations( TIME, DT, domain, solver->getDofManager(), localMatrix, localRhs );
    } );
  }

  void TestAssembleAccumulationTerms()
  {
    real64 const perturb = std::sqrt( EPS );
    real64 const tol = 1e-1; // 10% error margin

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    testNumericalJacobian( *solver, domain, perturb, tol,
                           [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
    {
      solver->wellSolver()->assembleAccumulationTerms( TIME, DT, domain, solver->getDofManager(), localMatrix, localRhs );
    } );
  }

  static real64 constexpr TIME = 0.0;
  static real64 constexpr DT = 1e4;
  static real64 constexpr EPS = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  SinglePhaseReservoirAndWells<> * solver;
};

real64 constexpr SinglePhaseReservoirSolverTest::TIME;
real64 constexpr SinglePhaseReservoirSolverTest::DT;
real64 constexpr SinglePhaseReservoirSolverTest::EPS;

/**
 * @brief Test SinglePhaseReservoirSolver with InternalWell generator
 *
 */
class SinglePhaseReservoirSolverInternalWellTest : public SinglePhaseReservoirSolverTest
{

public:

  SinglePhaseReservoirSolverInternalWellTest():
    SinglePhaseReservoirSolverTest()
  {}

protected:

  void SetUp() override
  {
    string const internalWells =
      R"(
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
        </InternalWell>)";

    string const xmlInput = PreXmlInput + internalWells + PostXmlInput;

    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str() );
    SinglePhaseReservoirSolverTest::SetUp();
  }
};


TEST_F( SinglePhaseReservoirSolverInternalWellTest, jacobianNumericalCheck_Perforation )
{
  TestAssembleCouplingTerms();
}

TEST_F( SinglePhaseReservoirSolverInternalWellTest, jacobianNumericalCheck_Flux )
{
  real64 const perturb = std::sqrt( EPS );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assembleFluxTerms( TIME, DT, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

TEST_F( SinglePhaseReservoirSolverInternalWellTest, jacobianNumericalCheck_PressureRel )
{
  TestAssemblePressureRelations();
}

TEST_F( SinglePhaseReservoirSolverInternalWellTest, jacobianNumericalCheck_Accum )
{
  TestAssembleAccumulationTerms();
}


/**
 * @brief Test SinglePhaseReservoirSolver with VTKWell generator
 *
 */
class SinglePhaseReservoirSolverVTKWellTest : public SinglePhaseReservoirSolverTest
{

public:

  SinglePhaseReservoirSolverVTKWellTest():
    SinglePhaseReservoirSolverTest()
  {}

protected:

  void SetUp() override
  {
    string const vtkWells =  GEOS_FMT( R"(
        <VTKWell name="well_producer1"
                      wellRegionName="wellRegion1"
                      wellControlsName="wellControls1"
                      file="{}"
                      radius="0.1"
                      numElementsPerSegment="1">
            <Perforation name="producer1_perf1"
                         distanceFromHead="1.45"/>
        </VTKWell>
        <VTKWell name="well_injector1"
                      wellRegionName="wellRegion2"
                      wellControlsName="wellControls2"
                      file="{}"
                      radius="0.1"
                      numElementsPerSegment="1">
            <Perforation name="injector1_perf1"
                         distanceFromHead="1.45"/>
        </VTKWell>)", testMeshDir + "/well1.vtk",
                                       testMeshDir + "/well2.vtk" );

    string const xmlInput = PreXmlInput + vtkWells + PostXmlInput;

    setupProblemFromXML( state.getProblemManager(), xmlInput.c_str());
    SinglePhaseReservoirSolverTest::SetUp();
  }
};


TEST_F( SinglePhaseReservoirSolverVTKWellTest, jacobianNumericalCheck_Perforation )
{
  TestAssembleCouplingTerms();
}

TEST_F( SinglePhaseReservoirSolverVTKWellTest, jacobianNumericalCheck_Flux )
{
  TestAssembleFluxTerms();
}

TEST_F( SinglePhaseReservoirSolverVTKWellTest, jacobianNumericalCheck_PressureRel )
{
  TestAssemblePressureRelations();
}

TEST_F( SinglePhaseReservoirSolverVTKWellTest, jacobianNumericalCheck_Accum )
{
  TestAssembleAccumulationTerms();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
