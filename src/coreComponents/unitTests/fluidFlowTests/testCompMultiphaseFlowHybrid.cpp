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

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "mainInterface/initialization.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
      <CompositionalMultiphaseHybridFVM name="compflow"
                                   logLevel="0"
                                   discretization="fluidHM"
                                   targetRegions="{Region}"
                                   temperature="297.15"
                                   useMass="1">

        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="2"/>
        <LinearSolverParameters solverType="gmres"
                                krylovTol="1.0e-10"/>
      </CompositionalMultiphaseHybridFVM>
    </Solvers>
    <Mesh>
      <InternalMesh name="mesh1"
                    elementTypes="{C3D8}"
                    xCoords="{0, 1}"
                    yCoords="{0, 1}"
                    zCoords="{0, 10}"
                    nx="{1}"
                    ny="{1}"
                    nz="{4}"
                    cellBlockNames="{cb1}"/>
    </Mesh>
    <NumericalMethods>
      <FiniteVolume>
        <HybridMimeticDiscretization name="fluidHM"
                                     innerProductType="beiraoDaVeigaLipnikovManzini"/>
      </FiniteVolume>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion name="Region" cellBlocks="{cb1}" materialList="{fluid1, rock, relperm}" />
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
                 objectPath="ElementRegions/Region/cb1"
                 fieldName="pressure"
                 functionName="initialPressureFunc"
                 scale="5e6"/>
      <FieldSpecification name="initialFacePressure"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="faceManager"
                 fieldName="facePressure"
                 functionName="initialFacePressureFunc"
                 scale="5e6"/>
      <FieldSpecification name="initialComposition_N2"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region/cb1"
                 fieldName="globalCompFraction"
                 component="0"
                 scale="0.099"/>
      <FieldSpecification name="initialComposition_C10"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region/cb1"
                 fieldName="globalCompFraction"
                 component="1"
                 scale="0.3"/>
      <FieldSpecification name="initialComposition_C20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region/cb1"
                 fieldName="globalCompFraction"
                 component="2"
                 scale="0.6"/>
      <FieldSpecification name="initialComposition_H20"
                 initialCondition="1"
                 setNames="{all}"
                 objectPath="ElementRegions/Region/cb1"
                 fieldName="globalCompFraction"
                 component="3"
                 scale="0.001"/>
    </FieldSpecifications>
    <Functions>
      <TableFunction name="initialPressureFunc"
                     inputVarNames="{elementCenter}"
                     coordinates="{0.0, 2.0, 4.0, 6.0, 8.0, 10.0 }"
                     values="{ 1.0, 0.5, 0.2, 3.0, 2.1, 1.0 }"/>
      <TableFunction name="initialFacePressureFunc"
                     inputVarNames="{faceCenter}"
                     coordinates="{0.0, 2.0, 4.0, 6.0, 8.0, 10.0 }"
                     values="{ 2.0, 0.1, 0.2, 2.1, 2.1, 0.1 }"/>
    </Functions>
  </Problem>
  )xml";

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseHybridFVM & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
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

  // fill jacobian FD for cell centered variables
  fillCellCenteredNumericalJacobian( solver,
                                     domain,
                                     false,
                                     perturbParameter,
                                     residual.toView(),
                                     residualOrig.toView(),
                                     jacobian.toView(),
                                     jacobianFD.toView(),
                                     assembleFunction );

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                       MeshLevel & mesh,
                                                                       arrayView1d< string const > const & )
  {

    FaceManager & faceManager = mesh.getFaceManager();

    // get the face-based pressure
    arrayView1d< real64 > const & facePres =
      faceManager.getField< fields::flow::facePressure >();
    facePres.move( hostMemorySpace, false );

    string const faceDofKey = dofManager.getKey( CompositionalMultiphaseHybridFVM::viewKeyStruct::faceDofFieldString() );

    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    faceDofNumber.move( hostMemorySpace );

    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    faceGhostRank.move( hostMemorySpace );

    for( localIndex iface = 0; iface < faceManager.size(); ++iface )
    {
      if( faceGhostRank[iface] >= 0 )
      {
        continue;
      }

      solver.resetStateToBeginningOfStep( domain );

      facePres.move( hostMemorySpace, true ); // to get the correct facePres after reset
      real64 const dFP = perturbParameter * ( facePres[iface] + perturbParameter );
      facePres[iface] += dFP;
#if defined(GEOS_USE_CUDA)
      facePres.move( parallelDeviceMemorySpace, false );
#endif


      residual.zero();
      jacobian.zero();
      assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

      fillNumericalJacobian( residual.toViewConst(),
                             residualOrig.toViewConst(),
                             faceDofNumber[iface],
                             dFP,
                             jacobianFD.toViewConstSizes() );
    }

    // assemble the analytical jacobian
    solver.resetStateToBeginningOfStep( domain );

    residual.zero();
    jacobian.zero();
    assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

    compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );

  } );
}

class CompositionalMultiphaseHybridFlowTest : public ::testing::Test
{
public:

  CompositionalMultiphaseHybridFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseHybridFVM >( "compflow" );

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
  CompositionalMultiphaseHybridFVM * solver;
};

real64 constexpr CompositionalMultiphaseHybridFlowTest::time;
real64 constexpr CompositionalMultiphaseHybridFlowTest::dt;
real64 constexpr CompositionalMultiphaseHybridFlowTest::eps;


TEST_F( CompositionalMultiphaseHybridFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 5e-3; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
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
