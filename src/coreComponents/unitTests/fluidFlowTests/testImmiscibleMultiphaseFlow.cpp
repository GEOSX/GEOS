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

#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "constitutive/fluid/twophasefluid/TwoPhaseFluid.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlow.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML
char const *xmlInput =
  R"xml(
<?xml version="1.0" ?>

<Problem>
  <Solvers
    gravityVector="{ 0.0, 0.0, -9.8 }">
    <ImmiscibleMultiphaseFlow
      name="FlowSolver"
      discretization="TPFA"
      targetRegions="{ Domain }"
      logLevel="4"
      writeLinearSystem="0"
      temperature="300"
      initialDt="0.01">
      <NonlinearSolverParameters
	newtonTol="1.0e-6"
        lineSearchAction="None"
        newtonMaxIter="10"/>
      <LinearSolverParameters
        directParallel="0"/>
    </ImmiscibleMultiphaseFlow>
  </Solvers>

  <Mesh>
    <InternalMesh
      name="mesh1"
      elementTypes="{ C3D8 }"
      xCoords="{ 0, 7 }"
      yCoords="{ 0, 3 }"
      zCoords="{ 0, 3 }"
      nx="{ 7 }"
      ny="{ 3 }"
      nz="{ 3 }"
      cellBlockNames="{ block1 }"/>
  </Mesh>

  <Geometry>
    <Box
      name="source"
      xMin="{ -0.01, -0.01, -0.01 }"
      xMax="{ 1.01, 3.01, 3.01 }"/>

    <Box
      name="sink"
      xMin="{ 5.99, -0.01, -0.01 }"
      xMax="{ 7.01, 3.01, 3.01 }"/>
  </Geometry>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation
        name="TPFA"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion
      name="Domain"
      cellBlocks="{ block1 }"
      materialList="{ fluid, rock, relperm, capPressure }"/>
  </ElementRegions>


    <Functions>
       <TableFunction
         name="densityTablePhase1"
         coordinates="{ 0.0}"
         values="{ 1 }" />
       <TableFunction
         name="densityTablePhase2"
         coordinates="{ 0.0 }"
         values="{ 1 }" />
       <TableFunction
         name="viscosityTablePhase1"
         coordinates="{ 0.0}"
         values="{ 1.0 }" />                  
       <TableFunction
         name="viscosityTablePhase2"
         coordinates="{ 0.0 }"
         values="{ 1.0 }" />
    </Functions>


  <Constitutive>
    <TwoPhaseFluid
           name="fluid"
           phaseNames="{water22, gas22}"
           densityTableNames="{densityTablePhase1, densityTablePhase2}"    
           viscosityTableNames="{viscosityTablePhase1, viscosityTablePhase2}" />

    <CompressibleSolidConstantPermeability
      name="rock"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm" />

    <NullModel
      name="nullSolid" />

    <PressurePorosity
      name="rockPorosity"
      defaultReferencePorosity="1.0"
      referencePressure="0.0"
      compressibility="0.0" />

    <ConstantPermeability
      name="rockPerm"
      permeabilityComponents="{ 1.0, 1.0, 1.0 }"/>

    <BrooksCoreyRelativePermeability
      name="relperm"
      phaseNames="{ water, oil }"
      phaseMinVolumeFraction="{ 0.0, 0.0 }"
      phaseRelPermExponent="{ 1.0, 1.0 }"
      phaseRelPermMaxValue="{ 1.0, 1.0 }" />

    <BrooksCoreyCapillaryPressure
      name="capPressure"
      phaseNames="{ water, gas }"
      phaseMinVolumeFraction="{ 0.0, 0.0 }"
      phaseCapPressureExponentInv="{ 4 , 1 }"
      phaseEntryPressure="{ 0.75, 0 }"
      capPressureEpsilon="1e-8" />      
  </Constitutive>

  <FieldSpecifications>
    <FieldSpecification
      name="Porosity"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/Domain/block1"
      fieldName="rockPorosity_referencePorosity"
	    scale="1.0"/>

    <FieldSpecification
      name="initialPressure"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/Domain/block1"
      fieldName="pressure"
      scale="0.0"/>

    <FieldSpecification
      name="initialSat1"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/Domain/block1"
      fieldName="phaseVolumeFraction"
      component="0"
      scale="0.3"/>

    <FieldSpecification
      name="initialSat2"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/Domain/block1"
      fieldName="phaseVolumeFraction"
      component="1"
      scale="0.7"/>

    <FieldSpecification
      name="sourceTermP"
      objectPath="ElementRegions/Domain/block1"
      fieldName="pressure"
      scale="1.0"
      setNames="{ source }"/>

    <FieldSpecification
      name="sinkTermP"
      objectPath="ElementRegions/Domain/block1"
      fieldName="pressure"
      scale="0.0"
      setNames="{ sink }"/>
  
    <FieldSpecification
      name="sourceTerm_S1"
      objectPath="ElementRegions/Domain/block1"
      fieldName="phaseVolumeFraction"
      component="0"
      scale="1"
      setNames="{ source }"/>

     <FieldSpecification
      name="sinkTerm_S1"
      objectPath="ElementRegions/Domain/block1"
      fieldName="phaseVolumeFraction"
      component="0"	     
      scale="0."
      setNames="{ sink }"/>

     <FieldSpecification
      name="sourceTerm_S2"
      objectPath="ElementRegions/Domain/block1"
      fieldName="phaseVolumeFraction"
      component="1"
      scale="0."
      setNames="{ source }"/>

     <FieldSpecification
      name="sinkTerm_S2"
      objectPath="ElementRegions/Domain/block1"
      fieldName="phaseVolumeFraction"
      component="1"
      scale="1."
      setNames="{ sink }"/>

  </FieldSpecifications>
</Problem>
  
  )xml";
// Sphinx end before input XML

template< typename LAMBDA >
void fillCellCenteredNumericalJacobian( ImmiscibleMultiphaseFlow & solver,
                                        DomainPartition & domain,
                                        real64 const perturbParameter,
                                        arrayView1d< real64 > residual,
                                        arrayView1d< real64 > residualOrig,
                                        CRSMatrixView< real64, globalIndex > jacobian,
                                        CRSMatrixView< real64, globalIndex > jacobianFD,
                                        LAMBDA assembleFunction )
{
  DofManager const & dofManager = solver.getDofManager();
  string const elemDofKey = dofManager.getKey( ImmiscibleMultiphaseFlow::viewKeyStruct::elemDofFieldString());

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );

      arrayView1d< real64 > const pres =
        subRegion.getField< fields::flow::pressure >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE >  const phaseVolumeFraction =
        subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();


      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( elemGhostRank[ei] >= 0 )
        {
          continue;
        }

        // Step 1: compute numerical derivatives wrt pressure
        solver.resetStateToBeginningOfStep( domain );
        pres.move( hostMemorySpace, true );

        real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
        pres[ei] += dP;

#if defined(GEOS_USE_CUDA)
        pres.move( parallelDeviceMemorySpace, false );
#endif
        solver.updateState( domain );

        residual.zero();
        jacobian.zero();
        assembleFunction( jacobian.toViewConstSizes(), residual.toView());

        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               elemDofNumber[ei],
                               dP,
                               jacobianFD.toViewConstSizes());


        // Step 2: compute numerical derivatives wrt saturation
        solver.resetStateToBeginningOfStep( domain );
        phaseVolumeFraction.move( hostMemorySpace, true );

        real64 const dS = perturbParameter * (phaseVolumeFraction[ei][0] + perturbParameter);
        phaseVolumeFraction[ei][0] += dS;
        phaseVolumeFraction[ei][1] = 1- phaseVolumeFraction[ei][0];

#if defined(GEOS_USE_CUDA)
        phaseVolumeFraction.move( parallelDeviceMemorySpace, false );
#endif

        solver.updateState( domain );

        residual.zero();
        jacobian.zero();
        assembleFunction( jacobian.toViewConstSizes(), residual.toView());

        fillNumericalJacobian( residual.toViewConst(),
                               residualOrig.toViewConst(),
                               elemDofNumber[ei]+1,
                               dS,
                               jacobianFD.toViewConstSizes());
      }
    } );
  } );
}

template< typename LAMBDA >
void testNumericalJacobian( ImmiscibleMultiphaseFlow & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  CRSMatrix< real64, globalIndex > const & jacobian = solver.getLocalMatrix();
  array1d< real64 > residual( jacobian.numRows());

  // assemble the analytical residual
  solver.resetStateToBeginningOfStep( domain );

  residual.zero();
  jacobian.zero();

  assembleFunction( jacobian.toViewConstSizes(), residual.toView());
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
  assembleFunction( jacobian.toViewConstSizes(), residual.toView());
  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

class ImmiscibleMultiphaseFlowTest : public ::testing::Test
{
public:
  ImmiscibleMultiphaseFlowTest(): state( std::make_unique< CommandLineOptions >( g_commandLineOptions ))
  {}

protected:
  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< ImmiscibleMultiphaseFlow >( "FlowSolver" );

    DomainPartition & domain = state.getProblemManager().getDomainPartition();

    solver->setupSystem( domain,
                         solver->getDofManager(),
                         solver->getLocalMatrix(),
                         solver->getSystemRhs(),
                         solver->getSystemSolution());

    solver->implicitStepSetup( time, dt, domain );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 1e4;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  ImmiscibleMultiphaseFlow *solver;
};

real64 constexpr ImmiscibleMultiphaseFlowTest::time;
real64 constexpr ImmiscibleMultiphaseFlowTest::dt;
real64 constexpr ImmiscibleMultiphaseFlowTest::eps;


TEST_F( ImmiscibleMultiphaseFlowTest, jacobianNumericalCheck_fluxTerm )
{
  real64 const perturb = 0.000001;
  real64 const tol = 1e-2; // 1% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&]( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
  {
    solver->assembleFluxTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}


TEST_F( ImmiscibleMultiphaseFlowTest, jacobianNumericalCheck_accumulationTerm )
{
  real64 const perturb = sqrt( eps );
  real64 const tol = 1e-2; // 1% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleAccumulationTerm( domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}


TEST_F( ImmiscibleMultiphaseFlowTest, jacobianNumericalCheck_full )
{
  real64 const perturb = 0.000001;
  real64 const tol = 1e-2; // 1% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleSystem( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}


int main( int argc, char * *argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
