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
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/ReactiveCompositionalMultiphaseOBLFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ReactiveCompositionalMultiphaseOBL.hpp"
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
      <ReactiveCompositionalMultiphaseOBL
                        name="compflow"
                        logLevel="1"
                        discretization="fluidTPFA"
                        targetRegions="{ region }"
                        componentNames="{N2, C10, C20}"
                        enableEnergyBalance="0"
                        maxCompFractionChange="1"
                        numComponents="3"
                        numPhases="2"
                        OBLOperatorsTableFile="obl_3comp_static.txt">

        <NonlinearSolverParameters newtonTol="1.0e-6"
                                   newtonMaxIter="2"/>
        <LinearSolverParameters solverType="gmres"
                                krylovTol="1.0e-10"/>
      </ReactiveCompositionalMultiphaseOBL>
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
      <CellElementRegion name="region" cellBlocks="{cb1}" materialList="{rock}" />
    </ElementRegions>
    <Constitutive>
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
    </FieldSpecifications>
    <Functions>
      <TableFunction name="initialPressureFunc"
                     inputVarNames="{elementCenter}"
                     coordinates="{0.0, 3.0}"
                     values="{1.0, 0.5}"/>
    </Functions>
  </Problem>
  )xml";

char const * oblInput =
  "3 28 \n"
  "2 1 300 \n"
  "2 1.0000000000000001e-09 0.99999999900000003 \n"
  "2 1.0000000000000001e-09 0.99999999900000003 \n"
  "3.3305578640016657e-08 3.3305578640016657e-08 33.305578573405498 0 0 0 6.6611157280033315e-08 6.6611157280033315e-08 66.611157146810996 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 600 0 0 1  \n"
  "1.2466496269485102e-08 12.466496257018607 0 2.4932992538970203e-07 249.32992514037213 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 200 0 0 0 1  \n"
  "4.5444217206792512 4.544421725223673e-09 0 90.888434413585017 9.088843450447346e-08 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 200 0 0 0 1  \n"
  "6.6607829710424458 6.6607829710424458 -6.6607829643816636 66.607829710424454 66.607829710424454 1.3321565955406458e-07 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 200 0 0 0 1  \n"
  "3.341515046368404e-08 3.341515046368404e-08 33.415150396853733 0 0 0 6.6810324640300631e-08 6.6810324640300631e-08 66.810324506679962 0 1.000299 0 0 0 0 0 0 0 0 0 0 0 0 0 601.79399999999998 0 0 1  \n"
  "1.6198820653678712e-08 16.198820637479891 0 3.2387957308122294e-07 323.87957275734334 0 0 0 0 1.000299 0 0 0 0 0 0 0 0 0 0 0 0 0 259.80000000000001 0 0 0 1  \n"
  "5.9049688731030807 5.9049688790080501e-09 0 118.06407630324695 1.1806407642131103e-07 0 0 0 0 1.000299 0 0 0 0 0 0 0 0 0 0 0 0 0 259.80000000000001 0 0 0 1  \n"
  "8.6549441341508739 8.6549441341508739 -8.6549441254959287 86.523570793841358 86.523570793841358 1.7304714176072988e-07 0 0 0 1.000299 0 0 0 0 0 0 0 0 0 0 0 0 0 259.80000000000001 0 0 0 1  \n";

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

template< int USD >
array1d< real64 > getPressureDerivative( arraySlice2d< real64 const, USD > const & input,
                                         localIndex N1 )
{
  array1d< real64 > output( N1 );

  for( localIndex i = 0; i < N1; ++i )
  {
    output( i ) = input( i, 0 ) * 1e-5; // Pa<--->bar conversion
  }

  return output;
}

template< int USD >
array2d< real64 > getCompFracDerivative( arraySlice2d< real64 const, USD > const & input,
                                         localIndex N1,
                                         localIndex N2 )
{
  array2d< real64 > output( N2 - 1, N1 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2 - 1; ++j )
    {
      output( j, i ) = input( i, j + 1 );
    }
  }

  return output;
}


// Sphinx end before input XML

void testOperatorsNumericalDerivatives( ReactiveCompositionalMultiphaseOBL & solver,
                                        DomainPartition & domain,
                                        real64 const perturbParameter,
                                        real64 const relTol )
{
  localIndex const NC = solver.numFluidComponents();
  localIndex const NOPS = solver.numOBLOperators();

  array1d< string > const operators( NOPS );

  // update component Fraction and check derivatives
  for( localIndex op = 0; op < NOPS; ++op )
  {
    operators[op] = std::to_string( op );
  }

  solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                         [&]( string const,
                                              MeshLevel & mesh,
                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions( regionNames,
                                               [&]( localIndex const,
                                                    ElementSubRegionBase & subRegion )
    {
      SCOPED_TRACE( subRegion.getParent().getParent().getName() + "/" + subRegion.getName() );

      arrayView1d< string const > const & components = solver.componentNames();

      arrayView1d< real64 > const & pres =
        subRegion.getField< fields::flow::pressure >();

      arrayView1d< real64 const > const & pres_n =
        subRegion.getField< fields::flow::pressure_n >();

      arrayView2d< real64, compflow::USD_COMP > const & compFrac =
        subRegion.getField< fields::flow::globalCompFraction >();

      arrayView2d< real64 const, compflow::USD_COMP > const & compFrac_n =
        subRegion.getField< fields::flow::globalCompFraction_n >();

      arrayView2d< real64 const, compflow::USD_OBL_VAL > const & OBLVals =
        subRegion.getField< fields::flow::OBLOperatorValues >();

      arrayView2d< real64 const, compflow::USD_OBL_VAL > const & OBLVals_n =
        subRegion.getField< fields::flow::OBLOperatorValues_n >();

      arrayView3d< real64 const, compflow::USD_OBL_DER > const & OBLDers =
        subRegion.getField< fields::flow::OBLOperatorDerivatives >();

      // reset the solver state to zero out variable updates
      solver.resetStateToBeginningOfStep( domain );

      // update pressure and check derivatives
      {
        // perturb pressure in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
        } );

        // recompute operators
        solver.updateOBLOperators( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          // get all pressure derivatives of every operator for current element in a 1d array
          auto dOp_dP = getPressureDerivative( OBLDers[ei].toSliceConst(), NOPS );

          real64 const delta = pres[ei] - pres_n[ei];
          checkDerivative( OBLVals[ei].toSliceConst(),
                           OBLVals_n[ei].toSliceConst(),
                           dOp_dP.toSliceConst(),
                           delta,
                           relTol,
                           "operator",
                           "Pres",
                           operators );
        } );
      }

      // update component fraction and check derivatives for all components except the last one (which  is not primary variable)
      for( localIndex jc = 0; jc < NC - 1; ++jc )
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver.resetStateToBeginningOfStep( domain );

        // perturb a single component fraction in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          real64 const dZ = perturbParameter * ( compFrac[ei][jc] + perturbParameter );
          compFrac[ei][jc] += dZ;
        } );

        // recompute operators
        solver.updateOBLOperators( subRegion );

        // check values in each cell
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
        {
          SCOPED_TRACE( "Element " + std::to_string( ei ) );

          // get NC-1 component fraction derivatives of every operator for current element in a 2d array (NC-1, NOPS)
          auto dOp_dZ = getCompFracDerivative( OBLDers[ei].toSliceConst(), NOPS, NC );
          string var = "compFrac[" + components[jc] + "]";

          real64 const delta = compFrac[ei][jc] - compFrac_n[ei][jc];
          checkDerivative( OBLVals[ei].toSliceConst(),
                           OBLVals_n[ei].toSliceConst(),
                           dOp_dZ[jc].toSliceConst(),
                           delta,
                           relTol,
                           "operator",
                           var,
                           operators );
        } );
      }
    } );
  } );
}

template< typename LAMBDA >
void testNumericalJacobian( ReactiveCompositionalMultiphaseOBL & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol,
                            LAMBDA assembleFunction )
{
  localIndex const NC = solver.numFluidComponents();

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

  string const dofKey = dofManager.getKey( ReactiveCompositionalMultiphaseOBL::viewKeyStruct::elemDofFieldString() );

  solver.forDiscretizationOnMeshTargets ( domain.getMeshBodies(),
                                          [&]( string const,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions( regionNames,
                                               [&]( localIndex const,
                                                    ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 > const & pres =
        subRegion.getField< fields::flow::pressure >();
      pres.move( hostMemorySpace, false );

      arrayView2d< real64, compflow::USD_COMP > const & compFrac =
        subRegion.getField< fields::flow::globalCompFraction >();
      compFrac.move( hostMemorySpace, false );

      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( elemGhostRank[ei] >= 0 )
        {
          continue;
        }


        {
          solver.resetStateToBeginningOfStep( domain );

          pres.move( hostMemorySpace, true );
          real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
          pres[ei] += dP;
#if defined(GEOS_USE_CUDA)
          pres.move( parallelDeviceMemorySpace, false );
#endif

          solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                                 [&]( string const,
                                                      MeshLevel & mesh2,
                                                      arrayView1d< string const > const & regionNames2 )
          {
            ElementRegionManager & elementRegionManager2 = mesh2.getElemManager();
            elementRegionManager2.forElementSubRegions( regionNames2,
                                                        [&]( localIndex const,
                                                             ElementSubRegionBase & subRegion2 )
            {
              solver.updateOBLOperators( subRegion2 );
            } );
          } );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 dofNumber[ei],
                                 dP,
                                 jacobianFD.toViewConstSizes() );
        }

        for( localIndex jc = 0; jc < NC - 1; ++jc )
        {
          solver.resetStateToBeginningOfStep( domain );

          compFrac.move( hostMemorySpace, true );
          compFrac[ei][jc] += perturbParameter;
#if defined(GEOS_USE_CUDA)
          compFrac.move( parallelDeviceMemorySpace, false );
#endif


          solver.forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                                 [&]( string const,
                                                      MeshLevel & mesh2,
                                                      arrayView1d< string const > const & regionNames2 )
          {
            ElementRegionManager & elementRegionManager2 = mesh2.getElemManager();
            elementRegionManager2.forElementSubRegions( regionNames2,
                                                        [&]( localIndex const,
                                                             ElementSubRegionBase & subRegion2 )
            {
              solver.updateOBLOperators( subRegion2 );
            } );
          } );

          residual.zero();
          jacobian.zero();
          assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

          fillNumericalJacobian( residual.toViewConst(),
                                 residualOrig.toViewConst(),
                                 dofNumber[ei] + jc + 1,
                                 perturbParameter,
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

class CompositionalMultiphaseFlowTest : public ::testing::Test
{
public:

  CompositionalMultiphaseFlowTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    writeTableToFile( "obl_3comp_static.txt", oblInput );
    setupProblemFromXML( state.getProblemManager(), xmlInput );
    removeFile( "obl_3comp_static.txt" );

    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< ReactiveCompositionalMultiphaseOBL >( "compflow" );

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
  ReactiveCompositionalMultiphaseOBL * solver;
};

real64 constexpr CompositionalMultiphaseFlowTest::time;
real64 constexpr CompositionalMultiphaseFlowTest::dt;
real64 constexpr CompositionalMultiphaseFlowTest::eps;

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_operators )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-4;

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testOperatorsNumericalDerivatives( *solver, domain, perturb, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-2;   // 1% error margin

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
TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_accumulation )
{
  real64 const perturb = sqrt( eps );
  real64 const tol = 1e-2;   // 1% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleAccumulationTerms( dt, domain, solver->getDofManager(), localMatrix, localRhs );
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
