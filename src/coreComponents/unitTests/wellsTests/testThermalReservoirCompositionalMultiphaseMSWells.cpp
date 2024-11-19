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
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::constitutive;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

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
char const * co2flash = "FlashModel CO2Solubility  1e6 7.5e7 5e5 299.15 369.15 10 0";
char const * pvtLiquid = "DensityFun PhillipsBrineDensity 1e6 7.5e7 5e5 299.15 369.15 10 0\n"
                         "ViscosityFun PhillipsBrineViscosity 0\n"
                         "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 299.15 369.15 10 0\n";

char const * pvtGas = "DensityFun SpanWagnerCO2Density 1e6 7.5e7 5e5 299.15 369.15 10\n"
                      "ViscosityFun FenghourCO2Viscosity 1e6 7.5e7 5e5 299.15 369.15 10\n"
                      "EnthalpyFun CO2Enthalpy 1e6 7.5e7 5e5 299.15 369.15 10\n";
char const * xmlInput =
  R"xml(

<?xml version="1.0" ?>
<Problem>
  <Solvers>
    <CompositionalMultiphaseReservoir
      name="reservoirSystem"
      flowSolverName="compflow"
      wellSolverName="compositionalMultiphaseWell"
      logLevel="4"
      initialDt="1e4"
      targetRegions="{ region, injwell }">
      <NonlinearSolverParameters
         logLevel="0"
        newtonTol="1.0e-5"
        lineSearchAction="None"
        newtonMaxIter="40"/>
      <LinearSolverParameters
        logLevel="4"
        directParallel="0"/>
    </CompositionalMultiphaseReservoir>

    <CompositionalMultiphaseFVM
      name="compflow"
      logLevel="4"
      discretization="fluidTPFA"
      temperature="368.15"
      useMass="1"
      isThermal="1"
      initialDt="1e4"
      targetRelativePressureChangeInTimeStep="1"
      targetRelativeTemperatureChangeInTimeStep="1"
      targetPhaseVolFractionChangeInTimeStep="1"      
      maxCompFractionChange="0.5"
      targetRegions="{ region }">
    </CompositionalMultiphaseFVM>

    <CompositionalMultiphaseWell
      name="compositionalMultiphaseWell"
      targetRegions="{ injwell }"
      isThermal="1"
      logLevel="1"
      initialDt="1e4"
      useMass="1">
      <WellControls
        name="WC_CO2_INJ"
        logLevel="2"
        type="injector"
        control="totalVolRate"
        referenceElevation="-0.01"
        targetBHP="1.45e7"
        enableCrossflow="0"
        useSurfaceConditions="1"
        surfacePressure="1.45e7"
        surfaceTemperature="300.15"
        targetTotalRate="0.001"
        injectionTemperature="300.15"
        injectionStream="{ 0.99, 0.01 }"/>
     </CompositionalMultiphaseWell>
  </Solvers>

  <Mesh>
    <InternalMesh
      name="mesh1"
      elementTypes="{ C3D8 }"
      xCoords="{ 0, 100 }"
      yCoords="{ 0, 100 }"
      zCoords="{ 0, 1 }"
      nx="{ 2 }"
      ny="{ 2 }"
      nz="{ 1 }"
      cellBlockNames="{ cb }">
      <InternalWell
        name="inj1"
        wellRegionName="injwell"
        wellControlsName="WC_CO2_INJ"
        logLevel="1"
        polylineNodeCoords="{ { 5.0, 5.0, 1.01 },
                              { 5.0, 5.0, -0.01 } }"
        polylineSegmentConn="{ { 0, 1 } }"
        radius="0.1"
        numElementsPerSegment="2">
        <Perforation
          name="injector1_perf1"
          distanceFromHead="0.75"/>

      </InternalWell>
    </InternalMesh>

  </Mesh>

  <Geometry>
    <Box
      name="sink"
      xMin="{ 89.99, 89.99, -0.01 }"
      xMax="{ 101.01, 101.01, 1.01 }"/>

    
  </Geometry>

  <Events
    maxTime="1.5e5">
    <PeriodicEvent
      name="outputs"
      timeFrequency="2.5e4"
      target="/Outputs/vtkOutput"/>

    <PeriodicEvent
      name="solverApplications"
      maxEventDt="2.5e4"
      target="/Solvers/coupledFlowAndWells"/>

    <PeriodicEvent
      name="restarts"
      timeFrequency="7.5e5"
      target="/Outputs/sidreRestart"/>
  </Events>

  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation
        name="fluidTPFA"/>
    </FiniteVolume>
  </NumericalMethods>

  <ElementRegions>
    <CellElementRegion
      name="region"
      cellBlocks="{ cb }"
      materialList="{ fluid, rock, relperm, thermalCond, diffusion }"/>
    <WellElementRegion
      name="injwell"
      materialList="{ fluid, relperm }"/>
  </ElementRegions>

  <Constitutive>

    <CompressibleSolidConstantPermeability
      name="rock"
      solidModelName="nullSolid"
      porosityModelName="rockPorosity"
      permeabilityModelName="rockPerm"
      solidInternalEnergyModelName="rockInternalEnergy"/>
    <NullModel
      name="nullSolid"/>
    <PressurePorosity
      name="rockPorosity"
      defaultReferencePorosity="0.2"
      referencePressure="0.0"
      compressibility="1.0e-9"/>
    <SolidInternalEnergy
      name="rockInternalEnergy"
      referenceVolumetricHeatCapacity="1.95e6"
      referenceTemperature="368.15"
      referenceInternalEnergy="0"/>
    <ConstantPermeability
      name="rockPerm"
      permeabilityComponents="{ 1.0e-13, 1.0e-13, 1.0e-13 }"/>

    <CO2BrinePhillipsThermalFluid
      name="fluid"
      logLevel="1"
      phaseNames="{ gas, water }"
      componentNames="{ co2, water }"
      componentMolarWeight="{ 44e-3, 18e-3 }"
      phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
      flashModelParaFile="co2flash.txt"/>

    <BrooksCoreyRelativePermeability
      name="relperm"
      phaseNames="{ gas, water }"
      phaseMinVolumeFraction="{ 0.0, 0.0 }"
      phaseRelPermExponent="{ 1.5, 1.5 }"
      phaseRelPermMaxValue="{ 0.9, 0.9 }"/>

    <MultiPhaseConstantThermalConductivity
      name="thermalCond"
      phaseNames="{ gas, water }"
      thermalConductivityComponents="{ 0.6, 0.6, 0.6 }"/>

    <ConstantDiffusion
      name="diffusion"
      phaseNames="{ gas, water }"
      defaultPhaseDiffusivityMultipliers="{ 20, 1 }"
      diffusivityComponents="{ 1e-9, 1e-9, 1e-9 }"/>
    
  </Constitutive>

  <FieldSpecifications>

    <FieldSpecification
      name="initialPressure"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="pressure"
      scale="9e6"/>
    <FieldSpecification
      name="initialTemperature"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="temperature"
      scale="368.15"/>
    <FieldSpecification
      name="initialComposition_co2"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="0"
      scale="0.005"/>
    <FieldSpecification
      name="initialComposition_water"
      initialCondition="1"
      setNames="{ all }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="1"
      scale="0.995"/>

    <FieldSpecification
      name="sinkPressure"
      setNames="{ sink }"       
      objectPath="ElementRegions/region/cb"
      fieldName="pressure"
      scale="7e6"/>
    <FieldSpecification
      name="sinkTemperature"
      setNames="{ sink }"       
      objectPath="ElementRegions/region/cb"
      fieldName="temperature"
      scale="368.15"/>
     <FieldSpecification
      name="sinkTermComposition_co2"
      setNames="{ sink }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="0"
      scale="0.005"/>
    <FieldSpecification
      name="sinkTermComposition_water"
      setNames="{ sink }"
      objectPath="ElementRegions/region/cb"
      fieldName="globalCompFraction"
      component="1"
      scale="0.995"/>

  </FieldSpecifications>

  <Outputs>
    <VTK
      name="vtkOutput"/>

    <Restart
      name="sidreRestart"/>
  </Outputs>
</Problem>
  )xml";



template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > & solver,
                            DomainPartition & domain,
                            real64 const perturbParameter,
                            real64 const relTol, bool diag_check,
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
  if( 1 )
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
              subRegion.getField< fields::flow::pressure >();
            pres.move( hostMemorySpace, false );

            arrayView2d< real64, compflow::USD_COMP > const & compDens =
              subRegion.getField< fields::flow::globalCompDensity >();
            compDens.move( hostMemorySpace, false );

            arrayView1d< real64 > const & temp =
              subRegion.getField< fields::flow::temperature >();
            temp.move( hostMemorySpace, false );

            // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei
            for( localIndex ei = 0; ei < subRegion.size(); ++ei )
            {
              if( ei !=0 )
                break;
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
                wellSolver.updateState( domain );
                residual.zero();
                jacobian.zero();
                assembleFunction( jacobian.toViewConstSizes(), residual.toView() );

                fillNumericalJacobian( residual.toViewConst(),
                                       residualOrig.toViewConst(),
                                       dofNumber[ei] + jc + 1,
                                       dRho,
                                       jacobianFD.toViewConstSizes() );
              }
              {
                solver.resetStateToBeginningOfStep( domain );

                // here is the perturbation in the pressure of the element
                real64 const dTemp = perturbParameter * (temp[ei] + perturbParameter);
                temp.move( hostMemorySpace, true );
                temp[ei] += dTemp;

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
                                       dofNumber[ei]+NC+1,
                                       dTemp,
                                       jacobianFD.toViewConstSizes() );
              }
            }
          } );
        }
      } );
      return;
    } );

  /////////////////////////////////////////////////
  // Step 2) Compute the terms in J_RW and J_WW //
  /////////////////////////////////////////////////

  // loop over the wells
  if( 1 )
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

        arrayView1d< real64 > const & wellElemTemperature =
          subRegion.getField< fields::well::temperature >();
        wellElemTemperature.move( hostMemorySpace, false );

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
          {
            solver.resetStateToBeginningOfStep( domain );
            residual.zero();
            jacobian.zero();
            if( diag_check || iwelem > 0 )
            {
              // here is the perturbation in the temperature of the well element
              real64 const dT = perturbParameter * ( wellElemTemperature[iwelem] + perturbParameter );
              wellElemTemperature.move( hostMemorySpace, true );
              wellElemTemperature[iwelem] += dT;

              // after perturbing, update the pressure-dependent quantities in the well
              wellSolver.updateState( domain );


              assembleFunction( jacobian.toViewConstSizes(), residual.toView() );
              fillNumericalJacobian( residual.toViewConst(),
                                     residualOrig.toViewConst(),
                                     wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1,
                                     dT,
                                     jacobianFD.toViewConstSizes() );
              if( iwelem == 1 )
              {
                real64 dRdX = 0.0;
                localIndex rowIndex = wellElemDofNumber[0] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1;;
                for( integer ider=0; ider< 3; ider++ )
                {
                  globalIndex colIndex = wellElemDofNumber[0]+ ider;
                  setNumericalJacobianValue( rowIndex, colIndex, dRdX, jacobianFD.toViewConstSizes() );
                }
                globalIndex colIndex = wellElemDofNumber[1]+3;
                setNumericalJacobianValue( rowIndex, colIndex, dRdX, jacobianFD.toViewConstSizes() );
              }
            }
            else
            {
              localIndex rowIndex = wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1;;
              globalIndex colIndex = wellElemDofNumber[iwelem] + compositionalMultiphaseWellKernels::ColOffset::DCOMP + NC+1;;
              setNumericalJacobianValue( rowIndex, colIndex, 1.0, jacobianFD.toViewConstSizes() );
            }

          }
        }


        // b) compute all the derivatives wrt to the connection in WELL elem
        // iwelem
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
  //printCompareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst());
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
    solver = &state.getProblemManager().getPhysicsSolverManager().getGroup< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > >( "reservoirSystem" );

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
  CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > * solver;
};

real64 constexpr CompositionalMultiphaseReservoirSolverTest::time;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::dt;
real64 constexpr CompositionalMultiphaseReservoirSolverTest::eps;

#if 0

TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Perforation )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain =  state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->assembleCouplingTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
    //solver->assembleSystem( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#endif

#if 0
TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol, false,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assembleSystem( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
    //solver->assembleSystem( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
    //solver->assembleCouplingTerms( time,dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}

#endif
#if 0

TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_flux )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol, true,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->computePerforationRates( time, dt, domain );
    solver->wellSolver()->assembleFluxTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
    solver->assembleCouplingTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}
#endif
#if 1

TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_Accum_Vol_Energy_Bal )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol, false,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assembleAccumulationTerms( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}
TEST_F( CompositionalMultiphaseReservoirSolverTest, jacobianNumericalCheck_PressureRel )
{
  real64 const perturb = std::sqrt( eps );
  real64 const tol = 1e-1; // 10% error margin

  DomainPartition & domain = state.getProblemManager().getDomainPartition();

  testNumericalJacobian( *solver, domain, perturb, tol, true,
                         [&] ( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
  {
    solver->wellSolver()->assemblePressureRelations( time, dt, domain, solver->getDofManager(), localMatrix, localRhs );
  } );
}
#endif
int main( int argc, char * * argv )
{
  writeTableToFile( "co2flash.txt", co2flash );
  writeTableToFile( "pvtliquid.txt", pvtLiquid );
  writeTableToFile( "pvtgas.txt", pvtGas );
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  removeFile( "co2flash.txt" );
  removeFile( "pvtliquid.txt" );
  removeFile( "pvtgas.txt" );

  return result;
}
