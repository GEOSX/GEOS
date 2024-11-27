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

// using some utility classes from the following unit test
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverBase.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/secondOrderEqn/isotropic/AcousticWaveEquationSEM.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// This unit test checks the interpolation done to extract seismic traces from a wavefield.
// It computes a seismogram at a receiver co-located with the source and compares it to the surrounding receivers.
char const * xmlInput =
  R"xml(
  <Problem>
    <Solvers>
      <AcousticSEM
        name="acousticSolver"
        cflFactor="0.25"
        discretization="FE1"
        targetRegions="{ Region }"
        sourceCoordinates="{ { 200, 500, 600 } }"
        timeSourceFrequency="1"
        receiverCoordinates="{ { 600.1, 700.1, 700.1 }, { 200.1, 500.1, 600.1 }, { 400.1, 650.1, 700.1 }, { 601.1, 701.1, 701.1 } }"
        outputSeismoTrace="0"
        dtSeismoTrace="0.005"/>
    </Solvers>
    <Mesh>
      <InternalMesh
        name="mesh"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 1000 }"
        yCoords="{ 0, 1200 }"
        zCoords="{ 0, 1500 }"
        nx="{ 40 }"
        ny="{ 50 }"
        nz="{ 60 }"
        cellBlockNames="{ cb }"/>
    </Mesh>
    <Events
      maxTime="2.5">
      <PeriodicEvent
        name="solverApplications"
        forceDt="0.005"
        targetExactStartStop="0"
        targetExactTimestep="0"
        target="/Solvers/acousticSolver"/>
      <PeriodicEvent
        name="waveFieldNp1Collection"
        timeFrequency="1"
        targetExactTimestep="0"
        target="/Tasks/waveFieldNp1Collection" />
      <PeriodicEvent
        name="waveFieldNCollection"
        timeFrequency="1"
        targetExactTimestep="0"
        target="/Tasks/waveFieldNCollection" />
      <PeriodicEvent
        name="waveFieldNm1Collection"
        timeFrequency="1"
        targetExactTimestep="0"
        target="/Tasks/waveFieldNm1Collection" />
      <PeriodicEvent
        name="vtk"
        timeFrequency="0.5"
        targetExactTimestep="0"
        target="/Outputs/vtkOutput"/>
    </Events>
    <NumericalMethods>
      <FiniteElements>
        <FiniteElementSpace
          name="FE1"
          order="1"
          formulation="SEM"/>
      </FiniteElements>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion
        name="Region"
        cellBlocks="{ cb }"
        materialList="{ nullModel }"/>
    </ElementRegions>
    <Constitutive>
      <NullModel
        name="nullModel"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification
        name="initialPressureN"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_n"
        scale="0.0"/>
      <FieldSpecification
        name="initialPressureNm1"
        initialCondition="1"
        setNames="{ all }"
        objectPath="nodeManager"
        fieldName="pressure_nm1"
        scale="0.0"/>
      <FieldSpecification
        name="cellVelocity"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="acousticVelocity"
        scale="2000"
        setNames="{ all }"/>
      <FieldSpecification
        name="cellDensity"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="acousticDensity"
        scale="2"
        setNames="{ all }"/>
      <FieldSpecification
        name="zposFreeSurface"
        objectPath="faceManager"
        fieldName="FreeSurface"
        scale="0.0"
        setNames="{ zpos }"/>
    </FieldSpecifications>
    <Tasks>
      <PackCollection
        name="waveFieldNp1Collection"
        objectPath="nodeManager"
        fieldName="pressure_np1"/>
      <PackCollection
        name="waveFieldNCollection"
        objectPath="nodeManager"
        fieldName="pressure_n"/>
      <PackCollection
        name="waveFieldNm1Collection"
        objectPath="nodeManager"
        fieldName="pressure_nm1"/>
    </Tasks>
    <Outputs>
      <VTK                                                                        
        name="vtkOutput"                                                          
        levelNames="{ FE1 }"                                                      
        plotLevel="3"/>                                                           
      <Restart                                                                    
        name="restartOutput"/>                                                    
    </Outputs>
  </Problem>
  )xml";

class AcousticWaveEquationSEMTest : public ::testing::Test
{
public:

  AcousticWaveEquationSEMTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 5e-3;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  AcousticWaveEquationSEM * propagator;
};

real64 constexpr AcousticWaveEquationSEMTest::time;
real64 constexpr AcousticWaveEquationSEMTest::dt;
real64 constexpr AcousticWaveEquationSEMTest::eps;

TEST_F( AcousticWaveEquationSEMTest, SeismoTrace )
{

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticWaveEquationSEM >( "acousticSolver" );

  // Check source term (sourceCoordinates and sourceValue)
  arrayView2d< real32 > const sourceValue = propagator->getReference< array2d< real32 > >( AcousticWaveEquationSEM::viewKeyStruct::sourceValueString() ).toView();
  // move it to CPU, if needed
  sourceValue.move( hostMemorySpace, false );
  /*for( int i = 0; i < 501; i++ )
  {
    std::cout << "Before explicit source value :" << sourceValue[i][0] << std::endl;
  }*/

  real64 time_n = time;
  std::cout << "Begin forward:" << time_n << std::endl;
  // run for 2.5s (500 steps)
  for( int i=0; i<500; i++ )
  {
    propagator->explicitStepForward( time_n, dt, i, domain, false );
    time_n += dt;
  }
  // cleanup (triggers calculation of the remaining seismograms data points)
  propagator->cleanup( 1.0, 500, 0, 0, domain );

  // retrieve seismo
  arrayView2d< real32 > const pReceivers = propagator->getReference< array2d< real32 > >( AcousticWaveEquationSEM::viewKeyStruct::pressureNp1AtReceiversString() ).toView();

  // move it to CPU, if needed
  pReceivers.move( hostMemorySpace, false );

  // check number of seismos and trace length
  ASSERT_EQ( pReceivers.size( 1 ), 5 );
  ASSERT_EQ( pReceivers.size( 0 ), 501 );

  /*----------Save receiver forward----------------------*/
  array2d< real32 > u_forward;
  u_forward.resize(501,1);
  array2d< real32 > rhs_forward;
  rhs_forward.resize(501,1);
  
  // save receiver value forward on u_forward.
  for( int i = 0; i < 501; i++ )
  {
    /*std::cout << "time: " << i*dt  << std::endl;
    std::cout << "pReceivers1 " << i << ":" << pReceivers[i][0] << std::endl;
    std::cout << "pReceivers2 " << i << ":" << pReceivers[i][1] << std::endl;
    std::cout << "pReceivers3  " << i << ":" << pReceivers[i][2] << std::endl;
    std::cout << "pReceivers4  " << i << ":" << pReceivers[i][3] << std::endl;*/
    u_forward[i][0] = pReceivers[i][0];
    pReceivers[i][0] = 0.;
    pReceivers[i][1] = 0.;
    pReceivers[i][2] = 0.;
    pReceivers[i][3] = 0.;
    rhs_forward[i][0] = sourceValue[i][0];
  }
  /*---------------------------------------------------*/

  std::cout << "Begin backward:" << time_n << std::endl;

  // modify source term
  /*std::cout << "Source value size(0):" << sourceValue.size(0) << std::endl;
  std::cout << "Source value size(1):" << sourceValue.size(1) << std::endl;
  for( int i = 0; i < sourceValue.size(0); i++ )
  {
    std::cout << "Source value " << i << ":" << sourceValue[i][0] << std::endl;
    std::cout << "u_forward " << i << ":" << u_forward[i][0] << std::endl;
  }*/
  ASSERT_EQ( sourceValue.size( 1 ), 1 );
  ASSERT_EQ( sourceValue.size( 0 ), 501 );

  //----------Switch source and receiver1 position for backward----------------------//
  //modify source coordinates (get receiver coordinate)
  // retrieve seismo
  arrayView2d< localIndex > const sNodeIds = propagator->getReference< array2d< localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::sourceNodeIdsString() ).toView();
  arrayView2d< real64 > const sConstants   = propagator->getReference< array2d< real64 > >( AcousticWaveEquationSEM::viewKeyStruct::sourceConstantsString() ).toView();
  arrayView1d< localIndex > const sBool = propagator->getReference< array1d<  localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::sourceIsAccessibleString() ).toView();
  arrayView2d< real64 > const sCoord = propagator->getReference< array2d< real64 > >( AcousticWaveEquationSEM::viewKeyStruct::sourceCoordinatesString() ).toView();
  arrayView2d< localIndex > const rNodeIds = propagator->getReference< array2d<  localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::receiverNodeIdsString() ).toView();
  arrayView2d< real64 > const rConstants   = propagator->getReference< array2d< real64 > >( AcousticWaveEquationSEM::viewKeyStruct::receiverConstantsString() ).toView();
  arrayView1d< localIndex > const rBool = propagator->getReference< array1d<  localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::receiverIsLocalString() ).toView();
  arrayView2d< real64 > const rCoord = propagator->getReference< array2d< real64 > >( AcousticWaveEquationSEM::viewKeyStruct::receiverCoordinatesString() ).toView();

  // move it to CPU, if needed
  sCoord.move( hostMemorySpace, false );
  sNodeIds.move( hostMemorySpace, false );
  sConstants.move( hostMemorySpace, false );
  sBool.move( hostMemorySpace, false );
  rCoord.move( hostMemorySpace, false );
  rNodeIds.move( hostMemorySpace, false );
  rConstants.move( hostMemorySpace, false );
  rBool.move( hostMemorySpace, false );

  //modify receiver coordinate (get state at original source coordinate)
  for( int i = 0; i < sNodeIds.size( 1 ); i++ )
  {
    localIndex tmp_ind;
    tmp_ind=rNodeIds[0][i];
    rNodeIds[0][i]=sNodeIds[0][i];
    sNodeIds[0][i]=tmp_ind;
    real64 tmp_double;
    tmp_double=rConstants[0][i];
    rConstants[0][i]=sConstants[0][i];
    sConstants[0][i]=tmp_double;
  }  
  sBool[0]=sBool[0]+rBool[0];
  rBool[0]=sBool[0]-rBool[0];
  sBool[0]=sBool[0]-rBool[0];
  for( int i = 0; i < 3; i++ )
  {
    real64 tmp_double;
    tmp_double=rCoord[0][i];
    rCoord[0][i]=sCoord[0][i];
    sCoord[0][i]=tmp_double;
  }  

  /*std::cout << "sCoord :" << sCoord[0][0] <<" "<< sCoord[0][1] <<" "<< sCoord[0][2] << std::endl;
  std::cout << "rCoord1 :" << rCoord[0][0] <<" "<< rCoord[0][1] <<" "<< rCoord[0][2] << std::endl;
  std::cout << "rCoord2 :" << rCoord[1][0] <<" "<< rCoord[1][1] <<" "<< rCoord[1][2] << std::endl;
  std::cout << "rCoord3 :" << rCoord[2][0] <<" "<< rCoord[2][1] <<" "<< rCoord[2][2] << std::endl;
  std::cout << "rCoord4 :" << rCoord[3][0] <<" "<< rCoord[3][1] <<" "<< rCoord[3][2] << std::endl;
  std::cout << "test sNodeIds :" << sNodeIds[0][0] <<" "<< sNodeIds[0][4] << std::endl;
  std::cout << "test rNodeIds1 :" << rNodeIds[0][0] <<" "<< rNodeIds[0][4] << std::endl;
  std::cout << "test rNodeIds2 :" << rNodeIds[1][0] <<" "<< rNodeIds[1][4] << std::endl;*/

  // move to GPU to be used by backward
  sNodeIds.move( parallelDeviceMemorySpace, false );
  sConstants.move( parallelDeviceMemorySpace, false );
  sBool.move( parallelDeviceMemorySpace, false );
  rNodeIds.move( parallelDeviceMemorySpace, false );
  rConstants.move( parallelDeviceMemorySpace, false );
  rBool.move( parallelDeviceMemorySpace, false );

  //reinit indexSeismoTrace
  localIndex * iSeismo_ptr = &propagator->getReference< localIndex >( AcousticWaveEquationSEM::viewKeyStruct::indexSeismoTraceString() );
  *iSeismo_ptr = pReceivers.size(0)-1;
  //reinit m_forward
  localIndex * forward_ptr = &propagator->getReference< localIndex >( AcousticWaveEquationSEM::viewKeyStruct::forwardString() );
  *forward_ptr = 0;
  
  state.getProblemManager().applyInitialConditions();

  /* //  NOT WORKING change source value amplitude (set 0.5*sourceValue[i])
  //sourceValue = propagator->getReference< array2d< real32 > >( AcousticWaveEquationSEM::viewKeyStruct::sourceValueString() ).toView();
  sourceValue.move( hostMemorySpace, false );
  for( int i = 0; i < sourceValue.size(0); i++ )
  {
    sourceValue[i][0] = 0.5 * sourceValue[i][0];
    std::cout << "Source value after init" << i << ":" << sourceValue[i][0] << std::endl;
  }
  sourceValue.move( parallelDeviceMemorySpace, false );*/

  /*---------------------------------------------------*/
  // run backward solver
  for( int i = 500; i > 0; i-- )
  {
    propagator->explicitStepBackward( time_n, dt, i, domain, false );
    time_n -= dt;
  }

  // move it to CPU, if needed
  pReceivers.move( hostMemorySpace, false );
  sourceValue.move( hostMemorySpace, false );

  std::cout << "pReceiver size(0):" << pReceivers.size(0) << std::endl;
  std::cout << "pReceiver size(1):" << pReceivers.size(1) << std::endl;


  /*----------Save receiver forward----------------------*/
  array2d< real32 > q_backward;
  q_backward.resize(501,1);

  real32 sum_ufb=0.;
  real32 sum_qff=0.;
  real32 sum_u2=0.;
  real32 sum_q2=0.;
  real32 sum_ff2=0.;
  real32 sum_fb2=0.;
  
  // fill backward field at receiver.
  for( int i=0; i<501; i++ )
  {
    /*std::cout << "back Source value " << i << ":" << sourceValue[i][0] << std::endl;
    std::cout << "back pReceivers1 " << i << ":" << pReceivers[i][0] << std::endl;
    std::cout << "back pReceivers2 " << i << ":" << pReceivers[i][1] << std::endl;
    std::cout << "back pReceivers3 " << i << ":" << pReceivers[i][2] << std::endl;
    std::cout << "back pReceivers4 " << i << ":" << pReceivers[i][3] << std::endl;*/
    q_backward[i][0] = pReceivers[i][0];
  }

  //check transitivity with sum
  for( int i=0; i<501; i++ )
  {
    sum_ufb += u_forward[i][0]*sourceValue[i][0];
    sum_qff += q_backward[i][0]*rhs_forward[i][0];

    sum_u2 += u_forward[i][0]*u_forward[i][0];
    sum_q2 += q_backward[i][0]*q_backward[i][0];
    sum_ff2 += rhs_forward[i][0]*rhs_forward[i][0];
    sum_fb2 += sourceValue[i][0]*sourceValue[i][0];
    /*std::cout << "sum evol sum_ufb:" << sum_ufb << " / sum_qff:" << sum_qff << std::endl;
    std::cout << "u_forward:" << u_forward[i][0] << " / q_backward:" << q_backward[i][0] << std::endl;
    std::cout << "ufb:" << u_forward[i][0]*sourceValue[i][0] << " / qff:" << q_backward[i][0]*rhs_forward[i][0] << std::endl;*/
  }

  // check ||<f,q> - <u,f'>||/max(||f||.||q||,||f'||.||u||) < 10^2or3 . epsilon_machine with f rhs direct and f' rhs backward
  std::cout << "<u,f'>:" << sum_ufb << " / <f,q>:" << sum_qff << std::endl;
  std::cout << "||<f,q> - <u,f'>||=" << std::abs(sum_ufb-sum_qff) << " / ||f||.||q||=" << std::sqrt(sum_q2*sum_ff2);  
  std::cout << " / ||f'||.||u||=" << std::sqrt(sum_fb2*sum_u2) << " / ||f||.||f'||=" << std::sqrt(sum_ff2*sum_fb2) << std::endl;
  real32 diffToCheck;
  diffToCheck=std::abs( sum_ufb-sum_qff ) / std::max( std::sqrt(sum_fb2*sum_ff2),std::sqrt(sum_q2*sum_ff2));
  std::cout << " Diff to compare with 100*1.e-8:" << diffToCheck << std::endl;
  ASSERT_TRUE( diffToCheck < 100*1.e-8 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
