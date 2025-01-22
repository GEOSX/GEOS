/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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
        sourceCoordinates="{ { 100, 250, 300 } }"
        timeSourceFrequency="1"
        receiverCoordinates="{ { 300.1, 350.1, 350.1 }, { 100.1, 250.1, 300.1 }, { 200.1, 325.1, 350.1 }, { 301.1, 351.1, 351.1 } }"
        outputSeismoTrace="0"
        dtSeismoTrace="0.005"/>
    </Solvers>
    <Mesh>
      <InternalMesh
        name="mesh"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 500 }"
        yCoords="{ 0, 600 }"
        zCoords="{ 0, 750 }"
        nx="{ 20 }"
        ny="{ 25 }"
        nz="{ 30 }"
        cellBlockNames="{ cb }"/>
    </Mesh>
    <Events
      maxTime=".25">
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
  array2d< real32 > rhsForward;
  rhsForward.resize( 51, 1 );
  real32 * ptrTimeSourceFrequency = &propagator->getReference< real32 >( AcousticWaveEquationSEM::viewKeyStruct::timeSourceFrequencyString() );
  real32 * ptrTimeSourceDelay = &propagator->getReference< real32 >( AcousticWaveEquationSEM::viewKeyStruct::timeSourceDelayString() );
  localIndex * ptrRickerOrder = &propagator->getReference< localIndex >( AcousticWaveEquationSEM::viewKeyStruct::rickerOrderString() );

  real64 time_n = time;
  std::cout << "Begin forward:" << time_n << std::endl;
  // run for 0.25s (100 steps)
  for( int i=0; i<50; i++ )
  {
    rhsForward[i][0]=WaveSolverUtils::evaluateRicker( time_n, *ptrTimeSourceFrequency, *ptrTimeSourceDelay, *ptrRickerOrder );
    propagator->explicitStepForward( time_n, dt, i, domain, false );
    time_n += dt;
  }
  // cleanup (triggers calculation of the remaining seismograms data points)
  propagator->cleanup( 1.0, 50, 0, 0, domain );

  // retrieve seismo
  arrayView2d< real32 > const pReceivers = propagator->getReference< array2d< real32 > >( AcousticWaveEquationSEM::viewKeyStruct::pressureNp1AtReceiversString() ).toView();

  // move it to CPU, if needed
  pReceivers.move( hostMemorySpace, false );

  // check number of seismos and trace length
  ASSERT_EQ( pReceivers.size( 1 ), 5 );
  ASSERT_EQ( pReceivers.size( 0 ), 51 );

  /*----------Save receiver forward----------------------*/
  array2d< real32 > uForward;
  uForward.resize( 51, 1 );

  // save receiver value forward on uForward.
  for( int i = 0; i < 51; i++ )
  {
    /*std::cout << "time: " << i*dt  << std::endl;
       std::cout << "pReceivers1 " << i << ":" << pReceivers[i][0] << std::endl;
       std::cout << "pReceivers2 " << i << ":" << pReceivers[i][1] << std::endl;
       std::cout << "pReceivers3 " << i << ":" << pReceivers[i][2] << std::endl;
       std::cout << "pReceivers4 " << i << ":" << pReceivers[i][3] << std::endl;
       std::cout << "rhsForward  " << i << ":" << rhsForward[i][0] << std::endl;*/
    uForward[i][0] = pReceivers[i][0];
    pReceivers[i][0] = 0.;
    pReceivers[i][1] = 0.;
    pReceivers[i][2] = 0.;
    pReceivers[i][3] = 0.;
  }

  ASSERT_EQ( rhsForward.size( 1 ), 1 );
  ASSERT_EQ( rhsForward.size( 0 ), 51 );

  arrayView2d< localIndex > const rNodeIds = propagator->getReference< array2d< localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::receiverNodeIdsString() ).toView();
  rNodeIds.move( hostMemorySpace, false );
  localIndex sNodesIdsAfterModif=rNodeIds[0][0];
  std::cout << "ref back sNodeIds[0][0]:" << sNodesIdsAfterModif << std::endl;

  /*---------------------------------------------------*/

  std::cout << "Begin backward:" << time_n << std::endl;

  //----------Switch source and receiver1 position for backward----------------------//
  arrayView2d< real64 > const sCoord = propagator->getReference< array2d< real64 > >( AcousticWaveEquationSEM::viewKeyStruct::sourceCoordinatesString() ).toView();
  arrayView2d< real64 > const rCoord = propagator->getReference< array2d< real64 > >( AcousticWaveEquationSEM::viewKeyStruct::receiverCoordinatesString() ).toView();

  for( int i = 0; i < 3; i++ )
  {
    real64 tmp_double;
    tmp_double=rCoord[0][i];
    rCoord[0][i]=sCoord[0][i];
    sCoord[0][i]=tmp_double;
  }

  sCoord.registerTouch( hostMemorySpace );
  rCoord.registerTouch( hostMemorySpace );

  std::cout << "sCoord  :" << sCoord[0][0] <<" "<< sCoord[0][1] <<" "<< sCoord[0][2] << std::endl;
  std::cout << "rCoord1 :" << rCoord[0][0] <<" "<< rCoord[0][1] <<" "<< rCoord[0][2] << std::endl;
  std::cout << "rCoord2 :" << rCoord[1][0] <<" "<< rCoord[1][1] <<" "<< rCoord[1][2] << std::endl;
  std::cout << "rCoord3 :" << rCoord[2][0] <<" "<< rCoord[2][1] <<" "<< rCoord[2][2] << std::endl;
  std::cout << "rCoord4 :" << rCoord[3][0] <<" "<< rCoord[3][1] <<" "<< rCoord[3][2] << std::endl;

  //change timeSourceFrequency
  std::cout << "timeSourceFrequency forward:" << *ptrTimeSourceFrequency << std::endl;
  real32 newTimeFreq=2;
  *ptrTimeSourceFrequency = newTimeFreq;
  std::cout << "timeSourceFrequency backward:" << *ptrTimeSourceFrequency << std::endl;

  //reinit m_indexSeismoTrace
  localIndex * ptrISeismo = &propagator->getReference< localIndex >( AcousticWaveEquationSEM::viewKeyStruct::indexSeismoTraceString() );
  *ptrISeismo = pReceivers.size( 0 )-1;
  //reinit m_forward
  localIndex * ptrForward = &propagator->getReference< localIndex >( AcousticWaveEquationSEM::viewKeyStruct::forwardString() );
  *ptrForward = 0;

  //"propagator->reinit()" not enough because state field not reinit to zero
  //propagator->reinit();
  state.getProblemManager().applyInitialConditions();

  array2d< real32 > rhsBackward;
  rhsBackward.resize( 51, 1 );

  arrayView2d< localIndex > const sNodeIds_new2 = propagator->getReference< array2d< localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::sourceNodeIdsString() ).toView();
  sNodeIds_new2.move( hostMemorySpace, false );
  std::cout << "sNodeIds[0][0] second get2:" << sNodeIds_new2[0][0] << std::endl;
  ASSERT_TRUE( sNodeIds_new2[0][0] == sNodesIdsAfterModif );

  /*---------------------------------------------------*/
  // run backward solver
  for( int i = 50; i > 0; i-- )
  {
    rhsBackward[i][0]=WaveSolverUtils::evaluateRicker( time_n, *ptrTimeSourceFrequency, *ptrTimeSourceDelay, *ptrRickerOrder );
    propagator->explicitStepBackward( time_n, dt, i, domain, false );
    time_n -= dt;
    //check source node in backward loop
    arrayView2d< localIndex > const sNodeIds_loop = propagator->getReference< array2d< localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::sourceNodeIdsString() ).toView();
    sNodeIds_loop.move( hostMemorySpace, false );
    ASSERT_TRUE( sNodeIds_loop[0][0] == sNodesIdsAfterModif );
  }

  // move it to CPU, if needed
  pReceivers.move( hostMemorySpace, false );

  localIndex mForward2 = propagator->getReference< localIndex >( AcousticWaveEquationSEM::viewKeyStruct::forwardString() );
  std::cout << "m_forward second get:" << mForward2 << std::endl;
  ASSERT_TRUE( mForward2 == 0 );

  arrayView2d< localIndex > const sNodeIds_new3 = propagator->getReference< array2d< localIndex > >( AcousticWaveEquationSEM::viewKeyStruct::sourceNodeIdsString() ).toView();
  sNodeIds_new3.move( hostMemorySpace, false );
  std::cout << "sNodeIds[0][0] get3:" << sNodeIds_new3[0][0] << std::endl;
  ASSERT_TRUE( sNodeIds_new3[0][0] == sNodesIdsAfterModif );

  real32 const timeSourceFrequency_new = propagator->getReference< real32 >( AcousticWaveEquationSEM::viewKeyStruct::timeSourceFrequencyString() );
  ASSERT_TRUE( std::abs( timeSourceFrequency_new - newTimeFreq ) < 1.e-8 );

  /*std::cout << "pReceiver size(0):" << pReceivers.size(0) << std::endl;
     std::cout << "pReceiver size(1):" << pReceivers.size(1) << std::endl;*/


  /*----------Save receiver backward----------------------*/
  array2d< real32 > qBackward;
  qBackward.resize( 51, 1 );

  real32 sum_ufb=0.;
  real32 sum_qff=0.;
  real32 sum_u2=0.;
  real32 sum_q2=0.;
  real32 sum_ff2=0.;
  real32 sum_fb2=0.;

  // fill backward field at receiver.
  for( int i=50; i > 0; i-- )
  {
    /*std::cout << "back time: " << i*dt  << std::endl;
       std::cout << "back pReceivers1 " << i << ":" << pReceivers[i][0] << std::endl;
       std::cout << "back pReceivers2 " << i << ":" << pReceivers[i][1] << std::endl;
       std::cout << "back pReceivers3 " << i << ":" << pReceivers[i][2] << std::endl;
       std::cout << "back pReceivers4 " << i << ":" << pReceivers[i][3] << std::endl;
       std::cout << "back rhsBackward " << i << ":" << rhsBackward[i][0] << std::endl;*/
    qBackward[i][0] = pReceivers[i][0];
  }

  //check transitivity with sum
  for( int i=0; i<51; i++ )
  {
    sum_ufb += uForward[i][0]*rhsBackward[i][0];
    sum_qff += qBackward[i][0]*rhsForward[i][0];

    sum_u2 += uForward[i][0]*uForward[i][0];
    sum_q2 += qBackward[i][0]*qBackward[i][0];
    sum_ff2 += rhsForward[i][0]*rhsForward[i][0];
    sum_fb2 += rhsBackward[i][0]*rhsBackward[i][0];
    /*std::cout << "sum evol sum_ufb:" << sum_ufb << " / sum_qff:" << sum_qff << std::endl;
       std::cout << "uForward:" << uForward[i][0] << " / qBackward:" << qBackward[i][0] << std::endl;
       std::cout << "ufb:" << uForward[i][0]*rhsBackward[i][0] << " / qff:" << qBackward[i][0]*rhsForward[i][0] << std::endl;*/
  }

  // check scalar products <u,f'> and <f,q> are non null
  ASSERT_TRUE( sum_ufb > 1.e-8 );
  ASSERT_TRUE( sum_qff > 1.e-8 );

  // check ||<f,q> - <u,f'>||/max(||f||.||q||,||f'||.||u||) < 10^1or2 x epsilon_machine with f rhs direct and f' rhs backward
  std::cout << "<u,f'>: " << sum_ufb << " / <f,q>: " << sum_qff << std::endl;
  std::cout << "||<f,q> - <u,f'>||=" << std::abs( sum_ufb-sum_qff ) << " / ||f||.||q||=" << std::sqrt( sum_q2*sum_ff2 );
  std::cout << " / ||f'||.||u||=" << std::sqrt( sum_fb2*sum_u2 ) << " / ||f||.||f'||=" << std::sqrt( sum_ff2*sum_fb2 ) << std::endl;
  real32 diffToCheck;
  diffToCheck=std::abs( sum_ufb-sum_qff ) / std::max( std::sqrt( sum_fb2*sum_u2 ), std::sqrt( sum_q2*sum_ff2 ));
  std::cout << " Diff to compare with 9e-3: " << diffToCheck << std::endl;
  ASSERT_TRUE( diffToCheck < 9e-3);
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
