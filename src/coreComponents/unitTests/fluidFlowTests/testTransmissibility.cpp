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


#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"
#include "unitTests/fluidFlowTests/testSingleFlowUtils.hpp"
#include "physicsSolvers/fluidFlow/StencilDataCollection.hpp"
#include "mainInterface/ProblemManager.hpp"


using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;


CommandLineOptions g_commandLineOptions;


/// Provide every common xml input for the transmissibility tests
constexpr string_view xmlInputCommon =
  R"xml(
<Problem>
  <Solvers gravityVector="{ 0.0, 0.0, -9.81 }">
    <SinglePhaseFVM name="singlePhaseFlow"
                    logLevel="1"
                    discretization="singlePhaseTPFA"
                    targetRegions="{Region1}">
    </SinglePhaseFVM>
  </Solvers>
  <NumericalMethods>
    <FiniteVolume>
      <TwoPointFluxApproximation name="singlePhaseTPFA"/>
    </FiniteVolume>
  </NumericalMethods>
  <ElementRegions>
    <CellElementRegion name="Region1"
                        cellBlocks="{cb1}"
                        materialList="{water, rock}"/>
  </ElementRegions>
  <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  defaultDensity="1000"
                                  defaultViscosity="0.001"
                                  compressibility="5e-10"/>
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

  <Events maxTime="3.0">
    <PeriodicEvent
      name="solverApplications"
      forceDt="10"
      target="/Solvers/singlePhaseFlow" />
    <SoloEvent
      name="cellToCellDataCollectionEvent"
      beginTime="1.0"
      target="/Tasks/cellToCellDataCollection" />
  </Events>

  <Tasks>
    <CellToCellDataCollection
      name="cellToCellDataCollection"
      flowSolverName="singlePhaseFlow"
      meshBody="mesh"
      logLevel="1" />
  </Tasks>

)xml";

/// Provide the ending of the xml input for the transmissibility tests
static string_view constexpr xmlInputEnd =
  R"xml(
</Problem>
)xml";

/// Path of the StencilDataCollection in the hierarchy
static string_view constexpr stencilDataCollectionPath = "/Tasks/cellToCellDataCollection";


/// a "stack" array to represent 3d floating point data (ie. coords)
using Float3 = std::array< real64, 3 >;
/// a "stack" array to represent 3d integer data (ie. cell count / axis)
using Int3 = std::array< integer, 3 >;

/// Enumeration of the 3D axis to take into account for a structured mesh.
enum class Axis : integer { X = 0, Y = 1, Z = 2 };

/**
 * @brief a map of the half transmissibilities, from the first to the second cell identified by the
 * key global id.
 */
using TransmissibilityMap = std::map< std::pair< globalIndex, globalIndex >, real64 >;

/**
 * @brief the parameters for a given test instance.
 */
struct TestParams
{
  Int3 m_cellCount;
  Float3 m_meshSize;
  Float3 m_cellDistance;

  constexpr TestParams( Int3 cellCount, Float3 meshSize ):
    m_cellCount( cellCount ),
    m_meshSize( meshSize ),
    m_cellDistance( { meshSize[0] / real64( cellCount[0] ),
                      meshSize[1] / real64( cellCount[1] ),
                      meshSize[2] / real64( cellCount[2] ) } )
  {}
};

static real64 constexpr g_transmissibilityTolerance = 1.0e-11;
static Float3 constexpr g_testPermeability = { 2.0e-16, 2.0e-16, 2.0e-16 };
static real64 constexpr g_testNetToGross = 1.0;


void verifyStencilOutputStructured( string_view xmlInput, TestParams const & params );


TEST( TransmissibilityTest, stencilOutputVerificationIso )
{
  static string_view constexpr meshInput =
    R"xml(
  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="{C3D8}"
                  xCoords="{0, 30}"
                  yCoords="{0, 30}"
                  zCoords="{0, 30}"
                  nx="{3}"
                  ny="{3}"
                  nz="{3}"
                  cellBlockNames="{cb1}">
    </InternalMesh>
  </Mesh>
)xml";
  std::ostringstream xmlInput;
  xmlInput << xmlInputCommon << meshInput << xmlInputEnd;

  static TestParams constexpr params {
    { 3, 3, 3 }, // cellCount
    { 30.0, 30.0, 30.0 }, // meshSize
  };

  verifyStencilOutputStructured( xmlInput.str(), params );
}

TEST( TransmissibilityTest, StencilOutputVerificationAniso )
{
  static string_view constexpr meshInput =
    R"xml(
  <Mesh>
    <InternalMesh name="mesh"
                  elementTypes="{C3D8}"
                  xCoords="{0, 70.0}"
                  yCoords="{0, 10.0}"
                  zCoords="{0, 54.321}"
                  nx="{3}"
                  ny="{4}"
                  nz="{5}"
                  cellBlockNames="{cb1}">
    </InternalMesh>
  </Mesh>
)xml";
  std::ostringstream xmlInput;
  xmlInput << xmlInputCommon << meshInput << xmlInputEnd;

  static TestParams constexpr params {
    { 3, 4, 5 }, // cellCount
    { 70.0, 10.0, 54.321 }, // meshSize
  };

  verifyStencilOutputStructured( xmlInput.str(), params );
}



/**
 * @return The theorical half transmissiblity (from A to B or B to A) within a structured mesh.
 * @param params The test parameters
 * @param axis The axis in which we want to compute the transmissibility (structured mesh).
 */
real64 computeTransmissiblityStructured( TestParams const & params, Axis axis )
{
  real64 const faceArea = axis == Axis::X ? params.m_cellDistance[1]*params.m_cellDistance[2]:
                          axis == Axis::Y ? params.m_cellDistance[0]*params.m_cellDistance[2]:
                          params.m_cellDistance[0]*params.m_cellDistance[1];
  real64 const halfDistance = 0.5 * params.m_cellDistance[integer( axis )];

  real64 const transmissibility = g_testPermeability[integer( axis )] * g_testNetToGross * ( faceArea / halfDistance );
  // return half transmissibility
  return g_testNetToGross * transmissibility * 0.5;
}

/**
 * @brief Verify the transmissibility data from the StencilDataCollection for each connection
 * within a structured mesh.
 * @param transmissibilities The transmissibility map
 * @param axis The axis in which we want to verify the transmissibility (structured mesh).
 * @param params The test parameters
 */
void verifyTransmissibilityDataStructured( TransmissibilityMap transmissibilities, Axis axis,
                                           TestParams const & params )
{
  integer const cellBIdOffset = axis == Axis::X ? 1:
                                axis == Axis::Y ? params.m_cellCount[0]:
                                params.m_cellCount[0]*params.m_cellCount[1];
  integer const endX = axis == Axis::X ? params.m_cellCount[0] - 1 : params.m_cellCount[0];
  integer const endY = axis == Axis::Y ? params.m_cellCount[1] - 1 : params.m_cellCount[1];
  integer const endZ = axis == Axis::Z ? params.m_cellCount[2] - 1 : params.m_cellCount[2];
  real64 const expectedAxisTransmissibility = computeTransmissiblityStructured( params, axis );

  for( integer z = 0; z < endZ; ++z )
  {
    for( integer y = 0; y < endY; ++y )
    {
      for( integer x = 0; x < endX; ++x )
      {
        integer const cellAId = x + params.m_cellCount[0] * (y + z * params.m_cellCount[1]);
        integer const cellBId = cellAId + cellBIdOffset;
        real64 const transmissibilityOutput = transmissibilities[std::make_pair( cellAId, cellBId )];
        real64 const transmissibilityTolerance = g_transmissibilityTolerance * std::abs( transmissibilityOutput );
        EXPECT_NEAR( transmissibilityOutput, expectedAxisTransmissibility,
                     transmissibilityTolerance )<< GEOS_FMT( "Transmissibility data from {} does not match with expectation.",
                                                             StencilDataCollection::catalogName() );
      }
    }
  }
}

/**
 * @brief Verify the source data consistency and return a transmissibility map for easy further data access.
 * @param stencilData A StencilDataCollection object which contains the data gathered for the stencil.
 * @param params The test parameters
 * @return a TransmissibilityMap object containing the stencil data.
 */
TransmissibilityMap getTransmissibilityMap( StencilDataCollection const & stencilData,
                                            TestParams const & params )
{
  using VK = StencilDataCollection::viewKeyStruct;
  auto const & cellAGlobalId = stencilData.getReference< array1d< globalIndex > >( VK::cellAGlobalIdString() );
  auto const & cellBGlobalId = stencilData.getReference< array1d< globalIndex > >( VK::cellBGlobalIdString() );
  auto const & transmissibilityAB = stencilData.getReference< array1d< real64 > >( VK::transmissibilityABString() );
  auto const & transmissibilityBA = stencilData.getReference< array1d< real64 > >( VK::transmissibilityBAString() );

  // verify data size
  integer const nx = params.m_cellCount[0];
  integer const ny = params.m_cellCount[1];
  integer const nz = params.m_cellCount[2];
  integer const arraySize = transmissibilityAB.size();
  EXPECT_EQ( arraySize, 3*nx*ny*nz - nx*ny - ny*nz - nx*nz );
  // verify that array size is always equal for all buffers
  EXPECT_EQ( cellAGlobalId.size(), arraySize );
  EXPECT_EQ( cellBGlobalId.size(), arraySize );
  EXPECT_EQ( transmissibilityBA.size(), arraySize );

  TransmissibilityMap transmissibilities;
  for( int i = 0; i < arraySize; ++i )
  {
    transmissibilities[std::make_pair( cellAGlobalId[i], cellBGlobalId[i] )] = transmissibilityAB[i];
  }
  return transmissibilities;
}


/**
 * @brief Launch a test to verify if:
 *        - data of the output of the stencil is consistent,
 *        - output transmissiblity values are conform to computed expectations.
 * @param xmlInput The XML input to launch the test on.
 * @param params The test parameters.
 */
void verifyStencilOutputStructured( string_view xmlInput, TestParams const & params )
{
  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );
  ProblemManager & problem = state.getProblemManager();
  setupProblemFromXML( problem, xmlInput.data() );

  EXPECT_FALSE( problem.runSimulation() ) << "Simulation exited early.";

  { // checking results
    StencilDataCollection const & stencilData =
      problem.getGroupByPath< StencilDataCollection >( string( stencilDataCollectionPath ) );
    TransmissibilityMap const transmissibilities = getTransmissibilityMap( stencilData, params );

    // let's check each couple of cell in the x direction
    verifyTransmissibilityDataStructured( transmissibilities, Axis::X, params );
    verifyTransmissibilityDataStructured( transmissibilities, Axis::Y, params );
    verifyTransmissibilityDataStructured( transmissibilities, Axis::Z, params );
  }
}



int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
