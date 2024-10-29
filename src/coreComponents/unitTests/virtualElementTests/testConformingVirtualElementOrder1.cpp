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

// Source includes
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "finiteElement/elementFormulations/ConformingVirtualElementOrder1.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/MeshManager.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

// TPL includes
#include "gtest/gtest.h"

using namespace geos;
using namespace finiteElement;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;
constexpr real64 absTol = geos::testing::DEFAULT_ABS_TOL;
constexpr real64 relTol = geos::testing::DEFAULT_REL_TOL*10;

template< typename VEM >
GEOS_HOST_DEVICE
static void checkIntegralMeanConsistency( FiniteElementBase const & feBase,
                                          typename VEM::StackVariables const & stack,
                                          real64 & sumBasisFunctions )
{
  static constexpr localIndex
    maxSupportPoints = VEM::maxSupportPoints;
  real64 basisFunctionsIntegralMean[maxSupportPoints];
  VEM::calcN( 0, stack, basisFunctionsIntegralMean );
  sumBasisFunctions = 0;
  for( localIndex iBasisFun = 0;
       iBasisFun < feBase.template numSupportPoints< VEM >( stack ); ++iBasisFun )
  {
    sumBasisFunctions += basisFunctionsIntegralMean[iBasisFun];
  }
}

template< typename VEM >
GEOS_HOST_DEVICE
static void
checkIntegralMeanDerivativesConsistency( FiniteElementBase const & feBase,
                                         typename VEM::StackVariables const & stack,
                                         real64 & sumXDerivatives,
                                         real64 & sumYDerivatives,
                                         real64 & sumZDerivatives )
{
  static constexpr localIndex
    maxSupportPoints = VEM::maxSupportPoints;
  real64 const dummy[VEM::numNodes][3] { { 0.0 } };
  localIndex const k = 0;
  for( localIndex q = 0; q < VEM::numQuadraturePoints; ++q )
  {
    real64 basisDerivativesIntegralMean[maxSupportPoints][3]{};
    feBase.template getGradN< VEM >( k, q, dummy, stack, basisDerivativesIntegralMean );
    sumXDerivatives = 0; sumYDerivatives = 0; sumZDerivatives = 0;
    for( localIndex iBasisFun = 0;
         iBasisFun < feBase.template numSupportPoints< VEM >( stack );
         ++iBasisFun )
    {
      sumXDerivatives += basisDerivativesIntegralMean[iBasisFun][0];
      sumYDerivatives += basisDerivativesIntegralMean[iBasisFun][1];
      sumZDerivatives += basisDerivativesIntegralMean[iBasisFun][2];
    }
  }
}

template< typename VEM >
GEOS_HOST_DEVICE
static void
checkStabilizationMatrixConsistency ( arrayView2d< real64 const,
                                                   nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                                      localIndex const & cellIndex,
                                      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const & cellToNodes,
                                      arrayView2d< real64 const > const & cellCenters,
                                      FiniteElementBase const & feBase,
                                      typename VEM::StackVariables const & stack,
                                      arraySlice1d< real64 > & stabTimeMonomialDofsNorm )
{
  static constexpr localIndex
    maxSupportPoints = VEM::maxSupportPoints;
  localIndex const numCellPoints = cellToNodes[cellIndex].size();

  real64 cellDiameter = 0;
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
  {
    for( localIndex numOthVertex = 0; numOthVertex < numVertex; ++numOthVertex )
    {
      real64 vertDiff[ 3 ]{};
      LvArray::tensorOps::copy< 3 >( vertDiff, nodesCoords[cellToNodes( cellIndex, numVertex )] );
      LvArray::tensorOps::subtract< 3 >( vertDiff,
                                         nodesCoords[cellToNodes( cellIndex, numOthVertex )] );
      real64 const candidateDiameter = LvArray::tensorOps::l2NormSquared< 3 >( vertDiff );
      if( cellDiameter < candidateDiameter )
      {
        cellDiameter = candidateDiameter;
      }
    }
  }
  cellDiameter = LvArray::math::sqrt< real64 >( cellDiameter );
  real64 const invCellDiameter = 1.0/cellDiameter;

  stackArray2d< real64, 3*VEM::numNodes > monomialVemDofs( 3, numCellPoints );
  for( localIndex numVertex = 0; numVertex < numCellPoints; ++numVertex )
  {
    for( localIndex pos = 0; pos < 3; ++pos )
      monomialVemDofs[ pos ][ numVertex ] = invCellDiameter*
                                            (nodesCoords( cellToNodes( cellIndex, numVertex ), pos ) - cellCenters( cellIndex, pos ));
  }

  stackArray1d< real64, VEM::numNodes > stabTimeMonomialDofs( numCellPoints );
  real64 stabilizationMatrix[maxSupportPoints][maxSupportPoints]{};
  feBase.template addGradGradStabilizationMatrix< VEM, 1, false >( stack,

                                                                   stabilizationMatrix );
  stabTimeMonomialDofsNorm( 0 ) = 0.0;
  for( localIndex i = 0; i < numCellPoints; ++i )
  {
    stabTimeMonomialDofs( i ) = 0.0;
    stabTimeMonomialDofsNorm( 0 ) = 0.0;
    for( localIndex j = 0; j < numCellPoints; ++j )
    {
      stabTimeMonomialDofs( i ) += stabilizationMatrix[ i ][ j ];
    }
    stabTimeMonomialDofsNorm( 0 ) += stabTimeMonomialDofs( i ) * stabTimeMonomialDofs( i );
  }

  for( localIndex monomInd = 0; monomInd < 3; ++monomInd )
  {
    stabTimeMonomialDofsNorm( monomInd+1 ) = 0;
    for( localIndex i = 0; i < numCellPoints; ++i )
    {
      stabTimeMonomialDofs( i ) = 0.0;
      for( localIndex j = 0; j < numCellPoints; ++j )
      {
        stabTimeMonomialDofs( i ) += stabilizationMatrix[ i ][ j ] * monomialVemDofs( monomInd, j );
      }
      stabTimeMonomialDofsNorm( monomInd+1 ) += stabTimeMonomialDofs( i ) * stabTimeMonomialDofs( i );
    }
  }

}

template< typename VEM >
GEOS_HOST_DEVICE
static void checkSumOfQuadratureWeights( typename VEM::StackVariables stack,
                                         real64 & sumOfQuadratureWeights )
{
  static constexpr localIndex
    maxSupportPoints = VEM::maxSupportPoints;
  sumOfQuadratureWeights = 0.0;
  real64 const dummy[maxSupportPoints][3]{};
  for( localIndex q = 0; q < VEM::numQuadraturePoints; ++q )
  {
    real64 const weight = VEM::transformedQuadratureWeight( q, dummy, stack );
    sumOfQuadratureWeights += weight;
  }
}

template< localIndex MAXCELLNODES, localIndex MAXFACENODES >
static void testCellsInMeshLevel( MeshLevel const & mesh )
{
  // Get managers.
  ElementRegionManager const & elementManager = mesh.getElemManager();
  CellElementRegion const & cellRegion =
    elementManager.getRegion< CellElementRegion >( 0 );
  CellElementSubRegion const & cellSubRegion =
    cellRegion.getSubRegion< CellElementSubRegion >( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();

  // Get geometric properties to be passed as inputs.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodesCoords =
    nodeManager.referencePosition();
  traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const & cellToNodeMap = cellSubRegion.nodeList();
  arrayView2d< real64 const > cellCenters = cellSubRegion.getElementCenter();
  arrayView1d< real64 const > cellVolumes = cellSubRegion.getElementVolume();

  // Allocate and fill a VEM::MeshData struct.
  using VEM = ConformingVirtualElementOrder1< MAXCELLNODES, MAXFACENODES >;
  typename VEM::template MeshData< CellElementSubRegion > meshData;
  FiniteElementBase::initialize< VEM >( nodeManager, edgeManager,
                                        faceManager, cellSubRegion,
                                        meshData );

  // Arrays that store quantities to be tested.
  localIndex const numCells = cellSubRegion.getElementVolume().size();
  array1d< real64 > sumBasisFunctions( numCells ), sumXDerivatives( numCells ),
  sumYDerivatives( numCells ), sumZDerivatives( numCells ), sumOfQuadWeights( numCells );
  array2d< real64 > stabTimeMonomialDofsNorm( numCells, 4 );

  // Get views of test arrays
  arrayView1d< real64 > sumBasisFunctionsView = sumBasisFunctions.toView(),
                        sumXDerivativesView = sumXDerivatives.toView(),
                        sumYDerivativesView = sumYDerivatives.toView(),
                        sumZDerivativesView = sumZDerivatives.toView(),
                        sumOfQuadWeightsView = sumOfQuadWeights.toView();
  arrayView2d< real64 > stabTimeMonomialDofsNormView = stabTimeMonomialDofsNorm.toView();

  // Loop over cells on the device.
  forAll< geos::parallelDevicePolicy< > >( numCells, [=] GEOS_HOST_DEVICE
                                             ( localIndex const cellIndex )
  {
    typename VEM::StackVariables stack;
    VEM virtualElement;
    virtualElement.template setup< VEM >( cellIndex, meshData, stack );

    checkIntegralMeanConsistency< VEM >( virtualElement, stack,
                                         sumBasisFunctionsView( cellIndex ) );
    checkIntegralMeanDerivativesConsistency< VEM >( virtualElement, stack,
                                                    sumXDerivativesView( cellIndex ),
                                                    sumYDerivativesView( cellIndex ),
                                                    sumZDerivativesView( cellIndex )
                                                    );
    arraySlice1d< real64 > thisCellStabTimeMonomialDofs = stabTimeMonomialDofsNormView[cellIndex];
    checkStabilizationMatrixConsistency< VEM >( nodesCoords, cellIndex, cellToNodeMap,
                                                cellCenters, virtualElement, stack,
                                                thisCellStabTimeMonomialDofs );
    checkSumOfQuadratureWeights< VEM >( stack, sumOfQuadWeightsView( cellIndex ) );
  } );

  // Perform checks.
  forAll< serialPolicy >( numCells, [=] ( localIndex const cellIndex )
  {
    checkRelativeError( sumXDerivativesView( cellIndex ), 0.0, relTol, absTol,
                        "Sum of integral means of x-derivatives of basis functions" );
    checkRelativeError( sumYDerivativesView( cellIndex ), 0.0, relTol, absTol,
                        "Sum of integral means of y-derivatives of basis functions" );
    checkRelativeError( sumZDerivativesView( cellIndex ), 0.0, relTol, absTol,
                        "Sum of integral means of z-derivatives of basis functions" );
    checkRelativeError( sumBasisFunctionsView( cellIndex ), 1.0, relTol, absTol,
                        "Sum of integral means of basis functions" );
    for( localIndex monomInd = 0; monomInd < 4; ++monomInd )
    {
      string const name = "Product of stabilization matrix and degrees of freedom of monomial " + std::to_string( monomInd );
      checkRelativeError( stabTimeMonomialDofsNormView( cellIndex, monomInd ), 0.0, relTol, absTol,
                          name );
    }
    checkRelativeError( sumOfQuadWeightsView( cellIndex ), cellVolumes( cellIndex ), relTol, absTol,
                        "Sum of quadrature weights" );
  } );
}

TEST( ConformingVirtualElementOrder1, hexahedra )
{
  string const inputStream=
    "<Problem>"
    "  <Mesh>"
    "    <InternalMesh"
    "      name=\"cube\""
    "      elementTypes=\"{C3D8}\""
    "      xCoords=\"{0.0, 1.0}\""
    "      yCoords=\"{0.0, 1.0}\""
    "      zCoords=\"{0.0, 1.0}\""
    "      nx=\"{1}\""
    "      ny=\"{1}\""
    "      nz=\"{1}\""
    "      cellBlockNames=\"{cb1}\""
    "    />"
    "  </Mesh>"
    "  <ElementRegions>"
    "    <CellElementRegion name=\"region1\" cellBlocks=\"{cb1}\""
    "                       materialList=\"{dummy_material}\" />"
    "  </ElementRegions>"
    "</Problem>";

  xmlWrapper::xmlDocument inputFile;
  xmlWrapper::xmlResult xmlResult = inputFile.loadString( inputStream );
  if( !xmlResult )
  {
    GEOS_LOG_RANK_0( "XML parsed with errors!" );
    GEOS_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }
  xmlWrapper::xmlNode xmlProblemNode = inputFile.getChild( dataRepository::keys::ProblemManager );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );

  ProblemManager & problemManager = state.getProblemManager();
  problemManager.processInputFileRecursive( inputFile, xmlProblemNode );

  // Open mesh levels
  DomainPartition & domain  = problemManager.getDomainPartition();
  MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys
                                                                        .meshManager );
  meshManager.generateMeshLevels( domain );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  ElementRegionManager & elementManager = mesh.getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
  elementManager.processInputFileRecursive( inputFile, topLevelNode );
  elementManager.postInputInitializationRecursive();
  problemManager.problemSetup();

  // Test computed projectors for all cells in MeshLevel
  testCellsInMeshLevel< 10, 6 >( mesh );
}

TEST( ConformingVirtualElementOrder1, wedges )
{
  string const inputStream=
    "<Problem>"
    "  <Mesh>"
    "    <InternalMesh"
    "      name=\"wedges\""
    "      elementTypes=\"{C3D6}\""
    "      xCoords=\"{0.0, 1.0}\""
    "      yCoords=\"{0.0, 1.0}\""
    "      zCoords=\"{0.0, 1.0}\""
    "      nx=\"{1}\""
    "      ny=\"{1}\""
    "      nz=\"{1}\""
    "      cellBlockNames=\"{cb1}\""
    "    />"
    "  </Mesh>"
    "  <ElementRegions>"
    "    <CellElementRegion name=\"region1\" cellBlocks=\"{cb1}\""
    "                       materialList=\"{dummy_material}\" />"
    "  </ElementRegions>"
    "</Problem>";
  xmlWrapper::xmlDocument inputFile;
  xmlWrapper::xmlResult xmlResult = inputFile.loadString( inputStream );
  if( !xmlResult )
  {
    GEOS_LOG_RANK_0( "XML parsed with errors!" );
    GEOS_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }
  xmlWrapper::xmlNode xmlProblemNode = inputFile.getChild( dataRepository::keys::ProblemManager );

  GeosxState state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) );

  ProblemManager & problemManager = state.getProblemManager();
  problemManager.processInputFileRecursive( inputFile, xmlProblemNode );

  // Open mesh levels
  DomainPartition & domain  = problemManager.getDomainPartition();
  MeshManager & meshManager = problemManager.getGroup< MeshManager >
                                ( problemManager.groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
  ElementRegionManager & elementManager = mesh.getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
  elementManager.processInputFileRecursive( inputFile, topLevelNode );
  elementManager.postInputInitializationRecursive();
  problemManager.problemSetup();

  // Test computed projectors for all cells in MeshLevel
  testCellsInMeshLevel< 8, 9 >( mesh );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
