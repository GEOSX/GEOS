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

#include "gtest/gtest.h"

#include "mainInterface/initialization.hpp"

#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/CellElementSubRegion.hpp"


using namespace geos;

constexpr double maxCoordInX = 1.0;
constexpr double maxCoordInY = 2.0;
constexpr double maxCoordInZ = 3.0;

constexpr localIndex numElemsInX = 10;
constexpr localIndex numElemsInY = 11;
constexpr localIndex numElemsInZ = 12;

constexpr localIndex numNodesInX = numElemsInX + 1;
constexpr localIndex numNodesInY = numElemsInY + 1;
constexpr localIndex numNodesInZ = numElemsInZ + 1;

constexpr double dx = maxCoordInX / numElemsInX;
constexpr double dy = maxCoordInY / numElemsInY;
constexpr double dz = maxCoordInZ / numElemsInZ;

constexpr localIndex node_dJ = numNodesInX;
constexpr localIndex node_dK = numNodesInX * numNodesInY;

constexpr localIndex elem_dJ = numElemsInX;
constexpr localIndex elem_dK = numElemsInX * numElemsInY;

constexpr localIndex minOrder = 1;
constexpr localIndex maxOrder = 5;

class MeshGenerationTest : public ::testing::Test
{
protected:

  void SetUp() override
  {
    DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
    MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

    m_nodeManager = &mesh.getNodeManager();
    m_faceManager = &mesh.getFaceManager();
    m_edgeManager = &mesh.getEdgeManager();

    ElementRegionManager & elemManager = mesh.getElemManager();
    ASSERT_EQ( elemManager.numRegions(), 1 );

    ElementRegionBase & elemRegion = elemManager.getRegion( 0 );
    ASSERT_EQ( elemRegion.numSubRegions(), 1 );

    m_subRegion = &elemRegion.getSubRegion< CellElementSubRegion >( 0 );
  }

  NodeManager * m_nodeManager{};
  FaceManager * m_faceManager{};
  EdgeManager * m_edgeManager{};
  CellElementSubRegion * m_subRegion{};

  static void SetUpTestCase()
  {
    string const inputStream = GEOS_FMT(
      "<Problem>"
      "  <Mesh>"
      "    <InternalMesh"
      "      name=\"mesh1\""
      "      elementTypes=\"{{C3D8}}\""
      "      xCoords=\"{{0,{}}}\""
      "      yCoords=\"{{0,{}}}\""
      "      zCoords=\"{{0,{}}}\""
      "      nx=\"{{{}}}\""
      "      ny=\"{{{}}}\""
      "      nz=\"{{{}}}\""
      "      cellBlockNames=\"{{cb1}}\"/>"
      "  </Mesh>"
      "  <ElementRegions>"
      "    <CellElementRegion name=\"region1\" cellBlocks=\"{{cb1}}\" materialList=\"{{}}\"/>"
      "  </ElementRegions>"
      "</Problem>",
      maxCoordInX, maxCoordInY, maxCoordInZ, numElemsInX, numElemsInY, numElemsInZ );

    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult xmlResult = xmlDocument.loadString( inputStream );
    ASSERT_TRUE( xmlResult );

    xmlWrapper::xmlNode xmlProblemNode = xmlDocument.getChild( dataRepository::keys::ProblemManager );
    ProblemManager & problemManager = getGlobalState().getProblemManager();
    problemManager.processInputFileRecursive( xmlDocument, xmlProblemNode );

    // Open mesh levels
    DomainPartition & domain = problemManager.getDomainPartition();
    MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys.meshManager );
    meshManager.generateMeshLevels( domain );

    ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager();
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
    elementManager.processInputFileRecursive( xmlDocument, topLevelNode );
    elementManager.postInputInitializationRecursive();

    problemManager.problemSetup();
    problemManager.applyInitialConditions();
  }
};

TEST_F( MeshGenerationTest, sizes )
{
  EXPECT_EQ( m_subRegion->numNodesPerElement(), 8 );
  EXPECT_EQ( m_subRegion->numFacesPerElement(), 6 );

  localIndex const numNodes = numNodesInX * numNodesInY * numNodesInZ;
  EXPECT_EQ( numNodes, m_nodeManager->size() );

  localIndex const numElements = numElemsInX * numElemsInY * numElemsInZ;
  EXPECT_EQ( numElements, m_subRegion->size() );

  localIndex const numFaces = numNodesInX * numElemsInY * numElemsInZ +
                              numElemsInX * numNodesInY * numElemsInZ +
                              numElemsInX * numElemsInY * numNodesInZ;
  EXPECT_EQ( numFaces, m_faceManager->size() );

  localIndex const numEdges = numNodesInX * numNodesInY * numElemsInZ +
                              numNodesInX * numElemsInY * numNodesInZ +
                              numElemsInX * numNodesInY * numNodesInZ;
  EXPECT_EQ( numEdges, m_edgeManager->size() );
}

TEST_F( MeshGenerationTest, nodePositions )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = m_nodeManager->referencePosition();

  localIndex nodeIndex = 0;
  for( localIndex k = 0; k < numNodesInZ; ++k )
  {
    for( localIndex j = 0; j < numNodesInY; ++j )
    {
      for( localIndex i = 0; i < numNodesInX; ++i )
      {
        EXPECT_DOUBLE_EQ( X( nodeIndex, 0 ), i * dx );
        EXPECT_DOUBLE_EQ( X( nodeIndex, 1 ), j * dy );
        EXPECT_DOUBLE_EQ( X( nodeIndex, 2 ), k * dz );
        ++nodeIndex;
      }
    }
  }
}

TEST_F( MeshGenerationTest, elementCentersAndVolumes )
{
  arrayView2d< real64 const > const centers = m_subRegion->getElementCenter();
  arrayView1d< real64 const > const volumes = m_subRegion->getElementVolume();

  constexpr double VOLUME = dx * dy * dz;

  localIndex elemID = 0;
  for( localIndex k = 0; k < numElemsInZ; ++k )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex i = 0; i < numElemsInX; ++i )
      {
        EXPECT_DOUBLE_EQ( centers[ elemID ][ 0 ], i * dx + dx / 2.0 );
        EXPECT_DOUBLE_EQ( centers[ elemID ][ 1 ], j * dy + dy / 2.0 );
        EXPECT_DOUBLE_EQ( centers[ elemID ][ 2 ], k * dz + dz / 2.0 );
        EXPECT_NEAR( volumes[ elemID ], VOLUME, 1e-13 );
        ++elemID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, elemToNodeMap )
{
  arrayView2d< localIndex const, cells::NODE_MAP_USD > const & nodeMap = m_subRegion->nodeList();
  GEOS_ERROR_IF_NE( nodeMap.size( 1 ), 8 );

  localIndex elemID = 0;
  for( localIndex k = 0; k < numElemsInZ; ++k )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex i = 0; i < numElemsInX; ++i )
      {
        localIndex const firstNodeID = i + j * node_dJ + k * node_dK;

        EXPECT_EQ( firstNodeID, nodeMap( elemID, 0 ) );
        EXPECT_EQ( firstNodeID + 1, nodeMap( elemID, 1 ) );
        EXPECT_EQ( firstNodeID + 1 + node_dJ, nodeMap( elemID, 3 ) );
        EXPECT_EQ( firstNodeID + node_dJ, nodeMap( elemID, 2 ) );

        EXPECT_EQ( firstNodeID + node_dK, nodeMap( elemID, 4 ) );
        EXPECT_EQ( firstNodeID + node_dK + 1, nodeMap( elemID, 5 ) );
        EXPECT_EQ( firstNodeID + node_dK + 1 + node_dJ, nodeMap( elemID, 7 ) );
        EXPECT_EQ( firstNodeID + node_dK + node_dJ, nodeMap( elemID, 6 ) );
        ++elemID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, nodeToElemMap )
{
  ArrayOfArraysView< localIndex const > const & nodeToElemMap = m_nodeManager->elementList().toViewConst();

  localIndex nodeIndex = 0;
  for( localIndex k = 0; k < numNodesInZ; ++k )
  {
    for( localIndex j = 0; j < numNodesInY; ++j )
    {
      for( localIndex i = 0; i < numNodesInX; ++i )
      {
        localIndex const elemID = i + j * elem_dJ + k * elem_dK;

        std::vector< localIndex > expectedElems;
        if( k < numElemsInZ )
        {
          if( i < numElemsInX && j < numElemsInY )
            expectedElems.push_back( elemID );
          if( i > 0 && j < numElemsInY )
            expectedElems.push_back( elemID - 1 );
          if( i > 0 && j > 0 )
            expectedElems.push_back( elemID - 1 - elem_dJ );
          if( i < numElemsInX && j > 0 )
            expectedElems.push_back( elemID - elem_dJ );
        }

        if( k > 0 )
        {
          if( i < numElemsInX && j < numElemsInY )
            expectedElems.push_back( elemID - elem_dK );
          if( i > 0 && j < numElemsInY )
            expectedElems.push_back( elemID - elem_dK - 1 );
          if( i > 0 && j > 0 )
            expectedElems.push_back( elemID - 1 - elem_dJ - elem_dK );
          if( i < numElemsInX && j > 0 )
            expectedElems.push_back( elemID - elem_dJ - elem_dK );
        }

        localIndex const numElems = expectedElems.size();
        ASSERT_EQ( numElems, nodeToElemMap.sizeOfArray( nodeIndex ) );

        localIndex const * const nodeElems = nodeToElemMap[ nodeIndex ];
        std::vector< localIndex > elems( nodeElems, nodeElems + numElems );

        std::sort( elems.begin(), elems.end() );
        std::sort( expectedElems.begin(), expectedElems.end() );

        for( localIndex a = 0; a < numElems; ++a )
        {
          EXPECT_EQ( elems[ a ], expectedElems[ a ] );
        }

        ++nodeIndex;
      }
    }
  }
}

TEST_F( MeshGenerationTest, faceNodeMaps )
{
  arrayView2d< localIndex const > const & elementToFaceMap = m_subRegion->faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = m_faceManager->nodeList().toViewConst();
  arrayView2d< localIndex const > const & faceToElementMap = m_faceManager->elementList().toViewConst();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = m_nodeManager->faceList().toViewConst();

  GEOS_ERROR_IF_NE( elementToFaceMap.size( 1 ), 6 );

  array1d< localIndex > faceNodesFromElem( 4 );
  array1d< localIndex > faceNodesFromFace( 4 );

  localIndex elemID = 0;
  for( localIndex k = 0; k < numElemsInZ; ++k )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex i = 0; i < numElemsInX; ++i )
      {
        for( localIndex f = 0; f < 6; ++f )
        {
          m_subRegion->getFaceNodes( elemID, f, faceNodesFromElem );
          ASSERT_EQ( faceNodesFromElem.size(), 4 );

          localIndex const faceIndex = elementToFaceMap( elemID, f );

          ASSERT_EQ( faceToNodeMap.sizeOfArray( faceIndex ), 4 );
          for( localIndex a = 0; a < 4; ++a )
          {
            faceNodesFromFace[ a ] = faceToNodeMap( faceIndex, a );
          }

          if( elemID != faceToElementMap( faceIndex, 0 ) )
          {
            EXPECT_EQ( elemID, faceToElementMap( faceIndex, 1 ) );
          }

          std::sort( faceNodesFromElem.begin(), faceNodesFromElem.end() );
          std::sort( faceNodesFromFace.begin(), faceNodesFromFace.end() );

          EXPECT_EQ( faceNodesFromElem[ 0 ], faceNodesFromFace[ 0 ] );
          EXPECT_EQ( faceNodesFromElem[ 1 ], faceNodesFromFace[ 1 ] );
          EXPECT_EQ( faceNodesFromElem[ 2 ], faceNodesFromFace[ 2 ] );
          EXPECT_EQ( faceNodesFromElem[ 3 ], faceNodesFromFace[ 3 ] );

          for( localIndex a = 0; a < 4; ++a )
          {
            EXPECT_TRUE( nodeToFaceMap.contains( faceNodesFromElem[ a ], faceIndex ) );
          }
        }
        ++elemID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, faceElementMaps )
{
  arrayView2d< localIndex const > const & elementToFaceMap = m_subRegion->faceList();
  arrayView2d< localIndex const > const & faceToElementMap = m_faceManager->elementList();

  GEOS_ERROR_IF_NE( elementToFaceMap.size( 1 ), 6 );

  localIndex const elemIDOffset[6] = { -elem_dJ, -elem_dK, -1, 1, elem_dJ, elem_dK };

  localIndex elemID = 0;
  for( localIndex k = 0; k < numElemsInZ; ++k )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex i = 0; i < numElemsInX; ++i )
      {
        for( localIndex f = 0; f < 6; ++f )
        {
          localIndex const faceIndex = elementToFaceMap( elemID, f );

          if( elemID == faceToElementMap( faceIndex, 0 ) )
          {
            bool external = false;
            if( f == 0 && j == 0 )
              external = true;
            if( f == 1 && k == 0 )
              external = true;
            if( f == 2 && i == 0 )
              external = true;
            if( f == 3 && i == numElemsInX - 1 )
              external = true;
            if( f == 4 && j == numElemsInY - 1 )
              external = true;
            if( f == 5 && k == numElemsInZ - 1 )
              external = true;

            if( external )
            {
              EXPECT_EQ( -1, faceToElementMap( faceIndex, 1 ) );
            }
            else
            {
              EXPECT_EQ( elemID + elemIDOffset[ f ], faceToElementMap( faceIndex, 1 ) );
            }
          }
          else
          {
            EXPECT_EQ( elemID, faceToElementMap( faceIndex, 1 ) );
            EXPECT_EQ( elemID + elemIDOffset[ f ], faceToElementMap( faceIndex, 0 ) );
          }
        }
        ++elemID;
      }
    }
  }
}

bool walkEdgesToFindNeighbor( localIndex const node0,
                              localIndex const node1,
                              arraySlice1d< localIndex const > const & nodeEdges,
                              arrayView2d< localIndex const > const & edgeToNodeMap )
{
  for( localIndex const edgeIndex : nodeEdges )
  {
    if( edgeToNodeMap[ edgeIndex ][ 0 ] == node0 )
    {
      if( edgeToNodeMap[ edgeIndex ][ 1 ] == node1 )
        return true;
    }
    else
    {
      EXPECT_EQ( edgeToNodeMap[ edgeIndex ][ 1 ], node0 );
      if( edgeToNodeMap[ edgeIndex ][ 0 ] == node1 )
        return true;
    }
  }

  return false;
}

TEST_F( MeshGenerationTest, edgeNodeMaps )
{
  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = m_nodeManager->edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = m_edgeManager->nodeList();

  GEOS_ERROR_IF_NE( edgeToNodeMap.size( 1 ), 2 );

  localIndex nodeIndex = 0;
  for( localIndex k = 0; k < numNodesInZ; ++k )
  {
    for( localIndex j = 0; j < numNodesInY; ++j )
    {
      for( localIndex i = 0; i < numNodesInX; ++i )
      {
        localIndex numEdges = 0;
        if( i != 0 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeIndex, nodeIndex - 1, nodeToEdgeMap[ nodeIndex ], edgeToNodeMap ) );
          ++numEdges;
        }
        if( i != numNodesInX - 1 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeIndex, nodeIndex + 1, nodeToEdgeMap[ nodeIndex ], edgeToNodeMap ) );
          ++numEdges;
        }
        if( j != 0 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeIndex, nodeIndex - node_dJ, nodeToEdgeMap[ nodeIndex ], edgeToNodeMap ) );
          ++numEdges;
        }
        if( j != numNodesInY - 1 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeIndex, nodeIndex + node_dJ, nodeToEdgeMap[ nodeIndex ], edgeToNodeMap ) );
          ++numEdges;
        }
        if( k != 0 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeIndex, nodeIndex - node_dK, nodeToEdgeMap[ nodeIndex ], edgeToNodeMap ) );
          ++numEdges;
        }
        if( k != numNodesInZ - 1 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeIndex, nodeIndex + node_dK, nodeToEdgeMap[ nodeIndex ], edgeToNodeMap ) );
          ++numEdges;
        }

        EXPECT_EQ( numEdges, nodeToEdgeMap.sizeOfSet( nodeIndex ) );
        ++nodeIndex;
      }
    }
  }
}

TEST_F( MeshGenerationTest, edgeFaceMaps )
{
  arrayView2d< localIndex const > const & elementToFaceMap = m_subRegion->faceList();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = m_faceManager->nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = m_faceManager->edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = m_edgeManager->nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = m_edgeManager->faceList().toViewConst();

  GEOS_ERROR_IF_NE( elementToFaceMap.size( 1 ), 6 );

  localIndex elemID = 0;
  for( localIndex k = 0; k < numElemsInZ; ++k )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex i = 0; i < numElemsInX; ++i )
      {
        for( localIndex f = 0; f < 6; ++f )
        {
          localIndex const faceIndex = elementToFaceMap( elemID, f );

          ASSERT_EQ( faceToNodeMap.sizeOfArray( faceIndex ), 4 );
          ASSERT_EQ( faceToEdgeMap.sizeOfArray( faceIndex ), 4 );

          for( localIndex a = 0; a < 4; ++a )
          {
            localIndex node0 = faceToNodeMap( faceIndex, a );
            localIndex node1 = faceToNodeMap( faceIndex, ( a + 1 ) % 4 );
            if( node0 > node1 )
              std::swap( node0, node1 );

            localIndex const edgeIndex = faceToEdgeMap( faceIndex, a );
            EXPECT_EQ( edgeToNodeMap( edgeIndex, 0 ), node0 );
            EXPECT_EQ( edgeToNodeMap( edgeIndex, 1 ), node1 );

            bool foundFace = false;
            for( localIndex const id : edgeToFaceMap[ edgeIndex ] )
            {
              if( id == faceIndex )
                foundFace = true;
            }

            EXPECT_TRUE( foundFace );
          }
        }
        ++elemID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, highOrderMapsSizes )
{
  ProblemManager & problemManager = getGlobalState().getProblemManager();
  DomainPartition & domain = problemManager.getDomainPartition();
  MeshBody & meshBody = domain.getMeshBody( 0 );
  MeshManager & meshManager = problemManager.getGroup< MeshManager >( problemManager.groupKeys.meshManager );
  meshManager.generateMeshes( domain );
  for( int order = minOrder; order < maxOrder; order++ )
  {
    MeshLevel & meshLevel = meshBody.createMeshLevel( MeshBody::groupStructKeys::baseDiscretizationString(), GEOS_FMT( "TestLevel{}", order ), order );
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    NodeManager & nodeManager = meshLevel.getNodeManager();
    FaceManager & faceManager = meshLevel.getFaceManager();
    EdgeManager & edgeManager = meshLevel.getEdgeManager();
    CellBlockManagerABC const & cellBlockManager = meshBody.getCellBlockManager();
    nodeManager.setGeometricalRelations( cellBlockManager, elemManager, false );
    edgeManager.setGeometricalRelations( cellBlockManager, false );
    faceManager.setGeometricalRelations( cellBlockManager, elemManager, nodeManager, false );

    ASSERT_EQ( elemManager.numRegions(), 1 );

    ElementRegionBase & elemRegion = elemManager.getRegion( 0 );
    ASSERT_EQ( elemRegion.numSubRegions(), 1 );

    CellElementSubRegion & subRegion = elemRegion.getSubRegion< CellElementSubRegion >( 0 );

    EXPECT_EQ( subRegion.numNodesPerElement(), pow( order + 1, 3 ) );

    localIndex const numVertices = numNodesInX * numNodesInY * numNodesInZ;
    localIndex const numEdges = numNodesInX * numNodesInY * numElemsInZ + numNodesInX * numElemsInY * numNodesInZ + numElemsInX * numNodesInY *numNodesInZ;
    localIndex const numFaces = numNodesInX * numElemsInY * numElemsInZ + numElemsInX * numElemsInY * numNodesInZ + numElemsInX * numNodesInY *numElemsInZ;
    localIndex const numElems = numElemsInX * numElemsInY * numElemsInZ;
    localIndex const numNodes = numVertices + numEdges * (order-1) + numFaces * pow((order-1), 2 ) + numElems * pow((order-1), 3 );

    EXPECT_EQ( numNodes, nodeManager.size() );

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & nodeMap = subRegion.nodeList();
    EXPECT_EQ( nodeMap.size( 1 ), pow( order+1, 3 ) );
    arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList();
    EXPECT_EQ( edgeToNodeMap.size( 1 ), order+1 );
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
    for( localIndex f = 0; f < faceManager.size(); ++f )
    {
      EXPECT_EQ( faceToNodeMap.sizeOfArray( f ), pow( order+1, 2 ) );
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
