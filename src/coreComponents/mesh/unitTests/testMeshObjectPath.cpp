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

/**
 * @file testMeshObjectPath.cpp
 */

#define MESH_OBJECT_PATH_PRIVATE_FUNCTION_UNIT_TESTING
#include "../MeshObjectPath.hpp"
#include "../MeshBody.hpp"
#include "dataRepository/Group.hpp"

#include <gtest/gtest.h>

namespace geos
{
using namespace dataRepository;

class TestMesh
{
public:


  static TestMesh & getTestMesh()
  {
    static TestMesh testMesh;
    return testMesh;
  }


  MeshObjectPath::permutationMapType const & pathPermutations() const
  {
    return m_pathPermutations;
  }

  Group const & meshBodies() const
  {
    return m_meshBodies;
  }

  Group & meshBodies()
  {
    return m_meshBodies;
  }

private:
  TestMesh():
    m_pathPermutations(),
    m_rootNode(),
    m_meshBodies( "meshBodies", m_rootNode )
  {
    createTestMesh();
  }

  void createTestMesh();

  MeshObjectPath::permutationMapType m_pathPermutations;
  conduit::Node m_rootNode;
  Group m_meshBodies;
};

void TestMesh::createTestMesh()
{

  m_pathPermutations["body0"]["level0"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body0"]["level0"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body0"]["level1"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body0"]["level1"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level0"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level0"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level1"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level1"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level1"]["region2"]={"subreg0", "subreg2"};
  m_pathPermutations["body3"]["level0"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body3"]["level0"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body3"]["level2"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body3"]["level2"]["region1"]={"subreg0", "subreg1"};


  for( auto const & meshBodyPair : m_pathPermutations )
  {
    MeshBody & meshBody = m_meshBodies.registerGroup< MeshBody >( meshBodyPair.first );
    for( auto const & meshLevelPair : meshBodyPair.second )
    {
      MeshLevel & level = meshBody.createMeshLevel( meshLevelPair.first );
      ElementRegionManager & elemRegMan = level.getElemManager();
      for( auto const & elemRegionPair : meshLevelPair.second )
      {
        elemRegMan.createChild( "CellElementRegion", elemRegionPair.first );
        CellElementRegion & elemRegion = elemRegMan.getRegion< CellElementRegion >( elemRegionPair.first );
        for( auto const & elemSubRegionName : elemRegionPair.second )
        {
          elemRegion.createElementSubRegion< CellElementSubRegion >( elemSubRegionName );
        }
      }
    }
  }
}

TEST( testMeshObjectPath, checkObjectTypeConsistency )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();

  {
    string const path = "nodeManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    EXPECT_TRUE( meshObjectPath.testCheckObjectTypeConsistency< NodeManager >() );
    EXPECT_FALSE( meshObjectPath.testCheckObjectTypeConsistency< EdgeManager >() );
    EXPECT_FALSE( meshObjectPath.testCheckObjectTypeConsistency< FaceManager >() );
    EXPECT_FALSE( meshObjectPath.testCheckObjectTypeConsistency< ElementRegionManager >() );
  }
}

TEST( testMeshObjectPath, fillPathTokens )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();

  map< string, std::vector< string > >
  entries =
  {
    { "ElementRegions", { "*", "*", "ElementRegions", "*", "*" } },
    { "nodeManager", { "*", "*", "nodeManager" } },
    { "edgeManager", { "*", "*", "edgeManager" } },
    { "faceManager", { "*", "*", "faceManager" } }
  };

  for( auto const & entry : entries )
  {
    string const & path = entry.first;
    std::vector< string > const & expectedPath = entry.second;
    size_t const pathSize = expectedPath.size();

    MeshObjectPath meshObjectPath( path, meshBodies );
    auto pathTokens = meshObjectPath.testFillPathTokens( path, meshBodies );

    EXPECT_TRUE( pathSize == pathTokens.size() );
    for( size_t a=0; a<pathTokens.size(); ++a )
    {
      EXPECT_TRUE( pathTokens[a] == expectedPath[a] );
    }
  }
}

TEST( testMeshObjectPath, meshObjectPathConstuction )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();


  std::map< string, MeshObjectPath::permutationMapType > pathsAndResults =
  {
    { "ElementRegions",
      { { "body0",
        { { "level0",
          { { "region0", {"subreg0", "subreg1"} },
            { "region1", {"subreg0", "subreg1"} }
          }
        },
          { "level1",
            { { "region0", {"subreg0", "subreg1"} },
              { "region1", {"subreg0", "subreg1"} }
            }
          }
        }
      },
        { "body1",
          { { "level0",
            { { "region0", {"subreg0", "subreg1"} },
              { "region1", {"subreg0", "subreg1"} }
            }
          },
            { "level1",
              { { "region0", {"subreg0", "subreg1"} },
                { "region1", {"subreg0", "subreg1"} },
                { "region2", {"subreg0", "subreg2"} }
              }
            }
          }
        },
        { "body3",
          { { "level0",
            { { "region0", {"subreg0", "subreg1"} },
              { "region1", {"subreg0", "subreg1"} }
            }
          },
            { "level2",
              { { "region0", {"subreg0", "subreg1"} },
                { "region1", {"subreg0", "subreg1"} }
              }
            }
          }
        }
      }
    },
    { "body0/level0/ElementRegions/{region0}",
      { { "body0",
        { { "level0",
          { { "region0", {"subreg0", "subreg1"} } }
        }
        }
      }
      }
    },
    { "body0/level0/ElementRegions/region0/subreg0",
      { { "body0",
        { { "level0",
          { { "region0", {"subreg0"} } }
        }
        }
      }
      }
    },

    { "{body0 body1}/*/nodeManager",
      { { "body0",
        { { "level0", {} },
          { "level1", {} }
        }
      },
        { "body1",
          { { "level0", {} },
            { "level1", {} }
          }
        }
      }
    },
    { "{body0 body3}/*/edgeManager",
      { { "body0",
        { { "level0", {} },
          { "level1", {} }
        }
      },
        { "body3",
          { { "level0", {} },
            { "level2", {} }
          }
        }
      }
    },
    { "{body0 body3}/level0/faceManager",
      { { "body0",
        { { "level0", {} } }
      },
        { "body3",
          { { "level0", {} }}
        }
      }
    }



  };

  for( auto const & entry : pathsAndResults )
  {
    string const & path = entry.first;
    MeshObjectPath::permutationMapType const & answer = entry.second;
    MeshObjectPath meshObjectPath( path, meshBodies );
    EXPECT_TRUE( meshObjectPath.pathPermutations() == answer );
  }
}


TEST( testMeshObjectPath, invalidMeshBody )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();
  {
    string const path = "*/level2/ElementRegions";
    ASSERT_THROW( MeshObjectPath meshObjectPath( path, meshBodies ), InputError );
  }
}


TEST( testMeshObjectPath, invalidMeshLevel )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();
  {
    string const path = "*/*/ElementRegions/{region2}";
    ASSERT_THROW( MeshObjectPath meshObjectPath( path, meshBodies ), InputError );
  }
}

TEST( testMeshObjectPath, invalidMeshRegion )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();
  {
    string const path = "*/*/ElementRegions/*/subreg2";
    ASSERT_THROW( MeshObjectPath meshObjectPath( path, meshBodies ), InputError );
  }
}


void checkSubRegionNames( std::vector< string > const & names )
{
  EXPECT_TRUE( names[0] == "subreg0" );
  EXPECT_TRUE( names[1] == "subreg1" );
  EXPECT_TRUE( names[2] == "subreg0" );
  EXPECT_TRUE( names[3] == "subreg1" );
  EXPECT_TRUE( names[4] == "subreg0" );
  EXPECT_TRUE( names[5] == "subreg1" );
  EXPECT_TRUE( names[6] == "subreg0" );
  EXPECT_TRUE( names[7] == "subreg1" );
  EXPECT_TRUE( names[8] == "subreg0" );
  EXPECT_TRUE( names[9] == "subreg1" );
  EXPECT_TRUE( names[10] == "subreg0" );
  EXPECT_TRUE( names[11] == "subreg1" );
  EXPECT_TRUE( names[12] == "subreg0" );
  EXPECT_TRUE( names[13] == "subreg1" );
  EXPECT_TRUE( names[14] == "subreg0" );
  EXPECT_TRUE( names[15] == "subreg1" );
  EXPECT_TRUE( names[16] == "subreg0" );
  EXPECT_TRUE( names[17] == "subreg2" );
  EXPECT_TRUE( names[18] == "subreg0" );
  EXPECT_TRUE( names[19] == "subreg1" );
  EXPECT_TRUE( names[20] == "subreg0" );
  EXPECT_TRUE( names[21] == "subreg1" );
  EXPECT_TRUE( names[22] == "subreg0" );
  EXPECT_TRUE( names[23] == "subreg1" );
  EXPECT_TRUE( names[24] == "subreg0" );
  EXPECT_TRUE( names[25] == "subreg1" );
}

void checkRegionNames( std::vector< string > const & names )
{
  EXPECT_TRUE( names[0] == "region0" );
  EXPECT_TRUE( names[1] == "region1" );
  EXPECT_TRUE( names[2] == "region0" );
  EXPECT_TRUE( names[3] == "region1" );
  EXPECT_TRUE( names[4] == "region0" );
  EXPECT_TRUE( names[5] == "region1" );
  EXPECT_TRUE( names[6] == "region0" );
  EXPECT_TRUE( names[7] == "region1" );
  EXPECT_TRUE( names[8] == "region2" );
  EXPECT_TRUE( names[9] == "region0" );
  EXPECT_TRUE( names[10] == "region1" );
  EXPECT_TRUE( names[11] == "region0" );
  EXPECT_TRUE( names[12] == "region1" );
}

TEST( testMeshObjectPath, forObjectsInPathFromMeshBodies )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group & meshBodies = testMesh.meshBodies();
  Group const & meshBodiesConst = meshBodies;
  string const path = "*/*/ElementRegions";
  MeshObjectPath meshObjectPath( path, meshBodiesConst );

  {
    std::vector< string > names;
    meshObjectPath.forObjectsInPath< CellElementSubRegion >( meshBodiesConst,
                                                             [&]( ElementSubRegionBase const & elemSubRegionBase )
    {
      names.push_back( elemSubRegionBase.getName() );
    } );
    checkSubRegionNames( names );
  }

  {
    std::vector< string > names;
    meshObjectPath.forObjectsInPath< CellElementSubRegion >( meshBodies,
                                                             [&]( ElementSubRegionBase & elemSubRegionBase )
    {
      names.push_back( elemSubRegionBase.getName() );
    } );
    checkSubRegionNames( names );
  }

  {
    std::vector< string > names;
    meshObjectPath.forObjectsInPath< CellElementRegion >( meshBodiesConst,
                                                          [&]( CellElementRegion const & elemRegionBase )
    {
      names.push_back( elemRegionBase.getName() );
    } );
    checkRegionNames( names );
  }
}

template< typename OBJECT_TYPE, typename CHECK_FUNC >
void testForObjectInPathsMeshLevel( Group & meshBodies,
                                    string const path,
                                    string const bodyName,
                                    string const levelName,
                                    CHECK_FUNC && checkFunc )
{
  MeshObjectPath meshObjectPath( path, meshBodies );
  MeshBody & meshBody = meshBodies.getGroup< MeshBody >( bodyName );
  MeshLevel & meshLevel = meshBody.getMeshLevel( levelName );
  std::vector< string > names;

  meshObjectPath.forObjectsInPath< OBJECT_TYPE >( meshLevel,
                                                  [&]( OBJECT_TYPE const & object )
  {
    names.push_back( object.getName() );
  } );

  checkFunc( names );
}

TEST( testMeshObjectPath, forObjectsInPath )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group & meshBodies = testMesh.meshBodies();

  testForObjectInPathsMeshLevel< CellElementSubRegion >( meshBodies,
                                                         "*/*/ElementRegions",
                                                         "body0",
                                                         "level0",
                                                         [&]( std::vector< string > const & names )
  {
    EXPECT_TRUE( names[0] == "subreg0" );
    EXPECT_TRUE( names[1] == "subreg1" );
    EXPECT_TRUE( names[2] == "subreg0" );
    EXPECT_TRUE( names[3] == "subreg1" );
  } );


  testForObjectInPathsMeshLevel< NodeManager >( meshBodies,
                                                "*/*/nodeManager",
                                                "body0",
                                                "level0",
                                                [&]( std::vector< string > const & names )
  {
    EXPECT_TRUE( names[0] == "nodeManager" );
  } );

  testForObjectInPathsMeshLevel< EdgeManager >( meshBodies,
                                                "*/*/edgeManager",
                                                "body0",
                                                "level0",
                                                [&]( std::vector< string > const & names )
  {
    EXPECT_TRUE( names[0] == "edgeManager" );
  } );

  testForObjectInPathsMeshLevel< FaceManager >( meshBodies,
                                                "*/*/faceManager",
                                                "body0",
                                                "level0",
                                                [&]( std::vector< string > const & names )
  {
    EXPECT_TRUE( names[0] == "faceManager" );
  } );



}

} /* namespace geos */
