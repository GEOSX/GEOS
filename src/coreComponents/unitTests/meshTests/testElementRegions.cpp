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

// TPL includes
#include <gtest/gtest.h>

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlockABC.hpp"

// special CMake-generated include
#include "tests/meshDirName.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::dataRepository;


CommandLineOptions g_commandLineOptions;

struct TestCase
{
  string name;
  bool isExpectedToPass = true;
  std::vector< string > stringsToMention;
  string_view xmlRegions;
};

class ElementRegionTestFixture : public ::testing::TestWithParam< TestCase >
{
public:
  ElementRegionTestFixture():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

  virtual ~ElementRegionTestFixture() = default;
private:
  GeosxState state;
};

TEST_P( ElementRegionTestFixture, testVTKImportRegionSyntaxes )
{
  TestCase const testCase = GetParam();

  string const pattern =
    R"xml(
      <Problem>
        <Mesh>
          <VTKMesh name="mesh" 
                   file="{}" />
        </Mesh>
        <ElementRegions>
          {}
        </ElementRegions>
      </Problem>
    )xml";
  string const xmlInput = GEOS_FMT( pattern,
                                    testMeshDir + "/box_hybrid_mesh.vtu",
                                    testCase.xmlRegions );

  ProblemManager & problem = getGlobalState().getProblemManager();
  problem.parseInputString( xmlInput );

  if( testCase.isExpectedToPass )
  {
    try
    {
      EXPECT_NO_FATAL_FAILURE( problem.problemSetup() ) << GEOS_FMT( "Test case '{}' did throw an error.",
                                                                     testCase.name );
    }
    catch( std::exception const & e )
    {
      GTEST_FAIL() << GEOS_FMT( "Test case '{}' did throw an exception ({}) but was not expected to :\n{}",
                                testCase.name, LvArray::system::demangle( typeid( e ).name() ), e.what() );
    }
  }
  else
  {
    try
    {
      EXPECT_NO_FATAL_FAILURE( problem.problemSetup() ) << GEOS_FMT( "Test case '{}' did throw an error.",
                                                                     testCase.name );
      GTEST_FAIL() << GEOS_FMT( "Test case '{}' did not thrown any exception but was expected to.",
                                testCase.name );
    }
    catch( InputError const & e )
    {
      string const expStr = e.what();
      for( auto const & str : testCase.stringsToMention )
      {
        bool isExceptionContainingStr = expStr.find( str ) != string::npos;
        EXPECT_TRUE( isExceptionContainingStr ) << GEOS_FMT( "Test case '{}' exception did not mention the string '{}'. Exception string:\n{}",
                                                             testCase.name, str, e.what() );
      }
    }
    catch( std::exception const & e )
    {
      GTEST_FAIL() << GEOS_FMT( "Test case '{}' did throw an exception ({}) but was not expected to :\n{}",
                                testCase.name, LvArray::system::demangle( typeid( e ).name() ), e.what() );
    }
    catch( ... )
    {
      GTEST_FAIL() << GEOS_FMT( "Test case '{}': Unexpected exception.",
                                testCase.name );
    }
  }
}

TestCase const vtkImportRegionSyntaxCases[] = {
  { // should not crash
    "regular cell-block list", true, { },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }"
                         cellBlocks="{ 1_hexahedra, 1_tetrahedra, 1_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                         cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
    )xml"
  },
  { // crash because of not existing primitive (2_pendecagonalPrism)
    "non existing 2_pendecagonalPrism", false, { "reservoir", "2_pendecagonalPrism" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" 
                         cellBlocks="{ 1_hexahedra, 1_tetrahedra, 1_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                         cellBlocks="{ 2_pendecagonalPrism, 2_hexahedra, 2_tetrahedra, 2_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
    )xml"
  },
  { // crash because of not unknown primitive name (helloWorld)
    "non existing helloWorld", false, { "overburden", "helloWorld" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" 
                          cellBlocks="{ helloWorld, 1_hexahedra, 1_tetrahedra, 1_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                          cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
    )xml"
  },
  { // crash because of lacking one primitive (1_hexahedra)
    "lacking 1_hexahedra", false, { "1_hexahedra" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }"
                          cellBlocks="{ 1_tetrahedra, 1_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                          cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
    )xml"
  },
  { // mentioning the same cell-blocks in multiple regions (1_hexahedra)
    "multiple 1_hexahedra", false, { "overburden", "reservoir", "1_hexahedra" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }"
                          cellBlocks="{ 1_hexahedra, 1_tetrahedra, 1_pyramids, 5_hexahedra, 5_tetrahedra, 5_pyramids }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                          cellBlocks="{ 1_hexahedra, 2_hexahedra, 2_tetrahedra, 2_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
    )xml"
  },
  { // should not crash
    "regular region attribute list", true, { },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" cellBlocks="{ 1, 5 }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlocks="{ 2, 6 }" />
    )xml"
  },
  { // mentioning the same region attribute in multiple cellBlocks (6 in overburden & reservoir)
    "multiple region 1", false, { "6", "region attribute", "overburden", "reservoir" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" cellBlocks="{ 1, 5, 6 }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlocks="{ 2, 6 }" />
    )xml"
  },
  { // forgetting a region attribute (5)
    "forget region 5", false, { "5", "region attribute" },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" cellBlocks="{ 1 }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlocks="{ 2, 6 }" />
    )xml"
  },
  { // should not crash
    "regular * wildcard", true, { },
    R"xml(
      <CellElementRegion name="everything" materialList="{ }" cellBlocks="{ * }" />
    )xml"
  },
  { // mentioning the same regions in multiple cellBlocks (because of "*")
    "* wildcard + region list", false, { "2", "6", "everything", "reservoir" },
    R"xml(
      <CellElementRegion name="everything" materialList="{ }" cellBlocks="{ * }" />
      <CellElementRegion name="reservoir" materialList="{ }" cellBlocks="{ 2, 6 }" />
    )xml"
  },
  { // should not crash
    "mixing selection methods on 2 regions", true, { },
    R"xml(
      <CellElementRegion name="overburden" materialList="{ }" cellBlocks="{ 1, 5 }" />
      <CellElementRegion name="reservoir" materialList="{ }"
                         cellBlocks="{ 2_hexahedra, 2_tetrahedra, 2_pyramids, 6_hexahedra, 6_tetrahedra, 6_pyramids }" />
    )xml"
  },
  { // should not crash
    "mixing all selection methods on 1 region", true, { },
    R"xml(
      <CellElementRegion name="everything" materialList="{ }"
                         cellBlocks="{ 1, 2_hexahedra, 2_tetrahedra, 2_pyramids, [5-6]_* }" />
    )xml"
  }
};
INSTANTIATE_TEST_SUITE_P( testElementRegions, ElementRegionTestFixture,
                          ::testing::ValuesIn( vtkImportRegionSyntaxCases ) );


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
