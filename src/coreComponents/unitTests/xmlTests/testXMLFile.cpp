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

// Source includes
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/FieldSpecificationBase.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/DataContext.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>
#include <algorithm>

using namespace geos;
using namespace geos::dataRepository;
using namespace xmlWrapper;

// Tests if the xml file parsing works with one file, one file with nested includes, and multiple files
TEST( testXML, testXMLFile )
{
  geos::ProblemManager & problemManager = geos::getGlobalState().getProblemManager();
  problemManager.parseCommandLineInput();
  problemManager.parseInputFile();

  // Check that we've read the full XML with all nested includes by inspecting boundary conditions
  EXPECT_TRUE( problemManager.getFieldSpecificationManager().hasGroup< FieldSpecificationBase >( "v0" ) );
}

/**
 * @brief Returns the xml nodes and attribute in hierarchy that GEOS is supposed to read in
 * the Group hierarchy.
 * @param document root xml document (which potentially contains includes).
 * @param targetNode the function will be called on this node and its children.
 * @param elementNames the list of xml elements names with the form "GroupName" or
 * "GroupName/WrapperName". The node name is supposed to be the name attribute value, it it exists,
 * or the tag name.
 */
void getElementsRecursive( xmlDocument const & document, xmlNode const & targetNode,
                           std::set< string > & elementNames )
{
  // Store here every node name that only exist in the xml (and not in the Group hierarchy).
  static const std::set< string > xmlOnlyNodes {
    "ElementRegions"
  };

  // The Group name will be the name attribute value, or the node tag name if the name attribute
  // doesn't exist.
  string const groupName = [&]() {
    xmlAttribute nameAtt = targetNode.attribute( "name" );
    return nameAtt ? string( nameAtt.value() ) : string( targetNode.name() );
  }();

  if( xmlOnlyNodes.find( groupName ) == xmlOnlyNodes.end() )
  {
    elementNames.emplace( groupName );

    for( xmlAttribute att : targetNode.attributes() )
    {
      if( !isFileMetadataAttribute( att.name() ) )
      {
        elementNames.emplace( groupName + '/' + att.name() );
      }
    }
  }

  for( xmlNode subNode : targetNode.children() )
  {
    getElementsRecursive( document, subNode, elementNames );
  }
}
/**
 * @brief Calls a lambda for the DataContext of this Group and the ones of its children Group and
 * Wrapper if they can be casted to a DataFileContext type.
 * @tparam LAMBDA the type of functor to call. The functor takes in parameter a DataContext of
 * TARGET_DC type.
 * @param func the functor to call for each DataFileContext of the group and its children
 */
template< typename FUNC >
void forAllDataFileContext( Group const & group, FUNC func )
{
  if( auto * groupCtx = dynamic_cast< DataFileContext const * >( &group.getDataContext() ) )
  {
    func( *groupCtx );
  }
  for( auto const & wrapperIterator : group.wrappers() )
  {
    auto * wrapperCtx = dynamic_cast< DataFileContext const * >( &wrapperIterator.second->getDataContext() );
    if( wrapperCtx )
    {
      func( *wrapperCtx );
    }
  }
  for( auto subGroup : group.getSubGroups() )
  {
    forAllDataFileContext( *subGroup.second, func );
  }
}
/**
 * @brief Verify if the DataFileContext data correctly locates from where the object comes from (in the
 * correct document/include, at the correct line and at the correct character offset).
 * @param document root xml document (which potentially contains includes).
 * @param fileContext the DataFileContext to verify.
 */
void verifyDataFileContext( DataFileContext const & fileContext,
                            xmlDocument const & document,
                            std::set< string > & verifications )
{
  verifications.emplace( fileContext.getTargetName() );

  string const & strToVerify = fileContext.getTypeName();
  string const & errInfo = "Verifying " + strToVerify + " in " + fileContext.toString();

  // verifying if all DataFileContext data have been found
  EXPECT_FALSE( fileContext.getFilePath().empty() ) << errInfo;
  EXPECT_FALSE( fileContext.getTargetName().empty() ) << errInfo;
  EXPECT_FALSE( fileContext.getTypeName().empty() ) << errInfo;
  EXPECT_NE( fileContext.getOffset(), xmlDocument::npos ) << errInfo;
  EXPECT_NE( fileContext.getLine(), xmlDocument::npos ) << errInfo;
  EXPECT_NE( fileContext.getOffsetInLine(), xmlDocument::npos ) << errInfo;

  // will crash if the source buffer specified by fileContext.getFilePath() isn't available
  string const * buffer = document.getOriginalBuffer( fileContext.getFilePath() );
  EXPECT_NE( buffer, nullptr ) << errInfo;

  size_t curLine = 1;
  bool lineFound = false;

  // Does fileContext.getOffset() locates the object?
  EXPECT_LT( fileContext.getOffset() + strToVerify.size(), buffer->size() ) << errInfo;
  EXPECT_EQ( strToVerify, buffer->substr( fileContext.getOffset(), strToVerify.size() ) ) << errInfo;

  // Were trying to reach the line return by DataFileContext::getLine()
  for( size_t offset=0;
       offset < buffer->size() && curLine < fileContext.getLine();
       ++offset )
  {
    if( (*buffer)[offset] == '\n' )
    {
      curLine++;
    }

    // the theorical line has been reach!
    if( curLine == fileContext.getLine() )
    {
      // Does fileContext.getLine() and fileContext.getOffsetInLine() locates the object?
      EXPECT_LT( offset + fileContext.getOffsetInLine() + strToVerify.size(),
                 buffer->size() ) << errInfo;
      EXPECT_EQ( strToVerify,
                 buffer->substr( offset + fileContext.getOffsetInLine(), strToVerify.size() ) ) << errInfo;
      lineFound = true;
    }
  }
  // does the fileContext line has been reached?
  EXPECT_TRUE( lineFound );
}

/**
 * @brief Returns the element that exists in setB but not in setA.
 */
std::set< string > getDifference( std::set< string > const & setA,
                                  std::set< string > const & setB )
{
  std::set< string > result;
  std::set_difference( setA.cbegin(), setA.cend(),
                       setB.cbegin(), setB.cend(),
                       std::inserter( result, result.begin() ) );
  return result;
}

// Tests
// - if the line information of each nodes and attributes can be all retrieved,
// - if the resulting Group & Wrapper hierarchy matches with the input xml documents and includes hierarchy.
TEST( testXML, testXMLFileLines )
{
  xmlDocument xmlDoc;
  ProblemManager & problemManager = getGlobalState().getProblemManager();

  {
    problemManager.parseCommandLineInput();
    Group & commandLine = problemManager.getGroup( problemManager.groupKeys.commandLine );
    string const & inputFileName = commandLine.getReference< string >( problemManager.viewKeys.inputFileName );
    xmlDoc.loadFile( inputFileName, true );
    problemManager.parseXMLDocument( xmlDoc );
  }

  GEOS_LOG( "Loaded files : " );
  for( auto const & buffer: xmlDoc.getOriginalBuffers() )
  {
    GEOS_LOG( "    " << buffer.first << " (" << buffer.second.size() << " chars)" );
  }

  std::set< string > expectedElements;
  getElementsRecursive( xmlDoc, xmlDoc.getFirstChild(), expectedElements );

  std::set< string > verifiedElements;
  forAllDataFileContext( problemManager, [&]( DataFileContext const & ctx )
  {
    verifyDataFileContext( ctx, xmlDoc, verifiedElements );
  } );

  std::set< string > const notFound = getDifference( expectedElements, verifiedElements );
  EXPECT_TRUE( notFound.empty() ) << "Info : There should not exist xml element that were not in "
                                     "the Group hierarchy.\nElements not found in Group hierarchy : {"
                                  << stringutilities::join( notFound, "," ) << "}";

  std::set< string > const notExpected = getDifference( verifiedElements, expectedElements );
  EXPECT_TRUE( notExpected.empty() ) << "Info : There should not exist an object in the Group "
                                        "hierarchy that contains a DataFileContext but which were not "
                                        "declared in the Xml.\nElements not found in XML hierarchy : {"
                                     << stringutilities::join( notExpected, "," ) << "}";
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  GeosxState state( basicSetup( argc, argv, true ) );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
