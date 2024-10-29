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
 * @file xmlWrapper.cpp
 */

#include <regex>

#include "xmlWrapper.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/MpiWrapper.hpp"
#include "dataRepository/KeyNames.hpp"

namespace geos
{
using namespace dataRepository;

namespace xmlWrapper
{

void validateString( string const & value, Regex const & regex )
{
  std::smatch m;
  bool inputValidated = std::regex_search( value, m, std::regex( regex.m_regexStr ) );
  if( !inputValidated || m.length() != ptrdiff_t( value.length() ) )
  {
    ptrdiff_t errorId = ( m.size()>0 && m.position( 0 )==0 ) ? m.length() : 0;
    GEOS_THROW( GEOS_FMT( "Input string validation failed at:\n"
                          "  \"{}\"\n"
                          "   {:>{}}\n"
                          "  Expected format: {}",
                          value, '^', errorId+1, regex.m_formatDescription ),
                InputError );
  }
}

void processInputException( std::exception const & ex, string const & targetAttributeName,
                            xmlNode const & targetNode,
                            xmlNodePos const & nodePos )
{
  xmlAttribute const attribute = targetNode.attribute( targetAttributeName.c_str() );
  string const inputStr = string( attribute.value() );
  std::ostringstream oss;
  string const exStr = ex.what();

  oss << "***** XML parsing error at node ";
  if( nodePos.isFound() )
  {
    xmlAttributePos const attPos = nodePos.getAttributeLine( targetAttributeName );
    string const & filePath = attPos.isFound() ? attPos.filePath : nodePos.filePath;
    int const line = attPos.isFound() ? attPos.line : nodePos.line;
    oss << "named " << targetNode.name() << ", attribute " << targetAttributeName
        << " (" << splitPath( filePath ).second << ", l." << line << ").";
  }
  else
  {
    oss << targetNode.path() << " (name='" << targetNode.attribute( "name" ).value() << "')/"
        << targetAttributeName;
  }
  oss << "\n***** Input value: '" << inputStr << '\'';
  oss << ( exStr[0]=='\n' ? exStr : "\n" + exStr );

  throw InputError( oss.str() );
}

template< typename T, int SIZE >
void stringToInputVariable( Tensor< T, SIZE > & target, string const & inputValue, Regex const & regex )
{
  validateString( inputValue, regex );

  std::istringstream ss( inputValue );
  auto const errorMsg = [&]( auto const & msg )
  {
    return GEOS_FMT( "{} for Tensor<{}> at position {} in input: \"{}\"", msg, SIZE, static_cast< int >( ss.tellg() ), inputValue );
  };

  // Read the head
  ss >> std::ws;
  GEOS_THROW_IF_NE_MSG( ss.peek(), '{', errorMsg( "Missing opening brace" ), InputError );
  ss.ignore(); // skip the opening brace

  // Read the values
  T value;
  int count = 0;
  while( ss >> value )
  {
    GEOS_THROW_IF_GE_MSG( count, SIZE, errorMsg( "Too many values" ), InputError );
    target[count++] = value;
    ss >> std::ws;
    if( count < SIZE )
    {
      GEOS_THROW_IF_NE_MSG( ss.peek(), ',', errorMsg( ss.peek() == '}' ? "Not enough values" : "Missing comma separator" ), InputError );
      ss.ignore(); // skip the comma
    }
  }
  GEOS_THROW_IF_LT_MSG( count, SIZE, errorMsg( "Not enough values" ), InputError );
  ss.clear();

  // Read the tail
  GEOS_THROW_IF_NE_MSG( ss.peek(), '}', errorMsg( "Missing closing brace" ), InputError );
  ss.ignore(); // skip the closing brace
  ss >> std::ws;
  GEOS_THROW_IF( ss.peek() != std::char_traits< char >::eof(), errorMsg( "Unparsed characters" ), InputError );
}

template void stringToInputVariable( Tensor< real32, 3 > & target, string const & inputValue, Regex const & regex );
template void stringToInputVariable( Tensor< real64, 3 > & target, string const & inputValue, Regex const & regex );
template void stringToInputVariable( Tensor< real64, 6 > & target, string const & inputValue, Regex const & regex );

/**
 * @brief Adds the filePath and character offset info on the node in filePathString
 * and charOffsetString attributes. This function allow to keep track of the source
 * filename & offset of each node.
 * @param targetNode the target node to add the informations on.
 * @param filePath the absolute path of the xml file containing the node.
 */
void addNodeFileInfo( xmlNode targetNode, string const & filePath )
{
  // we keep the file path and the character offset on each node so we keep track of these
  // informations, even if the nodes are manipulated within the xml hierarchy.
  targetNode.append_attribute( filePathString ).set_value( filePath.c_str() );
  targetNode.append_attribute( charOffsetString ).set_value( targetNode.offset_debug() );

  for( xmlNode subNode : targetNode.children() )
  {
    addNodeFileInfo( subNode, filePath );
  }
}
/**
 * @brief Returns true if the addNodeFileInfo() command has been called of the specified node.
 */
bool xmlDocument::hasNodeFileInfo() const
{ return !getFirstChild().attribute( filePathString ).empty(); }

void xmlDocument::addIncludedXML( xmlNode & targetNode, int const level )
{
  GEOS_THROW_IF( level > 100, "XML include level limit reached, please check input for include loops", InputError );

  string const currentFilePath = targetNode.attribute( filePathString ).value();

  // Schema currently allows a single unique <Included>, but a non-validating file may include multiple
  for( xmlNode includedNode : targetNode.children( includedListTag ) )
  {
    for( xmlNode fileNode : includedNode.children() )
    {
      // Extract the file name and construct full includedDirPath
      string const includedFilePath = [&]()
      {
        GEOS_THROW_IF_NE_MSG( string( fileNode.name() ), includedFileTag,
                              GEOS_FMT( "<{}> must only contain <{}> tags", includedListTag, includedFileTag ),
                              InputError );
        xmlAttribute const nameAttr = fileNode.attribute( "name" );
        string const fileName = nameAttr.value();
        GEOS_THROW_IF( !nameAttr || fileName.empty(),
                       GEOS_FMT( "<{}> tag must have a non-empty 'name' attribute", includedFileTag ),
                       InputError );
        return isAbsolutePath( fileName ) ? fileName : joinPath( splitPath( currentFilePath ).first, fileName );
      }();

      GEOS_LOG_RANK_0( "Included additionnal XML file: " << getAbsolutePath( includedFilePath ) );

      xmlDocument includedXmlDocument;
      xmlResult const result = includedXmlDocument.loadFile( includedFilePath, hasNodeFileInfo() );
      GEOS_THROW_IF( !result, GEOS_FMT( "Errors found while parsing included XML file {}\n"
                                        "Description: {}\nOffset: {}",
                                        includedFilePath, result.description(), result.offset ),
                     InputError );


      // All included files must contain a root node that must match the target node.
      // Currently, schema only allows <Included> tags at the top level (inside <Problem>).
      // We then proceed to merge each nested node from included file with the one in main.

      xmlNode includedRootNode = includedXmlDocument.getFirstChild();
      GEOS_THROW_IF_NE_MSG( string( includedRootNode.name() ), string( targetNode.name() ),
                            "Included document root does not match the including XML node", InputError );

      // Process potential includes in the included file to allow nesting
      addIncludedXML( includedRootNode, level + 1 );

      // Add each top level tag of imported document to current
      // This may result in repeated XML blocks, which will be implicitly merged when processed
      for( xmlNode importedNode : includedRootNode.children() )
      {
        targetNode.append_copy( importedNode );
      }

      m_originalBuffers[includedXmlDocument.getFilePath()] = includedXmlDocument.getOriginalBuffer();
    }
  }

  // Just in case, remove <Included> tags that have been processed
  while( targetNode.remove_child( includedListTag ) )
  {}
}

string buildMultipleInputXML( string_array const & inputFileList,
                              string const & outputDir )
{
  if( inputFileList.empty() )
  {
    return {};
  }
  if( inputFileList.size() == 1 )
  {
    return inputFileList[0];
  }

  string inputFileName = joinPath( outputDir, "composite_input.xml" );

  // Write the composite xml file on one rank
  if( MpiWrapper::commRank() == 0 )
  {
    xmlWrapper::xmlDocument compositeTree;
    xmlWrapper::xmlNode compositeRoot = compositeTree.appendChild( dataRepository::keys::ProblemManager );
    xmlWrapper::xmlNode includedRoot = compositeRoot.append_child( includedListTag );

    for( auto & fileName: inputFileList )
    {
      xmlWrapper::xmlNode fileNode = includedRoot.append_child( includedFileTag );
      fileNode.append_attribute( "name" ) = fileName.c_str();
    }

    compositeTree.saveFile( inputFileName );
  }

  // Everybody else has to wait before attempting to read
  MpiWrapper::barrier();

  return inputFileName;
}

bool isFileMetadataAttribute( string const & name )
{
  static const std::set< string > fileMetadataAttributes {
    "name", "xmlns:xsi", "xsi:noNamespaceSchemaLocation", xmlWrapper::filePathString, xmlWrapper::charOffsetString
  };
  return fileMetadataAttributes.find( name ) != fileMetadataAttributes.end();
}

constexpr size_t xmlDocument::npos;
size_t documentId=0;

xmlDocument::xmlDocument():
  pugiDocument(),
  m_rootFilePath( "CodeIncludedXML" + std::to_string( documentId++ ) )
{}

xmlResult xmlDocument::loadString( string_view content, bool loadNodeFileInfo )
{
  xmlResult result = pugiDocument.load_buffer( content.data(), content.size(),
                                               pugi::parse_default, pugi::encoding_auto );

  // keeping a copy of original buffer to allow line retrieval
  if( loadNodeFileInfo )
  {
    m_originalBuffers.clear();
    m_originalBuffers[m_rootFilePath] = content;

    addNodeFileInfo( getFirstChild(), m_rootFilePath );
  }

  return result;
}

xmlResult xmlDocument::loadFile( string const & path, bool loadNodeFileInfo )
{
  xmlResult result = pugiDocument.load_file( path.c_str(), pugi::parse_default, pugi::encoding_auto );
  m_rootFilePath = getAbsolutePath( path );

  // keeping a copy of original buffer to allow line retrieval
  if( loadNodeFileInfo )
  {
    std::ifstream t( path );
    std::stringstream buffer;
    buffer << t.rdbuf();

    m_originalBuffers.clear();
    m_originalBuffers[m_rootFilePath] = string( buffer.str() );

    addNodeFileInfo( getFirstChild(), getAbsolutePath( m_rootFilePath ) );
  }

  return result;
}

void xmlDocument::reset()
{ return pugiDocument.reset() ; }

xmlNode xmlDocument::appendChild( string const & name )
{ return pugiDocument.append_child( name.c_str() ); }

xmlNode xmlDocument::appendChild( xmlNodeType type )
{ return pugiDocument.append_child( type ); }

bool xmlDocument::saveFile( string const & path ) const
{ return pugiDocument.save_file( path.c_str() ); }

string const & xmlDocument::getFilePath() const
{ return m_rootFilePath; }

xmlNode xmlDocument::getFirstChild() const
{ return pugiDocument.first_child(); }

xmlNode xmlDocument::getChild( string const & name ) const
{ return pugiDocument.child( name.c_str() ); }

string const & xmlDocument::getOriginalBuffer() const
{ return m_originalBuffers.find( m_rootFilePath )->second; }

string const * xmlDocument::getOriginalBuffer( string const & filePath ) const
{
  map< string, string >::const_iterator it = m_originalBuffers.find( filePath );
  return it != m_originalBuffers.cend() ? &it->second : nullptr;
}

map< string, string > const & xmlDocument::getOriginalBuffers() const
{ return m_originalBuffers; }

xmlNodePos xmlDocument::getNodePosition( xmlNode const & node ) const
{
  size_t line = npos;
  size_t offsetInLine = npos;
  size_t offset = npos;
  xmlAttribute filePathAtt = node.attribute( filePathString );
  xmlAttribute charOffsetAtt = node.attribute( charOffsetString );
  string filePath;

  if( filePathAtt && charOffsetAtt )
  {
    filePath = string( filePathAtt.value() );
    offset = std::atoi( charOffsetAtt.value() );

    auto sourceBuffer = m_originalBuffers.find( filePath );
    string nodeName = node.name();

    if( sourceBuffer!=m_originalBuffers.end() &&
        offset > 0 && offset + nodeName.size() < sourceBuffer->second.size() &&
        sourceBuffer->second.substr( offset, nodeName.size() ) == nodeName )
    {
      line = std::count( sourceBuffer->second.begin(), sourceBuffer->second.begin() + offset, '\n' ) + 1;
      offsetInLine = offset - sourceBuffer->second.rfind( '\n', offset );
    }
    else
    {
      std::ostringstream oss;
      oss << "Node from file=" << filePath << ", path=" << node.path() << ", offset=" << offset;
      oss << " could not be found: ";
      oss << ( sourceBuffer == m_originalBuffers.end() ?
               ( filePathAtt.empty() ? "Source file not found." : "Source file info not loaded") :
               ( offset > 0 && offset + nodeName.size() < sourceBuffer->second.size() ?
                 "The offset doesn't lead to this node." :
                 "Offset out of bounds." )
               );
      oss << "\nSource files list:";
      for( auto const & buffer : m_originalBuffers )
        oss << "\n  " << buffer.first << " buffer size = " << std::to_string( buffer.second.size() );

      throw InputError( oss.str() );
    }
  }

  return xmlNodePos( *this, filePath, line, offsetInLine, offset );
}

xmlNodePos::xmlNodePos( xmlDocument const & document_, string const & filePath_, size_t line_,
                        size_t offsetInLine_, size_t offset_ ):
  xmlAttributePos( filePath_, line_, offsetInLine_, offset_ ),
  document( document_ )
{}

bool xmlNodePos::isFound() const
{ return line != xmlDocument::npos; }

size_t findTagEnd( string const & xmlBuffer, size_t const offset )
{
  bool outOfQuotes = true;
  auto bufferEnd = xmlBuffer.cend();
  for( auto it = xmlBuffer.cbegin() + offset; it != bufferEnd; ++it )
  {
    if( *it == '"' && *(it - 1) != '\\' )
    {
      outOfQuotes = !outOfQuotes;
    }
    else if( outOfQuotes && *it == '>' )
    {
      return it - xmlBuffer.cbegin();
    }
  }
  return string::npos;
}
size_t findAttribute( string const & attName, string const & xmlBuffer, size_t const tagBegin, size_t const tagEnd )
{
  if( !attName.empty())
  {
    size_t searchStart = tagBegin;
    try
    {
      std::smatch m;
      // As pugixml doesn't expose a way to get the attribute line or character offset, a regex
      // seem like a suficient solution for GEOS XMLs.
      // The regex search for the attribute name followed by an '=', eventually separated by spaces
      if( std::regex_search( xmlBuffer.cbegin() + searchStart, xmlBuffer.cbegin() + tagEnd,
                             m, std::regex( attName + "\\s*=\\s*\"" )))
      {
        size_t candidatePos = m.position() + searchStart;
        string previousString = xmlBuffer.substr( tagBegin, candidatePos - tagBegin );
        // We must be out of value surrounding quotes: the number of previous quotes '"' should be even
        // (ignoring the inner quotes preceded by '\\')
        size_t surroundingQuotesCount = 0;
        size_t quotePos = 0;
        while((quotePos = previousString.find( '"', quotePos + 1 )) != string::npos )
        {
          if( previousString[quotePos - 1] != '\\' )
            ++surroundingQuotesCount;
        }

        if(((surroundingQuotesCount % 1) == 0))
        {
          return candidatePos;
        }
        searchStart = candidatePos + attName.size();
      }
    }
    catch( std::regex_error const & )
    {}
  }
  return xmlDocument::npos;
}

xmlAttributePos xmlNodePos::getAttributeLine( string const & attName ) const
{
  string const * buffer = document.getOriginalBuffer( filePath );
  size_t attOffset = offset;
  size_t attLine = line;
  size_t attOffsetInLine = offsetInLine;

  if( isFound() && buffer != nullptr && offset < buffer->size() )
  {
    size_t tagEnd = findTagEnd( *buffer, offset );
    if( tagEnd != string::npos )
    {
      attOffset = findAttribute( attName, *buffer, offset, tagEnd );
      if( attOffset != string::npos )
      {
        attLine = line + std::count( buffer->cbegin() + offset, buffer->cbegin() + attOffset, '\n' );
        attOffsetInLine = attOffset - buffer->rfind( '\n', attOffset );
      }
    }
  }

  return xmlAttributePos( filePath, attLine, attOffsetInLine, attOffset );
}

xmlAttributePos::xmlAttributePos( string const & filePath_, size_t line_, size_t offsetInLine_,
                                  size_t offset_ ):
  filePath( filePath_ ),
  line( line_ ),
  offsetInLine( offsetInLine_ ),
  offset( offset_ )
{}

bool xmlAttributePos::isFound() const
{ return offset != xmlDocument::npos; }

string xmlAttributePos::toString() const
{
  if( line != xmlDocument::npos )
  {
    return splitPath( filePath ).second + ", l." + std::to_string( line );
  }
  else if( offset != xmlDocument::npos )
  {
    // line hasn't been found, we output the character offset.
    return splitPath( filePath ).second + ", offset " + std::to_string( offset );
  }
  else
  {
    // offset hasn't been found, filename is probably wrong too, we just output an error.
    return "Source file not found";
  }
}

} /* namespace xmlWrapper */

} /* namespace geos */
