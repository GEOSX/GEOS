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
 * @file xmlWrapper.hpp
 */

#ifndef GEOS_DATAREPOSITORY_XMLWRAPPER_HPP_
#define GEOS_DATAREPOSITORY_XMLWRAPPER_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "DefaultValue.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "LvArray/src/output.hpp"
#include "LvArray/src/input.hpp"
#include "codingUtilities/StringUtilities.hpp"

// TPL includes
#include <pugixml.hpp>

// System includes
#include <algorithm>
#include <sstream>

namespace geos
{


/**
 * Wraps/provides facilities to process entries from an xml file into the appropriate
 * data types. xmlWrapper provides some aliases that will wrap the underlying xml package being
 * used to extract data from the xml file, and a set of functions that facilitate
 * the parsing of string data from the xml into the variables that will hold those values in the
 * code.
 */
namespace xmlWrapper
{

/// Alias for the type of the result from an xml parse attempt.
using xmlResult = pugi::xml_parse_result;

/// Alias for the type of an xml node.
/// An xmlNode behave as a pointer object: passing it by value to a function which modify it
/// will modify the original xmlNode object.
using xmlNode = pugi::xml_node;

/// Alias for the type of an xml attribute.
/// An xmlAttribute behave as a pointer object: passing it by value to a function which modify it
/// will modify the original xmlAttribute object.
using xmlAttribute = pugi::xml_attribute;

/// Alias for the type variant of an xml node.
using xmlNodeType = pugi::xml_node_type;

class xmlDocument;

/**
 * @struct xmlAttributePos
 *
 * Stores the source file path, and position in the file of a xml attribute.
 */
struct xmlAttributePos
{
  /// Path of the file containing this element
  string filePath;
  /// Line where the element is defined. Start at 1.
  /// Default value is xmlDocument::npos
  size_t const line;
  /// Character offset of this element in the line that contains it (starting from 1)
  /// Equals to xmlDocument::npos if it couldn't be determined.
  size_t const offsetInLine;
  /// Character offset of this element in the file that contains it (starting from 0)
  /// Equals to xmlDocument::npos if it couldn't be determined.
  size_t const offset;

  /**
   * @brief Constructor of this struct.
   * @param filePath the path of the original xml file containing this attribute
   * @param line Line where the element is defined. Start at 1.
   * @param offsetInLine Character offset of this element in the line that contains it (starting from 1)
   * @param offset Character offset of this element in the file that contains it (starting from 0)
   */
  xmlAttributePos( string const & filePath, size_t line, size_t offsetInLine, size_t offset );
  /**
   * @return false if the position is undefined.
   */
  bool isFound() const;
  /**
   * @return a string containing the file name and position in the file.
   */
  string toString() const;
};

/**
 * @struct xmlNodePos
 *
 * Used to retrieve the position of dataRepository::Wrapper in XML
 */
struct xmlNodePos : xmlAttributePos
{
  /// Reference of the main xmlDocument that contains all original file buffers.
  xmlDocument const & document;

  /**
   * @brief Constructor of this struct.
   * @param document an xml document containing this node, or including a file which includes it
   * @param filePath the path of the original xml file containing this node
   * @param line Line where the node is defined. Start at 1.
   * @param offsetInLine Character offset of this node in the line that contains it (starting from 1)
   * @param offset Character offset of this node in the file that contains it (starting from 0)
   */
  xmlNodePos( xmlDocument const & document, string const & filePath, size_t line, size_t offsetInLine, size_t offset );
  /**
   * @return false if the position is undefined.
   */
  bool isFound() const;
  /**
   * @brief Compute the xmlAttributePos of an xml attribute
   * @param attName the name of the attribute to locate
   * @return an xmlAttributePos object that represents the position of the target node.
   */
  xmlAttributePos getAttributeLine( string const & attName ) const;
};

/**
 * @class xmlDocument
 *
 * Wrapper class for the type of xml document.
 * This class exists to intercept file / string loading methods, and to keep the loaded buffers,
 * in order to retrieve the source file and line of nodes and attributes.
 */
class xmlDocument
{
public:
  /// Error value for when an offset / line position is undefined.
  static constexpr size_t npos = string::npos;

  /**
   * @brief Construct an empty xmlDocument that waits to load something.
   */
  xmlDocument();
  /**
   * @brief non-copyable
   */
  xmlDocument( const xmlDocument & ) = delete;
  /**
   * @brief move constructor
   */
  xmlDocument( xmlDocument && ) = default;

  /**
   * @return the first child of this document (typically in GEOS, the Problem node)
   */
  xmlNode getFirstChild() const;
  /**
   * @return a child with the specified name
   * @param name the tag name of the node to find
   */
  xmlNode getChild( string const & name ) const;
  /**
   * @return the original file buffer loaded during the last load_X() call on this object.
   */
  string const & getOriginalBuffer() const;
  /**
   * @return If the specified file at the "filePath" is the loaded document of the instance,
   * or one of its includes, returns its original buffer as a string.
   * Returns nullptr if filePath is not a loaded and available document.
   * @param filePath the path of the file which buffer must be returned.
   */
  string const * getOriginalBuffer( string const & filePath ) const;
  /**
   * @return a map containing the original buffers of the document and its includes, indexed by file path.
   */
  map< string, string > const & getOriginalBuffers() const;
  /**
   * @return If loadFile() has been loaded, returns the path of the source file.
   * If another load method has been called, it returns a generated unique value.
   */
  string const & getFilePath() const;
  /**
   * @brief If the node file information are loaded, compute the position of a node.
   * @param node the node to locate
   * @return an xmlNodePos object that represents the position of the target node.
   * @throws an InputError if the node position is not found and source file is loaded.
   */
  xmlNodePos getNodePosition( xmlWrapper::xmlNode const & node ) const;

  /**
   * @brief Load document from zero-terminated string. No encoding conversions are applied.
   * Free any previously loaded xml tree.
   * Wrapper of pugi::xml_document::loadBuffer() method.
   * @param content the string containing the document content
   * @param loadNodeFileInfo Load the node source file info, allowing getNodePosition() to work.
   * @return an xmlResult object representing the parsing resulting status.
   */
  xmlResult loadString( string_view content, bool loadNodeFileInfo = false );

  /**
   * @brief Load document from file. Free any previously loaded xml tree.
   * Wrapper of pugi::xml_document::loadBuffer() method.
   * @param path the path of an xml file to load.
   * @param loadNodeFileInfo Load the node source file info, allowing getNodePosition() to work.
   * @return an xmlResult object representing the parsing resulting status.
   */
  xmlResult loadFile( string const & path, bool loadNodeFileInfo = false );

  /**
   * @brief Reset document
   */
  void reset( );

  /**
   * @brief Add a root element to the document
   * @param name the tag name of the node to add
   * @return the added node
   */
  xmlNode appendChild( string const & name );

  /**
   * @brief Add a root element to the document
   * @param type the type of the node to add to the root of the document.
   * As an exemple, node_declaration is useful to add the "<?xml ?>" node.
   * @return the added node
   */
  xmlNode appendChild( xmlNodeType type = xmlNodeType::node_element );

  /**
   * @brief Save the XML to a file
   * @param path the file path
   * @return true if the file has successfuly been saved
   * @return false otherwise
   */
  bool saveFile( string const & path ) const;

  /**
   * @brief Function to add xml nodes from included files.
   * @param targetNode the node for which to look for included children specifications
   * @param level include tree level used to detect circular includes
   *
   * This function looks for a "Included" node under the targetNode, loops over all subnodes under the "Included"
   * node, and then parses the file specified in those subnodes taking all the nodes in the file and adding them to
   * the targetNode.
   * Each found includes are added to xmlDocument::m_originalBuffers.
   */
  void addIncludedXML( xmlNode & targetNode, int level = 0 );

  /**
   * @return True if loadNodeFileInfo was true during the last load_X call.
   */
  bool hasNodeFileInfo() const;

private:
  /// original xml_document object that this class aims to wrap
  pugi::xml_document pugiDocument;

  /// Used to retrieve node positions as pugixml buffer is private and processed.
  map< string, string > m_originalBuffers;
  /// @see getFilePath()
  string m_rootFilePath;
};

/**
 * @brief constexpr variable to hold name for inserting the file path into the xml file.
 *
 * This is used because we would like the option to hold the file path in the xml structure.
 * The name is uglified with underscores to avoid collisions with real attribute names.
 */
constexpr char const filePathString[] = "__filePath__";

/**
 * @brief constexpr variable to hold node character offset from the start of the xml file.
 *
 * This is used because we would like the option to hold the offset in the xml structure.
 * The name is uglified with underscores to avoid collisions with real attribute names.
 */
constexpr char const charOffsetString[] = "__charOffset__";

/// XML tag name for included sections
constexpr char const includedListTag[] = "Included";

/// XML tag name for included files
constexpr char const includedFileTag[] = "File";

/**
 * @brief Function to handle multiple input xml files.
 * @param inputFileList the list of input xml files
 * @param outputDir the output directory to place the composite input file in
 * @return inputFilePath the input xml file name
 *
 * This function checks for multiple xml files, and will build
 * a new input xml file with an included block if neccesary
 */
string buildMultipleInputXML( string_array const & inputFileList,
                              string const & outputDir = {} );

/**
 * @return true if the attribute with the specified name declares metadata relative to the xml
 * @param name the name of an attribute
 */
bool isFileMetadataAttribute( string const & name );

/**
 * @name String to variable parsing.
 *
 * These functions take in @p value and parse that string based on the type of
 * @p target. The function implementation should provide sufficient error checking
 * in the case that @p value is formatted incorrectly for the type specified in @p target.
 */
///@{

/**
 * @throw An InputError if the string value could not be validated with the provided regular expression.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 * @param value The string to validate
 * @param regex The regular expression used for validating the string value.
 */
void validateString( string const & value, Regex const & regex );

/**
 * @brief Helper method to process an exception that has been thrown during xml parsing.
 * @param ex The caught exception.
 * @param targetAttributeName The name of the attribute in which an input exception occured.
 * @param targetNode The node from which this Group is interpreted.
 * @param nodePos the target node position.
 * @throw An InputError with the node & attribut info, including the nearest xml line possible.
 */
void processInputException( std::exception const & ex, string const & targetAttributeName,
                            xmlWrapper::xmlNode const & targetNode,
                            xmlWrapper::xmlNodePos const & nodePos );

/**
 * @brief Parse a string and fill a variable with the value(s) in the string.
 * @tparam T the type of variable fill with string value
 * @param[out] target the object to read values into
 * @param[in]  value  the string that contains the data to be parsed into target
 * @param[in]  regex  the regular expression used for validating the string value.
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 * @return
 */
template< typename T >
std::enable_if_t< traits::CanStreamInto< std::istringstream, T > >
stringToInputVariable( T & target, string const & value, Regex const & regex )
{
  validateString( value, regex );

  std::istringstream ss( value );
  ss >> target;
  GEOS_THROW_IF( ss.fail() || !ss.eof(),
                 "Error detected while parsing string \"" << value <<
                 "\" to type " << LvArray::system::demangleType< T >(),
                 InputError );
}

/**
 * @brief Parse a string and fill a R1Tensor with the value(s) in the string.
 * @param[out] target the object to read values into
 * @param[in]  value  the string that contains the data to be parsed into target
 * @param[in]  regex  the regular expression used for validating the string value.
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 */
template< typename T, int SIZE >
void
stringToInputVariable( Tensor< T, SIZE > & target, string const & value, Regex const & regex );

/**
 * @brief Parse a string and fill an Array with the value(s) in the string.
 * @tparam T    data type of the array
 * @tparam NDIM number of dimensions of the array
 * @tparam PERMUTATION the permutation of the array
 * @param[out] array the array to read values into
 * @param[in]  value the string that contains the data to be parsed into target
 * @param[in]  regex the regular expression used for validating the string value.
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 * @return
 */
template< typename T, int NDIM, typename PERMUTATION >
std::enable_if_t< traits::CanStreamInto< std::istringstream, T > >
stringToInputVariable( Array< T, NDIM, PERMUTATION > & array, string const & value, Regex const & regex )
{
  validateString( value, regex );

  LvArray::input::stringToArray( array, string( stringutilities::trimSpaces( value ) ) );
}

///@}

namespace internal
{

/// Defines a static constexpr bool canParseVariable that is true iff the template parameter T
/// is a valid argument to StringToInputVariable.
IS_VALID_EXPRESSION( canParseVariable, T, stringToInputVariable( std::declval< T & >(), string(), Regex() ) );

/**
 * @brief Set @p lhs equal to @p rhs.
 * @tparam T The type of @p lhs and @p rhs.
 * @param lhs The value to set to @p rhs.
 * @param rhs The value to set @p lhs to.
 */
template< typename T >
static void equate( T & lhs, T const & rhs )
{ lhs = rhs; }

/**
 * @brief Set the entries of @p lhs equal to @p rhs.
 * @tparam T The type of the values in @p lhs and @p rhs.
 * @tparam NDIM The dimension of @p lhs.
 * @tparam PERM The permutation of @p rhs.
 * @param lhs The array of value to set to @p rhs.
 * @param rhs The value to set @p lhs to.
 */
template< typename T, int NDIM, typename PERM >
static void equate( Array< T, NDIM, PERM > const & lhs, T const & rhs )
{ lhs.template setValues< serialPolicy >( rhs ); }

}   // namespace internal

/**
 * @name Attribute extraction from XML nodes.
 */
///@{

/**
 * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
 *        This SFINAE implementation is used if the value is not parsable.
 * @tparam T             the type of variable fill with xml attribute.
 * @tparam U             type of the default value for @p rval
 * @param[in] name       the name of the xml attribute to process
 * @return false
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 */
template< typename T, typename U >
std::enable_if_t< !internal::canParseVariable< T >, bool >
readAttributeAsType( T &, string const & name, Regex const &, xmlNode const &, U const & )
{
  GEOS_THROW( "Cannot parse key with name ("<<name<<") with the given type " << LvArray::system::demangleType< T >(), InputError );
}

/**
 * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
 * @tparam T             the type of variable fill with xml attribute.
 * @tparam T_DEF         type of the default value for @p rval
 * @param[out] rval      the variable to fill with value
 * @param[in] name       the name of the xml attribute to process
 * @param[in] regex      the regular expression used for validating the string value.
 * @param[in] targetNode the xml node that should contain the attribute
 * @param[in] defVal     default value of @p rval (or entries of @p rval, if it is an array)
 * @return true
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 */
template< typename T, typename T_DEF = T >
std::enable_if_t< internal::canParseVariable< T >, bool >
readAttributeAsType( T & rval,
                     string const & name,
                     Regex const & regex,
                     xmlNode const & targetNode,
                     T_DEF const & defVal )
{
  xmlAttribute const xmlatt = targetNode.attribute( name.c_str() );
  if( !xmlatt.empty() )
  {
    // parse the string/attribute into a value
    stringToInputVariable( rval, xmlatt.value(), regex );
    return true;
  }
  else
  {
    // set the value to the default value
    internal::equate( rval, defVal );
    return false;
  }
}

/**
 * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
 * @tparam T             the type of variable fill with xml attribute.
 * @param[out] rval      the variable to fill with value
 * @param[in] name       the name of the xml attribute to process
 * @param[in] regex      the regular expression used for validating the string value.
 * @param[in] targetNode the xml node that should contain the attribute
 * @param[in] required   whether or not the value is required
 * @return boolean value indicating whether the value was successfully read from XML.
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 */
template< typename T >
std::enable_if_t< internal::canParseVariable< T >, bool >
readAttributeAsType( T & rval,
                     string const & name,
                     Regex const & regex,
                     xmlNode const & targetNode,
                     bool const required )
{
  xmlAttribute const xmlatt = targetNode.attribute( name.c_str() );

  bool const success = !(xmlatt.empty() && required);

  if( success )
  {
    // parse the string/attribute into a value
    stringToInputVariable( rval, xmlatt.value(), regex );
  }
  return success;
}

/**
 * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
 * @tparam T             the type of variable fill with xml attribute.
 * @param[out] rval      the variable to fill with value
 * @param[in] name       the name of the xml attribute to process
 * @param[in] regex      the regular expression used for validating the string value.
 * @param[in] targetNode the xml node that should contain the attribute
 * @return boolean value indicating whether the value was successfully read from XML.
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 */
template< typename T >
std::enable_if_t< !dataRepository::DefaultValue< T >::has_default_value, bool >
readAttributeAsType( T & rval,
                     string const & name,
                     Regex const & regex,
                     xmlNode const & targetNode,
                     dataRepository::DefaultValue< T > const & )
{
  return readAttributeAsType( rval, name, regex, targetNode, false );
}

/**
 * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
 * @tparam T             the type of variable fill with xml attribute.
 * @param[out] rval      the variable to fill with value
 * @param[in] name       the name of the xml attribute to process
 * @param[in] regex      the regular expression used for validating the string value.
 * @param[in] targetNode the xml node that should contain the attribute
 * @param[in] defVal     default value of @p rval (or entries of @p rval, if it is an array)
 * @return boolean value indicating whether the value was successfully read from XML.
 * @throw An InputError if the string value could not be validated by the type regex or parsed to the destination type.
 * In the context of xml parsing, use processInputException to format this error with the proper xml line.
 */
template< typename T >
typename std::enable_if_t< dataRepository::DefaultValue< T >::has_default_value, bool >
readAttributeAsType( T & rval,
                     string const & name,
                     Regex const & regex,
                     xmlNode const & targetNode,
                     dataRepository::DefaultValue< T > const & defVal )
{
  return readAttributeAsType( rval, name, regex, targetNode, defVal.value );
}

///@}

}

} /* namespace geos */

#endif /*GEOS_DATAREPOSITORY_XMLWRAPPER_HPP_ */
