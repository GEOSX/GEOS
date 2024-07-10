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

/**
 * @file DataTypes.cpp
 */


#include "DataTypes.hpp"
#include "Logger.hpp"
#include "LvArray/src/system.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geos
{
#ifdef GEOSX_USE_MPI
MPI_Comm MPI_COMM_GEOSX;
#else
int MPI_COMM_GEOSX = 0;
#endif


Regex::Regex( string_view regexStr, string_view formatDescription ):
  m_regexStr( regexStr ),
  m_formatDescription( formatDescription )
{}

void printTypeSummary()
{
  GEOS_LOG_RANK_0( "real64 is alias of " <<LvArray::system::demangle( typeid(real64).name() ) );
  GEOS_LOG_RANK_0( "localIndex is alias of " <<LvArray::system::demangle( typeid(localIndex).name() ) );
  GEOS_LOG_RANK_0( "globalIndex is alias of "<<LvArray::system::demangle( typeid(globalIndex).name()) );
}

string rtTypes::getTypeName( std::type_index const key )
{
  static const std::unordered_map< std::type_index, string > type_names =
  {
    {std::type_index( typeid(integer)), "integer"},
    {std::type_index( typeid(real32)), "real32"},
    {std::type_index( typeid(real64)), "real64"},
    {std::type_index( typeid(localIndex)), "localIndex"},
    {std::type_index( typeid(globalIndex)), "globalIndex"},
    {std::type_index( typeid(R1Tensor)), "R1Tensor"},
    {std::type_index( typeid(R1Tensor32)), "R1Tensor32"},
    {std::type_index( typeid(R2SymTensor)), "R2SymTensor"},
    {std::type_index( typeid(integer_array)), "integer_array"},
    {std::type_index( typeid(real32_array)), "real32_array"},
    {std::type_index( typeid(real64_array)), "real64_array"},
    {std::type_index( typeid(localIndex_array)), "localIndex_array"},
    {std::type_index( typeid(globalIndex_array)), "globalIndex_array"},
    {std::type_index( typeid(integer_array2d)), "integer_array2d"},
    {std::type_index( typeid(real32_array2d)), "real32_array2d"},
    {std::type_index( typeid(real64_array2d)), "real64_array2d"},
    {std::type_index( typeid(localIndex_array2d)), "localIndex_array2d"},
    {std::type_index( typeid(globalIndex_array2d)), "globalIndex_array2d"},
    {std::type_index( typeid(integer_array3d)), "integer_array3d"},
    {std::type_index( typeid(real32_array3d)), "real32_array3d"},
    {std::type_index( typeid(real64_array3d)), "real64_array3d"},
    {std::type_index( typeid(localIndex_array3d)), "localIndex_array3d"},
    {std::type_index( typeid(globalIndex_array3d)), "globalIndex_array3d"},
    {std::type_index( typeid(real64_array4d)), "real64_array4d"},
    {std::type_index( typeid(string)), "string"},
    {std::type_index( typeid(Path)), "path"},
    {std::type_index( typeid(string_array)), "string_array"},
    {std::type_index( typeid(path_array)), "path_array"},
  };

  // If the data type is not defined here, return type_info.name()
  auto const iter = type_names.find( key );
  if( iter != type_names.end() )
  {
    return iter->second;
  }
  else
  {
    return LvArray::system::demangle( key.name());
  }
}

/**
 * @brief Recursive function to build array regexes.
 * @param subPattern pattern of the element to surround with braces and separate with commas
 * @param dimension 1 = bottom-level, 2 = array1d-level, 3 = array2d-level...
 * @param topLevelCall True if this is the first recursive call.
 * @return
 *
 * @note The sub pattern is the base object you are targeting.  It can either
 *       be a simple type or a lower-dimensional array. Sub-elements and
 *       axes are given as a comma-separated list enclosed in a curly brace.
 *       For example, a 2D string array would look like: {{"a", "b"}, {"c", "d"}}
 */
string constructArrayRegex( string_view subPattern, integer dimension, bool topLevelCall = true )
{
  string subPatternStr = dimension > 1 ?
                         constructArrayRegex( subPattern, dimension-1, false ) :
                         string( subPattern );
  // Add trailing space if is not already done
  if( !stringutilities::endsWith( subPatternStr, "\\s*" ) )
    subPatternStr+="\\s*";
  // Allow the bottom-level to be empty
  string const arrayRegex = dimension == 1 ?
                            "\\{\\s*((" + subPatternStr + ",\\s*)*" + subPatternStr + ")?\\}":
                            "\\{\\s*(" + subPatternStr + ",\\s*)*" + subPatternStr + "\\}";
  // accept spaces around surrounding braces at the top-level
  return topLevelCall ?
         "\\s*" + arrayRegex + "\\s*" :
         arrayRegex;
}
/**
 * @brief function to build array regexes.
 * @param subPattern pattern of the element to surround with braces and separate with commas
 * @param description description of the subPattern that starts by "Input value must "
 * @param dimension 1 = array1d, 2 = array2d...
 * @return
 */
Regex constructArrayRegex( string_view subPattern, string_view description, integer dimension )
{
  std::ostringstream arrayDesc;

  // Adapt the description so the form "Input value must be an int" is transformed to "Input value must be a 1d array. Each value must be an
  // int"
  {
    arrayDesc << "Input value must be a " << dimension << "d array (surrounded by ";
    if( dimension > 1 )
    {
      arrayDesc << dimension << " levels of ";
    }
    arrayDesc << "braces and separated by commas). Each value must ";

    // finish by the original description
    GEOS_ERROR_IF( !stringutilities::startsWith( description, "Input value must " ),
                   "Description \"" << description << "\" must start by \"Input value must \" to call constructArrayRegex() on it." );
    arrayDesc << description.substr( description.find( " must " ) );
  }

  return Regex( constructArrayRegex( subPattern, dimension ),
                arrayDesc.str() );
}

rtTypes::RegexMapType rtTypes::createBasicTypesRegexMap()
{
  // Define the component regexes:

  // Regex to match an unsigned int (123, etc.)
  // string_view const ru = "[\\d]+";// unused

  string_view const intDesc = "Input value must be a signed int (eg. -123, 455, +789, etc.)";
  string_view const intRegex = "[+-]?[\\d]+";

  // Explanation of parts:
  // [+-]?[\\d]*  matches an optional +/- at the beginning, any numbers preceding the decimal
  // ([\\d]\\.?|\\.[\\d]) matches the decimal region of the number (0, 1., 2.3, .4)
  // [\\d]*  matches any number of numbers following the decimal
  // ([eE][-+]?[\\d]+|\\s*)  matches an optional scientific notation number
  // Note: the xsd regex implementation does not allow an empty branch, so use allow whitespace at the end
  string_view const realDesc = "Input value must be a real number (eg. 1, .25, +2.3, -.4, 5.6e7, -8E-9, etc.)";
  string_view const realRegex = "[+-]?[\\d]*([\\d]\\.?|\\.[\\d])[\\d]*([eE][-+]?[\\d]+|\\s*)";

  string_view const R1Desc = "Input value must be a R1Tensor, an array of 3 real numbers surrounded by braces and separated by commas (eg.  \"{ 1, .25, +2.3}\", \"{ -.4, 5.6e7, -8E-9\", etc.) .";
  string const R1Regex = "\\s*\\{\\s*(" + string( realRegex ) + "\\s*,\\s*){2}" + string( realRegex ) + "\\s*\\}\\s*";
  string_view const R2Desc = "Input value must be a R2SymTensor, an array of 6 real numbers surrounded by braces and separated by commas (eg.  \"{ 1, .25, +2.3, -.4, 5.6e7, -8E-9\", etc.) .";
  string const R2Regex = "\\s*\\{\\s*(" + string( realRegex ) + "\\s*,\\s*){5}" + string( realRegex ) + "\\s*\\}\\s*";

  string_view const strDesc = "Input value must be a string that cannot be empty, contain any whitespaces nor the characters  , { }";
  string_view const strRegex = "[^,\\{\\}\\s]+\\s*";
  string_view const strEDesc = "Input value must be a string that cannot contain any whitespaces nor the characters  , { }";
  string_view const strERegex = "[^,\\{\\}\\s]*\\s*";

  string_view const pathDesc = "Input value must be a string that cannot be empty, contain any whitespaces nor the characters  * ? < > | : \" ";
  string_view const pathRegex = "[^*?<>\\|:\";,\\s]+\\s*";
  string_view const pathEDesc = "Input value must be a string that cannot contain any whitespaces nor the characters  * ? < > | : \" ";
  string_view const pathERegex = "[^*?<>\\|:\";,\\s]*\\s*";

  string_view const groupNameDesc = "Input value must be a string that cannot be empty and contains only upper/lower letters, digits, and the characters  . - _";
  string_view const groupNameRegex = "[a-zA-Z0-9.\\-_]+";
  string_view const groupNameRefDesc = "Input value must be a string that can contain only upper/lower letters, digits, and the characters  . - _ / *";
  string_view const groupNameRefRegex = "[a-zA-Z0-9.\\-_/*]*";


  // Build master list of regexes
  RegexMapType regexMap =
  {
    { "integer", Regex( intRegex, intDesc ) },
    { "localIndex", Regex( intRegex, intDesc ) },
    { "globalIndex", Regex( intRegex, intDesc ) },
    { "real32", Regex( realRegex, realDesc ) },
    { "real64", Regex( realRegex, realDesc ) },
    { "R1Tensor", Regex( R1Regex, R1Desc ) },
    { "R1Tensor32", Regex( R1Regex, R1Desc ) },
    { "R2SymTensor", Regex( R2Regex, R2Desc ) },
    { "integer_array", constructArrayRegex( intRegex, intDesc, 1 ) },
    { "localIndex_array", constructArrayRegex( intRegex, intDesc, 1 ) },
    { "globalIndex_array", constructArrayRegex( intRegex, intDesc, 1 ) },
    { "real32_array", constructArrayRegex( realRegex, realDesc, 1 ) },
    { "real64_array", constructArrayRegex( realRegex, realDesc, 1 ) },
    { "integer_array2d", constructArrayRegex( intRegex, intDesc, 2 ) },
    { "localIndex_array2d", constructArrayRegex( intRegex, intDesc, 2 ) },
    { "globalIndex_array2d", constructArrayRegex( intRegex, intDesc, 2 ) },
    { "real32_array2d", constructArrayRegex( realRegex, realDesc, 2 ) },
    { "real64_array2d", constructArrayRegex( realRegex, realDesc, 2 ) },
    { "integer_array3d", constructArrayRegex( intRegex, intDesc, 3 ) },
    { "localIndex_array3d", constructArrayRegex( intRegex, intDesc, 3 ) },
    { "globalIndex_array3d", constructArrayRegex( intRegex, intDesc, 3 ) },
    { "real32_array3d", constructArrayRegex( realRegex, realDesc, 3 ) },
    { "real64_array3d", constructArrayRegex( realRegex, realDesc, 3 ) },
    { "real64_array4d", constructArrayRegex( realRegex, realDesc, 4 ) },
    { "string", Regex( strERegex, strEDesc ) },
    { "path", Regex( pathERegex, pathEDesc ) },
    { "string_array", constructArrayRegex( strRegex, strDesc, 1 ) },
    { "path_array", constructArrayRegex( pathRegex, pathDesc, 1 ) },

    { string( CustomTypes::mapPair ), Regex( strERegex, strEDesc ) },
    { string( CustomTypes::plotLevel ), Regex( intRegex, intDesc ) },
    { string( CustomTypes::groupName ), Regex( groupNameRegex, groupNameDesc ) },
    { string( CustomTypes::groupNameRef ), Regex( groupNameRefRegex, groupNameRefDesc ) },
    { string( CustomTypes::groupNameRefArray ), constructArrayRegex( groupNameRefRegex, groupNameRefDesc, 1 ) }
  };
  return regexMap;
}



}
