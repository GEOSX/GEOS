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

/**
 * @file RTTypes.hpp
 *
 * This file contains various aliases and functions that provide operations regarding the
 * use of the runtime types.
 */

#ifndef GEOS_CODINGUTILITIES_RTTYPES_HPP
#define GEOS_CODINGUTILITIES_RTTYPES_HPP

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"

namespace geos
{

/**
 * @brief Perform a type cast of base to derived pointer.
 * @tparam NEW_TYPE      derived pointer type
 * @tparam EXISTING_TYPE base type
 * @param val            base pointer to cast
 * @return               pointer cast to derived type or @p nullptr
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE dynamicCast( EXISTING_TYPE * const val )
{
  static_assert( std::is_pointer< NEW_TYPE >::value, "NEW_TYPE must be a pointer." );
  return dynamic_cast< NEW_TYPE >( val );
}

/**
 * @brief Perform a type cast of base to derived reference.
 * @tparam NEW_TYPE      derived reference type
 * @tparam EXISTING_TYPE base type
 * @param val            base reference to cast
 * @return               reference cast to derived type or @p nullptr
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE dynamicCast( EXISTING_TYPE & val )
{
  static_assert( std::is_reference< NEW_TYPE >::value, "NEW_TYPE must be a reference." );

  using POINTER_TO_NEW_TYPE = std::remove_reference_t< NEW_TYPE > *;
  POINTER_TO_NEW_TYPE ptr = dynamicCast< POINTER_TO_NEW_TYPE >( &val );
  GEOS_ERROR_IF( ptr == nullptr, "Cast from " << LvArray::system::demangleType( val ) << " to " <<
                 LvArray::system::demangleType< NEW_TYPE >() << " failed." );

  return *ptr;
}

/**
 * @brief Print a short summary of a few select type aliases.
 */
void printTypeSummary();

/**
 * @brief The regular expression data for validating inputs. Use rtTypes to get the regex of a
 * type, and TypeRegex< T > to define a type regex.
 */
struct Regex
{
  /**
   * @brief the regular expression string.
   */
  string m_regexStr;
  /**
   * @brief the description of the expected format of the regular expression.
   */
  string m_formatDescription;
  /**
   * @brief Default constructor
   */
  Regex() {}
  /**
   * @param regexStr the regex string for validation (eg. "[\\d]+")
   * @param formatDescription the description of the expected format to be validated (eg. "Input value must be an integer.").
   */
  Regex( string_view regexStr, string_view formatDescription );
};

/**
 * @brief Extension point for custom types to provide a validation regexp to schema.
 * Do not use directly to obtain a type regex, rtTypes::getTypeRegex< T >() should be used instead.
 * @tparam T the type for which the regex is defined
 * @tparam ENABLE used to conditionally enable partial specializations
 *
 * Specializations should define the following method:
 * \code{cpp}
 *   static string get();
 * \endcode
 */
template< typename T, typename ENABLE = void >
struct TypeRegex
{
  /**
   * @brief Get the type's regex (default implementation returns nothing).
   * @return The Regex associated with T.
   */
  static Regex get() { return {}; }
};

/**
 * @brief Static class to manage the type selection of types at runtime and obtain the
 * regexes of these types. Typically, these types are 'xsd:simpleType' in the XSD.
 */
class rtTypes
{
public:

  /**
   * @brief the regex map type to store and find the regexes by the associated rtTypeName.
   */
  using RegexMapType = std::map< string, Regex >;

  /**
   * @brief Custom types are useful to customize the regexes of an existing type. The type name
   * can be one of the existing ones, or a totally new one (which can then be used in Wrapper::setRTTypename).
   */
  struct CustomTypes
  {
    /// @cond DO_NOT_DOCUMENT
    static constexpr string_view mapPair             = "mapPair";
    static constexpr string_view plotLevel           = "geos_dataRepository_PlotLevel";
    static constexpr string_view groupName           = "groupName";
    static constexpr string_view groupNameRef        = "groupNameRef";
    static constexpr string_view groupNameRefArray   = "groupNameRef_array";
    static constexpr string_view groupOfGroupNameRefArray   = "groupNameRef_array2d";
    /// @endcond
  };

  /**
   * @brief Convert a @p std::type_index to a string.
   * @param key the std::type_index of the type
   * @return a hard coded string that is related to the std::type_index
   */
  static string getTypeName( std::type_index const key );

  /**
   * @tparam T type we want the regex
   * @return the regex string for the default rtType of T to validate input values to this type.
   */
  template< typename T >
  static Regex const & getTypeRegex()
  { return getTypeRegex< T >( getTypeName( typeid( T ) ) ); }

  /**
   * @param typeName The rtType name of the type we want the regex (can be a CustomTypes).
   * @tparam T the type we want the regex. If none are available in createBasicTypesRegexMap(), one is
   * generated from TypeRegex< T >::get().
   * @return a regex string validating the type T.
   */
  template< typename T >
  static Regex const & getTypeRegex( string_view typeName )
  {
    RegexMapType & map = getTypeRegexMap();
    auto const it = map.find( string( typeName ) );
    if( it != map.end() )
    {
      return it->second;
    }
    else
    {
      return map.emplace( typeName, TypeRegex< T >::get() ).first->second;
    }
  }

  /**
   * @brief Construct the regexMap for all basic types (TypeRegex< T > extented types are not mentionned)
   * @return RegexMapType
   */
  static RegexMapType createBasicTypesRegexMap();

private:

  /**
   * @return A reference to the types regexes map
   */
  static RegexMapType & getTypeRegexMap()
  {
    static RegexMapType m = createBasicTypesRegexMap();
    return m;
  }

  /**
   * @brief Private constructor because of static class
   */
  rtTypes() {}

};

/**
 * @brief Utility class for querying type names at runtime.
 * @tparam T the target type
 *
 * This relies on LvArray's demangling facilities and simply
 * adds some convenience methods like getting the brief name.
 */
template< typename T >
struct TypeName
{
  /**
   * @brief @return Full name of the type.
   */
  static string full()
  {
    return ::LvArray::system::demangle( typeid( T ).name() );
  }

  /**
   * @brief @return brief name of the type (ignoring namespaces).
   */
  static string brief()
  {
    string const full_name = full();
    string::size_type const pos = full_name.find_last_of( "::" );
    return ( pos == string::npos ) ? full_name : full_name.substr( pos );
  }
};

}


#endif /* GEOS_CODINGUTILITIES_RTTYPES_HPP */
