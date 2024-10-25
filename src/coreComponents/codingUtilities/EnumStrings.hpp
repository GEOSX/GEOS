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
 * @file EnumStrings.hpp
 *
 * Collection of utilities to facilitate I/O of enumeration types.
 * Provides a macro definition that allows associating string names
 * with enumeration constants and a set of functions that make use
 * of these strings, like stream insertion/extraction operators.
 */

#ifndef GEOS_CODINGUTILITIES_ENUMSTRINGS_HPP
#define GEOS_CODINGUTILITIES_ENUMSTRINGS_HPP

#include "common/format/StringUtilities.hpp"
#include "codingUtilities/RTTypes.hpp"
#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"
#include "common/format/Format.hpp"

#include <iostream>
#include <type_traits>
#include <algorithm>

namespace geos
{

namespace internal
{
/**
 * @brief Simple compile-time variadic function that counts its arguments.
 * @tparam ARGS variadic pack of argument types
 * @return the number of arguments passed
 */
template< typename ... ARGS >
constexpr int countArgs( ARGS ... )
{
  return sizeof...( ARGS );
}
}

/**
 * @brief Associate a list of string names with enumeration values.
 * @param ENUM the enumeration type
 * @param ... list of names (C-string literals)
 *
 * Conditions (not enforced but won't work correctly if violated):
 *  - the macro must be called in the same namespace the enumeration type is defined in
 *  - the number and order of string arguments passed must match the enum values
 *  - enumeration constants must not have custom values assigned
 *
 * After the macro has been called, template instantiation EnumStrings<ENUM>
 * may be used to get access to strings at runtime. While not strictly necessary,
 * it is recommended that macro call immediately follows the enum definition
 * (or the class definition, if enum is defined inside a class).
 */
#define ENUM_STRINGS( ENUM, ... )                                     \
  inline auto const & getEnumStrings( ENUM const )                    \
  {                                                                   \
    static constexpr char const * ss[] { __VA_ARGS__ };               \
    return ss;                                                        \
  }                                                                   \
                                                                      \
  inline std::ostream & operator<<( std::ostream & os, ENUM const e ) \
  {                                                                   \
    os << EnumStrings< ENUM >::toString( e );                         \
    return os;                                                        \
  }                                                                   \
                                                                      \
  inline std::istream & operator>>( std::istream & is, ENUM & e )     \
  {                                                                   \
    string s; is >> s;                                                \
    e = EnumStrings< ENUM >::fromString( s );                         \
    return is;                                                        \
  }                                                                   \
                                                                      \
  inline string toString( ENUM const e )                              \
  {                                                                   \
    return EnumStrings< ENUM >::toString( e );                        \
  }                                                                   \
                                                                      \
  static_assert( std::is_enum< ENUM >::value, "Not an enumeration" )

/**
 * @brief Provides enum <-> string conversion facilities.
 * @tparam ENUM the enumeration type
 */
template< typename ENUM >
struct EnumStrings
{
  /// Alias for the enumeration type
  using enum_type = ENUM;

  /// Alias for enum's underlying fundamental type
  using base_type = std::underlying_type_t< ENUM >;

  /**
   * @brief @return An array of strings associated with enumeration.
   */
  static auto const & get()
  {
    return getEnumStrings( enum_type{} ); // invoke ADL
  }

  /**
   * @brief Get a list of valid options as a delimited string.
   * @param delim delimiter (defaults to single space)
   * @return the string containing all valid strings for this type
   */
  static string concat( string const & delim = " " )
  {
    auto const & strings = get();
    return stringutilities::join( std::begin( strings ), std::end( strings ), delim );
  }

  /**
   * @brief Convert enum to string.
   * @param e the enum value to convert
   * @return the corresponding string
   *
   * An error is raised if enum's numerical value is greater of equal than the number of strings.
   */
  static string toString( enum_type const & e )
  {
    auto const & strings = get();
    std::size_t size = std::distance( std::begin( strings ), std::end( strings ) );
    base_type const index = static_cast< base_type >( e );
    GEOS_THROW_IF( index >= LvArray::integerConversion< base_type >( size ),
                   "Invalid value " << index << " of type " << TypeName< ENUM >::brief() << ". Valid range is 0.." << size - 1,
                   InputError );
    return strings[ index ];
  }

  /**
   * @brief Convert string to enum
   * @param s the string to convert
   * @return the corresponding enum value
   */
  static enum_type fromString( string const & s )
  {
    auto const & strings = get();
    auto const it = std::find( std::begin( strings ), std::end( strings ), s );
    GEOS_THROW_IF( it == std::end( strings ),
                   "Invalid value '" << s << "' of type " << TypeName< enum_type >::brief() << ". Valid options are: " << concat( ", " ),
                   InputError );
    enum_type const e = static_cast< enum_type >( LvArray::integerConversion< base_type >( std::distance( std::begin( strings ), it ) ) );
    return e;
  }
};

namespace internal
{
IS_VALID_EXPRESSION( HasEnumStrings, ENUM, getEnumStrings( std::declval< ENUM >() ) );
}

/**
 * @brief Specialization of TypeRegex for enumeration types with strings attached (pun intended).
 * @tparam ENUM the type of enumeration
 */
template< typename ENUM >
struct TypeRegex< ENUM, std::enable_if_t< internal::HasEnumStrings< ENUM > > >
{
  /**
   * @brief @return Regex for validating enumeration inputs for @p ENUM type.
   */
  static Regex get()
  {
    return Regex( EnumStrings< ENUM >::concat( "|" ),
                  "Input value must be one of { " + EnumStrings< ENUM >::concat( ", " ) + " }." );
  }
};

} // namespace geos

// Formatter specialization for enums
template< typename Enum >
struct GEOS_FMT_NS::formatter< Enum, std::enable_if_t< std::is_enum< Enum >::value && geos::internal::HasEnumStrings< Enum >, char > >
  : GEOS_FMT_NS::formatter< std::string >
{
  template< typename FormatContext >
  auto format( Enum e, FormatContext & ctx ) const
  {
    return formatter< std::string >::format( toString( e ), ctx );
  }
};

// Formatter specialization for enums
template< typename Enum >
struct GEOS_FMT_NS::formatter< Enum, std::enable_if_t< std::is_enum< Enum >::value && !geos::internal::HasEnumStrings< Enum >, char > >
  : GEOS_FMT_NS::formatter< std::underlying_type_t< Enum > >
{
  template< typename FormatContext >
  auto format( Enum e, FormatContext & ctx ) const
  {
    return GEOS_FMT_NS::formatter< std::underlying_type_t< Enum > >::format( toUnderlying( e ), ctx );
  }
};

#endif //GEOS_CODINGUTILITIES_ENUMSTRINGS_HPP
