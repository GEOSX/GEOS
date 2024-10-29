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

#ifndef GEOS_COMMON_FORMAT_HPP_
#define GEOS_COMMON_FORMAT_HPP_

#include <type_traits>

#if __cplusplus < 202002L
#define GEOS_USE_FMT
#endif

#ifdef GEOS_USE_FMT
#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY
#endif
// Differentiate between standalone fmt path and umpire's fmt path
#include "../include/fmt/core.h"
#include "../include/fmt/chrono.h"
#include "../include/fmt/ranges.h"
#include "../include/fmt/xchar.h"
#define GEOS_FMT_NS fmt
#else // use C++20's <format>
#include <format>
#define GEOS_FMT_NS std
#endif

#ifdef GEOS_USE_FMT
/**
 * @brief fmtlib formatter for enum classes.
 * @tparam T The type of the object being formatted. This should be an
 *           enum class.
 */
template< typename T >
struct fmt::formatter< T, std::enable_if_t< std::is_enum< T >::value > >
{
  /**
   * @brief Parser for the fmtlib formatting library.
   * @param ctx The context provided by the fmtlib library, which includes
   *            the format string.
   * @return An iterator pointing to the end of the format string.
   */
  template< typename ParseContext >
  constexpr auto parse( ParseContext & ctx )
  {
    return ctx.end();
  }

  /**
   * @brief Formatter for the fmtlib formatting library.
   * @param value The enum class object to format.
   * @param ctx   The context provided by the fmtlib library, which includes
   *              the output iterator where the formatted string should be written.
   * @return An iterator pointing to the end of the formatted string.
   */
  template< typename FormatContext >
  auto format( const T & value, FormatContext & ctx ) const
  {
    return fmt::format_to( ctx.out(), "{}", static_cast< std::underlying_type_t< T > >( value ) );
  }
};
#endif

/**
 * @brief Interpolate arguments into a message format string.
 * @param msg the message format string, must be a constant expression
 */
#define GEOS_FMT( msg, ... ) GEOS_FMT_NS::format( msg, __VA_ARGS__ )

/**
 * @brief Interpolate arguments into a message format string and write into an output iterator.
 * @param iter the output iterator to write to
 * @param size maximum number of characters to write
 * @param msg the message format string, must be a constant expression
 * @note Ensures the output buffer is zero-terminated (std::format_to_n doesn't)
 */
#define GEOS_FMT_TO( iter, size, msg, ... ) *GEOS_FMT_NS::format_to_n( iter, size - 1, msg, __VA_ARGS__ ).out = '\0'

// The following workaround is needed to fix compilation with NVCC on some PowerPC machines.
// The issue causes the following assertion error message:
// "Cannot format an argument. To make type T formattable provide a formatter<T> specialization"
// The standard definition of the has_const_formatter check of fmt fails due to a compiler bug, see the issue below:
// https://github.com/fmtlib/fmt/issues/2746
// The workaround was originally implemented in fmt:
// https://github.com/fmtlib/fmt/commit/70de324aa801eaf52e94c402d526a3849797c620
// but later removed:
// https://github.com/fmtlib/fmt/commit/466e0650ec2d153d255a40ec230eb77d7f1c3334
// This workaround bypasse the check for a const formatter whenever the foramt context GEOS_FMT_NS::format_context is used
#ifdef GEOS_USE_FMT_CONST_FORMATTER_WORKAROUND
template< >
constexpr auto GEOS_FMT_NS::detail::has_const_formatter_impl< GEOS_FMT_NS::format_context >( ... ) -> bool
{
  return true;
}
#endif // End of the workaround for fmt compilation issues

/**
 * Evaluates at compile time if a fmt::formatter exists for a given type
 */
#if __cplusplus < 202002L
template< class T >
static constexpr bool has_formatter_v = fmt::has_formatter< fmt::remove_cvref_t< T >, fmt::format_context >();
#else
template< typename T >
concept has_formatter_v = requires ( T& v, std::format_context ctx )
{
  std::formatter< std::remove_cvref_t< T > >().format( v, ctx );
};
#endif

#endif //GEOS_COMMON_FORMAT_HPP_
