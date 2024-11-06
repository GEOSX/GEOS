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
 * @file LogPart.hpp
 */

#ifndef GEOS_COMMON_SECTION_HPP
#define GEOS_COMMON_SECTION_HPP

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"

namespace geos
{

/**
 * @brief Class for displaying different steps of simulation
 */
class LogPart
{
public:

  /**
   * @brief Construct a new LogPart
   * @param m_logPartTitle The logPart title
   */
  LogPart( string_view m_logPartTitle );

  /**
   * @brief Add a description to the logPart by concatening a description name and descriptions values.
   * @param name The description name
   * @param args Descriptions values to be concatened.
   * Descriptions values can be be any types and will be aligned
   */
  template< typename ... Args >
  void addDescription( string_view name, Args const & ... args );

  /**
   * @brief Add a description to the logPart
   * @param description The string value of the description
   */
  void addDescription( string_view description );

  /**
   * @brief Add a description to the end of the logPart by concatening a description name and descriptions values.
   * @param name The description name
   * @param args Descriptions values to be concatened.
   * Descriptions values can be be any types and will be aligned
   */
  template< typename ... Args >
  void addEndDescription( string_view name, Args const & ... args );

  /**
   * @brief Add a description to the end of the logPart
   * @param description The string value of the description
   */
  void addEndDescription( string_view description );

  /**
   * @brief Set the minimal width of a row
   * @param minWidth The minimal width of the table
   */
  void setMinWidth( integer const & minWidth );

  /**
   * @brief Draw the first part of the logPart. It include the title and optionnaly, the end description(s).
   * @param os An output stream (by default, std::cout)
   */
  void begin( std::ostream & os = std::cout ) const;

  /**
   * @brief Draw the last part of the logPart. It include the title
   * @param oss An output stream (by default, std::cout)
   */
  void end( std::ostream & oss = std::cout ) const;

private:

  /**
   * @brief Build a description from the name and description values
   * @param name The decription name
   * @param decriptionsValues The description values
   */
  void formatAndInsertDescriptions( std::vector< string > & formattedDescriptions,
                                    string_view name,
                                    std::vector< string > const & decriptionsValues );

  /**
   * @brief Constructs the string logPart title of the log part.
   * @param title The title to be set
   */
  string buildTitlePart( string_view title ) const;

  /**
   * @brief Constructs the string logPart descriptions of the log part.
   * @param descriptions The description to be formatted
   */
  string buildDescriptionPart( std::vector< string > const & formattedDescriptions ) const;

  /// Vector containing all descriptions
  std::vector< string > m_beginningDescs;
  /// title of logPart
  std::vector< string > m_endDescs;

  /// title of logPart
  string m_logPartTitle;
  /// Start title footer string
  string m_footerTitle;
  /// min width of logPart length
  integer m_rowMinWidth = 70;
  /// logPart length
  integer m_logPartWidth;

  /// description border margin
  static constexpr integer m_marginBorder = 2;
  /// numbers of character used as border
  static constexpr integer m_nbBorderChar = 2;
  /// character used for logPart construction
  static constexpr integer m_borderCharacter = '#';

  /// String containing horizontal border
  string m_horizontalBorder;
};

template< typename ... Args >
void LogPart::addDescription( string_view name, Args const &... args )
{
  std::vector< string > values;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) >, "Argument passed in addRow cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );
    values.push_back( value );
  } (), ...);

  formatAndInsertDescriptions( m_beginningDescs, name, values );
}

template< typename ... Args >
void LogPart::addEndDescription( string_view name, Args const &... args )
{
  std::vector< string > values;
  ( [&] {
    static_assert( has_formatter_v< decltype(args) >, "Argument passed in addRow cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );
    values.push_back( value );
  } (), ...);

  formatAndInsertDescriptions( m_endDescs, name, values );
}
}

#endif
