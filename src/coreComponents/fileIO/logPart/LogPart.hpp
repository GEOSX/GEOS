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
 * @brief Class for displaying section for different steps of simulation
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
   * @brief Add a description to the opening logPart by concatening a description name and descriptions values.
   * @param name The description name
   * @param args Descriptions values to be concatened.
   * @note Descriptions values can be be any formatted types. Values will be aligned altogether.
   */
  template< typename ... Args >
  void addDescription( string const & name, Args const & ... args );

  /**
   * @brief Add a description to the opening logPart
   * @param description The string value of the description
   */
  void addDescription( string const & description );

  /**
   * @brief Add a description to the closing logPart by concatening a description name and descriptions values.
   * @param name The description name
   * @param args Descriptions values to be concatened.
   * @note Descriptions values can be be any formatted types. Values will be aligned altogether.
   */
  template< typename ... Args >
  void addEndDescription( string const & name, Args const & ... args );

  /**
   * @brief Add a description to the closing logPart
   * @param description The string value of the description
   */
  void addEndDescription( string const & description );

  /**
   * @brief Set the minimal width of a row
   * @param minWidth The minimal width of the table
   */
  void setMinWidth( size_t const & minWidth );

  /**
   * @brief Set the minimal width of a row
   * @param minWidth The minimal width of the table
   */
  void setMaxWidth( size_t const & maxWidth );

  /**
   * @brief Draw the first part of the logPart. It include the title and optionnaly, the end description(s).
   * @param os An output stream (by default, std::cout)
   */
  void begin( std::ostream & os = std::cout );

  /**
   * @brief Draw the last part of the logPart. It include the title
   * @param oss An output stream (by default, std::cout)
   */
  void end( std::ostream & oss = std::cout );

private:

  struct Description
  {
    /// Title footer string
    string m_title;

    /// Name of the description, formatted to be : [Name] : [Values1]\n[Values2]
    std::vector< string > m_descriptionNames;
    /// Values in the description
    std::vector< std::vector< string > > m_descriptionsValues;
    /// Vector containing the descriptions formatted by formatDescriptions()
    std::vector< string > m_formattedDescriptionLines;

    /// logPart length
    size_t m_logPartWidth;
    /// logPart length
    size_t m_logPartMaxWidth;
    /// min width of logPart length
    size_t m_logPartMinWidth;
    /// min width of logPart length
    size_t m_logPartMaxNameWidth;
  };

  /**
   * @brief Add a description to a specific section
   * @param name The description name
   * @param args Descriptions values to be concatened.
   * @note Descriptions values can be be any formatted types. Values will be aligned altogether.
   */
  template< typename ... Args >
  void addDescriptionBySection( LogPart::Description & description, string const & name, Args const &... args );

  /**
   * @brief Computes the initial log width based on the provided description and names.
   * @param description A description object that contains the current log part
   * @param descriptionNames A vector containing all name used in the log Part
   */
  void computeInitialLogWidth( LogPart::Description & description,
                               std::vector< string > const & descriptionNames,
                               std::vector< std::vector< string > > m_descriptionsValues );
  /**
   * @brief Build a description from the name and description values
   * @param name The decription name
   * @param decriptionsValues The description values
   */
  void formatDescriptions( LogPart::Description & description );

  /**
   * @brief Constructs the string logPart title of the log part.
   */
  string buildTitlePart( LogPart::Description const & description );

  /**
   * @brief Constructs the string logPart descriptions of the log part.
   */
  string buildDescriptionPart( LogPart::Description const & description );

  /// prefix to append to the title of closing section
  string const m_prefixEndTitle = "End of ";

  Description m_startDesc = { "", {}, {}, {}, 50, SIZE_MAX, 50, 0 };
  Description m_endDesc  = { "", {}, {}, {}, 50, SIZE_MAX, 50, 0 };

  /// description border margin
  static constexpr size_t m_borderMargin = 2;
  /// numbers of character used as border
  static constexpr size_t m_nbBorderChar = 2;
  /// character used for logPart construction
  static constexpr size_t m_borderCharacter = '#';

  /// String containing horizontal border
  string m_horizontalBorder;
};

template< typename ... Args >
void LogPart::addDescriptionBySection( Description & description, string const & name, Args const &... args )
{
  std::vector< string > values;
  size_t maxValueSize = 0;

  ( [&] {
    static_assert( has_formatter_v< decltype(args) >,
                   "Argument passed cannot be converted to string" );
    string const value = GEOS_FMT( "{}", args );
    values.push_back( value );
    maxValueSize = std::max( maxValueSize, value.size());
  } (), ...);

  description.m_descriptionNames.push_back( name );
  description.m_descriptionsValues.push_back( values );
  description.m_logPartWidth = std::min( description.m_logPartWidth, maxValueSize + 6 );

}

template< typename ... Args >
void LogPart::addDescription( string const & name, Args const &... args )
{
  addDescriptionBySection( m_startDesc, name, args ... );
}

template< typename ... Args >
void LogPart::addEndDescription( string const & name, Args const &... args )
{
  addDescriptionBySection( m_endDesc, name, args ... );
}

}

#endif
