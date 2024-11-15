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
 * @file TableLayout.hpp
 */

#ifndef GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP
#define GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @brief Class for setup the table layout
 */
class TableLayout
{

public:

  /// Type of aligment for a column
  enum Alignment { right, left, center };

  /// Space to apply between all data and border
  enum MarginValue : integer
  {
    tiny = 0,
    small = 1,
    medium = 2,
    large = 3
  };

  /**
   * @brief Enumeration for table sections.
   */
  enum Section { header, values };

  /**
   * @brief Structure to set up each colum parameters.
   */
  struct ColumnParam
  {
    /// Name for a column
    string columnName;
    /// Alignment for a column. By default aligned to the right side
    Alignment alignment = Alignment::right;
    /// A boolean to display a colummn
    bool enabled = true;
    /// Vector containing substring column name delimited by "\n"
    std::vector< string > splitColumnNameLines;

    /**
     * @brief Construct a ColumnParam object with the specified name and alignment.
     * @param name The name of the column
     * @param align The alignment of the column
     */
    ColumnParam( std::string const & name, Alignment align )
      : columnName( name ), alignment( align )
    {}

    /**
     * @brief Construct a ColumnParam object with the specified name, alignment, and display flag.
     * @param name The name of the column
     * @param align The alignment of the column
     * @param display Flag indicating whether the column is enabled
     */
    ColumnParam( std::string const & name, Alignment align, bool display )
      : columnName( name ), alignment( align ), enabled( display )
    {}
  };

  /**
   * @brief Struct for a column.
   */
  struct Column
  {
    /// Structure who contains parameters for a column
    ColumnParam m_parameter;
    /// A vector containing all column values
    std::vector< string > m_columnValues;
    /// The largest string in the column
    string m_maxStringSize;
  };

  TableLayout() = default;

  /**
   * @brief Construct a new Table object, all values in the table are centered by default
   * @param columnNames The names of the columns
   * @param title The table name
   */
  TableLayout( std::vector< string > const & columnNames, string const & title = "" );

  /**
   * @brief Construct a new Table object by specifying value alignment and optionally their displays based to log levels
   * level
   * @param columnParameters List of structures to set up each colum parameters.
   * @param title The table name
   */
  TableLayout( std::vector< ColumnParam > const & columnParameters, string const & title = "" );

  /**
   * @return The columns vector
   */
  std::vector< Column > const & getColumns() const;

  /**
   * @return The table name
   */
  string_view getTitle() const;

  /**
   * @return The border margin, number of spaces at both left and right table sides
   */
  integer const & getBorderMargin() const;

  /**
   * @return The column margin, numbers of spaces separating both left and right side from each column content
   */
  integer const & getColumnMargin() const;

  /**
   * @return The margin title
   */
  integer const & getMarginTitle() const;

  /**
   * @brief Set the minimal margin width between cell content and borders.
   * @param marginValue The margin value
   */
  void setMargin( MarginValue marginValue );

private:

  std::vector< Column > m_columns;
  string m_tableTitle;
  integer m_borderMargin;
  integer m_columnMargin;
  integer m_titleMargin = 2;

};
}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP */
