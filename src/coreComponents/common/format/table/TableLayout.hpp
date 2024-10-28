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
#include <variant>

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
   * @brief Structure to set up values alignment for each colum.
   */
  struct ColumnAlignment
  {
    /// Alignment for column name. By default aligned to center
    Alignment headerAlignment = Alignment::center;
    /// Alignment for column values. By default aligned to right side
    Alignment valueAlignment = Alignment::right;
  };

  /**
   * @brief Structure to set up each colum.
   */
  struct Column
  {
    /// Name for a column
    string columnName;
    /// Alignment for  column name. By default aligned to center
    ColumnAlignment alignmentSettings;
    /// A boolean to display a colummn
    bool enabled = true;
    /// Vector containing sub columns name
    std::vector< string > subColumnNames;
    /// Vector containing substring column name delimited by "\n"
    std::vector< string > splitColumnNames = {""};

    /**
     * @brief Construct a Column object with the given parameters
     * @param name The Column name
     */
    Column( string_view name )
      : columnName( name )
    {}

    /**
     * @brief Construct a Column object with the given parameters
     * @param name The Column name
     * @param headerAlign The column name alignment
     */
    Column( string_view name, Alignment headerAlign )
      : columnName( name ), alignmentSettings( {headerAlign, Alignment::right} )
    {}

    /**
     * @brief Construct a Column object with the given parameters
     * @param name The Column name
     * @param subColumnNamesInit Vector containing subcolumn names
     */
    Column( string_view name, std::vector< string > && subColumnNamesInit )
      : columnName( name ), subColumnNames( subColumnNamesInit )
    {}

    /**
     * @brief Construct a Column object with the given parameters
     * @param name The Column name
     * @param headerAlign The column name alignment
     * @param display Flag indicating whether the column is enabled
     */
    Column( string_view name, Alignment headerAlign, bool display )
      : columnName( name ),
      alignmentSettings( {headerAlign, Alignment::right} ),
      enabled( display )
    {}

    /**
     * @brief Construct a Column object with the given parameters
     * @param name The Column name
     * @param headerAlign The column name alignment
     * @param display Flag indicating whether the column is enabled
     * @param subColumnNamesInit Vector containing subcolumn names
     */
    Column( string_view name, Alignment headerAlign, bool display, std::vector< string > && subColumnNamesInit )
      : columnName( name ),
      alignmentSettings( {headerAlign, Alignment::right} ),
      enabled( display ),
      subColumnNames( subColumnNamesInit )
    {}
  };

  /**
   * @brief Struct for a TableColumnData.
   * Each column contains its own parameters (such as name, alignment, etc.).
   */
  struct TableColumnData
  {
    /// Structure who contains parameters for a column
    Column column;
    /// A vector containing all the values of a column
    std::vector< string > columnValues;
    /// The largest string(s) in the column
    std::vector< string > maxStringSize;
    /// Vector containing all sub columns subdivison
    std::vector< TableColumnData > subColumn;

    /**
     * @brief Constructs a TableColumnData  with the given column.
     * @param col The parameters for the column.
     */
    TableColumnData ( Column col )
      : column( col )
    {}

    /**
     * @brief Constructs a TableColumnData with the given parameters.
     * @param col The parameters for the column.
     * @param subColumnInit The subcolumns contained in the colum
     */
    TableColumnData ( Column col, std::vector< TableColumnData > subColumnInit )
      : column( col ), subColumn( subColumnInit )
    {}

    /**
     * @brief Constructs a TableColumnData  with the given parameters.
     * @param col The parameters for the column.
     * @param columnValuesInit The values in the column.
     * @param maxStringSizeInit The largest string(s) in the column.
     * @param subColumnInit The sub-columns of the column.
     */
    TableColumnData ( Column col,
                      std::vector< string > columnValuesInit,
                      std::vector< string > maxStringSizeInit,
                      std::vector< TableColumnData > subColumnInit )
      : column( col ),
      columnValues( columnValuesInit ),
      maxStringSize( maxStringSizeInit ),
      subColumn( subColumnInit )
    {}

  };

  /// Alias for an initializer list of variants that can contain either a string or a layout column.
  using TableLayoutArgs = std::initializer_list< std::variant< string_view, TableLayout::Column > >;

  TableLayout() = default;

  /**
   * @brief Construct a new Table Layout object
   * @param title The table title
   * @param columns A vector containing all column initialized
   */
  TableLayout( string_view title,
               std::vector< TableLayout::Column > & columns )
  {
    setMargin( MarginValue::medium );
    setTitle( title );
    for( auto const & column :columns )
    {
      addToColumns( column );
    }
  }

  /**
   * @brief Construct a new Table Layout object
   * @param title The table title
   * @param args An initializer_list containing string / column
   */
  TableLayout( string_view title,
               TableLayoutArgs args )
  {
    setMargin( MarginValue::medium );
    setTitle( title );
    processArguments( args );
  }

  /**
   * @brief Construct a new Table Layout object
   * @param args An initializer_list containing string / column
   */

  TableLayout( TableLayoutArgs args )
  {
    setMargin( MarginValue::medium );
    processArguments( args );
  }

  /**
   * @brief Construct a new Table Layout object
   * @param title The table title
   * @param args An initializer_list containing string / column
   */
  TableLayout( string_view title,
               std::vector< string > args )
  {
    setMargin( MarginValue::medium );
    setTitle( title );
    addToColumns( args );
  }

  /**
   * @return The columns vector
   */
  std::vector< TableColumnData > const & getColumns() const;

  /**
   * @return The table name
   */
  string_view getTitle() const;

  /**
   * @param title The table title
   * @return The tableLayout reference
   */
  TableLayout & setTitle( string_view title );

  /**
   * @brief Remove the last return line a the end of the table
   * @return The tableLayout reference
   */
  TableLayout & disableLineWrap();

  /**
   * @brief Set the minimal margin width between cell content and borders.
   * @param marginValue The margin value
   * @return The tableLayout reference
   */
  TableLayout & setMargin( MarginValue marginValue );

  /**
   * @brief Set the column values alignment
   * @param alignment column values alignment
   * @return The tableLayout reference
   */
  TableLayout & setValuesAlignment( TableLayout::Alignment alignment );

  /**
   * @return whether we have a line return at the end of the table or not
   */
  bool isLineWrapEnabled() const;

  /**
   * @brief Remove all subcolumn in all columns
   * Can be used if we want to reuse a TableLayout without keep subcolumns
   */
  void removeSubColumn();

  /**
   * @return The border margin, number of spaces at both left and right table sides plus vertical character
   */
  integer const & getBorderMargin() const;

  /**
   * @return The column margin, numbers of spaces separating both left and right side from a vertical line
   */
  integer const & getColumnMargin() const;

  /**
   * @return The table margin value
   */
  integer const & getMarginValue() const;

  /**
   * @return The margin title
   */
  integer const & getMarginTitle() const;

private:

  /**
   * @brief Add a column to the table given an initializer_list of string & TableColumnData
   * @param args An initializer_list containing string / column
   */
  void processArguments( TableLayoutArgs args )
  {
    for( auto const & arg : args )
    {
      std::visit( [this]( auto const & value ) {
        addToColumns( value );
      }, arg );
    }
  }

  /**
   * @brief Recursively processes a variable number of arguments and adds them to the table data.
   * @tparam Ts The remaining arguments
   * @param args The remaining arguments to be processed
   */
  template< typename ... Ts >
  void processArguments( Ts &... args )
  {
    addToColumns( args ... );
  }

  /**
   * @brief Create and add columns to the columns vector given a string vector
   * @param columnNames The columns name
   */
  void addToColumns( std::vector< string > const & columnNames );

  /**
   * @brief Create and add a column to the columns vector given a string
   * @param columnName The column name
   */
  void addToColumns( string_view columnName );

/**
 * @brief Create and add a column to the columns vector given a Column
 * @param column Vector containing addition information on the column
 */
  void addToColumns( Column const & column );

  std::vector< TableColumnData > m_tableColumnsData;

  bool m_wrapLine = true;
  string m_tableTitle;
  integer m_borderMargin;
  integer m_columnMargin;
  integer m_marginValue;
  integer m_titleMargin = 2;

};
}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP */
