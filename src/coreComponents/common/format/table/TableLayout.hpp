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
  struct CellAlignment
  {
    /// Alignment for column name. By default aligned to center
    Alignment headerAlignment;
    /// Alignment for column values. By default aligned to right side
    Alignment valueAlignment;
  };

  struct Cell
  {
    char type;
    string value;
    /// Vector containing sub columns name
    std::vector< string > dividedValue;
    size_t nbRows;

    Cell()
      : type( ' ' ), value( "" ), nbRows( 1 )
    {}

    Cell( string_view val )
      : value( val )
    {}

    Cell( char t, string_view val )
      : type( t ), value( val )
    {}

    Cell( char t, string_view val, std::vector< std::string > const & subCols )
      : type( t ), value( val ), dividedValue( subCols )
    {}
  };

  /**
   * @brief Struct for a Column.
   * Each column contains its own parameters (such as name, alignment, etc.).
   */
  struct Column
  {
    /// Name for a column
    Cell columnName;
    /// A vector containing all the values of a column
    std::vector< Cell > cells;
    /// A boolean to display a colummn
    bool enabled;
    /// Vector containing all sub columns subdivison
    std::vector< Column > subColumn;
    CellAlignment cellAlignment;

    Column()
    {
      enabled = true;
      cellAlignment = CellAlignment{ Alignment::center, Alignment::right };
    }

    Column( Cell cell )
    {
      columnName = cell;
      enabled = true;
      cellAlignment = CellAlignment{ Alignment::center, Alignment::right };
    }


    Column & setName( string_view name )
    {
      columnName.value = name;
      return *this;
    }

    Column & setCells( std::vector< Cell > const & cellValues )
    {
      cells = cellValues;
      return *this;
    }

    Column & hide()
    {
      enabled = false;
      return *this;
    }

    Column & setMaxStringSize( size_t const size )
    {
      maxStringSize = size;
      return *this;
    }

    size_t getMaxStringSize() const
    {
      return maxStringSize;
    }

    Column & addSubColumns( std::vector< string > const & subColName )
    {
      std::vector< Column > subColumns;
      for( auto const & name : subColName )
      {
        Cell cell{'\x04', name};//TODO
        Column col{cell};
        subColumns.emplace_back( col );
      }
      subColumn = subColumns;
      return *this;
    }

    Column & setHeaderAlignment( Alignment headerAlignment )
    {
      cellAlignment.headerAlignment = headerAlignment;
      return *this;
    }

    Column & setValuesAlignment( Alignment valueAlignment )
    {
      cellAlignment.valueAlignment = valueAlignment;
      return *this;
    }

private:
    size_t nbHeaderRows;
    /// Vector of string containing the largest string for a column and its subColumns
    size_t maxStringSize;
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
  std::vector< Column > const & getColumns() const;

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
   * @brief Remove the return line at the end & begenning of the table
   * @return The tableLayout reference
   */
  TableLayout & disableLineBreak();

  /**
   * @brief Set the minimal margin width between cell content and borders.
   * @param marginValue The margin value
   * @return The tableLayout reference
   */
  TableLayout & setMargin( MarginValue marginValue );

  /**
   * @return check if the line break at the end & beginning is activated
   */
  bool isLineBreakEnabled() const;

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
   * @brief Add a column to the table given an initializer_list of string & Column
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

  std::vector< Column > m_tableColumnsData;

  bool m_wrapLine = true;
  string m_tableTitle;
  integer m_borderMargin;
  integer m_columnMargin;
  integer m_marginValue;
  integer m_titleMargin = 2;

};
}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP */
