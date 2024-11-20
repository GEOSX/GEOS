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
#include "TableTypes.hpp"
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
    Alignment headerAlignment = Alignment::center;
    /// Alignment for column values. By default aligned to right side
    Alignment valueAlignment = Alignment::right;
  };

  struct CellLayout
  {
    /// Vector containing sub values name
    std::vector< string > lines;
    CellType type;
    /// Cell alignment (left, right, center)
    Alignment alignment;

    /**
     * @brief Constructor to initialize a Cell with a specific type and value.
     * @param t The type of the cell.
     * @param val The value to be assigned to the cell.
     */
    CellLayout( CellType t, string_view val, Alignement alignment );

    void setAlignment( Alignment const align )
    {
      alignment = align;
    }

    void setCellSize( size_t const size )
    {
      cellSize = size;
    }
  };

  /**
   * @brief Struct for a Column.
   * Each column contains its own parameters (such as name, alignment, etc.).
   */
  class Column
  {
public:
    /// Name for a column
    CellLayout columnName;
    /// A boolean to display a colummn
    bool enabled;
    /// Vector containing all sub columns subdivison
    std::vector< Column > subColumn;
    /// struct containing alignment for the column (header and values)
    CellAlignment cellAlignment;

    Column()
    {
      enabled = true;
      columnName.value = "";
      columnName.type  = '\x03';
      columnName.alignment = Alignment::center;
    }

    Column( Cell cell )
    {
      columnName = cell;
      columnName.alignment = Alignment::center;
      enabled = true;
    }


    Column & setName( string_view name )
    {
      columnName.value = name;
      columnName.type = '\x03';
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

    void setMaxStringSize( size_t const size )
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
        Cell cell{'\x03', name};//TODO
        Column col{cell};
        subColumns.emplace_back( col );
      }
      subColumn = subColumns;
      return *this;
    }

    Column & setHeaderAlignment( Alignment headerAlignment )
    {
      cellAlignment.headerAlignment = headerAlignment;
      columnName.alignment = headerAlignment;
      return *this;
    }

    Column & setValuesAlignment( Alignment valueAlignment )
    {
      cellAlignment.valueAlignment = valueAlignment;
      return *this;
    }

    // std::vector< CellLayout > & getCells()
    // {
    //   return cells;
    // }

private:
    // /// A vector containing all the values of a column
    // std::vector< CellLayout > cells;
    // /// TODO DOCS
    // size_t nbHeaderRows;
    /// Vector of string containing the largest string for a column and its subColumns
    size_t maxStringSize; // TODO : Assigner cette stat
  };

  /**
   * @brief Allow to iterate among all deepest columns / sub columns.
   * An exemple of an iteration: A -> B.A -> B.B -> B.C -> C.A.A -> C.A.B -> C.B.A -> C.B.B -> D
   */
  class SubColumnIterator
  {
    SubColumnIterator() noexcept:
      m_currentColumn( m_spRoot )
      { }

    SubColumnIterator( Column const * columnPtr ) noexcept:
      m_currentColumn( columnPtr )
      { }

    SubColumnIterator & operator=( Column * columnPtr )
    {
      this->m_currentColumn= columnPtr;
      return *this;
    }

    // Prefix ++ overload
    SubColumnIterator & operator++()
    {
      // TODO!
      // if( m_currentColumn )
      //   m_currentColumn= m_currentColumn>pNext;
      return *this;
    }

    // Postfix ++ overload
    SubColumnIterator operator++( Column )
    {
      SubColumnIterator iterator = *this;
      ++*this;
      return iterator;
    }

    bool operator!=( SubColumnIterator const & iterator )
    {
      return m_currentColumn!= iterator.m_currentColumn
    }

    Column operator*()
    {
      return *m_currentColumn;
    }

private:
    Column const * m_currentColumn
  };

  struct Row
  {
    // maximum number of lines among the cells of a given row
    size_t maxLineCount; // TODO : Assigner cette stat
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
   * @brief
   */
  void setContainingSubColumn();

  /**
   * @return
   */
  bool isContainingSubColumn() const;

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
  bool m_containSubColumn = false;

  string m_tableTitle;

  integer m_borderMargin;
  integer m_columnMargin;
  integer m_marginValue;
  integer m_titleMargin = 2;

};
}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP */
