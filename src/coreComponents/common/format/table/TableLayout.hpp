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

    /**
     * @brief Set the Alignment object
     * @param align
     */
    void setAlignment( Alignment const align );

    /**
     * @brief Set the Cell Size object
     * @param size
     */
    void setCellSize( size_t const size );
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

    Column * m_parent;
    Column * m_next;

    Column();

    Column( Cell cell );

    /**
     * @brief Set the Name object
     * @param name
     * @return Column&
     */
    Column & setName( string_view name );

    /**
     * @brief Set the Cells object
     * @param cellValues
     * @return Column&
     */
    Column & setCells( std::vector< Cell > const & cellValues );

    /**
     * @brief
     * @return Column&
     */
    Column & hide();

    /**
     * @brief Set the Max String Size object
     * @param size
     */
    void setMaxStringSize( size_t const size );

    /**
     * @brief Get the Max String Size object
     * @return size_t
     */
    size_t getMaxStringSize() const;

    /**
     * @brief
     * @param subColName
     * @return Column&
     */
    Column & addSubColumns( std::vector< string > const & subColName );

    /**
     * @brief Set the Header Alignment object
     * @param headerAlignment
     * @return Column&
     */
    Column & setHeaderAlignment( Alignment headerAlignment );

    /**
     * @brief Set the Values Alignment object
     * @param valueAlignment
     * @return Column&
     */
    Column & setValuesAlignment( Alignment valueAlignment );

    /**
     * @brief
     * @param column
     */
    void updateMaxStringSize( Column * column );

private:
    size_t maxStringSize;
  };

  /**
   * @brief Allow to iterate among all deepest columns / sub columns.
   * An exemple of an iteration: A -> B.A -> B.B -> B.C -> C.A.A -> C.A.B -> C.B.A -> C.B.B -> D
   */
  class SubColumnIterator
  {

    SubColumnIterator( Column const * columnPtr ) noexcept:
      m_currentColumn( columnPtr )
    {}


    SubColumnIterator & operator=( Column * columnPtr )
    {
      this->m_currentColumn= columnPtr;
      return *this;
    }

    // Prefix ++ overload
    SubColumnIterator & operator++()
    {
      m_currentColumn++;
      if( m_currentColumn->m_next == nullptr ) GEOS_ERROR( "Column overflow" );
      if( m_currentColumn->m_parent == nullptr )
      {
        while( !m_currentColumn->subColumn.empty() )
        {
          m_currentColumn = m_currentColumn.subColumn.begin();
        }
      }
      else
      {
        if( m_currentColumn == m_currentColumn->m_parent->subColumn.end())
        {
          while( m_currentColumn->m_parent.empty() )
          {
            m_currentColumn = m_currentColumn->m_parent;
          }
        }
      }
      m_currentColumn = m_currentColumn->m_next;
      return *this;
    }

    // Postfix ++ overload
    SubColumnIterator operator++( Column )
    {
      SubColumnIterator iterator = *this;
      ++(*this);
      return iterator;
    }

    Column operator*()
    {
      return *m_currentColumn;
    }

    Column operator->()
    {
      return m_currentColumn;
    }

    friend bool operator== ( const SubColumnIterator & a, const SubColumnIterator & b )
    {
      return a.m_currentColumn == b.m_currentColumn;
    };
    friend bool operator!= ( const SubColumnIterator & a, const SubColumnIterator & b )
    {
      return a.m_currentColumn != b.m_currentColumn;
    };

private:
    Column const * m_currentColumn;
  };

  SubColumnIterator begin() { return SubColumnIterator( m_tableColumnsData.begin() );}
  SubColumnIterator end() { return SubColumnIterator( m_tableColumnsData.end() );}
  struct Row
  {
    // maximum number of lines among the cells of a given row
    size_t maxLineCount;   // TODO : Assigner cette stat
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

  size_t getMaxHeaderRow()
  {
    size_t depthMax=1;
    size_t maxLineCount=1;
    for( auto it = m_tableColumnsData.begin(), end = m_tableColumnsData.end(); it!=end; ++it )
    {
      if( !it->subColumn.empty()) depthMax++;
    }
    return depthMax;
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

  std::vector< Row > & getTrackerHeaderRows();

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
  // m_valueRows[0] = header then 1 line = 1 row.
  std::vector< Row > m_valueRows;
  std::vector< Row > m_headerRows;

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
