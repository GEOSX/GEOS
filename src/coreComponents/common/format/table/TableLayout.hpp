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
#include "common/logger/Logger.hpp"


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
   * @brief Structure to set up values m_alignment for each colum.
   */
  struct CellAlignment
  {
    /// Alignment for column name. By default aligned to center
    Alignment headerAlignment = Alignment::center;
    /// Alignment for column values. By default aligned to right side
    Alignment valueAlignment = Alignment::right;
  };

/**
 * @struct CellLayout
 * @brief Structure representing a cell in a table.
 * This structure contains information about the cell such as its type, alignment, maximum data length, and width.
 */
  struct CellLayout
  {
    /// vector containing the cell name separated by a '\n'.
    std::vector< string > m_lines;
    /// The type of the cell (Header,Value, Merge, ...).
    CellType m_cellType;
    /// The alignment of the cell (left, center, right).
    Alignment m_alignment;
    /// Maximum length of the data in the cell.
    size_t m_maxDataLength;

    /**
     * @brief Constructor to initialize a Cell with a specific type and value.
     */
    CellLayout();

    /**
     * @brief Constructor to initialize a cell given celltype, value and alignment.
     * @param cellType The type of the cell.
     * @param value The value to be assigned to the cell.
     * @param alignment The alignment of the cell (left, right, or center).
     */
    CellLayout( CellType cellType, string const & value, TableLayout::Alignment alignment );

    /**
     * @brief Sets the maximum size for the cell.
     * @param size The maximum size to set for the cell.
     */
    void setMaxCellSize( size_t size );//todo supp
  };

  /**
   * @class Column
   * @brief Class representing a column in a table layout.
   */
  class Column
  {
public:
    /// The name of the column.
    CellLayout m_columName;
    /// A vector containing all sub-columns in the column.
    std::vector< Column > m_subColumn;
    /// struct containing m_alignment for the column (header and values)
    CellAlignment m_cellAlignment;
    /// Pointer to the parent column (if any).
    Column * m_parent;
    /// Pointer to the next column (if any).
    Column * m_next;

    /// The width of the cell (e.g., for cell containing subColumns).
    size_t m_cellWidth = 0;

    /**
     * @brief Default constructor.
     * Initializes a column with default values.
     */
    Column();

    /**
     * @brief Constructor to initialize a column with a specific `CellLayout`.
     * @param cell The `CellLayout` object to initialize the column.
     *
     */
    Column( TableLayout::CellLayout cellLayout );

    /**
     * @brief Sets the name of the column.
     * @param name The name to set for the column.
     * @return The current column object.
     */
    Column & setName( string_view name );

    /**
     * @brief Hides the column.
     * @return The current column objec.
     */
    Column & hide();

    /**
     * @brief Sets the maximum string size for the column.
     * @param size The size to set as the maximum string length.
     */
    void setMaxStringSize( size_t const size );

    /**
     * @brief Gets the maximum string size of the column.
     * @return size_t The maximum string size of the column.
     */
    size_t getMaxStringSize() const;

    /**
     * @brief Adds multiple sub-columns to the column.
     * @param subColName A list of sub-column names to add.
     * @return The current column object
     */
    TableLayout::Column & addSubColumns( std::initializer_list< TableLayout::Column > subCol );

    /**
     * @brief Adds multiple sub-columns to the column.
     * @param subColName A list of sub-column names to add.
     * @return The current column object
     */
    TableLayout::Column & addSubColumns( std::initializer_list< string > subColName );

    /**
     * @brief Adds a single sub-column to the column.
     * @param subColName The name of the sub-column to add.
     * @return The current column object.
     */
    TableLayout::Column & addSubColumns( string const & subColName );

    /**
     * @brief Sets the header alignment for the column.
     * @param headerAlignment The alignment to set for the column header (left, right, or center).
     * @return The current column object
     */
    TableLayout::Column & setHeaderAlignment( Alignment headerAlignment );

    /**
     * @brief Sets the values alignment for the column.
     * @param headerAlignment The alignment to set for the column values (left, right, or center).
     * @return The current column object
     */
    TableLayout::Column & setValuesAlignment( Alignment valueAlignment );

    /**
     * @brief Checks if the column has any child columns.
     * @return bool True if the column has child columns, otherwise false.
     */
    bool hasChild() const;

    /**
     * @brief Checks if the column has a parent column.
     * @return bool True if the column has a parent, otherwise false.
     */
    bool hasParent() const;

private:
    /// The maximum string size of the column.
    size_t m_maxStringSize;
  };

  /**_est columns / sub columns.
   * An exemple of an iteration: A -> B.A -> B.B -> B.C -> C.A.A -> C.A.B -> C.B.A -> C.B.B -> D
   */
  class LeafIterator
  {
public:
    using ColumnType = Column;

    LeafIterator( ColumnType * columnPtr, size_t idxLayer ):
      m_currentColumn( columnPtr ), m_currentLayer( idxLayer )
    {}


    LeafIterator & operator=( Column * columnPtr )
    {
      this->m_currentColumn= columnPtr;
      return *this;
    }

    // Prefix ++ overload
    LeafIterator & operator++()
    {
      if( m_currentColumn->m_next != nullptr )
      {
        // chercher le dernier sous-enfant du suivant //todo
        m_currentColumn = m_currentColumn->m_next;
        while( m_currentColumn->hasChild() )
        {
          m_currentLayer++;
          m_currentColumn = &m_currentColumn->m_subColumn[0];
        }
      }
      else
      {
        bool const hasParent = (m_currentColumn->m_parent != nullptr);
        m_currentLayer -= size_t( hasParent );
        m_currentColumn = hasParent ? m_currentColumn->m_parent : nullptr;
      }
      return *this;
    }

    // Postfix ++ overload //todo
    // LeafIterator & operator++()
    // {
    //   LeafIterator iterator = *this;
    //   ++(*this);
    //   return iterator;
    // }

    ColumnType & operator*()
    {
      return *m_currentColumn;
    }

    ColumnType * operator->()
    {
      return m_currentColumn;
    }

    friend bool operator== ( LeafIterator const & a, LeafIterator const & b )
    {
      return a.m_currentColumn == b.m_currentColumn;
    };
    friend bool operator!= ( LeafIterator const & a, LeafIterator const & b )
    {
      return a.m_currentColumn != b.m_currentColumn;
    };

    size_t getCurrentLayer() const
    {
      return m_currentLayer;
    }

private:
    ColumnType * m_currentColumn;
    size_t m_currentLayer;
  };

  LeafIterator beginLeaf()
  {
    Column * startColumn = &(*m_tableColumnsData.begin());
    size_t idxLayer = 0;
    if( startColumn->hasChild() )
    {
      while( startColumn->hasChild() )
      {
        idxLayer++;
        startColumn = &startColumn->m_subColumn[0];
      }
    }
    return LeafIterator( startColumn, idxLayer );
  }

  LeafIterator endLeaf()
  {
    return LeafIterator( nullptr, 0 );
  }

  struct Row
  {
    // maximum number of lines among the cells of a given row
    size_t maxLineCount;   // TODO : Assigner cette stat
  };

  /// Alias for an initializer list of variants that can contain either a string or a layout column.
  using TableLayoutArgs = std::initializer_list< std::variant< string_view, TableLayout::Column > >;

  TableLayout() = default;

  // delete copy construct//todo
  //default move construct//tdo

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
    for( auto & column :columns )
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

  size_t getMaxHeaderRow() const;

  /**
   * @return The columns vector
   */
  std::vector< Column > & getColumns();

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

/**
 * @brief Get the Nb Rows object
 * @return std::vector< integer >&
 */
  std::vector< size_t > & getSublineInHeaderCounts();

/**
 * @brief Get the Nb Rows object
 * @return std::vector< integer >&
 */
  std::vector< size_t > & getNbSubDataLines();


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
   * @brief
   *
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
   * @param m_columName The column name
   */
  void addToColumns( string_view m_columName );

/**
 *
 * @brief Create and add a column to the columns vector given a Column
 * @param column Vector containing addition information on the column
 */
  void addToColumns( TableLayout::Column const & column );

  std::vector< Column > m_tableColumnsData;
  /// Contains the subdivision (line) counts for each line in header.
  std::vector< size_t > m_sublineHeaderCounts;
  /// Contains the subdivision (line) counts for each line in data.
  std::vector< size_t > m_sublineDataCounts;
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
