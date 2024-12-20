/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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
  struct ColumnAlignement
  {
    /// Alignment for column name. By default aligned to center
    Alignment headerAlignment = Alignment::center;
    /// Alignment for column values. By default aligned to right side
    Alignment valueAlignment = Alignment::right;
  };

/**
 * @struct CellLayout
 * @brief Structure grouping the cell information to display it in a table (content, type, alignment, ...).
 */
  struct CellLayout
  {
    /// vector containing each cell content, separated by lines.
    std::vector< string > m_lines;
    /// The type of the cell (Header,Value, Merge, ...).
    CellType m_cellType;
    /// The alignment of the cell (left, center, right).
    Alignment m_alignment;
    /// Maximum length of the data in the cell.
    size_t m_cellWidth;

    /**
     * @brief Constructor to initialize a Cell with a default settings.
     */
    CellLayout();

    /**
     * @brief Constructor to initialize a cell given celltype, value and alignment.
     * @param cellType The type of the cell.
     * @param value The value to be assigned to the cell.
     * @param alignment The alignment of the cell (left, right, or center).
     */
    CellLayout( CellType cellType, string const & value, TableLayout::Alignment alignment );
  };

  /**
   * @class Column
   * @brief Class representing a column in a table layout.
   */
  class Column
  {
public:
    /// The header cell layout.
    CellLayout m_header;
    /// A vector containing all sub-columns in the column.
    std::vector< Column > m_subColumn;
    /// struct containing m_alignment for the column (header and values)
    ColumnAlignement m_alignment;

    /**
     * @brief Default constructor.
     * Initializes a column with default values.
     */
    Column();

    /**
     * @brief Constructor to initialize a column with a specific `CellLayout`.
     * @param cellLayout The `CellLayout` object to initialize the column.
     *
     */
    Column( TableLayout::CellLayout cellLayout );

    /**
     * @brief Get the parent column.
     * @return Pointer to the parent column, or `nullptr` if no parent is set.
     */
    Column * getParent()
    { return m_parent; }

    /**
     * @brief Set the parent column.
     * @param parent Pointer to the parent column to set.
     */
    void setParent( Column * parent )
    { m_parent = parent; }

    /**
     * @brief GGet the next column in the layout.
     * @return  Pointer to the next column or `nullptr` if no next column exists.
     */
    Column * getNextCell()
    { return m_next; }

    /**
     * @brief Set the next column in the layout.
     * @param nextCell  The next column in the table layout.
     */
    void setNextCell( Column * nextCell )
    {  m_next = nextCell; }

    /**
     * @brief Sets the name of the column.
     * @param name The name to set for the column.
     * @return The current column object.
     */
    Column & setName( string_view name );

    /**
     * @brief Set the column visibility.
     * @param celltype Cell type to apply to hide the colmun
     * @return The current column .
     */
    Column & setVisibility( CellType celltype );

    /**
     * @brief Adds multiple sub-columns to the column.
     * @param subCol A list of sub-column names to add.
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
     * @param valueAlignment The alignment to set for the column values (left, right, or center).
     * @return The current column object
     */
    TableLayout::Column & setValuesAlignment( Alignment valueAlignment );

    /**
     * @return number of times we will divide the current cell
     */
    size_t getNumberCellMerge()
    { return m_headerMergeCount; }

    /**
     * @brief Increment number of times we will divide the current cell
     * @param value number of division to add
     */
    void incrementMergeHeaderCount( size_t value )
    { m_headerMergeCount+= value;}

    /**
     * @brief Decremente number of times we will divide the current cell
     */
    void decrementMergeHeaderCount()
    { m_headerMergeCount--; }

    /**
     * @brief Checks if the column has any child columns.
     * @return bool True if the column has child columns, otherwise false.
     */
    bool hasChild() const
    { return !this->m_subColumn.empty(); }

    /**
     * @brief Checks if the column has a parent column.
     * @return bool True if the column has a parent, otherwise false.
     */
    bool hasParent() const
    { return this->m_parent != nullptr; }

private:
    /// Pointer to the parent cell (if any).
    Column * m_parent;
    /// Pointer to the next cell (if any).
    Column * m_next;
    /// The width of the cell (e.g., for cell containing subColumns).
    size_t m_headerMergeCount  = 0;
  };

  /**
   * @brief Iterator to loop over all columns, starting by the deepest sub columns,
   * then to their parents, then to their siblings.
   */
  class DeepFirstIterator
  {
public:
    ///alias for column
    using ColumnType = Column;

    /**
     * @brief Construct a new Leaf Iterator object
     * @param columnPtr The first deepest column of vector
     * @param idxLayer the layer associated with the column
     */
    DeepFirstIterator( ColumnType * columnPtr, size_t idxLayer ):
      m_currentColumn( columnPtr ), m_currentLayer( idxLayer )
    {}

    /**
     * @brief Copy assignment operator
     * @param[in] columnPtr Coulmn  to copy
     * @return Leaf iterator
     */
    DeepFirstIterator & operator=( Column * columnPtr )
    {
      this->m_currentColumn= columnPtr;
      return *this;
    }

    /**
     * @brief Prefix ++ overload
     * @return Leaf iterator
     */
    DeepFirstIterator & operator++();

    /**
     * @brief Postfix ++ overload
     * @return Leaf iterator
     */
    DeepFirstIterator operator++( int );

    /**
     * @brief Dereference operator.
     * @return Reference to the current Column object pointed to by the iterator.
     */
    ColumnType & operator*()
    { return *m_currentColumn; }

    /**
     * @brief Arrow operator.
     * @return Pointer to the current Column object.
     */
    ColumnType * operator->()
    { return m_currentColumn; }

    /**
     * @brief Equality comparison operator.
     * @param a The first iterator.
     * @param b The second iterator.
     * @return True if both iterators point to the same column; false otherwise.
     */
    friend bool operator== ( DeepFirstIterator const & a, DeepFirstIterator const & b )
    { return a.m_currentColumn == b.m_currentColumn; };
    /**
     * @brief Inequality comparison operator.
     * @param a The first iterator.
     * @param b The second iterator.
     * @return True if the iterators point to different columns; false otherwise.
     */
    friend bool operator!= ( DeepFirstIterator const & a, DeepFirstIterator const & b )
    { return a.m_currentColumn != b.m_currentColumn; };

    /**
     * @brief Gets the current layer (depth) of the iterator.
     * @return The current layer (depth) of the iterator.
     */
    size_t getCurrentLayer() const
    { return m_currentLayer; }

    /**
     * @brief Check if the current cell belong the last column
     * @return true
     * @return false
     */
    bool isLastColumn();

private:
    /// Pointer to the current column
    ColumnType * m_currentColumn;
    /// The current depth of the iterator
    size_t m_currentLayer;
  };

  /**
   * @return Return an itarator pointing on the first leaf of the first columns vector
   * Example on 2 column with Column A : 2 layer and Column B : 3 layers
   * A.A -> A-B -> A-C -> A -> B-A-A -> B-A-B -> B-A -> B-B-A -> B-B-B -> B-B -> B
   */
  DeepFirstIterator beginDeepFirst();

  /**
   * @return Return a end itarator
   * This iterator is initialized with a null pointer
   * representing a position after the last valid element
   */
  DeepFirstIterator endDeepFirst()
  {
    return DeepFirstIterator( nullptr, 0 );
  }

  /// Alias for an initializer list of variants that can contain either a string or a layout column.
  using TableLayoutArgs = std::initializer_list< std::variant< string_view, TableLayout::Column > >;

  TableLayout() = default;

  /**
   * @brief Construct a new Table Layout object
   * @param title The table title
   * @param columns A vector containing all column initialized
   */
  TableLayout( string_view title,
               std::vector< TableLayout::Column > const & columns )
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
   * @brief Get the max depth of a column
   * @return The max column depth
   */
  size_t getMaxDepth() const;

  /**
   * @return The columns vector
   */
  std::vector< Column > & getColumns()
  { return m_tableColumnsData; }

  /**
   * @return The columns vector
   */
  std::vector< Column > const & getColumns() const
  { return m_tableColumnsData; }

  /**
   * @return The table name
   */
  string_view getTitle() const
  { return m_tableTitle; }

  /**
   * @param title The table title
   * @return The tableLayout reference
   */
  TableLayout & setTitle( string_view title );

  /**
   * @brief Remove the return line at the end & begenning of the table
   * @param value Value to desactivate or not wrapLine at the end
   * @return The tableLayout reference
   */
  TableLayout & disableLineBreak( bool value );

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
   * @return The border margin,
   * number of spaces at both left and right table sides plus vertical character
   */
  integer const & getBorderMargin() const
  { return m_borderMargin; }

  /**
   * @return The column margin,
   * numbers of spaces separating both left and right side from a vertical line
   */
  integer const & getColumnMargin() const
  { return m_columnMargin; }

  /**
   * @return The table margin value
   */
  integer const & getMarginValue() const
  { return m_marginValue; }

  /**
   * @return The margin title
   */
  integer const & getMarginTitle() const
  { return m_titleMargin; }

/**
 * @brief Get the Nb Rows object
 * @return std::vector< integer >&
 */
  std::vector< size_t > & getSublineInHeaderCounts()
  { return m_sublineHeaderCounts; }

/**
 * @brief Get the Nb Rows object
 * @return std::vector< integer >&
 */
  std::vector< size_t > & getNbSubDataLines()
  { return m_sublineDataCounts; }

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
   * @param m_header The column name
   */
  void addToColumns( string_view m_header );

/**
 *
 * @brief Create and add a column to the columns vector given a Column
 * @param column Vector containing addition information on the column
 */
  void addToColumns( TableLayout::Column const & column );

  /// Contains the columns layout
  std::vector< Column > m_tableColumnsData;
  /// Contains the subdivision (line) counts for each line in header.
  std::vector< size_t > m_sublineHeaderCounts;
  /// Contains the subdivision (line) counts for each line in data.
  std::vector< size_t > m_sublineDataCounts;
  bool m_wrapLine = true;

  string m_tableTitle;

  integer m_borderMargin;
  integer m_columnMargin;
  integer m_marginValue;
  integer m_titleMargin = 2;

};
}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLELAYOUT_HPP */
