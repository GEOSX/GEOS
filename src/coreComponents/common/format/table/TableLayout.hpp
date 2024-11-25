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
    CellType cellType;
    /// Cell alignment (left, right, center)
    Alignment alignment;
    size_t maxDataLength;

    /**
     * @brief Constructor to initialize a Cell with a specific type and value.
     */
    CellLayout();

    /**
     * @brief Constructor to initialize a Cell with a specific type and value.
     * @param t The type of the cell.
     * @param val The value to be assigned to the cell.
     */
    CellLayout( CellType t, string const & val, TableLayout::Alignment alignment );
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

    Column( TableLayout::CellLayout cellLayout );

    /**
     * @brief Set the Name object
     * @param name
     * @return Column&
     */
    Column & setName( string_view name );

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
    void compareAndSetMaxStringSize( Column * column, std::vector< size_t > & subColumnsLength );


    /**
     * @brief
     * @return Column&
     */
    bool hasChild() const;

    /**
     * @brief
     * @return Column&
     */
    bool hasParent() const;

private:
    size_t maxStringSize;
  };

  /**_est columns / sub columns.
   * An exemple of an iteration: A -> B.A -> B.B -> B.C -> C.A.A -> C.A.B -> C.B.A -> C.B.B -> D
   */

  class LeafIterator
  {
public:
    LeafIterator( Column const * columnPtr ):
      m_currentColumn( columnPtr )
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
        // chercher le dernier sous-enfant du suivant
        m_currentColumn = m_currentColumn->m_next;
        while( m_currentColumn->hasChild() ) m_currentColumn = &m_currentColumn->subColumn[0];
      }
      else
      {
        m_currentColumn = m_currentColumn->m_parent != nullptr ? m_currentColumn->m_parent : nullptr;
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

    Column const & operator*()
    {
      return *m_currentColumn;
    }

    Column const * operator->()
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

private:
    Column const * m_currentColumn;
  };

  LeafIterator beginLeaf()
  {
    Column const * startColumn = &(*m_tableColumnsData.begin());
    if( startColumn->hasChild() )
    {
      while( startColumn->hasChild() ) startColumn = &startColumn->subColumn[0];
    }
    return LeafIterator( startColumn );
  }

  LeafIterator endLeaf() { return LeafIterator( &(*m_tableColumnsData.end()) );}

  //

  class RootIterator
  {
public:
    RootIterator( Column const * columnPtr ) noexcept:
      m_currentColumn( columnPtr ), m_checkpointColumn( 0, nullptr )
    {}


    RootIterator & operator=( Column * columnPtr )
    {
      this->m_currentColumn= columnPtr;
      return *this;
    }

    RootIterator & operator++()
    {
      if( m_currentColumn->hasChild())
      {
        this->setCheckpoint( m_currentColumn );
        m_currentColumn = &m_currentColumn->subColumn[0];
      }
      else
      {
        if( m_currentColumn->m_next != nullptr )
        {
          m_currentColumn = m_currentColumn->m_next;
        }
        else
        {
          while( this->hasCheckPoint() && this->getLastCheckPoint()->m_next == nullptr )
          {
            m_currentColumn = this->getLastCheckPoint();
            m_checkpointColumn.pop_back();
          }
          if( this->hasCheckPoint() )
          {
            m_currentColumn = this->getLastCheckPoint()->m_next;
            m_checkpointColumn.pop_back();
          }
        }
      }
      return *this;
    }

    // Postfix ++ overload //todo
    // RootIterator operator++( Column )
    // {
    //   RootIterator iterator = *this;
    //   ++(*this);
    //   return iterator;
    // }

    Column const & operator*()
    {
      return *m_currentColumn;
    }

    Column const * operator->()
    {
      return m_currentColumn;
    }

    friend bool operator== ( RootIterator const & a, RootIterator const & b )
    {
      return a.m_currentColumn == b.m_currentColumn;
    };
    friend bool operator!= ( RootIterator const & a, RootIterator const & b )
    {
      return a.m_currentColumn != b.m_currentColumn;
    };

    Column const * getLastCheckPoint() const
    {
      if( !m_checkpointColumn.empty())
        return m_checkpointColumn.back();
      return nullptr;
    }

    RootIterator & setCheckpoint( Column const * column )
    {
      this->m_checkpointColumn.push_back( column );
      return *this;
    }

    bool hasCheckPoint()
    {
      return m_checkpointColumn.empty();
    }

private:
    Column const * m_currentColumn;
    std::vector< Column const * > m_checkpointColumn;
  };

  RootIterator beginRoot()
  {return RootIterator( &(*m_tableColumnsData.begin()) );}

  RootIterator endRoot()
  {
    Column const * startColumn = &(*m_tableColumnsData.end());
    if( startColumn->hasChild() )
    {
      while( startColumn->hasChild() )
      {
        startColumn = &startColumn->subColumn[startColumn->subColumn.size() - 1];
      }


    }
    return RootIterator( &(*m_tableColumnsData.end()) );
  }


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
      std::visit( [this]( auto & value ) {
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
  void addToColumns( Column & column );

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
