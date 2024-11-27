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
 * @file TableData.hpp
 */
#include "TableLayout.hpp"
#include <numeric>

namespace geos
{

void TableLayout::addToColumns( std::vector< string > const & columnNames )
{
  for( auto const & columnName : columnNames )
  {
    addToColumns( columnName );
  }
}

void TableLayout::addToColumns( string_view columnName )
{
  TableLayout::Column column = TableLayout::Column().setName( columnName );
  m_tableColumnsData.push_back( column );
}

void TableLayout::addToColumns( TableLayout::Column const & column )
{
  m_tableColumnsData.push_back( column );
}

TableLayout & TableLayout::setTitle( string_view title )
{
  m_tableTitle = title;
  return *this;
}

TableLayout & TableLayout::disableLineBreak()
{
  m_wrapLine = false;
  return *this;
}

TableLayout & TableLayout::setMargin( MarginValue marginValue )
{
  m_marginValue = marginValue;
  m_borderMargin = marginValue + 1; // margin + border character
  m_columnMargin = integer( marginValue ) * 2 + 1;

  return *this;
}

bool TableLayout::isLineBreakEnabled() const
{
  return m_wrapLine;
}

void TableLayout::setContainingSubColumn()
{
  m_containSubColumn = true;
}

bool TableLayout::isContainingSubColumn() const
{
  return m_containSubColumn;
}

void TableLayout::removeSubColumn()
{
  for( auto & column : m_tableColumnsData )
  {
    if( !column.subColumn.empty() )
    {
      column.subColumn.clear();
    }
  }
}

size_t TableLayout::getMaxHeaderRow() const
{
  size_t depthMax = 0;
  size_t currDepth = 0;
  for( auto const & column : m_tableColumnsData )
  {
    currDepth = 0;
    TableLayout::Column const * currColumn = &column;
    while( !currColumn->subColumn.empty())
    {
      currColumn = &currColumn->subColumn[0];
      currDepth++;
    }
    depthMax = std::max( currDepth, depthMax );
  }
  return depthMax;
}

std::vector< TableLayout::Column > const & TableLayout::getColumns() const
{
  return m_tableColumnsData;
}

string_view TableLayout::getTitle() const
{
  return m_tableTitle;
}

integer const & TableLayout::getBorderMargin() const
{
  return m_borderMargin;
}

integer const & TableLayout::getColumnMargin() const
{
  return m_columnMargin;
}

integer const & TableLayout::getMarginValue() const
{
  return m_marginValue;
}

integer const & TableLayout::getMarginTitle() const
{
  return m_titleMargin;
}

void divideCell( std::vector< string > & lines, string const & value )
{
  std::istringstream strStream( value );
  std::string line;
  while( getline( strStream, line, '\n' ))
  {
    lines.push_back( line );
  }
}

TableLayout::CellLayout::CellLayout():
  lines( {} ), cellType( CellType::Header ), alignment( TableLayout::Alignment::center )
{}

TableLayout::CellLayout::CellLayout( CellType type, string const & cellValue, TableLayout::Alignment cellAlignment ):
  cellType( type ),
  alignment( cellAlignment )
{
  divideCell( lines, cellValue );
  maxDataLength = std::max_element( lines.begin(), lines.end(), []( const auto & a, const auto & b )
  {
    return a.length() < b.length();
  } )->length();
}

void TableLayout::CellLayout::setMaxCellSize( size_t size )
{
  maxDataLength = size;
}


//
// COLUMN
//

TableLayout::Column::Column():
  m_parent( nullptr ), m_next( nullptr )
{
  enabled = true;
  columnName.lines = {};
  columnName.cellType  = CellType::Header;
  columnName.alignment = Alignment::center;
}

TableLayout::Column::Column( TableLayout::CellLayout cell ):
  m_parent( nullptr ), m_next( nullptr )
{
  columnName = cell;
  enabled = true;
}


TableLayout::Column & TableLayout::Column::setName( string_view name )
{
  columnName.lines.push_back( std::string( name ) );
  columnName.cellType = CellType::Header;
  return *this;
}

TableLayout::Column & TableLayout::Column::hide()
{
  enabled = false;
  return *this;
}

void TableLayout::Column::setMaxStringSize( size_t const size )
{
  maxStringSize = size;
}

size_t TableLayout::Column::getMaxStringSize() const
{
  return maxStringSize;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::initializer_list< string > subColName )
{
  std::vector< TableLayout::Column > subColumns;
  for( auto const & name : subColName )
  {
    TableLayout::CellLayout cell{CellType::Header, name, TableLayout::Alignment::center};
    TableLayout::Column col{cell};
    subColumns.emplace_back( col );
  }
  subColumn = subColumns;
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( string const & subColName )
{
  TableLayout::CellLayout cell{CellType::Header, subColName, TableLayout::Alignment::center};
  TableLayout::Column col{cell};
  this->subColumn.push_back( col );
  return *this;
}

TableLayout::Column & TableLayout::Column::setHeaderAlignment( Alignment headerAlignment )
{
  cellAlignment.headerAlignment = headerAlignment;
  columnName.alignment = headerAlignment;
  return *this;
}

TableLayout::Column & TableLayout::Column::setValuesAlignment( Alignment valueAlignment )
{
  cellAlignment.valueAlignment = valueAlignment;
  return *this;
}

void TableLayout::Column::compareAndSetMaxStringSize( TableLayout::Column * currentColumn,
                                                      std::vector< size_t > & subColumnsLength )//todo renaùme
{
  auto sumSubColumnsLen = std::reduce( subColumnsLength.begin(), subColumnsLength.end());

  size_t columnNameLinesMaxLength = std::max_element(
    currentColumn->columnName.lines.begin(),
    currentColumn->columnName.lines.end(),
    []( const auto & a, const auto & b ) {
    return a.length() < b.length();
  } )->length();

  size_t maxSize = std::max( sumSubColumnsLen, columnNameLinesMaxLength );

  currentColumn->columnName.maxDataLength = maxSize;
  currentColumn->setMaxStringSize( maxSize );
}

bool TableLayout::Column::hasChild() const
{
  return !this->subColumn.empty();
}
bool TableLayout::Column::hasParent() const
{
  return this->m_parent != nullptr;
}

}
