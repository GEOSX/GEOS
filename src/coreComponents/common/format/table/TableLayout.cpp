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

namespace geos
{

void TableLayout::addToColumns( const std::vector< string > & columnNames )
{
  for( const auto & columnName : columnNames )
  {
    addToColumns( columnName );
  }
}

void TableLayout::addToColumns( string_view columnName )
{
  Column column = TableLayout::Column().setName( columnName );
  if( !m_tableColumnsData.empty( ))
  {
    m_tableColumnsData.end()->next = &column;
  }
  m_tableColumnsData.push_back( column );
}

void TableLayout::addToColumns( Column const & column )
{
  if( !m_tableColumnsData.empty( ))
  {
    m_tableColumnsData.end()->next = &column;
  }
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

integer & TableLayout::getTrackerHeaderRows() const
{
  return m_headerRows;
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

//
// CellLayout
//


TableLayout::CellLayout::CellLayout( CellType type, string_view cellValue, Alignment cellAlignment ):
  cellType( type ),
  alignment( cellAlignment )
{
  divideCell( lines, cellValue );
}

TableLayout::CellLayout::CellLayout( CellType type, string_view cellValue ):
  cellType( type )
{
  divideCell( lines, cellValue );
}

void TableLayout::CellLayout: setAlignment( Alignment const align )
{
  alignment = align;
}

void TableLayout::setCellSize( size_t const size )
{
  cellSize = size;
}

//
// COLUMN
//

TableLayout::Column()
  : m_parent( nullptr ), m_m_next( nullptr )
{
  enabled = true;
  columnName.value = "";
  columnName.type  = '\x03';
  columnName.alignment = Alignment::center;
}

TableLayout::Column( Cell cell )
{
  columnName = cell;
  columnName.alignment = Alignment::center;
  enabled = true;
}


TableLayout::Column & TableLayout::Column::setName( string_view name )
{
  columnName.value = name;
  columnName.type = '\x03';
  return *this;
}

TableLayout::Column & TableLayout::Column::setCells( std::vector< Cell > const & cellValues )
{
  cells = cellValues;
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
  return *this;
}

size_t TableLayout::Column::getMaxStringSize() const
{
  return maxStringSize;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::vector< string > const & subColName )
{
  std::vector< TableLayout::Column > subColumns;
  for( auto const & name : subColName )
  {
    Cell cell{'\x03', name};    //TODO
    TableLayout::Column col{cell};
    subColumns.parent = this;
    if( !subColumns.empty())
    {
      subColumns.end()->next = &col;
    }
    subColumns.emplace_back( col );
  }
  subColumn = subColumns;
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

void TableLayout::Column::updateMaxStringSizeForColumn( Column * currentColumn,
                                                        std::vector< size_t > & subColumnsLength )//todo renaùme
{
  auto sumSubColumnsLen = std::reduce(
    subColumnsLength.begin(),
    subColumnsLength.end()).length();
  currentColumn->setMaxStringSize( std::max(
                                     sumSubColumnsLen,
                                     getMaxStringLen( currentColumn->columnName.lines )));
}

}
