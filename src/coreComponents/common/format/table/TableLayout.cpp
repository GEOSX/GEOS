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


CellLayout::CellLayout( CellType type, string_view cellValue, Alignment cellAlignment ):
  cellType( type ),
  alignment( cellAlignment )
{
  divideCell( lines, cellValue );
}

CellLayout::CellLayout( CellType type, string_view cellValue ):
  cellType( type )
{
  divideCell( lines, cellValue );
}

}
