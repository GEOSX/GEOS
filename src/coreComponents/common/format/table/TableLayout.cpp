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
 * @file TableData.hpp
 */
#include "TableLayout.hpp"
#include <numeric>

namespace geos
{

void TableLayout::addToColumns( std::vector< string > const & columnNames )
{
  for( auto const & m_header : columnNames )
  {
    addToColumns( m_header );
  }
}

void TableLayout::addToColumns( string_view m_header )
{
  TableLayout::Column column = TableLayout::Column().setName( m_header );
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
  m_borderMargin = marginValue;
  m_columnMargin = integer( marginValue ) * 2 + 1;

  return *this;
}

bool TableLayout::isLineBreakEnabled() const
{ return m_wrapLine; }

size_t TableLayout::getMaxDepth() const
{
  size_t depthMax = 1;
  size_t currDepth = 1;
  for( auto const & column : m_tableColumnsData )
  {
    currDepth = 1;
    TableLayout::Column const * currColumn = &column;
    while( !currColumn->m_subColumn.empty())
    {
      currColumn = &currColumn->m_subColumn[0];
      currDepth++;
    }
    depthMax = std::max( currDepth, depthMax );
  }
  return depthMax;
}

std::vector< TableLayout::Column > & TableLayout::getColumns()
{ return m_tableColumnsData; }

std::vector< TableLayout::Column > const & TableLayout::getColumns() const
{ return m_tableColumnsData; }

string_view TableLayout::getTitle() const
{ return m_tableTitle; }

integer const & TableLayout::getBorderMargin() const
{ return m_borderMargin; }

integer const & TableLayout::getColumnMargin() const
{ return m_columnMargin; }

integer const & TableLayout::getMarginValue() const
{ return m_marginValue; }

integer const & TableLayout::getMarginTitle() const
{ return m_titleMargin; }

std::vector< size_t > & TableLayout::getSublineInHeaderCounts()
{ return m_sublineHeaderCounts; }

std::vector< size_t > & TableLayout::getNbSubDataLines()
{ return m_sublineDataCounts; }

void divideCell( std::vector< string > & lines, string const & value )
{
  std::istringstream strStream( value );
  std::string line;
  lines.clear();
  while( getline( strStream, line, '\n' ))
  {
    lines.push_back( line );
  }

  if( line.empty())
  {
    lines.push_back( "" );
  }
}

TableLayout::CellLayout::CellLayout():
  m_lines( {""} ),
  m_cellType( CellType::Header ),
  m_alignment( TableLayout::Alignment::center ),
  m_cellWidth( 0 )
{}

TableLayout::CellLayout::CellLayout( CellType type, string const & cellValue, TableLayout::Alignment alignment ):
  m_cellType( type ),
  m_alignment( alignment )
{
  divideCell( m_lines, cellValue );
  if( !m_lines.empty())
  {
    m_cellWidth = std::max_element( m_lines.begin(), m_lines.end(), []( const auto & a, const auto & b )
    {
      return a.length() < b.length();
    } )->length();
  }
  else
  {
    m_cellWidth = 0;
  }
}

TableLayout::Column::Column():
  m_parent( nullptr ), m_next( nullptr )
{
  m_header.m_lines = {};
  m_header.m_cellType  = CellType::Header;
  m_header.m_alignment = Alignment::center;
}

TableLayout::Column::Column( TableLayout::CellLayout cell ):
  m_parent( nullptr ), m_next( nullptr )
{ m_header = cell; }


TableLayout::Column & TableLayout::Column::setName( string_view name )
{
  m_header.m_lines.push_back( std::string( name ) );
  divideCell( m_header.m_lines, m_header.m_lines[0] );
  m_header.m_cellType = CellType::Header;
  return *this;
}

TableLayout::Column & TableLayout::Column::setVisibility( CellType celltype )
{
  m_header.m_cellType = celltype;
  return *this;
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
  m_subColumn = subColumns;
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( std::initializer_list< TableLayout::Column > subCol )
{
  m_subColumn = subCol;
  return *this;
}

TableLayout::Column & TableLayout::Column::addSubColumns( string const & subColName )
{
  TableLayout::CellLayout cell{CellType::Header, subColName, TableLayout::Alignment::center};
  TableLayout::Column col{cell};
  this->m_subColumn.push_back( col );
  return *this;
}

TableLayout::Column & TableLayout::Column::setHeaderAlignment( Alignment headerAlignment )
{
  m_alignment.headerAlignment = headerAlignment;
  m_header.m_alignment = headerAlignment;
  return *this;
}

TableLayout::Column & TableLayout::Column::setValuesAlignment( Alignment valueAlignment )
{
  m_alignment.valueAlignment = valueAlignment;
  return *this;
}

size_t TableLayout::Column::getNumberCellMerge()
{ return m_headerMergeCount; }

void TableLayout::Column::incrementMergeHeaderCount()
{ m_headerMergeCount++;}

void TableLayout::Column::decrementMergeHeaderCount()
{ m_headerMergeCount--; }

bool TableLayout::Column::hasChild() const
{ return !this->m_subColumn.empty(); }
bool TableLayout::Column::hasParent() const
{ return this->m_parent != nullptr; }

}
