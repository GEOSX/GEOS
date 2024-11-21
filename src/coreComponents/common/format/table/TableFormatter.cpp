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
 * @file TableFormatter.cpp
 */

#include "TableFormatter.hpp"
#include <numeric>
#include "common/format/StringUtilities.hpp"
#include "common/logger/Logger.hpp"
#include "TableFormatter.hpp"

namespace geos
{

TableFormatter::TableFormatter( TableLayout const & tableLayout ):
  m_tableLayout( tableLayout )
{}

///////////////////////////////////////////////////////////////////////
////// CSV Formatter implementation
///////////////////////////////////////////////////////////////////////

TableCSVFormatter::TableCSVFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{
  m_tableLayout = tableLayout;
}

string TableCSVFormatter::headerToString() const
{
  std::stringstream oss;
  static constexpr string_view separator = ",";

  for( std::size_t idxColumn = 0; idxColumn < m_tableLayout.getColumns().size(); ++idxColumn )
  {
    oss << m_tableLayout.getColumns()[idxColumn].columnName.value;
    if( idxColumn < m_tableLayout.getColumns().size() - 1 )
    {
      oss << separator;
    }
  }
  oss << "\n";
  return oss.str();
}

string TableCSVFormatter::dataToString( TableData const & tableData ) const
{

  RowsCellInput const rowsValues( tableData.getTableDataRows() );
  std::ostringstream oss;
  for( const auto & row : rowsValues )
  {
    std::vector< string > rowConverted;
    for( const auto & item : row )
    {
      rowConverted.push_back( item.value );
    }
    oss << stringutilities::join( rowConverted.cbegin(), rowConverted.cend(), "," ) << "\n";
  }
  return oss.str();
}

template<>
string TableCSVFormatter::toString< TableData >( TableData const & tableData ) const
{
  return headerToString() + dataToString( tableData );
}

///////////////////////////////////////////////////////////////////////
////// Log Formatter implementation
///////////////////////////////////////////////////////////////////////

/**
 * @brief Build cell given an alignment, a value and spaces
 * @param alignment The aligment of cell value
 * @param value The cell value
 * @param spaces The number of spaces in the cell
 * @return A formated cell
 */
string buildCell( TableLayout::Alignment const alignment, string_view value, size_t const spaces )
{
  switch( alignment )
  {
    case TableLayout::right:   return GEOS_FMT( "{:>{}}", value, spaces );
    case TableLayout::left:    return GEOS_FMT( "{:<{}}", value, spaces );
    case TableLayout::center:  return GEOS_FMT( "{:^{}}", value, spaces );
    default:                   return GEOS_FMT( "{:>{}}", value, spaces );
  }
}

/**
 * @brief Detect columns who are not displayed from TableLayout and therefore modify columns / tableDataRows vectors
 * @param columns Vector built in TableLayout containing all columns with their parameters
 * @param tableDataRows Vector built in TableData containing all rows values
 */
void updateVisibleColumns( std::vector< TableLayout::Column > & columns,
                           RowsCellInput & tableDataRows )
{
  integer idxColumn = 0;
  for( auto iterColumn = columns.begin(); iterColumn != columns.end(); )
  {
    if( !iterColumn->enabled )
    {
      iterColumn = columns.erase( iterColumn );
      for( auto & row : tableDataRows )
      {
        row.erase( row.begin() + idxColumn );
      }
    }
    else
    {
      ++iterColumn;
      ++idxColumn;
    }
  }
}

TableTextFormatter::TableTextFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{}

string TableTextFormatter::toString() const
{
  std::ostringstream tableOutput;
  std::vector< TableLayout::Column > columns( m_tableLayout.getColumns());
  TableData tableData;
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( columns, tableData, sectionSeparatingLine, topSeparator );
  tableOutput << '\n';
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  outputHeader( columns, tableOutput, sectionSeparatingLine );

  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::vector< TableLayout::Column > columns( m_tableLayout.getColumns());
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( columns, tableData, sectionSeparatingLine, topSeparator );
  std::ostringstream tableOutput;
  outputTable( tableOutput, columns, sectionSeparatingLine, topSeparator );
  return tableOutput.str();
}

void TableTextFormatter::prepareAndBuildTable( std::vector< TableLayout::Column > & columns,
                                               TableData const & tableData,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{
  RowsCellInput tableDataRows( tableData.getTableDataRows());
  if( !tableDataRows.empty())
  {
    updateVisibleColumns( columns, tableDataRows );
    computeHeaderRows( columns, m_tableLayout.getTrackerHeaderRows());
    populateColumnsFromTableData( columns, tableDataRows );
  }
  // dividesCells( columns );
  for( auto & column : columns )
  {
    findAndSetLongestColumnString( column );
  }

  computeAndBuildTableSeparator( columns, sectionSeparatingLine, topSeparator );
}

void TableTextFormatter::outputTable( std::ostringstream & tableOutput,
                                      std::vector< TableLayout::Column > & columns,
                                      string_view sectionSeparatingLine,
                                      string_view topSeparator ) const
{
  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  outputHeader( columns, tableOutput, sectionSeparatingLine );
  outputValues( columns, tableOutput, sectionSeparatingLine );
  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

void TableTextFormatter::computeHeaderRows( std::vector< TableLayout::Column > & columns,
                                            std::vector< TableLayout::Row > & headersRows ) const
{
  size_t maxLines = 1;
  for( auto & column : columns )
  {
    if( !column.subColumn.empty())
    {
      computeHeaderRows( column.subColumn, headersRows );
    }
  }
  maxLines = std::max( maxLines, column.columnName.lines.size());
  headersRows.insert( v.begin(), maxLines );
}

void TableTextFormatter::populateColumnsFromTableData( std::vector< TableLayout::Column > & columns,
                                                       RowsCellLayout & rowsCellsLayout,
                                                       RowsCellInput const & rowsCellsInput ) const
{
  // let's reserve the layout cells buffer
  RowsCellLayout initRowCellsLayout
  {
    std::vector< TableLayout::CellLayout >( rowsCellsInput.size() ),
    rowsCellsInput[0].size()
  };

  // to insert directly the values in each columns, we fill with the transposed rowsCells (row major->column major)
  for( size_t idxRow = 0; idxRow < rowsCellsInput.size(); ++idxRow )
  {
    GEOS_ERROR_IF( columns.size() != rowsCellsInput[idxRow].size(), "Dest matrix must have the number of rows equal to the number of columns in the" \
                                                                    "source matrix" );
    size_t maxLineCount = 0;

    for( size_t idxCol = 0; idxCol < rowsCellsInput[idxRow].size(); ++idxCol )
    {
      TableData::CellData const & cell = rowsCellsInput[idxRow][idxCol];

      TableLayout::CellAlignment const cellAlignement = columns[idxCol].cellAlignment;
      TableLayout::Alignment const alignement = cell.type == CellType::Header ?
                                                cellAlignement.headerAlignment :
                                                cellAlignement.valueAlignment;

      TableLayout::CellLayout & layoutCell = initRowCellsLayout[idxRow][idxCol];
      layoutCell = TableLayout::CellLayout( cell.type, cell.value, alignement );

      maxLineCount = max( maxLineCount, layoutCell.lines.size() );
    }
    m_tableLayout.m_valueRows.push_back( { maxLineCount } );
  }
}

void TableTextFormatter::findAndSetLongestColumnString( TableLayout & tableLayout,
                                                        RowsCellLayout & rowsCellsLayout ) const
{
  auto getMaxStringLen = [&]( std::vector< string > & lines )
  {
    if( lines.empty())
      return 0;
    return *std::max_element( lines.begin(),
                              lines.end(),
                              []( const auto & a, const auto & b )
    {
      return a.length() < b.length();
    } );
  };

  std::vector< size_t > subColumnsLength;
  for( auto it = tableLayout.begin(); it != tableLayout.end(); ++it )
  {
    Column * currentColumn = it;

    //1. process cells first
    size_t maxDataLength = 1;
    size_t idxLine = currentColumn - columnsLayout.begin();
    for( auto idxColumnData = 0; idxColumnData< rowsCellsLayout[idxLine].size(); idxColumnData++ )
    {
      TableLayout::CellLayout & cell = rowsCellsLayout[idxColumnData][idxLine];
      maxDataLength = std::max( maxDataLength, getMaxStringLen( cell.lines ));
    }
    currentColumn->setMaxStringSize( std::max(
                                       getMaxStringLen( currentColumn->columnName.lines ),
                                       maxDataLength )
                                     );

    //2. then header
    if( currentColumn.parent != nullptr )
    {
      // subcolumns case, adding B.A/B.B/...
      subColumnsLength.push_back( currentColumn->getMaxStringSize());

      if( currentColumn->next == nullptr )
      {
        while( currentColumn.parent != nullptr )
        {
          currentColumn->updateMaxStringSize( currentColumn.parent, subColumnsLength );
          currentColumn = currentColumn->parent;
        }
      }
      // detect 2 differents column
      while( currentColumn.parent != nullptr && currentColumn.parent != currentColumn->next.parent )
      {
        currentColumn = currentColumn->parent;
        currentColumn->updateMaxStringSize( currentColumn, subColumnsLength );
      }
      //clear if we detect that we change column
      if( currentColumn.parent == nullptr ||  currentColumn.parent != currentColumn->next.parent )
      {
        subColumnsLength.clear();
      }
    }
  }
}

void TableTextFormatter::computeAndBuildTableSeparator( TableLayout & tableLayout,
                                                        string & sectionSeparatingLine,
                                                        string & topSeparator ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );

  size_t const numColumns = tableLayout.getColumns().size() - 1;
  size_t const spacingBetweenColumns = numColumns * (size_t) columnMargin;
  size_t const margins = (size_t) borderMargin * 2;
  size_t sectionlineLength = spacingBetweenColumns + margins;
  for( auto const & column : tableLayout.getColumns() )
  {
    sectionlineLength += column.getMaxStringSize();
  }

  size_t maxTopLineLength =  tableTitle.length() + margins;
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );
  if( sectionlineLength < maxTopLineLength )
  {
    size_t const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( tableLayout, extraCharacters );
    sectionlineLength = maxTopLineLength;
  }

  sectionSeparatingLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
  integer const topSeparatorLength = maxTopLineLength - 2;   // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::increaseColumnsSize( TableLayout & tableLayout,
                                              size_t const extraCharacters ) const
{
  size_t const extraCharactersPerColumn = std::floor( (extraCharacters) / tableLayout.getColumns().size() );
  size_t const overflowCharacters = extraCharacters - (extraCharactersPerColumn *  tableLayout.getColumns().size() );
  for( auto it = tableLayout.begin(); it != tableLayout.end(); ++it )
  {
    Column * currentColumn = it;
    if( currentColumn.parent != nullptr )
    {
      if( currentColumn->next == nullptr )
      {
        while( currentColumn.parent != nullptr )
        {
          integer const newMaxStringSize = extraCharactersPerColumn +
                                           currentColumn.getMaxStringSize() +
                                           overflowCharacters;
          currentColumn->setMaxStringSize( newMaxStringSize );
          currentColumn = currentColumn->parent;
        }
      }

      while( currentColumn.parent != nullptr && currentColumn.parent != currentColumn->next.parent )
      {
        integer const newMaxStringSize = extraCharactersPerColumn +
                                         currentColumn.getMaxStringSize();
        currentColumn->setMaxStringSize( newMaxStringSize );
        currentColumn = currentColumn->parent;
      }
    }

    integer const newMaxStringSize =  currentColumn.getMaxStringSize() + extraCharactersPerColumn;
    newMaxStringSize += currentColumn == tableLayout.end() ?  overflowCharacters : 0;
    currentColumn->setMaxStringSize( extraCharactersPerColumn + currentColumn.getMaxStringSize() );
  }
}

void TableTextFormatter::outputTitleRow( std::ostringstream & tableOutput,
                                         string_view topSeparator ) const
{
  string const tableTitle = string( m_tableLayout.getTitle());
  if( !tableTitle.empty() )
  {
    tableOutput << GEOS_FMT( "{}\n", topSeparator );
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, m_tableLayout.getBorderMargin());
    tableOutput << buildCell( TableLayout::Alignment::center,
                              tableTitle,
                              (topSeparator.length() - (m_tableLayout.getBorderMargin() *  2)));
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, m_tableLayout.getBorderMargin() );
  }
}

void TableTextFormatter::CellFormatterStrategy::formatCellCommon( std::ostringstream & tableOutput, TableLayout::Column const & column,
                                                                  TableLayout const & tableLayout, TableLayout::CellLayout const & cell, size_t const idxRowCell,
                                                                  bool isFirstColumn, bool isNotLastColumn,
                                                                  string const & cellChar )
{
  const size_t cellSize = column.getMaxStringSize();

  if( isFirstColumn )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, cellChar.front() );
  }

  tableOutput << buildCell( cell.alignment, cell.dividedValues[idxRowCell], cellSize );

  if( isNotLastColumn )
  {
    tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, tableLayout.getColumnMargin());
  }
  else
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
    tableOutput << m_verticalLine << "\n";
  }
}


void TableTextFormatter::formatCell( std::ostringstream & tableOutput, TableLayout::Column const & column,
                                     TableLayout const & tableLayout, TableLayout::CellLayout const & cell,
                                     bool isFirstColumn, bool isNotLastColumn )
{
  for( auto const & line : cell.lines )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
    tableOutput << buildCell( cell.alignment,
                              line,
                              column.getMaxStringSize() );
    tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, tableLayout.getColumnMargin());
    tableOutput << m_verticalLine << "\n";
  }
}

void TableTextFormatter::outputHeader( std::vector< TableLayout::Column > const & columns,
                                       std::ostringstream & tableOutput,
                                       string_view sectionSeparatingLine ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  string const spaces =  string( columnMargin, ' ' );
  size_t nbRows =  columns[0].columnName.nbRows;


  SubColumnIterator * currentColumn = tableLayout.begin();
  for( int i = 0; i< nbRow; i++ )// ou pas juste incrmeenter une variable
  {
    currentColumn++;
    while( currentColumn.next != nullptr )
    {
      m_currentColumn++;
      currentColumn = currentColumn.next;
    }

  }

  for( auto it = tableLayout.begin(); it != tableLayout.end(); ++it )
  {
    std::vector< TableLayout::CellLayout > temp;

    TableLayout::CellLayout emptyCell = TableLayout::CellLayout( currentCoumn.type, "", currentCoumn.alignement );
    emptyCell.setMaxStringSize( currentCoumn.getMaxStringSize );
    temp.push_back( emptyCell );
  }

  for( )
  {
    formatCell( t )
  }
}

void TableTextFormatter::outputValues( std::vector< TableLayout::Column > & columns,
                                       RowsCellLayout & rowsCellsLayout,
                                       std::ostringstream & tableOutput,
                                       string_view sectionSeparatingLine ) const
{
  for( auto it = tableLayout.begin(); it != tableLayout.end(); ++it )
  {
    Column * currentColumn = it;
    size_t idxLine = currentColumn - columnsLayout.begin();
    for( auto idxColumnData = 0; idxColumnData< rowsCellsLayout[idxLine].size(); idxColumnData++ )
    {
      formatCell( rowsCellsLayout[idxColumnData][idxLine]; );
    }
  }
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
}

}
