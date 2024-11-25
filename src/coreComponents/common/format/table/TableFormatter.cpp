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

string TableCSVFormatter::headerToString() const //todo
{
   std::stringstream oss;
  // static constexpr string_view separator = ",";

  // for( std::size_t idxColumn = 0; idxColumn < m_tableLayout.getColumns().size(); ++idxColumn )
  // {
  //   oss << m_tableLayout.getColumns()[idxColumn].columnName.value;
  //   if( idxColumn < m_tableLayout.getColumns().size() - 1 )
  //   {
  //     oss << separator;
  //   }
  // }
  // oss << "\n";
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

TableTextFormatter::TableTextFormatter( TableLayout const & tableLayout ):
  TableFormatter( tableLayout )
{}

string TableTextFormatter::toString() const
{
  std::ostringstream tableOutput;
  TableData tableData;
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( tableData, sectionSeparatingLine, topSeparator );
  tableOutput << '\n';
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  //outputHeader( tableOutput, sectionSeparatingLine ); todo

  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  RowsCellLayout cellsHeaderLayout;
  RowsCellLayout cellsDataLayout;
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( tableData,
                        cellsHeaderLayout, cellsDataLayout,
                        sectionSeparatingLine, topSeparator );
  std::ostringstream tableOutput;
  outputTable( tableOutput, sectionSeparatingLine, topSeparator );
  return tableOutput.str();
}

void TableTextFormatter::prepareAndBuildTable( TableData const & tableData,
                                               RowsCellLayout & cellsDataLayout,
                                               RowsCellLayout & cellsHeaderLayout,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{
  if( !inputDataRows.empty())
  {
    updateVisibleColumns( inputDataRows );
    populateDataCellsLayout( cellsDataLayout, tableData );
  }

  gridifyHeaders();

  for( auto & column : m_tableLayout.getColumns() )
  {
    findAndSetLongestColumnString( column, cellsDataLayout );
  }

  computeAndBuildTableSeparator( sectionSeparatingLine, topSeparator );

  populateHeaderCellsLayout( cellsHeaderLayout );
}

void TableTextFormatter::outputTable( std::ostringstream & tableOutput,
                                      RowsCellLayout const & cellsHeader,
                                      RowsCellLayout const & cellsData,
                                      string_view sectionSeparatingLine,
                                      string_view topSeparator ) const
{
  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputLines( headerLine, tableOutput, sectionSeparatingLine );
  outputLines( cellsData, tableOutput, sectionSeparatingLine );

  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

/**
 * @brief Detect columns who are not displayed from TableLayout and therefore modify columns / inputDataRows vectors
 * @param inputDataRows Vector built in TableData containing all rows values
 */
void TableTextFormatter::updateVisibleColumns( RowsCellInput & inputDataRows ) const
{
  integer idxColumn = 0;
  auto columns =  m_tableLayout.getColumns();
  for( auto iterColumn = columns.begin(); iterColumn !=  columns.end(); )
  {
    if( !iterColumn->enabled )
    {
      iterColumn = columns.erase( iterColumn );
      for( auto & row : inputDataRows )
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

void TableTextFormatter::gridifyHeaders()
{
  size_t countNbLayer = 1;
  size_t nbHeaderDepth = m_tableLayout.getMaxHeaderRow();
  for( auto it = m_tableLayout.beginTop(); it != m_tableLayout.endTop(), ++it )
  {
    if( !it->hasChild() )
    {
      for( size_t idxLayer = countNbLayer; idxLayer < nbHeaderDepth; idxLayer++ )
      {
        Column const * column = it;
        column.addSubColumn( "" );
      }
      countNbLayer = 1;
    }
    else
    {
      countNbLayer++;
    }
  }
}

void TableTextFormatter::populateDataCellsLayout( RowsCellLayout & cellsDataLayout,
                                                  TableData const & tableData ) const
{
  RowsCellInput inputDataRows( tableData.getTableDataRows());
  // let's reserve the layout cells buffer
  cellsDataLayout =
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
      TableLayout::Alignment const alignement = cell.cellType == CellType::Header ?
                                                cellAlignement.headerAlignment :
                                                cellAlignement.valueAlignment;

      cellsDataLayout[idxRow][idxCol] = TableLayout::CellLayout( cell.cellType, cell.value, alignement );
      maxLineCount = max( maxLineCount, layoutCell.lines.size() );
    }
    m_tableLayout.m_valueRows.push_back( { maxLineCount } );
  }
}

void TableTextFormatter::findAndSetLongestColumnString( RowsCellLayout & cellsDataLayout ) const
{
  auto getMaxStringLen = [&]( std::vector< string > & lines )
  {
    if( lines.empty())
      return 0;
    return *std::max_element( lines.begin(), lines.end(), []( const auto & a, const auto & b )
    {
      return a.length() < b.length();
    } );
  };

  for( auto it = m_tableLayout.beginLeaf(); it != m_tableLayout.endLeaf(); ++it )
  {
    size_t maxSize = 1;
    Column * currentColumn = it;
    if( !it->hasChild() )
    {
      //1. currentColumn = max(subcolumn, dataColumn)
      size_t idxLine = currentColumn - columnsLayout.begin();
      for( auto idxColumnData = 0; idxColumnData< cellsDataLayout[idxLine].size(); idxColumnData++ )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[idxColumnData][idxLine];
        cell.maxDataLength = std::max( cell.maxDataLength, getMaxStringLen( cell.lines ));
        maxSize = std::max( cell.maxDataLength, maxSize );
      }

      currentColumn->setMaxStringSize( std::max( getMaxStringLen( currentColumn->columnName.lines ),
                                                 maxSize ));
    }
    else //2. and max(parent, subcolumn)
    {
      if( !it->hasParent())
      {
        size_t maxSize = std::max( getMaxStringLen( currentColumn->columnName.lines ),
                                   maxSize );
        currentColumn->columnName.maxDataLength = maxSize;
        currentColumn->setMaxStringSize( maxSize );
      }
      else
      {
        std::vector< size_t > subColumnLength;
        for( auto const & subColumn : currentColumn->subColumn )
        {
          subColumnLength.push_back( subcolumn->getMaxStringSize() );
        }
        currentColumn->compareAndSetMaxStringSize( currentColumn, subColumnsLength );
        currentColumn->columnName.maxDataLength = currentColumn.getMaxStringSize()
      }
    }
  }
}

void TableTextFormatter::computeAndBuildTableSeparator( string & sectionSeparatingLine,
                                                        string & topSeparator ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );
  size_t const numColumns = m_tableLayout.getColumns().size() - 1;

  size_t const spacingBetweenColumns = numColumns * (size_t) columnMargin;
  size_t const margins = (size_t) borderMargin * 2;
  size_t sectionlineLength = spacingBetweenColumns + margins;
  for( auto const & column : m_tableLayout.getColumns() )
  {
    sectionlineLength += column.getMaxStringSize();
  }

  size_t maxTopLineLength =  tableTitle.length() + margins;
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );
  if( sectionlineLength < maxTopLineLength )
  {
    size_t const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( m_tableLayout, extraCharacters );
    sectionlineLength = maxTopLineLength;
  }

  sectionSeparatingLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
  integer const topSeparatorLength = maxTopLineLength - 2;   // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::populateHeaderCellsLayout(
  RowsCellLayout & cellsHeaderLayout ) const
{
  for( const auto column : columns )
  {
    size_t idxRow = 0;
    for( auto it = column.beginRoot(); it< column.endRoot(); ++it )
    {
      cellsHeaderLayout[idxRow].emplace_back( it->columnName );
      idxRow = it->hasChildren() ? idxRow++; idxRow--;
    }
  }
}

void TableTextFormatter::increaseColumnsSize( size_t const extraCharacters ) const
{
  size_t const columnsCount = m_tableLayout.getColumns().size();
  size_t const extraCharactersPerColumn = std::floor( extraCharacters / columnsCount );
  size_t const overflowCharacters = extraCharacters - (extraCharactersPerColumn * columnsCount);

  for( auto it = m_tableLayout.beginLeaf(); it != m_tableLayout.endTop(); ++it )
  {
    Column * currentColumn = it;

    size_t divider = currentColumn->hasParent() ? currentColumn->parent.subColumn.size() : columnsCount;
    size_t extraCharactersForCurrentColumn = std::floor( extraCharacters / divider );
    size_t overflowForCurrentColumn = (currentColumn->next == nullptr) ? overflowCharacters : 0;

    integer const newMaxStringSize = currentColumn->getMaxStringSize() +
                                     extraCharactersForCurrentColumn +
                                     overflowForCurrentColumn;
    currentColumn->setMaxStringSize( newMaxStringSize );
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

// void TableTextFormatter::CellFormatterStrategy::formatCellCommon( std::ostringstream & tableOutput, TableLayout::Column const & column,
//                                                                   TableLayout const & tableLayout, TableLayout::CellLayout const & cell,
// size_t const idxRowCell,
//                                                                   bool isFirstColumn, bool isNotLastColumn,
//                                                                   string const & cellChar )
// {
//   const size_t cellSize = column.getMaxStringSize();

//   if( isFirstColumn )
//   {
//     tableOutput << string( tableLayout.getBorderMargin() - 1, cellChar.front() );
//   }

//   tableOutput << buildCell( cell.alignment, cell.dividedValues[idxRowCell], cellSize );

//   if( isNotLastColumn )
//   {
//     tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, tableLayout.getColumnMargin());
//   }
//   else
//   {
//     tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
//     tableOutput << m_verticalLine << "\n";
//   }
}


void TableTextFormatter::formatCell( std::ostringstream & tableOutput,
                                     TableLayout::CellLayout const & cell ) const
{
  for( auto const & line : cell.lines )
  {
    tableOutput << string( m_tableLayout.getBorderMargin() - 1, ' ' );
    tableOutput << buildCell( cell.alignment,
                              line,
                              cell.maxDataLength );
    tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, m_tableLayout.getColumnMargin());
    tableOutput << m_verticalLine << "\n";
  }
}

void TableTextFormatter::outputLines( RowsCellLayout & cellsLayout,
                                      std::ostringstream & tableOutput,
                                      string_view sectionSeparatingLine ) const
{
  for( auto const line : cellsLayout )
  {
    for( auto const cell : line )
    {
      formatCell( tableOutput, cell );
    }
  }
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
}
