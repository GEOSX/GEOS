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

  std::vector< std::vector< TableData::CellData > > const rowsValues( tableData.getTableDataRows() );
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
                           std::vector< std::vector< TableData::CellData > > & tableDataRows )
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
  std::vector< std::vector< TableData::CellData > > tableDataRows( tableData.getTableDataRows());
  if( !tableDataRows.empty())
  {
    updateVisibleColumns( columns, tableDataRows );
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

void populateSubColumnsFromTableData( TableLayout::Column & column,
                                      std::vector< std::vector< TableLayout::CellLayout > > const & valuesByColumn,
                                      size_t & idxColumn )
{
  size_t idxSubColumn = idxColumn;
  for( auto & subColumn : column.subColumn )
  {
    subColumn.cells = valuesByColumn[idxSubColumn++];
  }

  size_t nbSubColumns = column.subColumn.size();
  auto subColumnStartIter = valuesByColumn.begin() + idxColumn;
  std::vector< std::vector< TableLayout::CellLayout > > valuesBySubColumn( subColumnStartIter,
                                                                           subColumnStartIter + nbSubColumns );
  for( const auto & columnValues : valuesBySubColumn )
  { // add all subcolumn values in parent column
    column.cells.insert( column.cells.end(),
                         columnValues.begin(),
                         columnValues.end() );
  }
  idxColumn += nbSubColumns;
}

void TableTextFormatter::
  populateColumnsFromTableData( std::vector< TableLayout::Column > & columns,
                                std::vector< std::vector< TableLayout::CellLayout > > rowsCellsLayout,
                                std::vector< std::vector< TableData::CellData > > const & rowsCellsInput ) const
{
  size_t currentColumn = 0;
  bool containSubColumn=  false;

  // let's reserve the layout cells buffer
  std::vector< std::vector< TableLayout::CellLayout > > rowCellsLayout
  {
    std::vector< TableLayout::CellLayout >( tableDataRows.size() ),
    tableDataRows[0].size()
  };

  // to insert directly the values in each columns, we fill with the transposed rowsCells (row major->column major)
  for( size_t idxRow = 0; idxRow < rowsCells.size(); ++idxRow )
  {
    GEOS_ERROR_IF( columns.size() != rowsCells[idxRow].size(), "Dest matrix must have the number of rows equal to the number of columns in the" \
                                                               "source matrix" );
    size_t maxLineCount = 0;

    for( size_t idxCol = 0; idxCol < rowsCells[idxRow].size(); ++idxCol )
    {
      TableData::CellData const & cell = rowsCells[idxRow][idxCol];

      TableLayout::CellAlignment const cellAlignement = columns[idxCol].cellAlignment;
      TableLayout::Alignment const alignement = cell.type == CellType::Header ?
                                                cellAlignement.headerAlignment :
                                                cellAlignement.valueAlignment;

      TableLayout::CellLayout & layoutCell = rowCellsLayout[idxRow][idxCol];
      layoutCell = TableLayout::CellLayout( cell.type, cell.value, alignement );
      maxLineCount = max( maxLineCount, layoutCell.lines.size() );
    }

    m_rows.push_back( maxLineCount );
  }

  for( auto & column : columns )
  {
    if( column.subColumn.empty())
    {
      // column.cells = valuesByColumn[currentColumn++];
      for( auto & cell : column.cells )  //todo
      {
        // cell.alignment = column.cellAlignment.valueAlignment;
      }
    }
    else
    {
      containSubColumn = true;
      populateSubColumnsFromTableData( column, valuesByColumn, currentColumn );
    }
  }

  if( containSubColumn )
  {
    for( auto & column : columns )
    {
      if( column.subColumn.empty())
      {
        column.subColumn = std::vector< TableLayout::Column >{TableLayout::Column()};
      }
    }
  }
}


void TableTextFormatter::dividesCells( std::vector< TableLayout::Column > & columns ) const
{
  size_t maxNbRow = 1;
  size_t maxCellRow = 1;

  auto splitValue = []( const std::string & value, std::vector< std::string > & dividedValue ) {
      std::istringstream ss( value );
      std::string subValue;
      while( getline( ss, subValue, '\n' ))
      {
        dividedValue.push_back( subValue );
      }
    };

  for( auto & column : columns )
  {
    splitValue( column.columnName.dividedValues );
    maxNbRow = std::max( maxNbRow, column.columnName.dividedValues.size());

    for( auto & cell : column.cells )
    {
      splitValue( cell.value, cell.dividedValues );
      maxCellRow = std::max( maxCellRow, cell.dividedValues.size());
    }

    if( !column.subColumn.empty())
    {
      dividesCells( column.subColumn );
    }
  }

  // set the same size for all cells
  for( auto & column : columns )
  {

    if( column.columnName.dividedValues.size() < maxNbRow )
    {
      column.columnName.dividedValues.resize( maxNbRow, " " );
    }
    column.columnName.nbRows = maxNbRow;

    for( auto & cell : column.cells )
    {
      if( cell.dividedValues.size() < maxCellRow )
      {
        cell.dividedValues.resize( maxCellRow, " " );

      }
      cell.nbRows = maxCellRow;
    }
  }
}

void TableTextFormatter::findAndSetLongestColumnString( TableLayout::Column & column ) const
{
  {   // header case
    auto maxCellIt = *std::max_element( column.columnName.dividedValues.begin(),
                                        column.columnName.dividedValues.end(),
                                        []( const auto & a, const auto & b )
      {
        return a.length() < b.length();
      } );
    column.setMaxStringSize( maxCellIt.length());
  }

  {   // cells case
    if( column.subColumn.empty() && !column.cells.empty())
    {
      auto maxCellIt = std::max_element( column.cells.begin(),
                                         column.cells.end(),
                                         []( const auto & a, const auto & b )
        {
          return a.value.length() < b.value.length();
        } );

      if( column.getMaxStringSize() < maxCellIt->value.length() )
      {
        column.setMaxStringSize( maxCellIt->value.length() );
      }
    }
  }

  {   // subcolumn values case
    size_t totalLengthSubColumns = 0;
    if( !column.subColumn.empty() )
    {
      for( auto & subColumn : column.subColumn )
      {
        findAndSetLongestColumnString( subColumn );  //todo
        totalLengthSubColumns += subColumn.getMaxStringSize();
        if( &subColumn != &column.subColumn[0] )
        {
          totalLengthSubColumns += m_tableLayout.getColumnMargin();
        }
      }
    }

    if( totalLengthSubColumns > column.getMaxStringSize() )
    {
      column.setMaxStringSize( totalLengthSubColumns );
    }
  }

}

void TableTextFormatter::computeAndBuildTableSeparator( std::vector< TableLayout::Column > & columns,
                                                        string & sectionSeparatingLine,
                                                        string & topSeparator ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );

  size_t const numColumns = columns.size() - 1;
  size_t const spacingBetweenColumns = numColumns * (size_t) columnMargin;
  size_t const margins = (size_t) borderMargin * 2;
  size_t sectionlineLength = spacingBetweenColumns + margins;
  for( auto const & column : columns )
  {
    sectionlineLength += column.getMaxStringSize();
  }

  size_t maxTopLineLength =  tableTitle.length() + margins;
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );
  if( sectionlineLength < maxTopLineLength )
  {
    size_t const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( columns, extraCharacters );
    sectionlineLength = maxTopLineLength;
  }

  sectionSeparatingLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
  integer const topSeparatorLength = maxTopLineLength - 2;   // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::increaseColumnsSize( std::vector< TableLayout::Column > & columns,
                                              size_t const extraCharacters ) const
{
  size_t const extraCharactersPerColumn = std::floor( (extraCharacters) / columns.size() );
  size_t const overflowCharacters = extraCharacters - (extraCharactersPerColumn *  columns.size() );
  for( auto & column : columns )
  {
    if( !column.subColumn.empty())
    {
      increaseColumnsSize( column.subColumn, extraCharactersPerColumn );
    }

    size_t const cellSize = column.getMaxStringSize();
    integer const newMaxStringSize = &column == &columns[0] ?
                                     extraCharactersPerColumn + cellSize + overflowCharacters :
                                     extraCharactersPerColumn + cellSize;
    column.setMaxStringSize( newMaxStringSize );
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

void TableTextFormatter::MergingCell::formatCell( std::ostringstream & tableOutput, TableLayout::Column const & column,
                                                  TableLayout const & tableLayout, TableLayout::CellLayout const & cell, size_t const idxRowCell,
                                                  bool isFirstColumn, bool isNotLastColumn )
{
  formatCellCommon( tableOutput, column, tableLayout, cell, idxRowCell, isFirstColumn, isNotLastColumn, cell.value );
  tableOutput << string( tableLayout.getColumnMargin(), cell.value.front() );
}

void TableTextFormatter::SeparatingCell::formatCell( std::ostringstream & tableOutput, TableLayout::Column const & column,
                                                     TableLayout const & tableLayout, TableLayout::CellLayout const & cell, size_t const idxRowCell,
                                                     bool isFirstColumn, bool isNotLastColumn )
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

void TableTextFormatter::outputSubValues( std::vector< TableLayout::Column > const & columns,
                                          std::ostringstream & tableOutput,
                                          size_t idxRow ) const
{
  for( auto const & column : columns )
  {
    const auto & cell = column.cells.at( idxRow );
    tableOutput << buildCell( cell.alignment, cell.value, column.getMaxStringSize());

    if( &column < &columns.back())
    {
      tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, m_tableLayout.getColumnMargin());
    }
  }
}

void TableTextFormatter::outputHeader( std::vector< TableLayout::Column > const & columns,
                                       std::ostringstream & tableOutput,
                                       string_view sectionSeparatingLine ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  string const spaces =  string( columnMargin, ' ' );
  size_t nbRows =  columns[0].columnName.nbRows;
  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << m_verticalLine;

    for( auto & column : columns )
    {
      bool isNotLastColumn = &column < &columns.back();
      bool isFirstColumn = &column == &columns.front();
      auto formatter = createCellFormatter( static_cast< TableData::CellType >(column.columnName.type));
      formatter->formatCell( tableOutput, column, m_tableLayout, column.columnName, idxRow, isFirstColumn, isNotLastColumn );
    }
  }

  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  if( columns[0].subColumn.size() > 0 )
  {
    std::vector< TableLayout::Column > rowSubColumns;
    for( auto const & column : columns )
    {
      rowSubColumns.insert( rowSubColumns.end(), column.subColumn.begin(), column.subColumn.end());
    }
    outputHeader( rowSubColumns, tableOutput, sectionSeparatingLine );
  }
}

void TableTextFormatter::outputValues( std::vector< TableLayout::Column > const & columns,
                                       std::ostringstream & tableOutput,
                                       string_view sectionSeparatingLine ) const
{
  size_t const nbRows = columns[0].cells.size();
  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << m_verticalLine;

    for( auto & column : columns )
    {
      bool isNotLastColumn = &column < &columns.back();
      bool isFirstColumn = &column == &columns.front();
      if( column.subColumn.size() > 1 )
      {
        outputSubValues( column.subColumn, tableOutput, idxRow );

        if( isNotLastColumn )
        {
          tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, m_tableLayout.getColumnMargin() );
        }
      }
      else
      {
        size_t nbCellRow = column.cells[idxRow].nbRows;
        for( size_t idxCellRow = 0; idxCellRow < nbCellRow; ++idxCellRow )
        {
          auto formatter = createCellFormatter( static_cast< TableData::CellType >(column.cells[idxRow].type));
          formatter->formatCell( tableOutput, column, m_tableLayout, column.cells[idxRow], idxCellRow,
                                 isFirstColumn, isNotLastColumn );
        }
      }
    }
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }
}

}
