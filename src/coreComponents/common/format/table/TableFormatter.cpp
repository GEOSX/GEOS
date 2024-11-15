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

  std::vector< std::vector< TableData::DataType > > const rowsValues( tableData.getTableDataRows() );
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

void transpose( std::vector< std::vector< TableLayout::Cell > > & dest,
                std::vector< std::vector< TableData::DataType > > const & source )
{

  for( size_t idxRow = 0; idxRow < source.size(); ++idxRow )
  {
    GEOS_ERROR_IF( dest.size() != source[idxRow].size(), "Dest matrix must have the number of rows equal to the number of columns in the" \
                                                         "source matrix" );
    for( size_t idxCol = 0; idxCol < source[idxRow].size(); ++idxCol )
    {
      dest[idxCol][idxRow].type = (char) source[idxRow][idxCol].type;
      dest[idxCol][idxRow].value = source[idxRow][idxCol].value;
    }
  }
}

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
                           std::vector< std::vector< TableData::DataType > > & tableDataRows )
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
  outputHeaderSectionRows( columns, tableOutput, sectionSeparatingLine );

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
  std::vector< std::vector< TableData::DataType > > tableDataRows( tableData.getTableDataRows());
  if( !tableDataRows.empty())
  {
    updateVisibleColumns( columns, tableDataRows );
    populateColumnsFromTableData( columns, tableDataRows );
  }
  splitAndMergeColumnHeaders( columns );
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
  outputHeaderSectionRows( columns, tableOutput, sectionSeparatingLine );
  outputValuesSectionRows( columns, tableOutput, sectionSeparatingLine );
  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

void populateSubColumnsFromTableData( TableLayout::Column & column,
                                      std::vector< std::vector< TableLayout::Cell > > const & valuesByColumn,
                                      size_t & idxColumn )
{
  size_t idxSubColumn = idxColumn;
  for( auto & subColumn : column.subColumn )
  {
    subColumn.cells = valuesByColumn[idxSubColumn++];
  }

  size_t nbSubColumns = column.subColumn.size();
  auto subColumnStartIter = valuesByColumn.begin() + idxColumn;
  std::vector< std::vector< TableLayout::Cell > > valuesBySubColumn( subColumnStartIter,
                                                                     subColumnStartIter + nbSubColumns );
  for( const auto & columnValues : valuesBySubColumn )
  { // add all subcolumn values in parent column
    column.cells.insert( column.cells.end(),
                         columnValues.begin(),
                         columnValues.end() );
  }
  idxColumn += nbSubColumns;
}

void TableTextFormatter::populateColumnsFromTableData( std::vector< TableLayout::Column > & columns,
                                                       std::vector< std::vector< TableData::DataType > > const & tableDataRows ) const
{
  size_t currentColumn = 0;
  std::vector< std::vector< TableLayout::Cell > > valuesByColumn( tableDataRows[0].size(),
                                                                  std::vector< TableLayout::Cell >( tableDataRows.size()));

  // to insert directly the values in each columns, we fill with the transposed tableDataRows (row major->column major)
  transpose( valuesByColumn, tableDataRows );

  for( auto & column : columns )
  {
    if( column.subColumn.empty())
    {
      column.cells = valuesByColumn[currentColumn++];
    }
    else
    {
      populateSubColumnsFromTableData( column, valuesByColumn, currentColumn );
    }
  }
}


void TableTextFormatter::splitAndMergeColumnHeaders( std::vector< TableLayout::Column > & columns ) const
{
  size_t maxNbRow = 1;
  for( auto & column : columns )
  {
    std::vector< string > & headerDivided = column.columnName.dividedValue;
    std::istringstream ss( column.columnName.value );
    string subColumnNames;
    while( getline( ss, subColumnNames, '\n' ))
    {
      headerDivided.push_back( subColumnNames );
    }

    if( headerDivided.size() > maxNbRow )
    {
      maxNbRow = headerDivided.size();
    }

    if( !column.subColumn.empty())
    {
      splitAndMergeColumnHeaders( column.subColumn );
    }
  }
  for( auto & column : columns )
  {
    std::vector< string > & headerDivided = column.columnName.dividedValue;

    if( headerDivided.size() < maxNbRow )
    {
      headerDivided.resize( maxNbRow, " " );

    }
    column.columnName.nbRows = maxNbRow;
  }

}

void TableTextFormatter::findAndSetLongestColumnString( TableLayout::Column & column ) const
{
  { // header case
    auto maxCellIt = *std::max_element( column.columnName.dividedValue.begin(),
                                        column.columnName.dividedValue.end(),
                                        []( const auto & a, const auto & b )
    {
      return a.length() < b.length();
    } );
    column.setMaxStringSize( maxCellIt.length());
  }

  {  // cells case
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

  { // subcolumn values case
    size_t totalLengthSubColumns = 0;
    if( !column.subColumn.empty() )
    {
      for( auto & subColumn : column.subColumn )
      {
        findAndSetLongestColumnString( subColumn );
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
  integer const topSeparatorLength = maxTopLineLength - 2; // Adjust for border characters
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

void TableTextFormatter::outputHeaderSectionRows( std::vector< TableLayout::Column > const & columns,
                                                  std::ostringstream & tableOutput,
                                                  string_view sectionSeparatingLine ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const spaces =  string( columnMargin, ' ' );
  bool containSubColumn = false;
  size_t nbRows =  columns[0].columnName.nbRows;
  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, borderMargin );

    for( auto const & column : columns )
    {
      bool isNotLastColumn = &column < &columns.back();
      //bool isFirstColumn = &column == &columns.front();
      string const cell = column.columnName.dividedValue.size() == 0 ? ""
                          : column.columnName.dividedValue.at( idxRow );
      size_t const cellSize = column.getMaxStringSize();
      tableOutput << buildCell( column.cellAlignment.headerAlignment,
                                cell,
                                cellSize );
      // Add space between columns
      if( isNotLastColumn )
      {
        tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
      }

      if( !column.subColumn.empty())
      {
        containSubColumn = true; // TOODODOD
      }

    }
    // Append right border with line return
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, borderMargin );
  }


  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  // Check and build subrow header
  if( containSubColumn )
  {
    std::vector< TableLayout::Column > rowSubColumns;

    for( auto const & column : columns )
    {
      if( column.subColumn.empty())
      {
        rowSubColumns.push_back( {TableLayout::Column().setMaxStringSize( column.getMaxStringSize() )} );
      }
      else
      {
        rowSubColumns.insert( rowSubColumns.end(),
                              column.subColumn.begin(),
                              column.subColumn.end());
      }
    }
    outputHeaderSectionRows( rowSubColumns, tableOutput, sectionSeparatingLine );

  }
}

void TableTextFormatter::ValueCell::formatCell( std::ostringstream & tableOutput,
                                                TableLayout::Column const & column,
                                                TableLayout const & tableLayout,
                                                TableLayout::Cell & cell, bool isFirstColumn, bool isNotLastColumn )
{
  const size_t cellSize = column.getMaxStringSize();
  if( isFirstColumn )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
  }

  tableOutput << buildCell( column.cellAlignment.valueAlignment, cell.value, cellSize );

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

void TableTextFormatter::HeaderCell::formatCell( std::ostringstream & tableOutput,
                                                 TableLayout::Column const & column,
                                                 TableLayout const & tableLayout,
                                                 TableLayout::Cell & cell, bool isFirstColumn, bool isNotLastColumn )
{
  const size_t cellSize = column.getMaxStringSize();
  if( isFirstColumn )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
  }

  tableOutput << buildCell( column.cellAlignment.headerAlignment, cell.value, cellSize );

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

void TableTextFormatter::MergingCell::formatCell( std::ostringstream & tableOutput,
                                                  TableLayout::Column const & column,
                                                  TableLayout const & tableLayout,
                                                  TableLayout::Cell & cell, bool isFirstColumn, bool isNotLastColumn )
{
  const size_t cellSize = column.getMaxStringSize();
  if( isFirstColumn )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
  }
  cell.value = " ";
  tableOutput << buildCell( column.cellAlignment.headerAlignment, cell.value, cellSize );
  tableOutput << string( tableLayout.getColumnMargin(), ' ' );
  if( !isNotLastColumn )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
    tableOutput << m_verticalLine << "\n";
  }

}

void TableTextFormatter::SeparatingCell::formatCell( std::ostringstream & tableOutput,
                                                     TableLayout::Column const & column,
                                                     TableLayout const & tableLayout,
                                                     TableLayout::Cell & cell, bool isFirstColumn, bool isNotLastColumn )
{
  const size_t cellSize = column.getMaxStringSize();
  if( isFirstColumn )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, '-' );
  }

  if( isNotLastColumn )
  {
    cell.value = string( cellSize + tableLayout.getColumnMargin(), '-' );
    tableOutput << buildCell( column.cellAlignment.headerAlignment, cell.value, cellSize );
  }
  else
  {
    cell.value = string( cellSize + tableLayout.getBorderMargin(), '-' );
    tableOutput << buildCell( column.cellAlignment.headerAlignment, cell.value, cellSize );
    tableOutput << '\n';
  }
}


void TableTextFormatter::outputCell( std::ostringstream & tableOutput,
                                     TableLayout::Column const & column,
                                     size_t const idxRow,
                                     bool isFirstColumn,
                                     bool isNotLastColumn ) const
{
  TableLayout::Cell cell = column.cells.at( idxRow );

  std::unique_ptr< TableTextFormatter::CellFormatterStrategy > formatter;

  if( cell.type == (char) TableData::CellType::MERGE )
  {
    formatter = std::make_unique< TableTextFormatter::MergingCell >();
  }
  else if( cell.type == (char)TableData::CellType::SEPARATOR )
  {
    formatter = std::make_unique< TableTextFormatter::SeparatingCell >();
  }
  else if( cell.type == (char) TableData::CellType::Header )
  {
    formatter = std::make_unique< TableTextFormatter::ValueCell >();
  }
  else if( cell.type == (char) TableData::CellType::Value )
  {
    formatter = std::make_unique< TableTextFormatter::ValueCell >();
  }

  formatter->formatCell( tableOutput, column, m_tableLayout, cell, isFirstColumn, isNotLastColumn );
}

void TableTextFormatter::outputSubSection( std::vector< TableLayout::Column > const & columns,
                                           std::ostringstream & tableOutput,
                                           size_t const idxRow ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  for( auto const & column : columns )
  {
    string const cell = column.cells.at( idxRow ).value;
    size_t const cellSize =  column.getMaxStringSize();
    tableOutput << buildCell( column.cellAlignment.valueAlignment,
                              cell,
                              cellSize );
    if( &column < &columns.back() )
    {
      tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
    }
  }
}

void TableTextFormatter::outputValuesSectionRows( std::vector< TableLayout::Column > const & columns,
                                                  std::ostringstream & tableOutput,
                                                  string_view sectionSeparatingLine ) const
{
  size_t const nbRows = columns[0].cells.size();
  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << m_verticalLine;

    for( auto const & column : columns )
    {
      bool isNotLastColumn = &column < &columns.back();
      bool isFirstColumn = &column == &columns.front();
      if( !column.subColumn.empty())
      {
        outputSubSection( column.subColumn, tableOutput, idxRow );

        if( isNotLastColumn )
        {
          tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, m_tableLayout.getColumnMargin() );
        }
      }
      else
      {
        outputCell( tableOutput, column, idxRow, isFirstColumn, isNotLastColumn );
      }
    }
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }
}

}
