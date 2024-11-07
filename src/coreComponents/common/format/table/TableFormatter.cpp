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
    oss << m_tableLayout.getColumns()[idxColumn].column.columnName;
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

  std::vector< std::vector< string > > const rowsValues( tableData.getTableDataRows() );
  std::ostringstream oss;

  for( const auto & row : rowsValues )
  {
    oss << stringutilities::join( row.cbegin(), row.cend(), "," ) << "\n";
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

void transpose( std::vector< std::vector< string > > & dest,
                std::vector< std::vector< string > > const & source )
{

  for( size_t idxRow = 0; idxRow < source.size(); ++idxRow )
  {
    GEOS_ERROR_IF( dest.size() != source[idxRow].size(), "Dest matrix must have the number of rows equal to the number of columns in the" \
                                                         "source matrix" );
    for( size_t idxCol = 0; idxCol < source[idxRow].size(); ++idxCol )
    {
      dest[idxCol][idxRow] = source[idxRow][idxCol];
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
 * @brief Detect tableColumnsData who are not displayed from TableLayout and therefore modify tableColumnsData / tableDataRows vectors
 * @param tableColumnsData Vector built in TableLayout containing all tableColumnsData with their parameters
 * @param tableDataRows Vector built in TableData containing all rows values
 */
void updateVisibleColumns( std::vector< TableLayout::ColumnStructure > & tableColumnsData,
                           std::vector< std::vector< string > > & tableDataRows )
{
  integer idxColumn = 0;
  for( auto iterColumn = tableColumnsData.begin(); iterColumn != tableColumnsData.end(); )
  {
    if( !iterColumn->column.enabled )
    {
      iterColumn = tableColumnsData.erase( iterColumn );
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
  std::vector< TableLayout::ColumnStructure > tableColumnsData( m_tableLayout.getColumns());
  TableData tableData;
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( tableColumnsData, tableData, sectionSeparatingLine, topSeparator );

  tableOutput << '\n';
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  outputHeaderSectionRows( tableColumnsData, tableOutput, sectionSeparatingLine );

  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::vector< TableLayout::ColumnStructure > tableColumnsData( m_tableLayout.getColumns());
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( tableColumnsData, tableData, sectionSeparatingLine, topSeparator );

  std::ostringstream tableOutput;
  outputTable( tableOutput, tableColumnsData, sectionSeparatingLine, topSeparator );

  return tableOutput.str();
}

void TableTextFormatter::prepareAndBuildTable( std::vector< TableLayout::ColumnStructure > & tableColumnsData,
                                               TableData const & tableData,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{
  std::vector< std::vector< string > > tableDataRows( tableData.getTableDataRows());
  if( !tableDataRows.empty())
  {
    updateVisibleColumns( tableColumnsData, tableDataRows );
    populateColumnsFromTableData( tableColumnsData, tableDataRows );
  }

  std::vector< std::vector< string > > allDividedHeaderParts;
  splitAndMergeColumnHeaders( tableColumnsData, allDividedHeaderParts );

  for( auto & tableColumnData : tableColumnsData )
  {
    findAndSetLongestColumnString( tableColumnData );
  }

  computeAndBuildTableSeparator( tableColumnsData, sectionSeparatingLine, topSeparator );
}

void TableTextFormatter::outputTable( std::ostringstream & tableOutput,
                                      std::vector< TableLayout::ColumnStructure > & tableColumnsData,
                                      string_view sectionSeparatingLine,
                                      string_view topSeparator ) const
{
  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputHeaderSectionRows( tableColumnsData, tableOutput, sectionSeparatingLine );

  outputValuesSectionRows( tableColumnsData, tableOutput, sectionSeparatingLine );
  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

void populateSubColumnsFromTableData( TableLayout::ColumnStructure & tableColumnData,
                                      std::vector< std::vector< string > > const & valuesByColumn,
                                      size_t & idxColumn )
{
  size_t idxSubColumn = idxColumn;
  for( auto & subColumn : tableColumnData.subColumn )
  {
    subColumn.columnValues = valuesByColumn[idxSubColumn++];
  }

  size_t nbSubColumns = tableColumnData.column.subColumnNames.size();
  auto subColumnStartIter = valuesByColumn.begin() + idxColumn;
  std::vector< std::vector< string > > valuesBySubColumn( subColumnStartIter,
                                                          subColumnStartIter + nbSubColumns );
  for( const auto & columnValues : valuesBySubColumn )
  { // add all subcolumn values in parent tableColumnData
    tableColumnData.columnValues.insert( tableColumnData.columnValues.end(),
                                         columnValues.begin(),
                                         columnValues.end() );
  }
  idxColumn += nbSubColumns;
}

void TableTextFormatter::populateColumnsFromTableData( std::vector< TableLayout::ColumnStructure > & tableColumnsData,
                                                       std::vector< std::vector< string > > const & tableDataRows ) const
{
  size_t currentColumn = 0;
  std::vector< std::vector< string > > valuesByColumn( tableDataRows[0].size(),
                                                       std::vector< string >( tableDataRows.size()));

  // to insert directly the values in each columns, we fill with the transposed tableDataRows (row major->column major)
  transpose( valuesByColumn, tableDataRows );

  for( auto & tableColumnData : tableColumnsData )
  {
    if( tableColumnData.subColumn.empty())
    {
      tableColumnData.columnValues = valuesByColumn[currentColumn++];
    }
    else
    {
      populateSubColumnsFromTableData( tableColumnData, valuesByColumn, currentColumn );
    }
  }
}


void TableTextFormatter::splitAndMergeColumnHeaders( std::vector< TableLayout::ColumnStructure > & tableColumnsData,
                                                     std::vector< std::vector< string > > & allDividedHeaderParts ) const
{
  for( auto & tableColumnData : tableColumnsData )
  {
    std::vector< string > dividedHeaderParts;
    std::istringstream ss( tableColumnData.column.columnName );
    string subColumnNames;
    while( getline( ss, subColumnNames, '\n' ))
    {
      dividedHeaderParts.push_back( subColumnNames );
    }

    allDividedHeaderParts.push_back( dividedHeaderParts );

    if( !tableColumnData.subColumn.empty())
    {
      std::vector< std::vector< string > > dividedSubColHeaders;
      splitAndMergeColumnHeaders( tableColumnData.subColumn, dividedSubColHeaders );
    }
  }

  size_t nbHeaderRows = std::max_element( allDividedHeaderParts.begin(), allDividedHeaderParts.end(),
                                          []( auto const & v1, auto const & v2 )
  {
    return v1.size() < v2.size();
  } )->size();

  for( auto & headerDivided : allDividedHeaderParts )
  {
    if( headerDivided.size() < nbHeaderRows )
    {
      headerDivided.resize( nbHeaderRows, " " );
    }
    tableColumnsData[&headerDivided - &allDividedHeaderParts[0]].column.splitColumnNames = headerDivided;
  }

}

void TableTextFormatter::findAndSetLongestColumnString( TableLayout::ColumnStructure & tableColumnData ) const
{
  size_t & maxStringSizeColumn = tableColumnData.maxStringSize;

  { // header case
    auto const maxStringSizeHeader = *std::max_element( tableColumnData.column.splitColumnNames.begin(),
                                                        tableColumnData.column.splitColumnNames.end(),
                                                        []( const auto & a, const auto & b )
    {
      return a.size() < b.size();
    } );
    maxStringSizeColumn = maxStringSizeHeader.length();
  }

  {  // cells case
    if( tableColumnData.subColumn.empty() && !tableColumnData.columnValues.empty())
    {
      auto const maxStringSizeCell = *std::max_element( tableColumnData.columnValues.begin(),
                                                        tableColumnData.columnValues.end(),
                                                        []( const auto & a, const auto & b )
      {
        return a.size() < b.size();
      } );
      if( maxStringSizeColumn < maxStringSizeCell.length())
      {
        maxStringSizeColumn = maxStringSizeCell.length();
      }
    }
  }

  { // subcolumn values case
    size_t totalLengthSubColumns = 0;
    if( !tableColumnData.subColumn.empty() )
    {
      for( auto & subColumn : tableColumnData.subColumn )
      {
        findAndSetLongestColumnString( subColumn );
        totalLengthSubColumns += subColumn.maxStringSize;
        if( &subColumn != &tableColumnData.subColumn[0] )
        {
          totalLengthSubColumns += m_tableLayout.getColumnMargin();
        }
      }
    }

    if( totalLengthSubColumns > tableColumnData.maxStringSize )
    {
      tableColumnData.maxStringSize = totalLengthSubColumns;
    }
  }

}

void TableTextFormatter::computeAndBuildTableSeparator( std::vector< TableLayout::ColumnStructure > & tableColumnsData,
                                                        string & sectionSeparatingLine,
                                                        string & topSeparator ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );

  size_t const numColumns = tableColumnsData.size() - 1;
  size_t const spacingBetweenColumns = numColumns * (size_t) columnMargin;
  size_t const margins = (size_t) borderMargin * 2;
  size_t sectionlineLength = spacingBetweenColumns + margins;
  for( auto const & tableColumnData : tableColumnsData )
  {
    sectionlineLength += tableColumnData.maxStringSize;
  }

  size_t maxTopLineLength =  tableTitle.length() + margins;
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );
  if( sectionlineLength < maxTopLineLength )
  {
    size_t const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( tableColumnsData, extraCharacters );
    sectionlineLength = maxTopLineLength;
  }

  sectionSeparatingLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
  integer const topSeparatorLength = maxTopLineLength - 2; // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::increaseColumnsSize( std::vector< TableLayout::ColumnStructure > & tableColumnsData,
                                              size_t const extraCharacters ) const
{
  size_t const extraCharactersPerColumn = std::floor( (extraCharacters) / tableColumnsData.size() );
  size_t const overflowCharacters = extraCharacters - (extraCharactersPerColumn *  tableColumnsData.size() );
  for( auto & tableColumnData : tableColumnsData )
  {
    if( !tableColumnData.subColumn.empty())
    {
      increaseColumnsSize( tableColumnData.subColumn, extraCharactersPerColumn );
    }

    size_t const cellSize = tableColumnData.maxStringSize;
    integer const newMaxStringSize = &tableColumnData == &tableColumnsData[0] ?
                                     extraCharactersPerColumn + cellSize + overflowCharacters :
                                     extraCharactersPerColumn + cellSize;
    tableColumnData.maxStringSize = newMaxStringSize;
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

void TableTextFormatter::outputHeaderSectionRows( std::vector< TableLayout::ColumnStructure > const & tableColumnsData,
                                                  std::ostringstream & tableOutput,
                                                  string_view sectionSeparatingLine ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const spaces =  string( columnMargin, ' ' );
  bool containSubColumn = false;
  size_t nbRows = tableColumnsData[0].column.splitColumnNames.size();
  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, borderMargin );

    for( std::size_t idxColumn = 0; idxColumn < tableColumnsData.size(); ++idxColumn )
    {
      auto const & tableColumnData = tableColumnsData[idxColumn];
      string const cell = tableColumnData.column.splitColumnNames.at( idxRow );
      size_t const cellSize = tableColumnData.maxStringSize;
      tableOutput << buildCell( tableColumnData.column.alignmentSettings.headerAlignment,
                                cell,
                                cellSize );

      // Add space between tableColumnsData
      if( idxColumn < tableColumnsData.size() - 1 )
      {
        tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
      }

      if( !tableColumnData.subColumn.empty())
      {
        containSubColumn = true;
      }
    }
    // Append right border with line return
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, borderMargin );
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }

  // Check and build subrow header
  if( containSubColumn )
  {
    std::vector< TableLayout::ColumnStructure > rowSubColumns;

    for( auto const & tableColumnData : tableColumnsData )
    {
      if( tableColumnData.subColumn.empty())
      {
        rowSubColumns.push_back( {TableLayout::Column{""}, {}, tableColumnData.maxStringSize, {}} );
      }
      else
      {
        rowSubColumns.insert( rowSubColumns.end(),
                              tableColumnData.subColumn.begin(),
                              tableColumnData.subColumn.end());
      }
    }
    outputHeaderSectionRows( rowSubColumns, tableOutput, sectionSeparatingLine );

  }
}

void TableTextFormatter::outputSubSection( std::vector< TableLayout::ColumnStructure > const & tableColumnsData,
                                           std::ostringstream & tableOutput,
                                           integer idxRow ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  for( auto const & tableColumnData : tableColumnsData )
  {
    string const cell = tableColumnData.columnValues[idxRow];
    size_t const cellSize =  tableColumnData.maxStringSize;
    tableOutput << buildCell( tableColumnData.column.alignmentSettings.valueAlignment,
                              cell,
                              cellSize );
    if( &tableColumnData < &tableColumnsData.back() )
    {
      tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
    }
  }
}

void TableTextFormatter::outputValuesSectionRows( std::vector< TableLayout::ColumnStructure > const & tableColumnsData,
                                                  std::ostringstream & tableOutput,
                                                  string_view sectionSeparatingLine ) const
{
  size_t const nbRows = tableColumnsData[0].columnValues.size();
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const spaces =  string( columnMargin - 1, ' ' );

  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, borderMargin );

    for( std::size_t idxColumn = 0; idxColumn < tableColumnsData.size(); ++idxColumn )
    {
      auto const & tableColumnData = tableColumnsData[idxColumn];

      if( !tableColumnData.subColumn.empty())
      {
        outputSubSection( tableColumnData.subColumn, tableOutput, idxRow );
      }
      else
      {
        string const cell = tableColumnData.columnValues.at( idxRow );
        size_t const cellSize =  tableColumnData.maxStringSize;
        tableOutput << buildCell( tableColumnData.column.alignmentSettings.valueAlignment,
                                  cell,
                                  cellSize );
      }

      if( idxColumn < tableColumnsData.size() - 1 )
      {
        tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
      }

    }

    // Append right border
    tableOutput << GEOS_FMT( "{:>{}}", m_verticalLine, borderMargin );


    tableOutput << "\n";

  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }
}

}
