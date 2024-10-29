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

/**
 * @brief distributes spaces equally from a string vector
 * @param vec The vector where we adding spaces
 * @param totalSpaces Total spaces to be distributed
 */
void distributeSpaces( std::vector< string > & vec, int totalSpaces )
{
  int numElements = vec.size();
  int baseSpaces = totalSpaces / numElements;
  int extraSpaces = totalSpaces % numElements;

  for( int i = 0; i < numElements; ++i )
  {
    vec[i] += string( baseSpaces, ' ' );

    if( i < extraSpaces )
    {
      vec[i] += ' ';
    }
  }
}

void transpose( std::vector< std::vector< string > > & dest,
                std::vector< std::vector< string > > const & source )
{
  for( size_t idxRow = 0; idxRow < source.size(); ++idxRow )
  {
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
void updateVisibleColumns( std::vector< TableLayout::TableColumnData > & tableColumnsData,
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
  std::vector< TableLayout::TableColumnData > tableColumnsData( m_tableLayout.getColumns());
  size_t nbHeaderRows = 0;
  TableData tableData;
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( tableColumnsData, tableData, nbHeaderRows, sectionSeparatingLine, topSeparator );

  tableOutput << '\n';
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  outputHeaderSectionRows( tableColumnsData, tableOutput, nbHeaderRows, sectionSeparatingLine );

  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::vector< TableLayout::TableColumnData > tableColumnsData( m_tableLayout.getColumns());
  size_t nbHeaderRows = 0;
  string sectionSeparatingLine;
  string topSeparator;

  prepareAndBuildTable( tableColumnsData, tableData, nbHeaderRows, sectionSeparatingLine, topSeparator );

  std::ostringstream tableOutput;
  outputTable( tableOutput, tableColumnsData, nbHeaderRows, sectionSeparatingLine, topSeparator );

  return tableOutput.str();
}

void TableTextFormatter::prepareAndBuildTable( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                                               TableData const & tableData,
                                               size_t & nbHeaderRows,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{
  std::vector< std::vector< string > > tableDataRows( tableData.getTableDataRows());
  if( !tableDataRows.empty())
  {
    updateVisibleColumns( tableColumnsData, tableDataRows );
    populateColumnsFromTableData( tableColumnsData, tableDataRows, false );
  }

  std::vector< std::vector< string > > splitHeaders;
  splitAndMergeColumnHeaders( tableColumnsData, nbHeaderRows, splitHeaders );

  for( auto & tableColumnData : tableColumnsData )
  {
    findAndSetLongestColumnString( tableColumnData, tableColumnData.maxStringSize, 0 );
  }

  computeTableWidth( tableColumnsData );
  buildTableSeparators( tableColumnsData, sectionSeparatingLine, topSeparator );
}

void TableTextFormatter::outputTable( std::ostringstream & tableOutput,
                                      std::vector< TableLayout::TableColumnData > & tableColumnsData,
                                      size_t & nbHeaderRows,
                                      string_view sectionSeparatingLine,
                                      string_view topSeparator ) const
{
  integer const nbValuesRows = tableColumnsData[0].columnValues.size();

  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputHeaderSectionRows( tableColumnsData, tableOutput, nbHeaderRows, sectionSeparatingLine );

  outputValuesSectionRows( tableColumnsData, tableOutput, nbValuesRows, sectionSeparatingLine );
  if( m_tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

void TableTextFormatter::populateColumnsFromTableData( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                                                       std::vector< std::vector< string > > const & tableDataRows,
                                                       bool isSubColumn ) const
{
  size_t currentColumn = 0;
  std::vector< std::vector< string > > tColumnsValues( tableDataRows[0].size(),
                                                       std::vector< string >( tableDataRows.size()));
  if( !isSubColumn )
  {
    transpose( tColumnsValues, tableDataRows );
  }

  for( auto & tableColumnData : tableColumnsData )
  {
    if( tableColumnData.subColumn.empty())
    {
      tableColumnData.columnValues = !isSubColumn ?
                                     tColumnsValues[currentColumn++] : tableDataRows[currentColumn++];
    }
    else
    {
      size_t nbSubColumns = tableColumnData.column.subColumnNames.size();
      auto subColumnStartIdx = tColumnsValues.begin() + currentColumn;
      std::vector< std::vector< string > > subColumnValues( subColumnStartIdx,
                                                            subColumnStartIdx + nbSubColumns );
      populateColumnsFromTableData( tableColumnData.subColumn, subColumnValues, true );

      std::vector< std::vector< string > > tSubColumnValues( subColumnValues[0].size(),
                                                             std::vector< string >( subColumnValues.size()) );
      transpose( tSubColumnValues, subColumnValues );
      for( const auto & columnValues : tSubColumnValues )
      { // add all subcolumn values in parent tableColumnData
        tableColumnData.columnValues.insert( tableColumnData.columnValues.end(),
                                             columnValues.begin(),
                                             columnValues.end() );
      }
      currentColumn += subColumnValues.size();
    }
  }
}

void TableTextFormatter::splitAndMergeColumnHeaders( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                                                     size_t & nbHeaderRows,
                                                     std::vector< std::vector< string > > & splitHeaders ) const
{
  splitHeaders.reserve( tableColumnsData.size() );

  for( auto & tableColumnData : tableColumnsData )
  {
    std::vector< string > splitHeaderParts;
    std::istringstream ss( tableColumnData.column.columnName );
    string subColumnNames;
    while( getline( ss, subColumnNames, '\n' ))
    {
      splitHeaderParts.push_back( subColumnNames );
    }

    splitHeaders.push_back( splitHeaderParts );

    if( !tableColumnData.subColumn.empty())
    {
      std::vector< std::vector< string > > splitSubColHeaders;
      size_t nbHeaderSubColRows = 0;
      splitAndMergeColumnHeaders( tableColumnData.subColumn, nbHeaderSubColRows, splitSubColHeaders );
    }
  }

  nbHeaderRows = std::max_element( splitHeaders.begin(), splitHeaders.end(),
                                   []( auto const & v1, auto const & v2 )
  {
    return v1.size() < v2.size();
  } )->size();

  for( auto & headerParts : splitHeaders )
  {
    if( headerParts.size() < nbHeaderRows )
    {
      headerParts.resize( nbHeaderRows, " " );
    }
    tableColumnsData[&headerParts - &splitHeaders[0]].column.splitColumnNames = headerParts;
  }

}

void TableTextFormatter::findAndSetLongestColumnString( TableLayout::TableColumnData & tableColumnData,
                                                        std::vector< string > & maxStringSize,
                                                        integer const idxMaxString ) const
{
  string maxStringColumn;
  { // header case
    auto const maxStringSizeHeader = *std::max_element( tableColumnData.column.splitColumnNames.begin(),
                                                        tableColumnData.column.splitColumnNames.end(),
                                                        []( const auto & a, const auto & b )
    {
      return a.size() < b.size();
    } );

    maxStringColumn = maxStringSizeHeader;
    maxStringSize.push_back( maxStringSizeHeader );
  }

  {  // values case
    if( tableColumnData.subColumn.empty() && !tableColumnData.columnValues.empty())
    {
      auto const maxStringSizeCell = *std::max_element( tableColumnData.columnValues.begin(),
                                                        tableColumnData.columnValues.end(),
                                                        []( const auto & a, const auto & b )
      {
        return a.size() < b.size();
      } );

      if( maxStringColumn.length() < maxStringSizeCell.length())
      {
        maxStringColumn = maxStringSizeCell;
      }
    }
  }

  { // Update max string size if necessary
    if( maxStringSize[idxMaxString].length() < maxStringColumn.length() )
    {
      maxStringSize[idxMaxString] = maxStringColumn;
    }
  }

  { // subcolumn values case
    if( !tableColumnData.subColumn.empty() )
    {
      tableColumnData.maxStringSize.clear();
      for( size_t idxSubColumn = 0; idxSubColumn < tableColumnData.subColumn.size(); ++idxSubColumn )
      {
        findAndSetLongestColumnString( tableColumnData.subColumn[idxSubColumn],
                                       tableColumnData.maxStringSize,
                                       idxSubColumn );
      }
    }
  }

  if( tableColumnData.maxStringSize.empty() )
  {
    tableColumnData.maxStringSize.push_back( maxStringColumn );
  }
}

void TableTextFormatter::computeTableWidth( std::vector< TableLayout::TableColumnData > & tableColumnsData ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );

  string::size_type numColumns = tableColumnsData.size() - 1;
  string::size_type spacing = numColumns * columnMargin;
  string::size_type margins = borderMargin * 2;

  string::size_type const sectionLengthWithSpacing = spacing + margins;

  string::size_type sectionlineLength = sectionLengthWithSpacing;
  string const spaces =  string( m_tableLayout.getColumnMargin(), ' ' );

  { // Compute total length of all tableColumnsData with margins
    sectionlineLength += std::accumulate( tableColumnsData.begin(), tableColumnsData.end(), 0,
                                          [&]( auto sum, auto & tableColumnData ) -> auto
    { // take into account subColumn
      string sumOfString = stringutilities::join( tableColumnData.maxStringSize, spaces );
      return static_cast< decltype(sum) >(sum + sumOfString.length());
    } );
  }

  string::size_type maxTopLineLength =  tableTitle.length() + margins;
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );
  if( sectionlineLength < maxTopLineLength )
  {
    real64 const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( tableColumnsData, extraCharacters );
  }
}

void TableTextFormatter::increaseColumnsSize( std::vector< TableLayout::TableColumnData > & tableColumnsData,
                                              real64 const extraCharacters ) const
{
  real64 const extraCharactersPerColumn = std::floor( (extraCharacters) / tableColumnsData.size() );
  integer overflowCharacters = extraCharacters - (extraCharactersPerColumn *  tableColumnsData.size() );
  for( std::size_t idxColumn = 0; idxColumn < tableColumnsData.size(); ++idxColumn )
  {
    if( !tableColumnsData[idxColumn].subColumn.empty())
    {
      distributeSpaces( tableColumnsData[idxColumn].maxStringSize, (int)extraCharactersPerColumn );
      increaseColumnsSize( tableColumnsData[idxColumn].subColumn, extraCharactersPerColumn );
    }
    else
    {
      string & cell = tableColumnsData[idxColumn].maxStringSize[0];
      integer newMaxStringSize = idxColumn == 0 ?
                                 extraCharactersPerColumn + cell.size() + overflowCharacters :
                                 extraCharactersPerColumn + cell.size();
      cell = GEOS_FMT( "{:>{}}", cell, newMaxStringSize );
    }
  }
}

void TableTextFormatter::buildTableSeparators( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{
  std::vector< string > maxStringsPerColumn;
  for( auto const & tableColumnData : tableColumnsData )
  {
    std::for_each( tableColumnData.maxStringSize.begin(), tableColumnData.maxStringSize.end(),
                   [&] ( string maxString ) {
      maxStringsPerColumn.push_back( string( maxString.length(), m_horizontalLine ) );
    } );
  }

  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const patternBetweenColumns = GEOS_FMT( "{:-^{}}", m_horizontalLine, columnMargin );
  string const leftBorder = GEOS_FMT( "{:-<{}}", m_horizontalLine, borderMargin );
  string const rightBorder = GEOS_FMT( "{:-<{}}", m_horizontalLine, borderMargin );

  string const columnJoin = stringutilities::join( maxStringsPerColumn, patternBetweenColumns );
  std::ostringstream oss;
  oss << leftBorder << columnJoin << rightBorder;
  sectionSeparatingLine = oss.str();

  integer const topSeparatorLength = sectionSeparatingLine.size() - 2; // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::outputSubSection( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
                                           std::ostringstream & tableOutput,
                                           integer idxRow ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  for( size_t idxCol = 0; idxCol< tableColumnsData.size(); ++idxCol )
  {
    tableOutput << buildCell( tableColumnsData[idxCol].column.alignmentSettings.valueAlignment,
                              tableColumnsData[idxCol].columnValues[idxRow],
                              tableColumnsData[idxCol].maxStringSize[0].length() );
    if( idxCol < tableColumnsData.size() - 1 )
    {
      tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
    }
  }
}

void TableTextFormatter::outputValuesSectionRows( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
                                                  std::ostringstream & tableOutput,
                                                  size_t const nbRows,
                                                  string_view sectionSeparatingLine ) const
{
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
        string const cellSize = stringutilities::join( tableColumnData.maxStringSize, spaces );
        tableOutput << buildCell( tableColumnData.column.alignmentSettings.valueAlignment,
                                  cell,
                                  cellSize.length());
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

void TableTextFormatter::outputHeaderSectionRows( std::vector< TableLayout::TableColumnData > const & tableColumnsData,
                                                  std::ostringstream & tableOutput,
                                                  size_t const nbRows,
                                                  string_view sectionSeparatingLine ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const spaces =  string( columnMargin, ' ' );
  bool containSubColumn = false;

  for( size_t idxRow = 0; idxRow < nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, borderMargin );

    for( std::size_t idxColumn = 0; idxColumn < tableColumnsData.size(); ++idxColumn )
    {
      auto const & tableColumnData = tableColumnsData[idxColumn];
      string cell = tableColumnData.column.splitColumnNames.at( idxRow );
      string cellSize =  stringutilities::join( tableColumnData.maxStringSize, spaces );
      tableOutput << buildCell( tableColumnData.column.alignmentSettings.headerAlignment,
                                cell,
                                cellSize.length());

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
    std::vector< TableLayout::TableColumnData > rowSubColumns;

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
    outputHeaderSectionRows( rowSubColumns, tableOutput, 1, sectionSeparatingLine );
  }
}

}
