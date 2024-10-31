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

/**
 * @brief Given spaces number, distribute it over all string vector
 * @param vec The vector where we adding spaces
 * @param totalSpaces Total spaces to distribute over all strings
 */
void distributeSpaces( std::vector< string > & vec, int totalSpaces )
{
  int numElements = vec.size();
  auto dv = std::div( totalSpaces, numElements );
  int baseSpaces = dv.quot;
  int extraSpaces = dv.rem;

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
    GEOS_ERROR_IF( dest.size() != source[idxRow].size(), "Dest matrix must have the number of rows equal to the number of columns in the source matrix" );
    GEOS_ERROR_IF( dest[idxRow].size() != source.size(), "Dest matrix must have the number of columns equal to the number of rows in the source matrix." );
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

  std::vector< std::vector< string > > splitHeaders;
  splitAndMergeColumnHeaders( tableColumnsData, splitHeaders );

  for( auto & tableColumnData : tableColumnsData )
  {
    findAndSetLongestColumnString( tableColumnData, tableColumnData.maxStringSize, 0 );
  }

  computeTableWidth( tableColumnsData );
  buildTableSeparators( tableColumnsData, sectionSeparatingLine, topSeparator );
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
  for( auto & subColumnData : tableColumnData.subColumn )
  {
    subColumnData.columnValues = valuesByColumn[idxSubColumn++];
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
                                                     std::vector< std::vector< string > > & splitHeaders ) const
{
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
      splitAndMergeColumnHeaders( tableColumnData.subColumn, splitSubColHeaders );
    }
  }

  size_t nbHeaderRows = std::max_element( splitHeaders.begin(), splitHeaders.end(),
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

void TableTextFormatter::findAndSetLongestColumnString( TableLayout::ColumnStructure & tableColumnData,
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

void TableTextFormatter::computeTableWidth( std::vector< TableLayout::ColumnStructure > & tableColumnsData ) const
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

void TableTextFormatter::increaseColumnsSize( std::vector< TableLayout::ColumnStructure > & tableColumnsData,
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

void TableTextFormatter::buildTableSeparators( std::vector< TableLayout::ColumnStructure > const & tableColumnsData,
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

}
