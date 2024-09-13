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
    oss << m_tableLayout.getColumns()[idxColumn].m_parameter.columnName;
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
  std::vector< std::vector< string > > const rowsValues = tableData.getTableDataRows();
  std::ostringstream oss;

  for( const auto & row : rowsValues )
  {
    oss << stringutilities::join( row.cbegin(), row.cend(), "," ) << "\n";
  }

  return oss.str();
}

string TableCSVFormatter::toString( TableData const & tableData ) const
{
  return headerToString() + dataToString( tableData );
}

///////////////////////////////////////////////////////////////////////
////// Log Formatter implementation
///////////////////////////////////////////////////////////////////////

/**
 * @brief Build  cell given an alignment, a value and spaces
 * @param alignment The aligment of cell value
 * @param value The cell value
 * @param spaces The number of spaces in the cell
 * @return A formated cell
 */
string buildCell( TableLayout::Alignment const alignment, string_view value, integer const spaces )
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
void formatColumnsFromLayout( std::vector< TableLayout::Column > & columns,
                              std::vector< std::vector< string > > & tableDataRows )
{
  integer idxColumn = 0;
  for( auto iterColumn = columns.begin(); iterColumn!=columns.end(); )
  {
    if( !iterColumn->m_parameter.enabled )
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

void TableTextFormatter::fillTableColumnsFromRows( std::vector< TableLayout::Column > & columns,
                                                   std::vector< std::vector< string > > & rows ) const
{
  for( size_t idxRow = 0; idxRow < rows.size(); idxRow++ )
  {
    if( rows[idxRow].size() < columns.size() )
    {
      rows[idxRow].resize( columns.size(), " " );
    }

    for( size_t idxColumn = 0; idxColumn < columns.size(); idxColumn++ )
    {
      if( m_tableLayout.getColumns()[idxColumn].m_parameter.enabled )
      {
        columns[idxColumn].m_columnValues.push_back( rows[idxRow][idxColumn] );
      }
    }
  }
}

string TableTextFormatter::layoutToString() const
{
  std::ostringstream tableOutput;
  string sectionSeparatingLine;
  std::vector< TableLayout::Column > columns = m_tableLayout.getColumns();

  outputLayout( tableOutput, columns, {}, sectionSeparatingLine );
  return tableOutput.str();
}

string TableTextFormatter::toString( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  string sectionSeparatingLine;
  std::vector< TableLayout::Column > columns         = m_tableLayout.getColumns();
  std::vector< std::vector< string > > tableDataRows = tableData.getTableDataRows();
  std::vector< string > const & msgTableError        = tableData.getErrorMsgs();
  integer const nbValuesRows                         = tableDataRows.size();

  formatColumnsFromLayout( columns, tableDataRows );
  fillTableColumnsFromRows( columns, tableDataRows );

  outputLayout( tableOutput, columns, msgTableError, sectionSeparatingLine );
  outputSectionRows( columns, sectionSeparatingLine, tableOutput, nbValuesRows, TableLayout::Section::values );
  tableOutput << '\n';

  return tableOutput.str();
}

void TableTextFormatter::outputLayout( std::ostringstream & tableOutput,
                                       std::vector< TableLayout::Column > & columns,
                                       std::vector< string > const & msgTableError,
                                       string & sectionSeparatingLine ) const
{
  string topSeparator;
  size_t nbHeaderRows = 0;
  std::vector< std::vector< string > > splitHeaders;
  string const tableTitle = string( m_tableLayout.getTitle());

  splitAndSetColumnNames( columns, nbHeaderRows, splitHeaders );
  findAndSetMaxStringSize( columns );

  computeTableWidth( columns, msgTableError );
  buildTableSeparators( columns, topSeparator, sectionSeparatingLine );

  tableOutput << '\n';
  outputTopRows( tableOutput, {tableTitle}, topSeparator, TableLayout::Alignment::center );
  outputTopRows( tableOutput, msgTableError, topSeparator, TableLayout::Alignment::left );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputSectionRows( columns, sectionSeparatingLine, tableOutput, nbHeaderRows, TableLayout::Section::header );
}

void TableTextFormatter::splitAndSetColumnNames( std::vector< TableLayout::Column > & columns,
                                                 size_t & nbHeaderRows,
                                                 std::vector< std::vector< string > > & splitHeaders ) const
{

  splitHeaders.reserve( columns.size() );
  for( auto const & column : columns )
  {
    std::vector< string > splitHeaderParts;
    std::istringstream ss( column.m_parameter.columnName );
    string subColumnNames;

    while( getline( ss, subColumnNames, '\n' ))
    {
      splitHeaderParts.push_back( subColumnNames );
    }

    splitHeaders.push_back( std::move( splitHeaderParts ) );
  }

  nbHeaderRows = std::max_element( splitHeaders.begin(), splitHeaders.end(),
                                   []( auto const & v1, auto const & v2 ) { return v1.size() < v2.size(); } )->size();

  for( auto & headerParts : splitHeaders )
  {
    if( headerParts.size() < nbHeaderRows )
    {
      headerParts.resize( nbHeaderRows, " " );
    }
    columns[&headerParts - &splitHeaders[0]].m_parameter.splitColumnNameLines = headerParts;
  }

}

void TableTextFormatter::findAndSetMaxStringSize( std::vector< TableLayout::Column > & columns ) const
{
  for( auto & column : columns )
  {
    auto const maxStringSizeHeader = *std::max_element( column.m_parameter.splitColumnNameLines.begin(),
                                                        column.m_parameter.splitColumnNameLines.end(),
                                                        []( const auto & a, const auto & b ) {return a.size() < b.size();} );
    column.m_maxStringSize = maxStringSizeHeader;

    for( auto const & cell : column.m_columnValues )
    {
      if( column.m_maxStringSize.length() < cell.length())
      {
        column.m_maxStringSize = cell;
      }
    }
  }
}

void TableTextFormatter::increaseColumnsSize( std::vector< TableLayout::Column > & columns,
                                              integer const extraCharacters ) const
{
  integer const extraCharactersPerColumn = std::ceil( extraCharacters / columns.size() );
  for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {
    integer newMaxStringSize = extraCharactersPerColumn + columns[idxColumn].m_maxStringSize.size();
    if( idxColumn == columns.size() - 1 )
    {
      newMaxStringSize += m_tableLayout.getColumnMargin();
    }

    columns[idxColumn].m_maxStringSize = GEOS_FMT( "{:>{}}",
                                                   columns[idxColumn].m_maxStringSize,
                                                   newMaxStringSize );
  }
}

void computeTableErrorLength( std::vector< string > const & msgTableError, string::size_type & msgTableErrorLength )
{
  if( !msgTableError.empty() )
  {
    auto maxStringSize = *(std::max_element( msgTableError.begin(), msgTableError.end(),
                                             []( const auto & a, const auto & b ) {
      return a.size() < b.size();
    } ));

    msgTableErrorLength += maxStringSize.size() + 1; // max string size + line return at the end
  }
}

void computeTableSectionLength( std::vector< TableLayout::Column > & columns, string::size_type & sectionlineLength )
{
  sectionlineLength += std::accumulate( columns.begin(), columns.end(), 0,
                                        []( auto sum, const auto & column )
  { return sum + column.m_maxStringSize.length();} );
}

void TableTextFormatter::computeTableWidth( std::vector< TableLayout::Column > & columns,
                                            std::vector< string > const & msgTableError ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();
  string const tableTitle = string( m_tableLayout.getTitle() );

  string::size_type sectionLengthWithSpacing = ( ( columns.size() - 1 ) *  columnMargin ) + (borderMargin * 2);
  string::size_type sectionlineLength = sectionLengthWithSpacing;
  string::size_type maxTopLineLength =  tableTitle.length();
  string::size_type msgTableErrorLength = borderMargin;

  computeTableErrorLength( msgTableError, msgTableErrorLength );
  computeTableSectionLength( columns, sectionlineLength );

  maxTopLineLength = std::max( maxTopLineLength, msgTableErrorLength );
  if( sectionlineLength < maxTopLineLength )
  {
    integer const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( columns, extraCharacters );
  }
}


void TableTextFormatter::buildTableSeparators( std::vector< TableLayout::Column > const & columns,
                                               string & topSeparator,
                                               string & sectionSeparatingLine ) const
{
  { // section separator line
    integer const columnMargin = m_tableLayout.getColumnMargin();
    integer const borderMargin = m_tableLayout.getBorderMargin();

    std::vector< string > maxStringPerColumn;
    for( auto const & column : columns )
    {
      maxStringPerColumn.push_back( string( column.m_maxStringSize.length(), m_horizontalLine ) );
    }

    string const patternBetweenColumns = GEOS_FMT( "{:-^{}}", m_horizontalLine, columnMargin );

    std::string const leftBorder = GEOS_FMT( "{}{:-<{}}", m_horizontalLine, "", borderMargin );
    std::string const rightBorder = GEOS_FMT( "{}{:-<{}}", m_horizontalLine, "", borderMargin );
    std::string const columnJoin = stringutilities::join( maxStringPerColumn, patternBetweenColumns );

    std::ostringstream oss;
    oss << leftBorder << columnJoin << rightBorder;
    sectionSeparatingLine = oss.str();
  }

  { // top line separator
    // -2 because we can have '+' to the extremity of top separator
    integer const topSeparatorLength = sectionSeparatingLine.size() - 2;
    topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
  }
}

void TableTextFormatter::outputTopRows( std::ostringstream & tableOutput,
                                        std::vector< string > const & msg,
                                        string_view topSeparator,
                                        TableLayout::Alignment alignment ) const
{
  if( msg.size() != 0 && msg[0] != "" )
  {
    tableOutput << GEOS_FMT( "{}\n", topSeparator );
    for( std::string const & str : msg )
    {
      tableOutput << m_verticalLine << string( m_tableLayout.getBorderMargin(), ' ' );
      tableOutput << buildCell( alignment, str, (topSeparator.length() - 6));
      tableOutput << string( m_tableLayout.getBorderMargin(), ' ' ) << "|\n";
    }
  }
}

void TableTextFormatter::outputSectionRows( std::vector< TableLayout::Column > const & columns,
                                            string_view sectionSeparatingLine,
                                            std::ostringstream & tableOutput,
                                            integer const nbRows,
                                            TableLayout::Section const section ) const
{
  integer const columnMargin = m_tableLayout.getColumnMargin();
  integer const borderMargin = m_tableLayout.getBorderMargin();

  for( integer idxRow = 0; idxRow< nbRows; ++idxRow )
  {
    // Append the left border
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, 1 +  borderMargin );
    for( std::size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
    {
      TableLayout::Column const currentColumn = columns[idxColumn];
      auto const & columnContent = section == TableLayout::Section::header ?
                                   columns[idxColumn].m_parameter.splitColumnNameLines :
                                   columns[idxColumn].m_columnValues;
      string cell = columnContent.at( idxRow );
      integer const cellSize = currentColumn.m_maxStringSize.length();

      tableOutput << buildCell( currentColumn.m_parameter.alignment, cell, cellSize );

      // Add space between column
      if( idxColumn < columns.size() - 1 )
      {
        tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, columnMargin );
      }

    }

    // Append right border with line return
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, borderMargin + 1 );
  }

  if( nbRows != 0 )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }
}
}
