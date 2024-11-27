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
#include "TableTypes.hpp"

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
  std::ostringstream tableOutput;//todo
  // TableData tableData;
  // RowsCellLayout cellsHeaderLayout;
  // RowsCellLayout cellsDataLayout;
  // string sectionSeparatingLine;
  // string topSeparator;

  // initalizeTableLayout( tableData, cellsHeaderLayout, cellsDataLayout,
  //                       sectionSeparatingLine, topSeparator );
  // tableOutput << '\n';
  // outputTitleRow( tableOutput, topSeparator );
  // tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  // outputTable( tableOutput, cellsHeaderLayout, cellsDataLayout,
  //              sectionSeparatingLine, topSeparator );
  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  RowsCellLayout cellsHeaderLayout;
  RowsCellLayout cellsDataLayout;
  string sectionSeparatingLine;
  string topSeparator;

  TableLayout tableLayout = m_tableLayout;
  initalizeTableLayout( tableLayout, tableData, cellsHeaderLayout, cellsDataLayout,
                        sectionSeparatingLine, topSeparator );
  outputTable( tableLayout, tableOutput, cellsHeaderLayout, cellsDataLayout,
               sectionSeparatingLine, topSeparator );
  return tableOutput.str();
}

void TableTextFormatter::initalizeTableLayout( TableLayout & tableLayout,
                                               TableData const & tableData,
                                               RowsCellLayout & cellsDataLayout,
                                               RowsCellLayout & cellsHeaderLayout,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{
  RowsCellInput inputDataValues( tableData.getTableDataRows());
  if( !cellsDataLayout.empty())
  {
    updateVisibleColumns( tableLayout, inputDataValues );
    populateDataCellsLayout( tableLayout, cellsDataLayout, inputDataValues );
  }

  findAndSetLongestColumnString( tableLayout, cellsDataLayout );

  computeAndBuildTableSeparator( tableLayout, sectionSeparatingLine, topSeparator );

  cellsHeaderLayout = { tableLayout.getColumns().size(),
                        std::vector< TableLayout::CellLayout >( cellsDataLayout.size() )
  };

  gridifyHeaders( tableLayout, cellsHeaderLayout );

}

void TableTextFormatter::outputTable( TableLayout & tableLayout,
                                      std::ostringstream & tableOutput,
                                      RowsCellLayout const & cellsHeader,
                                      RowsCellLayout const & cellsData,
                                      string_view sectionSeparatingLine,
                                      string_view topSeparator ) const
{
  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableLayout, tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputLines( tableLayout, cellsHeader, tableOutput, sectionSeparatingLine );
  outputLines( tableLayout, cellsData, tableOutput, sectionSeparatingLine );

  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

/**
 * @brief Detect columns who are not displayed from TableLayout and therefore modify columns / inputDataValues vectors
 * @param inputDataValues Vector built in TableData containing all rows values
 */
void TableTextFormatter::updateVisibleColumns( TableLayout & tableLayout,
                                               RowsCellInput & inputDataValues ) const
{
  integer idxColumn = 0;
  auto columns =  tableLayout.getColumns();
  for( auto iterColumn = columns.begin(); iterColumn !=  columns.end(); )
  {
    if( !iterColumn->enabled )
    {
      iterColumn = columns.erase( iterColumn );
      for( auto & row : inputDataValues )
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

void TableTextFormatter::gridifyHeaders( TableLayout & tableLayout,
                                         RowsCellLayout & cellsDataLayout ) const
{
  size_t nbLayers = tableLayout.getMaxHeaderRow();

  for( auto it = tableLayout.beginLeaf(); it !=  tableLayout.endLeaf(); ++it )
  {
    size_t currLayer = it.getCurrentLayer();
    if( !it->hasChild()  )
    {
      if( nbLayers != currLayer )
      {
        for( size_t i = currLayer; i< nbLayers; i++ )
        {
          TableLayout::CellLayout cell{CellType::Header, "", TableLayout::Alignment::center};
          cellsDataLayout[i].push_back( cell );
        }
      }
    }

    TableLayout::CellLayout currentCell = it->columnName;
    if( it->hasParent() )
    {
      currentCell.cellWidth++;
    }
    cellsDataLayout[currLayer].push_back( currentCell );

    for( size_t idxColumn = 1; idxColumn <= currentCell.cellWidth; idxColumn++ )
    {
      TableLayout::CellLayout mergingCell{ CellType::Merge, "", TableLayout::Alignment::center };
      cellsDataLayout[currLayer].push_back( mergingCell );
    }
  }
}

void TableTextFormatter::populateDataCellsLayout( TableLayout & tableLayout,
                                                  RowsCellLayout & cellsDataLayout,
                                                  RowsCellInput & inputDataValues ) const
{
  // let's reserve the layout cells buffer
  cellsDataLayout = { std::vector< TableLayout::CellLayout >( inputDataValues.size() )};

  // to insert directly the values in each columns, we fill with the transposed rowsCells (row major->column major)
  for( size_t idxRow = 0; idxRow < inputDataValues.size(); ++idxRow )
  {
    size_t maxLineCount = 0;

    for( size_t idxCol = 0; idxCol < inputDataValues[idxRow].size(); ++idxCol )
    {
      TableData::CellData const & cell = inputDataValues[idxRow][idxCol];
      TableLayout::CellAlignment const cellAlignement = tableLayout.getColumns()[idxCol].cellAlignment;
      TableLayout::Alignment const alignement = cell.type == CellType::Header ?
                                                cellAlignement.headerAlignment :
                                                cellAlignement.valueAlignment;

      cellsDataLayout[idxRow][idxCol] = TableLayout::CellLayout( cell.type, cell.value, alignement );
      maxLineCount = std::max( maxLineCount, cellsDataLayout[idxRow][idxCol].lines.size() );
    }
    // tableLayout.m_valueRows.push_back( { maxLineCount } ); //todo
  }
}

void TableTextFormatter::findAndSetLongestColumnString( TableLayout & tableLayout,
                                                        RowsCellLayout & cellsDataLayout ) const
{
  auto getMaxStringLen = []( std::vector< std::string > const & lines ) -> size_t
  {
    auto maxLineLength = std::max_element( lines.begin(), lines.end(),
                                           []( const std::string & a, const std::string & b ) {
      return a.length() < b.length();
    } );

    return (maxLineLength != lines.end()) ? maxLineLength->length() : 0;
  };

  for( auto it = tableLayout.beginLeaf(); it != tableLayout.endLeaf(); ++it )
  {
    size_t maxSize = 1;
    TableLayout::Column * currentColumn = &(*it);
    if( !it->hasChild() )
    {
      //1. currentColumn = max(subcolumn, dataColumn)
      size_t idxLine = it.getCurrentLayer();
      for( size_t idxColumnData = 0; idxColumnData< cellsDataLayout[idxLine].size(); idxColumnData++ )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[idxColumnData][idxLine];
        cell.maxDataLength = std::max( cell.maxDataLength, getMaxStringLen( cell.lines ));
        maxSize = std::max( cell.maxDataLength, maxSize );
      }

      currentColumn->setMaxStringSize( std::max( getMaxStringLen( currentColumn->columnName.lines ),
                                                 maxSize ));
    }
    else     //2. and max(parent, subcolumn)
    {
      if( !it->hasParent())
      {
        maxSize = std::max( getMaxStringLen( currentColumn->columnName.lines ),
                            maxSize );
        currentColumn->columnName.maxDataLength = maxSize;
        currentColumn->setMaxStringSize( maxSize );
      }
      else
      {
        std::vector< size_t > subColumnsLength;
        for( auto const & subColumn : currentColumn->subColumn )
        {
          subColumnsLength.push_back( subColumn.getMaxStringSize() );
        }
        currentColumn->compareAndSetMaxStringSize( currentColumn, subColumnsLength );
        currentColumn->columnName.maxDataLength = currentColumn->getMaxStringSize();
      }
    }
  }
}

void TableTextFormatter::computeAndBuildTableSeparator( TableLayout & tableLayout,
                                                        string & sectionSeparatingLine,
                                                        string & topSeparator ) const
{
  integer const columnMargin = tableLayout.getColumnMargin();
  integer const borderMargin = tableLayout.getBorderMargin();
  string const tableTitle = string( tableLayout.getTitle() );
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
  integer const topSeparatorLength = maxTopLineLength - 2;     // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::increaseColumnsSize( TableLayout & tableLayout,
                                              size_t const extraCharacters ) const
{
  size_t const columnsCount = tableLayout.getColumns().size();
  size_t const extraCharactersPerColumn = std::floor( extraCharacters / columnsCount );
  size_t const overflowCharacters = extraCharacters - (extraCharactersPerColumn * columnsCount);

  for( auto it = tableLayout.beginRoot(); it != tableLayout.endRoot(); ++it )
  {
    TableLayout::Column * currentColumn = &(*it);

    size_t divider = currentColumn->hasParent() ? currentColumn->m_parent->subColumn.size() : columnsCount;
    size_t extraCharactersForCurrentColumn = std::floor( extraCharacters / divider );
    size_t overflowForCurrentColumn = (currentColumn->m_next == nullptr) ? overflowCharacters : 0;

    integer const newMaxStringSize = currentColumn->getMaxStringSize() +
                                     extraCharactersForCurrentColumn +
                                     overflowForCurrentColumn;
    currentColumn->setMaxStringSize( newMaxStringSize );
  }
}

void TableTextFormatter::outputTitleRow( TableLayout & tableLayout,
                                         std::ostringstream & tableOutput,
                                         string_view topSeparator ) const
{
  string const tableTitle = string( tableLayout.getTitle());
  if( !tableTitle.empty() )
  {
    tableOutput << GEOS_FMT( "{}\n", topSeparator );
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin());
    tableOutput << buildCell( TableLayout::Alignment::center,
                              tableTitle,
                              (topSeparator.length() - (tableLayout.getBorderMargin() *  2)));
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, tableLayout.getBorderMargin() );
  }
}

void TableTextFormatter::formatCell( TableLayout & tableLayout,
                                     std::ostringstream & tableOutput,
                                     TableLayout::CellLayout const & cell ) const
{
  for( auto const & line : cell.lines )
  {
    tableOutput << string( tableLayout.getBorderMargin() - 1, ' ' );
    tableOutput << buildCell( cell.alignment,
                              line,
                              cell.maxDataLength );
    tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, tableLayout.getColumnMargin());
    tableOutput << m_verticalLine << "\n";
  }
}

void TableTextFormatter::outputLines( TableLayout & tableLayout,
                                      RowsCellLayout const & cellsLayout,
                                      std::ostringstream & tableOutput,
                                      string_view sectionSeparatingLine ) const
{
  for( auto const & line : cellsLayout )
  {
    for( auto const & cell : line )
    {
      formatCell( tableLayout, tableOutput, cell );
    }
  }
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
}
}
