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
  //   oss << m_tableLayout.getColumns()[idxColumn].m_columName.value;
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
 * @brief Build cell given an m_alignment, a value and spaces
 * @param m_alignment The aligment of cell value
 * @param value The cell value
 * @param spaces The number of spaces in the cell
 * @return A formated cell
 */
string buildCell( TableLayout::Alignment const m_alignment, string_view value, size_t const spaces )
{
  switch( m_alignment )
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
  CellLayoutRows cellsHeaderLayout;
  CellLayoutRows cellsDataLayout;
  string sectionSeparatingLine;
  string topSeparator;

  TableLayout tableLayout = m_tableLayout;
  initalizeTableLayout( tableLayout, tableData, cellsHeaderLayout, cellsDataLayout,
                        sectionSeparatingLine, topSeparator );
  outputTable( tableLayout, tableOutput, cellsHeaderLayout, cellsDataLayout,
               sectionSeparatingLine, topSeparator );
  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  CellLayoutRows cellsHeaderLayout;
  CellLayoutRows cellsDataLayout;
  string sectionSeparatingLine;
  string topSeparator;

  TableLayout tableLayout = m_tableLayout;
  initalizeTableLayout( tableLayout, tableData, cellsHeaderLayout, cellsDataLayout,
                        sectionSeparatingLine, topSeparator );
  outputTable( tableLayout, tableOutput, cellsHeaderLayout, cellsDataLayout,
               sectionSeparatingLine, topSeparator );
  std::cout <<  tableOutput.str()<< std::endl;
  return tableOutput.str();
}

void TableTextFormatter::initalizeTableLayout( TableLayout & tableLayout,
                                               TableData const & tableData,
                                               CellLayoutRows & cellsHeaderLayout,
                                               CellLayoutRows & cellsDataLayout,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{
  setLinks( tableLayout.getColumns()  );

  RowsCellInput inputDataValues( tableData.getTableDataRows());

  if( inputDataValues.size() > 0 )
  {
    populateDataCellsLayout( tableLayout, cellsDataLayout, inputDataValues );
  }

  updateColumnMaxLength( tableLayout, cellsDataLayout );

  calculateTableSeparators( tableLayout, cellsDataLayout, sectionSeparatingLine, topSeparator );

  cellsHeaderLayout.resize( tableLayout.getMaxHeaderRow() );

  gridifyHeaders( tableLayout, cellsHeaderLayout );
}

void TableTextFormatter::outputTable( TableLayout & tableLayout,
                                      std::ostringstream & tableOutput,
                                      CellLayoutRows const & cellsHeader,
                                      CellLayoutRows const & cellsData,
                                      string_view sectionSeparatingLine,
                                      string_view topSeparator ) const
{
  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableLayout, tableOutput, topSeparator );
  tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );

  outputLines( tableLayout, cellsHeader, tableOutput, tableLayout.getSublineInHeaderCounts(), CellType::Header, sectionSeparatingLine );
  outputLines( tableLayout, cellsData, tableOutput, tableLayout.getNbSubDataLines(), CellType::Value, sectionSeparatingLine );
  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

void TableTextFormatter::setLinks( std::vector< TableLayout::Column > & columns ) const
{
  for( size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {
    if( idxColumn < columns.size() - 1 )
    {
      columns[idxColumn].m_next = &columns[idxColumn + 1];
    }

    if( !columns[idxColumn].m_subColumn.empty())
    {
      for( auto & subCol : columns[idxColumn].m_subColumn )
      {
        subCol.m_parent = &columns[idxColumn];
      }

      setLinks( columns[idxColumn].m_subColumn );
    }
  }
}

void TableTextFormatter::populateDataCellsLayout( TableLayout & tableLayout,
                                                  CellLayoutRows & cellsDataLayout,
                                                  RowsCellInput & inputDataValues ) const
{
  // let's reserve the layout cells buffer
  cellsDataLayout = {
    inputDataValues.size(),
    std::vector< TableLayout::CellLayout >( inputDataValues[0].size(), TableLayout::CellLayout())
  };

  std::vector< size_t > & subDataLineCounts = tableLayout.getNbSubDataLines();
  for( size_t idxRow = 0; idxRow < inputDataValues.size(); ++idxRow )
  {
    size_t maxLinesPerRow  = 0;
    size_t idxColumn = 0;
    for( auto it = tableLayout.beginLeaf(); it != tableLayout.endLeaf(); ++it )
    {
      if( !it->hasChild())
      {
        TableData::CellData & cell = inputDataValues[idxRow][idxColumn];
        TableLayout::CellAlignment const cellAlignement = it->m_cellAlignment;
        TableLayout::Alignment const alignement = cell.type == CellType::Header ?
                                                  cellAlignement.headerAlignment :
                                                  cellAlignement.valueAlignment;

        if( it->m_columName.m_cellType == CellType::Hidden )
        {
          cell.type = it->m_columName.m_cellType;
        }
        cellsDataLayout[idxRow][idxColumn] = TableLayout::CellLayout( cell.type, cell.value, alignement );
        maxLinesPerRow  = std::max( maxLinesPerRow, cellsDataLayout[idxRow][idxColumn].m_lines.size() );

        idxColumn++;

      }
    }
    subDataLineCounts.push_back( { maxLinesPerRow } );
  }
}

void TableTextFormatter::updateColumnMaxLength( TableLayout & tableLayout,
                                                CellLayoutRows & cellsDataLayout ) const
{
  auto getMaxLineLength = []( std::vector< std::string > const & lines ) -> size_t
  {
    auto maxLineLength = std::max_element( lines.begin(), lines.end(),
                                           []( const std::string & a, const std::string & b ) {
      return a.length() < b.length();
    } );
    return (maxLineLength != lines.end()) ? maxLineLength->length() : 0;
  };

  size_t idxColumn = 0;
  for( auto it = tableLayout.beginLeaf(); it != tableLayout.endLeaf(); ++it )
  {
    size_t maxColumnSize = 1;
    TableLayout::Column * currentColumn = &(*it);
    if( !it->hasChild())
    {
      // find max cell length between cell data column and current cell layout
      size_t cellLayoutLength = getMaxLineLength( currentColumn->m_columName.m_lines );
      for( size_t rowIdx = 0; rowIdx < cellsDataLayout.size(); ++rowIdx )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[rowIdx][idxColumn];
        size_t const cellDataLength = getMaxLineLength( cell.m_lines );
        cell.m_maxDataLength = std::max( cellLayoutLength, cellDataLength );
        maxColumnSize = std::max( cell.m_maxDataLength, maxColumnSize );
      }

      // update all current data cell, route by column
      for( size_t rowIdx = 0; rowIdx < cellsDataLayout.size(); ++rowIdx )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[rowIdx][idxColumn];
        bool isPreviousCellMerged = (idxColumn != 0) && (cellsDataLayout[rowIdx][idxColumn - 1].m_cellType == CellType::MergeNext);
        if( isPreviousCellMerged )
        {
          size_t const previousCellMaxLength = cellsDataLayout[rowIdx][idxColumn - 1].m_maxDataLength;
          size_t const additionalMargin = tableLayout.getColumnMargin();
          cell.m_maxDataLength = maxColumnSize + previousCellMaxLength + additionalMargin;
          cellsDataLayout[rowIdx][idxColumn - 1].m_maxDataLength = 0;
        }
        else
        {
          cell.m_maxDataLength = maxColumnSize;
        }

        if( cellsDataLayout[rowIdx][idxColumn].m_cellType == CellType::Separator )
        {
          size_t separatorLength = cell.m_maxDataLength;//
          separatorLength += idxColumn ==  cellsDataLayout[0].size() - 1 ? 
                             tableLayout.getBorderMargin() * 2 + 2:
                             tableLayout.getColumnMargin();

          cell.m_lines[0] = string( separatorLength, '-' );
        }
      }
      currentColumn->setMaxStringSize( maxColumnSize );
      ++idxColumn;
    }
    else
    { // Handle columns with subcolumns
      size_t subColumnsLength = 0;
      for( auto const & m_subColumn : currentColumn->m_subColumn )
      {
        subColumnsLength += m_subColumn.getMaxStringSize();
      }

      subColumnsLength += (currentColumn->m_subColumn.size() - 1) * (size_t)tableLayout.getColumnMargin();
      currentColumn->setMaxStringSize( std::max( subColumnsLength, currentColumn->getMaxStringSize() ) );
    }
  }
}

void TableTextFormatter::calculateTableSeparators( TableLayout & tableLayout,
                                                   CellLayoutRows & cellsDataLayout,
                                                   string & sectionSeparatingLine,
                                                   string & topSeparator ) const
{
  string const tableTitle = string( tableLayout.getTitle() );
  size_t const margins = (size_t) tableLayout.getBorderMargin() * 2;
  size_t const horizontalBar = 2;

  size_t sectionlineLength = 0;
  size_t nbColumns = 0;

  for( auto const & column : tableLayout.getColumns() )
  {
    if( column.m_columName.m_cellType != CellType::Hidden )
    {
      sectionlineLength += column.getMaxStringSize();
      nbColumns++;
    }
  }

  size_t const spacingBetweenColumns = (nbColumns - 1) * (size_t) tableLayout.getColumnMargin();
  sectionlineLength += spacingBetweenColumns + margins + horizontalBar;

  size_t maxTopLineLength =  tableTitle.length() + margins + horizontalBar;
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );

  if( sectionlineLength < maxTopLineLength )
  {
    size_t const paddingCharacters  = maxTopLineLength - sectionlineLength;
    adjustColumnWidths( tableLayout, cellsDataLayout, nbColumns, paddingCharacters );
    sectionlineLength = maxTopLineLength;
  }

  sectionSeparatingLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
  integer const topSeparatorLength = maxTopLineLength - 2; // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::adjustColumnWidths( TableLayout & tableLayout,
                                             CellLayoutRows & cellsDataLayout,
                                             size_t const nbColumns,
                                             size_t const paddingCharacters ) const
{
  for( auto it = tableLayout.beginRoot(); it != tableLayout.endRoot(); ++it )
  {
    size_t remainingPaddingForLastColumn = paddingCharacters % nbColumns;
    size_t paddingPerColumn = std::floor( paddingCharacters / nbColumns );

    TableLayout::Column *currentColumn = &(*it);

    if( it.getCurrentLayer() != 0 )
    {
      std::vector< size_t > nbDivider;
      TableLayout::Column *columnHierarchy = currentColumn;
      size_t lastParentPaddingRemaining = 0;
      while( columnHierarchy->m_parent != nullptr )
      {
        nbDivider.push_back( columnHierarchy->m_parent->m_subColumn.size());
        columnHierarchy = columnHierarchy->m_parent;
      }

      if( columnHierarchy->m_next == nullptr )
      {
        lastParentPaddingRemaining = remainingPaddingForLastColumn;
      }

      for( auto subColumnDivider = nbDivider.rbegin(); subColumnDivider != nbDivider.rend(); ++subColumnDivider )
      {
        remainingPaddingForLastColumn = (paddingPerColumn % (*subColumnDivider) + lastParentPaddingRemaining);
        paddingPerColumn = std::floor( paddingPerColumn / (*subColumnDivider));
      }
    }
    TableLayout::Column * nextVisibleColumn = currentColumn;
    bool isLastVisibleColumn = false;
    while( nextVisibleColumn->m_next != nullptr && nextVisibleColumn->m_next->m_columName.m_cellType == CellType::Hidden )
    {
      nextVisibleColumn = nextVisibleColumn->m_next;
    }
    if( nextVisibleColumn->m_next == nullptr )
    {
      isLastVisibleColumn = true;
    }

    size_t const additionalPadding = (isLastVisibleColumn) ? remainingPaddingForLastColumn : 0;
    currentColumn->setMaxStringSize( currentColumn->getMaxStringSize() + paddingPerColumn + additionalPadding );
  }
  //set new max string size
  size_t idxColumn = 0;
  for( auto it = tableLayout.beginLeaf(); it !=  tableLayout.endLeaf(); ++it )
  {
    if( !it->hasChild() )
    {
      for( size_t rowIdx = 0; rowIdx< cellsDataLayout.size(); rowIdx++ )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[rowIdx][idxColumn];
        cell.m_maxDataLength = it->getMaxStringSize();
      }
      idxColumn++;
    }
  }
}

void TableTextFormatter::gridifyHeaders( TableLayout & tableLayout,
                                         CellLayoutRows & cellsHeaderLayout ) const
{
  size_t const headerLayersCount = tableLayout.getMaxHeaderRow();
  std::vector< size_t > & sublineHeaderCounts = tableLayout.getSublineInHeaderCounts();
  sublineHeaderCounts.resize( headerLayersCount, 1 );
  for( auto it = tableLayout.beginLeaf(); it !=  tableLayout.endLeaf(); ++it )
  {
    size_t const currentLayer = it.getCurrentLayer();
    TableLayout::CellLayout currentCell = it->m_columName;
    if( !it->hasChild() && headerLayersCount - 1 != currentLayer )
    {
      for( size_t idxRow = currentLayer; idxRow< headerLayersCount - 1; idxRow++ )
      {
        TableLayout::CellLayout emptyCell{CellType::Header, "", TableLayout::Alignment::center};
        emptyCell.m_maxDataLength = it->getMaxStringSize();
        cellsHeaderLayout[idxRow + 1].push_back( emptyCell );
      }
    }

    currentCell.m_maxDataLength = it->getMaxStringSize();
    if( it->hasParent() )
    {
      it->m_parent->m_columName.m_cellWidth++;
    }

    sublineHeaderCounts[currentLayer] = std::max( sublineHeaderCounts[currentLayer], currentCell.m_lines.size() );
    cellsHeaderLayout[currentLayer].push_back( currentCell );

    //subdivision of parent cell = current cell + subColumns - 1
    if( currentCell.m_cellWidth > 1 )
    {
      currentCell.m_cellWidth--;
    }

    for( size_t idxColumn = 0; idxColumn < currentCell.m_cellWidth; idxColumn++ )
    {
      TableLayout::CellLayout mergingCell{ CellType::MergeNext, "", TableLayout::Alignment::center };
      cellsHeaderLayout[currentLayer].push_back( mergingCell );
    }
  }

  size_t idxLayer = 0;
  for( auto & lines: cellsHeaderLayout )
  {
    size_t nbLinesInLayer =  sublineHeaderCounts[idxLayer];
    if( nbLinesInLayer != 1 )
    {
      for( auto & cell : lines )
      {
        cell.m_lines.resize( nbLinesInLayer, "" );
      }
    }
    idxLayer++;
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
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
    tableOutput << buildCell( TableLayout::Alignment::center,
                              tableTitle,
                              topSeparator.length() - (tableLayout.getBorderMargin() * 2) - 2 );
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, tableLayout.getBorderMargin() + 1 );
  }
}

void TableTextFormatter::formatCell( TableLayout & tableLayout,
                                     std::ostringstream & tableOutput,
                                     TableLayout::CellLayout const & cell,
                                     size_t idxLine ) const
{
  GEOS_UNUSED_VAR( tableLayout );
  tableOutput << buildCell( cell.m_alignment,
                            cell.m_lines[idxLine],
                            cell.m_maxDataLength );
}

void TableTextFormatter::outputLines( TableLayout & tableLayout,
                                      CellLayoutRows const & cellsLayout,
                                      std::ostringstream & tableOutput,
                                      std::vector< size_t > const & nbLinesRow,
                                      CellType sectionType,
                                      string_view sectionSeparatingLine ) const
{
  size_t idxLine = 0;
  for( auto const & line : cellsLayout )
  {
    for( size_t idxSubLine = 0; idxSubLine < nbLinesRow[idxLine]; idxSubLine++ )
    {

      for( auto const & cell : line )
      {

        if( &cell == &(line.front()) && cell.m_cellType != CellType::Separator )
        {
          tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
        }

        if( cell.m_cellType != CellType::Hidden && cell.m_cellType != CellType::MergeNext )
        {
          formatCell( tableLayout, tableOutput, cell, idxSubLine );
        }

        if( cell.m_cellType == CellType::MergeNext && (&cell == &(line.back())))
        {
          formatCell( tableLayout, tableOutput, cell, idxSubLine );
          tableOutput << GEOS_FMT( "{:>{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
        }

        if( cell.m_cellType  == CellType::Header || cell.m_cellType  == CellType::Value )
        {
          if( &cell == &(line.back()) )
          {
            tableOutput << GEOS_FMT( "{:>{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
          }
          else
          {
            tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, tableLayout.getColumnMargin());
          }
        }

      }
      tableOutput << "\n";
    }

    if( sectionType == CellType::Header )
    {
      tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
    }
    idxLine++;
  }

  if( sectionType == CellType::Value && !cellsLayout.empty())
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }

}
}
