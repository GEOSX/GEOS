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
  string result;
  static constexpr string_view separator = ",";

  size_t total_size = 0;
  for( const auto & column : m_tableLayout.getColumns())
  {
    for( const auto & str : column.m_columName.m_lines )
    {
      total_size += str.size();
    }
    total_size += 1;
  }
  result.reserve( total_size );

  for( std::size_t idxColumn = 0; idxColumn < m_tableLayout.getColumns().size(); ++idxColumn )
  {
    std::ostringstream strValue;
    for( auto const & str :  m_tableLayout.getColumns()[idxColumn].m_columName.m_lines )
    {
      result.append( str );
    }

    if( idxColumn < m_tableLayout.getColumns().size() - 1 )
    {
      result.append( separator );
    }
  }
  result.append( "\n" );
  return result;
}

string TableCSVFormatter::dataToString( TableData const & tableData ) const
{

  RowsCellInput const rowsValues( tableData.getTableDataRows() );
  string result;
  size_t total_size = 0;
  for( const auto & row : rowsValues )
  {
    for( const auto & item : row )
    {
      total_size += item.value.size();
    }
    total_size += row.size() - 1;
    total_size += 1;
  }

  result.reserve( total_size );

  for( const auto & row : rowsValues )
  {
    std::vector< string > rowConverted;
    for( const auto & item : row )
    {
      rowConverted.push_back( item.value );
    }
    result.append( stringutilities::join( rowConverted.cbegin(), rowConverted.cend(), "," ));
    result.append( "\n" );
  }

  return result;
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
  string separatorLine;

  TableLayout tableLayout = m_tableLayout;
  initalizeTableLayout( tableLayout, tableData, cellsHeaderLayout, cellsDataLayout, separatorLine );
  outputTable( tableLayout, tableOutput, cellsHeaderLayout, cellsDataLayout, separatorLine );
  return tableOutput.str();
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  CellLayoutRows cellsHeaderLayout;
  CellLayoutRows cellsDataLayout;
  string separatorLine;

  TableLayout tableLayout = m_tableLayout;
  initalizeTableLayout( tableLayout, tableData,
                        cellsHeaderLayout, cellsDataLayout,
                        separatorLine );
  outputTable( tableLayout, tableOutput,
               cellsHeaderLayout, cellsDataLayout,
               separatorLine );
  return tableOutput.str();
}

void TableTextFormatter::initalizeTableLayout( TableLayout & tableLayout,
                                               TableData const & tableData,
                                               CellLayoutRows & cellsHeaderLayout,
                                               CellLayoutRows & cellsDataLayout,
                                               string & separatorLine ) const
{
  setLinks( tableLayout.getColumns()  );

  populateHeaderCellsLayout( tableLayout, cellsHeaderLayout );

  RowsCellInput inputDataValues( tableData.getTableDataRows());
  if( inputDataValues.size() > 0 )
  {
    populateDataCellsLayout( tableLayout, cellsDataLayout, inputDataValues );
  }

  updateColumnMaxLength( tableLayout, cellsHeaderLayout, cellsDataLayout );

  calculateTableSeparators( tableLayout, cellsHeaderLayout, cellsDataLayout, separatorLine );
}

void TableTextFormatter::outputTable( TableLayout & tableLayout,
                                      std::ostringstream & tableOutput,
                                      CellLayoutRows const & cellsHeader,
                                      CellLayoutRows const & cellsData,
                                      string_view separatorLine ) const
{
  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableLayout, tableOutput, separatorLine );
  tableOutput << GEOS_FMT( "{}\n", separatorLine );
  outputLines( tableLayout, cellsHeader, tableOutput, tableLayout.getSublineInHeaderCounts(),
               CellType::Header, separatorLine );
  outputLines( tableLayout, cellsData, tableOutput, tableLayout.getNbSubDataLines(),
               CellType::Value, separatorLine );
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

void TableTextFormatter::populateHeaderCellsLayout( TableLayout & tableLayout,
                                                    CellLayoutRows & cellsHeaderLayout ) const
{
  cellsHeaderLayout.resize( tableLayout.getMaxHeaderRow() );
  size_t const headerLayersCount = tableLayout.getMaxHeaderRow();
  std::vector< size_t > & sublineHeaderCounts = tableLayout.getSublineInHeaderCounts();
  sublineHeaderCounts.resize( headerLayersCount, 1 );

  for( auto it = tableLayout.beginLeaf(); it !=  tableLayout.endLeaf(); ++it )
  {
    size_t const currentLayer = it.getCurrentLayer();
    TableLayout::CellLayout currentCell = it->m_columName;

    if( !it->hasChild() && headerLayersCount - 1 != currentLayer )
    {
      for( size_t idxRow = currentLayer; idxRow < headerLayersCount - 1; idxRow++ )
      {
        TableLayout::CellLayout emptyCell{CellType::Header, "", TableLayout::Alignment::center};
        emptyCell.m_maxDataLength = it->getMaxStringSize();
        cellsHeaderLayout[idxRow + 1].push_back( emptyCell );
      }
    }
    currentCell.m_maxDataLength = it->getMaxStringSize();

    if( it->hasParent() )
    {
      it->m_parent->m_cellWidth += it->m_cellWidth == 0 ? 1 : it->m_cellWidth;
    }
    if( it->m_cellWidth > 1 )
    {
      it->m_cellWidth--;
    }

    sublineHeaderCounts[currentLayer] = std::max( sublineHeaderCounts[currentLayer],
                                                  currentCell.m_lines.size() );
    for( size_t idxColumn = 0; idxColumn < it->m_cellWidth; idxColumn++ )
    {
      TableLayout::CellLayout mergingCell{ CellType::MergeNext, "", TableLayout::Alignment::center };
      cellsHeaderLayout[currentLayer].push_back( mergingCell );
    }

    cellsHeaderLayout[currentLayer].push_back( currentCell );
  }

  size_t idxLayer = 0;
  for( auto & lines: cellsHeaderLayout )
  {
    size_t nbLinesInLayer = sublineHeaderCounts[idxLayer];

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

void TableTextFormatter::populateDataCellsLayout( TableLayout & tableLayout,
                                                  CellLayoutRows & cellsDataLayout,
                                                  RowsCellInput & inputDataValues ) const
{
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
        if( it->m_columName.m_cellType == CellType::Separator )
        {
          cell.value = m_horizontalLine;
        }
        cellsDataLayout[idxRow][idxColumn] = TableLayout::CellLayout( cell.type, cell.value, alignement );
        maxLinesPerRow  = std::max( maxLinesPerRow, cellsDataLayout[idxRow][idxColumn].m_lines.size() );
        idxColumn++;
      }
    }
    subDataLineCounts.push_back( { maxLinesPerRow } );
  }

  size_t idxLayer = 0;
  for( auto & lines: cellsDataLayout )
  {
    size_t nbLinesInLayer =  subDataLineCounts[idxLayer];
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

void TableTextFormatter::updateColumnMaxLength( TableLayout & tableLayout,
                                                CellLayoutRows & cellsHeaderLayout,
                                                CellLayoutRows & cellsDataLayout ) const
{
  auto getMaxLineLength = []( std::vector< std::string > const & lines ) -> size_t
  {
    auto maxLineLength = std::max_element( lines.begin(), lines.end(),
                                           []( std::string const & a, std::string const & b ) {
      return a.length() < b.length();
    } );
    return (maxLineLength != lines.end()) ? maxLineLength->length() : 0;
  };

  auto updateCellMaxLength = [&tableLayout]( TableLayout::CellLayout & cell, size_t maxColumnSize,
                                             TableLayout::CellLayout * previousCell = nullptr )
  {
    if( previousCell && previousCell->m_cellType == CellType::MergeNext )
    {
      size_t const previousCellMaxLength = previousCell->m_maxDataLength;
      size_t const additionalMargin = tableLayout.getColumnMargin();
      cell.m_maxDataLength = maxColumnSize + previousCellMaxLength + additionalMargin;
      previousCell->m_maxDataLength = 0;
    }
    else
    {
      cell.m_maxDataLength = maxColumnSize;
    }
  };

  size_t const numColumns = cellsHeaderLayout[0].size();
  //each idx per row
  std::vector< size_t > accMaxStringColumn( cellsDataLayout.size(), 0 );
  for( size_t idxColumn = 0; idxColumn < numColumns; ++idxColumn )
  {
    size_t maxColumnSize = 1;

    // init header column max length
    for( size_t rowIdx = 0; rowIdx < cellsDataLayout.size(); ++rowIdx )
    {
      size_t const cellDataLength = getMaxLineLength( cellsDataLayout[rowIdx][idxColumn].m_lines );
      maxColumnSize = std::max( {maxColumnSize, cellDataLength} );
    }

    for( size_t rowIdx = 0; rowIdx < cellsHeaderLayout.size(); ++rowIdx )
    {
      size_t const cellHeaderLength = getMaxLineLength( cellsHeaderLayout[rowIdx][idxColumn].m_lines );
      maxColumnSize = std::max( {maxColumnSize, cellHeaderLength} );
      cellsHeaderLayout[rowIdx][idxColumn].m_maxDataLength = maxColumnSize;
    }

    // update maxColumnSize for data cell
    for( size_t rowIdx = 0; rowIdx < cellsDataLayout.size(); ++rowIdx )
    {
      TableLayout::CellLayout & dataCell = cellsDataLayout[rowIdx][idxColumn];
      TableLayout::CellLayout * previousDataCell = (idxColumn > 0) ? &cellsDataLayout[rowIdx][idxColumn - 1] : nullptr;

      if( dataCell.m_cellType == CellType::MergeNext )
      {
        accMaxStringColumn[rowIdx] += cellsHeaderLayout[0][idxColumn].m_maxDataLength + tableLayout.getColumnMargin();
      }

      if( idxColumn > 0 &&
          previousDataCell->m_cellType == CellType::MergeNext && dataCell.m_cellType != CellType::MergeNext )
      {
        size_t sumOfMergingCell = accMaxStringColumn[rowIdx] + cellsHeaderLayout[0][idxColumn].m_maxDataLength;
        if( sumOfMergingCell <  dataCell.m_maxDataLength )
        {
          maxColumnSize += dataCell.m_maxDataLength - sumOfMergingCell;
          for( size_t rowIdx2 = 0; rowIdx2 < rowIdx; rowIdx2++ )
          {
            TableLayout::CellLayout * previousDataCellTemp = (idxColumn > 0) ? &cellsDataLayout[rowIdx2][idxColumn - 1] : nullptr;
            if( previousDataCellTemp )
            {
              previousDataCellTemp->m_maxDataLength = cellsHeaderLayout[0][idxColumn - 1].m_maxDataLength;
            }
            updateCellMaxLength( cellsDataLayout[rowIdx2][idxColumn], maxColumnSize, previousDataCellTemp );
            if( cellsDataLayout[rowIdx2][idxColumn].m_cellType == CellType::Separator )
            {
              size_t const additionnalMargin = (idxColumn == numColumns - 1) ?
                                               tableLayout.getBorderMargin() * 2 + 2 :
                                               tableLayout.getColumnMargin();
              cellsDataLayout[rowIdx2][idxColumn].m_lines[0]  = std::string( maxColumnSize + additionnalMargin, '-' );
            }
          }
        }

        accMaxStringColumn[rowIdx] = 0;
      }
      else
      {
        size_t const cellDataLength = getMaxLineLength( cellsDataLayout[rowIdx][idxColumn].m_lines );
        maxColumnSize = std::max( {maxColumnSize, cellDataLength} );
      }

      updateCellMaxLength( dataCell, maxColumnSize, previousDataCell );

      if( dataCell.m_cellType == CellType::Separator )
      {
        size_t separatorLength = dataCell.m_maxDataLength;
        separatorLength += (idxColumn == numColumns - 1) ?
                           tableLayout.getBorderMargin() * 2 + 2 :
                           tableLayout.getColumnMargin();

        dataCell.m_lines[0] = std::string( separatorLength, '-' );
      }
    }

    // last pass for updating cells header
    for( size_t rowIdx = 0; rowIdx < cellsHeaderLayout.size(); ++rowIdx )
    {
      TableLayout::CellLayout & headerCell = cellsHeaderLayout[rowIdx][idxColumn];
      TableLayout::CellLayout * previousHeaderCell = (idxColumn > 0) ?
                                                     &cellsHeaderLayout[rowIdx][idxColumn - 1]:
                                                     nullptr;

      updateCellMaxLength( headerCell, maxColumnSize, previousHeaderCell );
    }
  }
}

void TableTextFormatter::calculateTableSeparators( TableLayout & tableLayout,
                                                   CellLayoutRows & cellsHeaderLayout,
                                                   CellLayoutRows & cellsDataLayout,
                                                   string & separatorLine ) const
{
  std::string const tableTitle = std::string( tableLayout.getTitle() );
  size_t const margins = (size_t) tableLayout.getBorderMargin() * 2;
  size_t const horizontalBar = 2;

  size_t sectionlineLength = 0;
  size_t nbColumns = 0;
  size_t nbHiddenColumns = 0;

  for( auto const & column : cellsHeaderLayout[0] )
  {
    if( column.m_cellType == CellType::Hidden )
    {
      nbHiddenColumns++;
    }
    if( column.m_cellType != CellType::Hidden && column.m_cellType != CellType::MergeNext )
    {
      sectionlineLength += column.m_maxDataLength;
      nbColumns++;
    }
  }

  size_t const spacingBetweenColumns = (nbColumns - 1) * (size_t) tableLayout.getColumnMargin();
  sectionlineLength += spacingBetweenColumns + margins + horizontalBar;

  size_t titleRowLength = tableTitle.length() + margins + horizontalBar;
  size_t maxLength = std::max( {titleRowLength, sectionlineLength} );
  if( sectionlineLength < maxLength )
  {
    size_t const paddingCharacters = maxLength - sectionlineLength;
    adjustColumnWidths( cellsHeaderLayout, nbHiddenColumns, paddingCharacters );
    adjustColumnWidths( cellsDataLayout, nbHiddenColumns, paddingCharacters );
    sectionlineLength = maxLength;
  }

  separatorLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
}

void TableTextFormatter::adjustColumnWidths( CellLayoutRows & cells,
                                             size_t nbHiddenColumns,
                                             size_t const paddingCharacters ) const
{
  size_t const numRows = cells.size();
  size_t const nbColumns = cells[0].size();

  size_t remainingPaddingForLastColumn = paddingCharacters % (nbColumns - nbHiddenColumns);
  size_t paddingPerColumn = std::floor( paddingCharacters / (nbColumns - nbHiddenColumns));
  for( size_t idxRow = 0; idxRow < numRows; ++idxRow )
  {
    for( size_t idxColumn = 0; idxColumn < nbColumns; ++idxColumn )
    {
      auto & currentCell = cells[idxRow][idxColumn];

      if( currentCell.m_cellType != CellType::Hidden )
      {
        size_t nextIdxColumn = idxColumn + 1;

        while( nextIdxColumn < nbColumns && cells[idxRow][nextIdxColumn].m_cellType == CellType::Hidden )
        {
          nextIdxColumn++;
        }
        bool isLastVisibleColumn = nextIdxColumn == nbColumns;

        if( idxColumn > 0 && cells[idxRow][idxColumn - 1].m_cellType == CellType::MergeNext )
        {
          auto * previousCell = &cells[idxRow][idxColumn - 1];
          currentCell.m_maxDataLength += previousCell->m_maxDataLength;
          previousCell->m_maxDataLength = 0;
        }

        size_t const additionalPadding = (isLastVisibleColumn || idxColumn == nbColumns - 1) ?
                                         remainingPaddingForLastColumn :
                                         0;

        currentCell.m_maxDataLength += paddingPerColumn + additionalPadding;
      }
    }
  }

}

void TableTextFormatter::outputTitleRow( TableLayout & tableLayout,
                                         std::ostringstream & tableOutput,
                                         string_view separatorLine ) const
{
  string const tableTitle = string( tableLayout.getTitle());
  if( !tableTitle.empty() )
  {
    tableOutput << GEOS_FMT( "{}\n", separatorLine );
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
    tableOutput << buildCell( TableLayout::Alignment::center,
                              tableTitle,
                              separatorLine.length() - (tableLayout.getBorderMargin() * 2) - 2 );
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
                                      string_view separatorLine ) const
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
      tableOutput <<std::endl;
    }

    if( sectionType == CellType::Header )
    {
      tableOutput << GEOS_FMT( "{}\n", separatorLine );
    }
    idxLine++;
  }

  if( sectionType == CellType::Value && !cellsLayout.empty())
  {
    tableOutput << GEOS_FMT( "{}\n", separatorLine );
  }

}
}
