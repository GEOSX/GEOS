/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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

static constexpr string_view separator = ",";
string TableCSVFormatter::headerToString() const
{
  string result;

  size_t total_size = 0;
  for( const auto & column : m_tableLayout.getColumns())
  {
    for( const auto & str : column.m_header.m_lines )
    {
      total_size += str.size();
    }
    total_size += separator.size();
  }
  result.reserve( total_size );

  for( std::size_t idxColumn = 0; idxColumn < m_tableLayout.getColumns().size(); ++idxColumn )
  {
    std::ostringstream strValue;
    for( auto const & str :  m_tableLayout.getColumns()[idxColumn].m_header.m_lines )
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
    total_size += row.size();
  }

  result.reserve( total_size );

  for( const auto & row : rowsValues )
  {
    std::vector< string > rowConverted;
    for( const auto & item : row )
    {
      std::istringstream strStream( item.value );
      string line;
      bool detectNewLine = false;
      while( getline( strStream, line, '\n' ))
      {
        rowConverted.push_back( line );
        detectNewLine = true;
      }

      if( !detectNewLine )
        rowConverted.push_back( item.value );
    }
    result.append( stringutilities::join( rowConverted.cbegin(), rowConverted.cend(), separator ));
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
  TableData tableData;
  return toString( tableData );
}

template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const
{
  std::ostringstream tableOutput;
  CellLayoutRows cellsHeaderLayout;
  CellLayoutRows cellsDataLayout;
  string separatorLine;

  TableLayout tableLayout = m_tableLayout;
  size_t nbVisibleColumn = 0;

  initalizeTableLayout( tableLayout, tableData,
                        cellsHeaderLayout, cellsDataLayout,
                        separatorLine,
                        nbVisibleColumn );
  outputTable( tableLayout, tableOutput,
               cellsHeaderLayout, cellsDataLayout,
               separatorLine,
               nbVisibleColumn );
  return tableOutput.str();
}

void TableTextFormatter::initalizeTableLayout( TableLayout & tableLayout,
                                               TableData const & tableData,
                                               CellLayoutRows & cellsHeaderLayout,
                                               CellLayoutRows & cellsDataLayout,
                                               string & separatorLine,
                                               size_t & nbVisibleColumn ) const
{
  setLinks( tableLayout.getColumns()  );

  populateHeaderCellsLayout( tableLayout, cellsHeaderLayout );

  RowsCellInput inputDataValues( tableData.getTableDataRows());
  if( inputDataValues.size() > 0 )
  {
    populateDataCellsLayout( tableLayout, cellsDataLayout, inputDataValues );
  }

  updateColumnMaxLength( tableLayout, cellsHeaderLayout, cellsDataLayout );

  adjustTableWidth( tableLayout, cellsHeaderLayout, cellsDataLayout, separatorLine, nbVisibleColumn );
}

void TableTextFormatter::outputTable( TableLayout & tableLayout,
                                      std::ostringstream & tableOutput,
                                      CellLayoutRows const & cellsHeader,
                                      CellLayoutRows const & cellsData,
                                      string_view separatorLine,
                                      size_t & nbVisibleColumn ) const
{
  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
  outputTitleRow( tableLayout, tableOutput, separatorLine );
  tableOutput << GEOS_FMT( "{}\n", separatorLine );
  outputLines( tableLayout, cellsHeader, tableOutput, tableLayout.getSublineInHeaderCounts(),
               CellType::Header, separatorLine, nbVisibleColumn );
  outputLines( tableLayout, cellsData, tableOutput, tableLayout.getNbSubDataLines(),
               CellType::Value, separatorLine, nbVisibleColumn );
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
      columns[idxColumn].setNextCell( &columns[idxColumn + 1] );
    }

    if( !columns[idxColumn].m_subColumn.empty())
    {
      for( auto & subCol : columns[idxColumn].m_subColumn )
      {
        subCol.setParent( &columns[idxColumn] );
      }

      setLinks( columns[idxColumn].m_subColumn );
    }
  }
}

void TableTextFormatter::populateHeaderCellsLayout( TableLayout & tableLayout,
                                                    CellLayoutRows & cellsHeaderLayout ) const
{
  cellsHeaderLayout.resize( tableLayout.getMaxDepth() );
  size_t const headerLayersCount = tableLayout.getMaxDepth();
  std::vector< size_t > & sublineHeaderCounts = tableLayout.getSublineInHeaderCounts();
  sublineHeaderCounts.resize( headerLayersCount, 1 );

  for( auto it = tableLayout.beginDeepFirst(); it !=  tableLayout.endDeepFirst(); ++it )
  {
    size_t const currentLayer = it.getCurrentLayer();
    TableLayout::CellLayout currentCell = it->m_header;

    if( !it->hasChild() && headerLayersCount - 1 != currentLayer )
    {
      for( size_t idxRow = currentLayer; idxRow < headerLayersCount - 1; idxRow++ )
      {
        TableLayout::CellLayout emptyCell{CellType::Header, "", TableLayout::Alignment::center};
        emptyCell.m_cellWidth = it->m_header.m_cellWidth;
        cellsHeaderLayout[idxRow + 1].push_back( emptyCell );
      }
    }
    currentCell.m_cellWidth = it->m_header.m_cellWidth;

    if( it->hasParent() )
    {
      if( it->getNumberCellMerge() == 0 )
      {
        it->getParent()->incrementMergeHeaderCount( 1 );
      }
      else
      {
        it->getParent()->incrementMergeHeaderCount( it->getNumberCellMerge() );
      }
    }

    if( it->getNumberCellMerge()  > 1 )
    {
      it->decrementMergeHeaderCount();
    }

    sublineHeaderCounts[currentLayer] = std::max( sublineHeaderCounts[currentLayer],
                                                  currentCell.m_lines.size() );

    for( size_t idxColumn = 0; idxColumn < it->getNumberCellMerge(); idxColumn++ )
    {
      TableLayout::CellLayout mergingCell{ CellType::MergeNext, "", TableLayout::Alignment::center };
      cellsHeaderLayout[currentLayer].push_back( mergingCell );
    }

    cellsHeaderLayout[currentLayer].push_back( currentCell );
  }

  size_t idxLayer = 0;
  for( auto & lines: cellsHeaderLayout )
  {
    size_t const nbLinesInLayer = sublineHeaderCounts[idxLayer];

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
    for( auto it = tableLayout.beginDeepFirst(); it != tableLayout.endDeepFirst(); ++it )
    {
      if( !it->hasChild())
      {
        TableData::CellData & cell = inputDataValues[idxRow][idxColumn];
        TableLayout::ColumnAlignement const cellAlignement = it->m_alignment;
        TableLayout::Alignment const alignement = cell.type == CellType::Header ?
                                                  cellAlignement.headerAlignment :
                                                  cellAlignement.valueAlignment;

        if( it->m_header.m_cellType == CellType::Hidden ||
            it->m_header.m_cellType == CellType::Disabled )
        {
          cell.type = it->m_header.m_cellType;
        }
        if( it->m_header.m_cellType == CellType::Separator )
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
      size_t const previousCellMaxLength = previousCell->m_cellWidth;
      size_t const additionalMargin = tableLayout.getColumnMargin();
      cell.m_cellWidth = maxColumnSize + previousCellMaxLength + additionalMargin;
      previousCell->m_cellWidth = 0;
    }
    else
    {
      cell.m_cellWidth = maxColumnSize;
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
      if( idxColumn == 0 ||
          (idxColumn > 0 && cellsDataLayout[rowIdx][idxColumn - 1].m_cellType != CellType::MergeNext))
      {
        maxColumnSize = std::max( maxColumnSize, cellDataLength );
      }
    }

    for( size_t rowIdx = 0; rowIdx < cellsHeaderLayout.size(); ++rowIdx )
    {
      size_t const cellHeaderLength = getMaxLineLength( cellsHeaderLayout[rowIdx][idxColumn].m_lines );
      maxColumnSize = std::max( {maxColumnSize, cellHeaderLength} );
      cellsHeaderLayout[rowIdx][idxColumn].m_cellWidth = maxColumnSize;
    }

    // update maxColumnSize for data cell
    for( size_t rowIdx = 0; rowIdx < cellsDataLayout.size(); ++rowIdx )
    {
      TableLayout::CellLayout & dataCell = cellsDataLayout[rowIdx][idxColumn];
      TableLayout::CellLayout * previousDataCell = (idxColumn > 0) ?
                                                   &cellsDataLayout[rowIdx][idxColumn - 1] :
                                                   nullptr;

      if( dataCell.m_cellType == CellType::MergeNext )
      {
        accMaxStringColumn[rowIdx] += cellsHeaderLayout[0][idxColumn].m_cellWidth + tableLayout.getColumnMargin();
      }

      if( idxColumn > 0 &&
          previousDataCell->m_cellType == CellType::MergeNext && dataCell.m_cellType != CellType::MergeNext )
      {
        // root header cells know the maximum string size in the column
        size_t const sumOfMergingCell = accMaxStringColumn[rowIdx] + cellsHeaderLayout[0][idxColumn].m_cellWidth;
        if( sumOfMergingCell <  dataCell.m_cellWidth )
        {
          maxColumnSize -= dataCell.m_cellWidth - sumOfMergingCell;
          for( size_t rowIdx2 = 0; rowIdx2 < rowIdx; rowIdx2++ )
          {
            updateCellMaxLength( cellsDataLayout[rowIdx2][idxColumn], maxColumnSize );
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
        size_t separatorLength = dataCell.m_cellWidth;
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

void TableTextFormatter::adjustTableWidth( TableLayout & tableLayout,
                                           CellLayoutRows & cellsHeaderLayout,
                                           CellLayoutRows & cellsDataLayout,
                                           string & separatorLine,
                                           size_t & nbVisibleColumn ) const
{
  std::string const tableTitle = std::string( tableLayout.getTitle() );
  size_t const margins = (size_t) tableLayout.getBorderMargin() * 2;
  size_t const horizontalBar = 2;

  size_t sectionlineLength = 0;
  size_t nbHiddenColumns = 0;
  size_t nbParentColumn = 0;

  for( auto const & column : cellsHeaderLayout[0] )
  {
    std::cout << column.m_lines[0] << std::endl;
    if( column.m_cellType == CellType::Hidden ||  column.m_cellType == CellType::Disabled )
    {
      nbHiddenColumns++;
    }
    if( column.m_cellType != CellType::Disabled )
    {
      nbVisibleColumn++;
    }
    if( column.m_cellType == CellType::Value || column.m_cellType == CellType::Header )
    {
      sectionlineLength += column.m_cellWidth;
      //parent column has the max string size
      nbParentColumn++;
    }
  }


  size_t const spacingBetweenColumns = nbParentColumn == 0 ? (size_t) tableLayout.getColumnMargin():
                                       (nbParentColumn - 1) * (size_t) tableLayout.getColumnMargin();

  sectionlineLength += spacingBetweenColumns + margins + horizontalBar;

  size_t const titleRowLength = tableTitle.length() + margins + horizontalBar;
  size_t const maxLength = std::max( {titleRowLength, sectionlineLength} );
  if( sectionlineLength < maxLength )
  {
    size_t const paddingCharacters = maxLength - sectionlineLength;
    adjustColumnWidth( cellsHeaderLayout, nbHiddenColumns, paddingCharacters );
    adjustColumnWidth( cellsDataLayout, nbHiddenColumns, paddingCharacters );
    sectionlineLength = maxLength;
  }

  separatorLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
}

void TableTextFormatter::adjustColumnWidth( CellLayoutRows & cells,
                                            size_t const nbHiddenColumns,
                                            size_t const paddingCharacters ) const
{
  size_t const numRows = cells.size();
  size_t const nbColumns = cells[0].size();
  size_t const remainingPaddingForLastColumn = paddingCharacters % (nbColumns - nbHiddenColumns);
  size_t const paddingPerColumn = paddingCharacters / (nbColumns - nbHiddenColumns);
  for( size_t idxRow = 0; idxRow < numRows; ++idxRow )
  {
    for( size_t idxColumn = 0; idxColumn < nbColumns; ++idxColumn )
    {
      auto & currentCell = cells[idxRow][idxColumn];

      if( currentCell.m_cellType != CellType::Hidden &&
          currentCell.m_cellType != CellType::Disabled )
      {
        size_t nextIdxColumn = idxColumn + 1;
        while( nextIdxColumn < nbColumns &&
               ( cells[idxRow][nextIdxColumn].m_cellType == CellType::Hidden ||
                 cells[idxRow][nextIdxColumn].m_cellType == CellType::Disabled ) )
        {
          nextIdxColumn++;
        }
        bool isLastVisibleColumn = nextIdxColumn == nbColumns;

        if( idxColumn > 0 && cells[idxRow][idxColumn - 1].m_cellType == CellType::MergeNext )
        {
          auto * previousCell = &cells[idxRow][idxColumn - 1];
          currentCell.m_cellWidth += previousCell->m_cellWidth;
          previousCell->m_cellWidth = 0;
        }

        size_t const additionalPadding = (isLastVisibleColumn || idxColumn == nbColumns - 1) ?
                                         remainingPaddingForLastColumn :
                                         0;

        currentCell.m_cellWidth += paddingPerColumn + additionalPadding;
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
                            cell.m_cellWidth );
}

void TableTextFormatter::outputLines( TableLayout & tableLayout,
                                      CellLayoutRows const & cellsLayout,
                                      std::ostringstream & tableOutput,
                                      std::vector< size_t > const & nbLinesRow,
                                      CellType sectionType,
                                      string_view separatorLine,
                                      size_t & nbVisibleColumn ) const
{
  size_t idxLine = 0;
  for( auto const & row : cellsLayout )
  {
    for( size_t idxSubLine = 0; idxSubLine < nbLinesRow[idxLine]; idxSubLine++ )
    {
      size_t idxColumn = 0;
      size_t nbVisibleColumnTemp = nbVisibleColumn;
      while( nbVisibleColumnTemp > 0 )
      {
        auto & cell = row[idxColumn];
        switch( cell.m_cellType )
        {
          case CellType::Value:
          case CellType::Header:
            if( &cell == &(row.front()) )
            {
              tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
            }

            formatCell( tableLayout, tableOutput, cell, idxSubLine );

            if( &cell == &(row.back()) || nbVisibleColumnTemp == 1 )
            {
              tableOutput << GEOS_FMT( "{:>{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
            }
            else
            {
              tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, tableLayout.getColumnMargin());
            }
            nbVisibleColumnTemp--;
            break;
          case CellType::MergeNext:
            if( &cell == &(row.front()) )
            {
              tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
            }
            if( &cell == &(row.back()) )
            {
              formatCell( tableLayout, tableOutput, cell, idxSubLine );
              tableOutput << GEOS_FMT( "{:>{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
            }
            nbVisibleColumnTemp--;
            break;
          case CellType::Separator:
            formatCell( tableLayout, tableOutput, cell, idxSubLine );
            nbVisibleColumnTemp--;
            break;
          case CellType::Disabled:
            if( &cell == &(row.back()) )
            {
              tableOutput << GEOS_FMT( "{:>{}}", m_verticalLine, tableLayout.getBorderMargin() + 1 );
            }
          default:
            break;
        }
        idxColumn++;
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
    tableOutput << separatorLine;
  }
}
}
