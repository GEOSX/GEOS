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
                                               RowsCellLayout & cellsHeaderLayout,
                                               RowsCellLayout & cellsDataLayout,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{

  setLinks( tableLayout.getColumns()  );

  RowsCellInput inputDataValues( tableData.getTableDataRows());

  if( inputDataValues.size() > 0 )
  {
    populateDataCellsLayout( tableLayout, cellsDataLayout, inputDataValues );
  }

  std::cout << "initalizeTableLayout verif ";
  for( auto & aa : cellsDataLayout )
  {
    for( auto & bb : aa )
    {
      std::cout << " - " <<bb.m_lines[0] << " - ";
    }
    std::cout << std::endl;
  }

  updateColumnMaxLength( tableLayout, cellsDataLayout );

  calculateTableSeparators( tableLayout, cellsDataLayout, sectionSeparatingLine, topSeparator );

  cellsHeaderLayout.resize( tableLayout.getMaxHeaderRow() );

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
  std::cout << "last verif ";
  for( auto & aa : cellsData )
  {
    for( auto & bb : aa )
    {
      std::cout << " - " <<bb.m_lines[0] << " - ";
    }
    std::cout << std::endl;
  }

  outputLines( tableLayout, cellsHeader, tableOutput, tableLayout.getNbSubHeaderLines(), CellType::Header, sectionSeparatingLine );
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
                                                  RowsCellLayout & cellsDataLayout,
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
                                                RowsCellLayout & cellsDataLayout ) const
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
    if( !it->hasChild() )
    {
      // Find the maximum size between leaf and data column
      size_t columnNameLength = getMaxLineLength( currentColumn->m_columName.m_lines );
      for( size_t rowIdx = 0; rowIdx < cellsDataLayout.size(); ++rowIdx )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[rowIdx][idxColumn];
        size_t cellDataLength = getMaxLineLength( cell.m_lines );
        cell.m_maxDataLength = std::max( columnNameLength, cellDataLength );
        maxColumnSize = std::max( cell.m_maxDataLength, maxColumnSize );
      }
      // Update all data cells
      for( size_t rowIdx = 0; rowIdx< cellsDataLayout.size(); rowIdx++ )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[rowIdx][idxColumn];
        cell.m_maxDataLength = maxColumnSize;
      }

      currentColumn->setMaxStringSize( std::max(
                                         getMaxLineLength( currentColumn->m_columName.m_lines ),
                                         maxColumnSize ));
      idxColumn++;
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
                                                   RowsCellLayout & cellsDataLayout,
                                                   string & sectionSeparatingLine,
                                                   string & topSeparator ) const
{
  string const tableTitle = string( tableLayout.getTitle() );
  size_t const margins = (size_t) tableLayout.getBorderMargin() * 2;

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
  sectionlineLength += spacingBetweenColumns + margins;

  size_t maxTopLineLength =  tableTitle.length() + margins;
  maxTopLineLength = std::max( {maxTopLineLength, sectionlineLength} );

  if( sectionlineLength < maxTopLineLength )
  {
    size_t const paddingCharacters  = maxTopLineLength - sectionlineLength;
    adjustColumnWidths( tableLayout, cellsDataLayout, nbColumns, paddingCharacters );
    sectionlineLength = maxTopLineLength;
  }

  sectionSeparatingLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
  integer const topSeparatorLength = maxTopLineLength - 2;     // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::adjustColumnWidths( TableLayout & tableLayout,
                                             RowsCellLayout & cellsDataLayout,
                                             size_t const nbColumns,
                                             size_t const paddingCharacters ) const
{
  std::cout << "--adjustColumnWidths-- " << paddingCharacters << std::endl;
  size_t const paddingPerColumn  = std::floor( paddingCharacters  / nbColumns );
  size_t const remainingPadding  = paddingCharacters  - (paddingPerColumn  * nbColumns);
  std::cout << "--nbColumns-- " <<  nbColumns << "paddingPerColumn  "<< paddingPerColumn  << std::endl;
  size_t idxColumn = 0;
  for( auto it = tableLayout.beginRoot(); it != tableLayout.endRoot(); ++it )
  {
    TableLayout::Column * currentColumn = &(*it);

    size_t divider = currentColumn->hasParent() ? currentColumn->m_parent->m_subColumn.size() : nbColumns;
    size_t extraPaddingForCurrentColumn  = std::floor( paddingCharacters  / divider );
    size_t remainingPaddingForCurrentColumn = (currentColumn->m_next == nullptr) ? remainingPadding  : 0;
    std::cout << "--remainingPaddingForCurrentColumn-- " << std::endl;

    integer const newMaxStringSize = currentColumn->getMaxStringSize() +
                                     extraPaddingForCurrentColumn  +
                                     remainingPaddingForCurrentColumn;
    currentColumn->setMaxStringSize( newMaxStringSize );

    for( size_t idxRow  = 0; idxRow < cellsDataLayout.size(); idxRow++ )
    {
      TableLayout::CellLayout & cell = cellsDataLayout[idxRow ][idxColumn];
      cell.m_maxDataLength = currentColumn->getMaxStringSize();
    }
    idxColumn++;
  }
}

void TableTextFormatter::gridifyHeaders( TableLayout & tableLayout,
                                         RowsCellLayout & cellsHeaderLayout ) const
{
  size_t nbHeaderLayers = tableLayout.getMaxHeaderRow();
  std::vector< size_t > & subHeaderRowCount = tableLayout.getNbSubHeaderLines();

  subHeaderRowCount.resize( nbHeaderLayers, 1 );
  std::cout << "nbHeaderLayers "<< nbHeaderLayers << std::endl;
  for( auto it = tableLayout.beginLeaf(); it !=  tableLayout.endLeaf(); ++it )
  {
    size_t currentLayer = it.getCurrentLayer();
    TableLayout::CellLayout currentCell = it->m_columName;
    std::cout << "name : "<< currentCell.m_lines[0] << " currentLayer "<< currentLayer << std::endl;
    if( !it->hasChild()  )
    {
      if( nbHeaderLayers - 1 != currentLayer )//todo doc
      {
        for( size_t i = currentLayer; i< nbHeaderLayers - 1; i++ )
        {
          TableLayout::CellLayout emptyCell{CellType::Header, "", TableLayout::Alignment::center};//todo
          emptyCell.m_maxDataLength = it->getMaxStringSize();
          cellsHeaderLayout[i + 1].push_back( emptyCell );//todo doc
        }
      }
    }

    currentCell.m_maxDataLength = it->getMaxStringSize();
    if( it->hasParent() )
    {
      it->m_parent->m_columName.m_cellWidth++;
    }

    subHeaderRowCount[currentLayer] = std::max( subHeaderRowCount[currentLayer], currentCell.m_lines.size() );
    cellsHeaderLayout[currentLayer].push_back( currentCell );
    std::cout<< " length . "<< currentCell.m_maxDataLength  << " m_cellWidth " << currentCell.m_cellWidth << std::endl;
    for( size_t idxColumn = 0; idxColumn < currentCell.m_cellWidth; idxColumn++ )
    {
      TableLayout::CellLayout mergingCell{ CellType::Merge, "", TableLayout::Alignment::center };
      cellsHeaderLayout[currentLayer].push_back( mergingCell );
    }
  }

  //fonction resize;
  size_t idxLayer = 0;
  for( auto & lines: cellsHeaderLayout )
  {
    size_t nbLinesInLayer =  subHeaderRowCount[idxLayer];
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
    tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin());
    tableOutput << buildCell( TableLayout::Alignment::center,
                              tableTitle,
                              (topSeparator.length() - (tableLayout.getBorderMargin() *  2)));
    tableOutput << GEOS_FMT( "{:>{}}\n", m_verticalLine, tableLayout.getBorderMargin() );
  }
}

void TableTextFormatter::formatCell( TableLayout & tableLayout,
                                     std::ostringstream & tableOutput,
                                     TableLayout::CellLayout const & cell,
                                     size_t idxLine ) const
{
  tableOutput << buildCell( cell.m_alignment,
                            cell.m_lines[idxLine],
                            cell.m_maxDataLength );
  if( cell.m_cellType != CellType::Merge )
  {
    tableOutput << GEOS_FMT( "{:^{}}", m_verticalLine, tableLayout.getColumnMargin());
  }
}

void TableTextFormatter::outputLines( TableLayout & tableLayout,
                                      RowsCellLayout const & cellsLayout,
                                      std::ostringstream & tableOutput,
                                      std::vector< size_t > const & nbLinesRow,
                                      CellType sectionType,
                                      string_view sectionSeparatingLine ) const
{
  size_t idxLine = 0;
  for( auto const & line : cellsLayout )
  {
    std::cout << "size line " << line.size() << std::endl;
    for( size_t idxSubLine = 0; idxSubLine < nbLinesRow[idxLine]; idxSubLine++ )
    {
      tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin());
      for( auto const & cell : line )
      {
        if( cell.m_cellType != CellType::Hidden )
        {
          formatCell( tableLayout, tableOutput, cell, idxSubLine );
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
  if( sectionType == CellType::Value )
  {
    tableOutput << GEOS_FMT( "{}\n", sectionSeparatingLine );
  }
  std::cout << tableOutput.str() <<std::endl;
}
}
