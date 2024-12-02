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
                                               RowsCellLayout & cellsHeaderLayout,
                                               RowsCellLayout & cellsDataLayout,
                                               string & sectionSeparatingLine,
                                               string & topSeparator ) const
{

  setLinks( tableLayout.getColumns()  );
  for( auto it = tableLayout.beginLeaf(), end = tableLayout.endLeaf(); it!=end; ++it )
  {}
  RowsCellInput inputDataValues( tableData.getTableDataRows());
  if( inputDataValues.size() > 0 )
  {
    populateDataCellsLayout( tableLayout, cellsDataLayout, inputDataValues );
  }

  findAndSetLongestColumnString( tableLayout, cellsDataLayout );

  computeAndBuildTableSeparator( tableLayout, cellsDataLayout, sectionSeparatingLine, topSeparator );

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

  outputLines( tableLayout, cellsHeader, tableOutput, tableLayout.getNbSubHeaderLines(), CellType::Header, sectionSeparatingLine );
  outputLines( tableLayout, cellsData, tableOutput, tableLayout.getNbSubDataLines(), CellType::Value, sectionSeparatingLine );

  if( tableLayout.isLineBreakEnabled())
  {
    tableOutput << '\n';
  }
}

/**
 * @brief
 *
 * @param tableLayout
 */
void TableTextFormatter::setLinks( std::vector< TableLayout::Column > & columns ) const
{
  for( size_t idxColumn = 0; idxColumn < columns.size(); ++idxColumn )
  {
    if( idxColumn < columns.size() - 1 )
    {
      columns[idxColumn].m_next = &columns[idxColumn + 1];
    }

    if( !columns[idxColumn].subColumn.empty())
    {
      for( auto & subCol : columns[idxColumn].subColumn )
      {
        subCol.m_parent = &columns[idxColumn];
      }

      setLinks( columns[idxColumn].subColumn );
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
  std::vector< size_t > & nbSubDataLines = tableLayout.getNbSubDataLines();
  for( size_t idxRow = 0; idxRow < inputDataValues.size(); ++idxRow )
  {
    size_t maxLineCount = 0;
    size_t idxCol = 0;
    for( auto it = tableLayout.beginLeaf(); it != tableLayout.endLeaf(); ++it )
    {
      if( !it->hasChild())
      {
        TableData::CellData & cell = inputDataValues[idxRow][idxCol];
        TableLayout::CellAlignment const cellAlignement = it->cellAlignment;
        TableLayout::Alignment const alignement = cell.type == CellType::Header ?
                                                  cellAlignement.headerAlignment :
                                                  cellAlignement.valueAlignment;
        if( it->columnName.cellType == CellType::Hidden )
        {
          cell.type = it->columnName.cellType;
        }

        cellsDataLayout[idxRow][idxCol] = TableLayout::CellLayout( cell.type, cell.value, alignement );
        maxLineCount = std::max( maxLineCount, cellsDataLayout[idxRow][idxCol].lines.size() );

        idxCol++;

      }
    }
    nbSubDataLines.push_back( { maxLineCount } );
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

  size_t idxColumn = 0;
  for( auto it = tableLayout.beginLeaf(); it != tableLayout.endLeaf(); ++it )
  {
    size_t maxSize = 1;
    TableLayout::Column * currentColumn = &(*it);
    if( !it->hasChild() )
    {
      //1. currentColumn = max(subcolumn, dataColumn)
      for( size_t idxColumnData = 0; idxColumnData< cellsDataLayout.size(); idxColumnData++ )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[idxColumnData][idxColumn];
        cell.maxDataLength = std::max( getMaxStringLen( currentColumn->columnName.lines ),
                                       getMaxStringLen( cell.lines ));
        maxSize = std::max( cell.maxDataLength, maxSize );
      }

      for( size_t idxColumnData = 0; idxColumnData< cellsDataLayout.size(); idxColumnData++ )
      {
        TableLayout::CellLayout & cell = cellsDataLayout[idxColumnData][idxColumn];
        cell.maxDataLength = maxSize;
      }

      currentColumn->setMaxStringSize( std::max( getMaxStringLen( currentColumn->columnName.lines ),
                                                 maxSize )); //todo on le garde ?
      idxColumn++;
    }
    else      //2. and max(parent, subcolumn)
    {

      size_t subColumnsLength = 0;
      for( auto const & subColumn : currentColumn->subColumn )
      {
        subColumnsLength += subColumn.getMaxStringSize();
      }

      subColumnsLength += (currentColumn->subColumn.size() - 1) * (size_t)tableLayout.getColumnMargin();
      currentColumn->setMaxStringSize( std::max( subColumnsLength, currentColumn->getMaxStringSize() ) );
    }
  }
}

void TableTextFormatter::computeAndBuildTableSeparator( TableLayout & tableLayout,
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
    if( column.columnName.cellType != CellType::Hidden )
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
    size_t const extraCharacters = maxTopLineLength - sectionlineLength;
    increaseColumnsSize( tableLayout, cellsDataLayout, nbColumns, extraCharacters );
    sectionlineLength = maxTopLineLength;
  }

  sectionSeparatingLine = GEOS_FMT( "{:-^{}}", m_horizontalLine, sectionlineLength );
  integer const topSeparatorLength = maxTopLineLength - 2;     // Adjust for border characters
  topSeparator = GEOS_FMT( "{}{:-<{}}{}", m_horizontalLine, "", topSeparatorLength, m_horizontalLine );
}

void TableTextFormatter::increaseColumnsSize( TableLayout & tableLayout,
                                              RowsCellLayout & cellsDataLayout,
                                              size_t const nbColumns,
                                              size_t const extraCharacters ) const
{
  std::cout << "--increaseColumnsSize-- " << extraCharacters<< std::endl;
  size_t const extraCharactersPerColumn = std::floor( extraCharacters / nbColumns );
  size_t const overflowCharacters = extraCharacters - (extraCharactersPerColumn * nbColumns);
  std::cout << "--nbColumns-- " <<  nbColumns << "extraCharactersPerColumn "<< extraCharactersPerColumn << std::endl;
  size_t idxColumn = 0;
  for( auto it = tableLayout.beginRoot(); it != tableLayout.endRoot(); ++it )
  {
    TableLayout::Column * currentColumn = &(*it);

    size_t divider = currentColumn->hasParent() ? currentColumn->m_parent->subColumn.size() : nbColumns;
    size_t extraCharactersForCurrentColumn = std::floor( extraCharacters / divider );
    size_t overflowForCurrentColumn = (currentColumn->m_next == nullptr) ? overflowCharacters : 0;
    std::cout << "--overflowForCurrentColumn-- " << std::endl;

    integer const newMaxStringSize = currentColumn->getMaxStringSize() +
                                     extraCharactersForCurrentColumn +
                                     overflowForCurrentColumn;
    currentColumn->setMaxStringSize( newMaxStringSize );

    for( size_t idxColumnData = 0; idxColumnData< cellsDataLayout.size(); idxColumnData++ )
    {
      TableLayout::CellLayout & cell = cellsDataLayout[idxColumnData][idxColumn];
      cell.maxDataLength = currentColumn->getMaxStringSize();
    }
    idxColumn++;
  }
}

void TableTextFormatter::gridifyHeaders( TableLayout & tableLayout,
                                         RowsCellLayout & cellsHeaderLayout ) const
{
  size_t nbLayers = tableLayout.getMaxHeaderRow();
  std::vector< size_t > & nbRows = tableLayout.getNbSubHeaderLines();
  nbRows.resize( nbLayers, 1 );
  std::cout << "nbLayers "<< nbLayers << std::endl;
  for( auto it = tableLayout.beginLeaf(); it !=  tableLayout.endLeaf(); ++it )
  {
    size_t currLayer = it.getCurrentLayer();
    TableLayout::CellLayout currentCell = it->columnName;
    std::cout << "name :"<< currentCell.lines[0] << "currLayer "<< currLayer << std::endl;
    if( !it->hasChild()  )
    {
      if( nbLayers - 1 != currLayer )//todo doc
      {
        for( size_t i = currLayer; i< nbLayers - 1; i++ )
        {
          TableLayout::CellLayout cell{CellType::Header, "", TableLayout::Alignment::center};//todo
          cell.maxDataLength = it->getMaxStringSize();
          cellsHeaderLayout[i + 1].push_back( cell );//todo doc
        }
      }
    }


    currentCell.maxDataLength = it->getMaxStringSize();
    if( it->hasParent() )
    {
      it->m_parent->columnName.cellWidth++;
    }
    nbRows[currLayer] = std::max( nbRows[currLayer], currentCell.lines.size() );
    cellsHeaderLayout[currLayer].push_back( currentCell );
    std::cout<< " length . "<< currentCell.maxDataLength  << " cellWidth " << currentCell.cellWidth << std::endl;
    for( size_t idxColumn = 0; idxColumn < currentCell.cellWidth; idxColumn++ )
    {
      TableLayout::CellLayout mergingCell{ CellType::Merge, "", TableLayout::Alignment::center };
      cellsHeaderLayout[currLayer].push_back( mergingCell );
    }
  }

  //fonction resize;
  size_t idxLayer = 0;
  for( auto & lines: cellsHeaderLayout )
  {
    size_t nbLineLayer =  nbRows[idxLayer];
    if( nbLineLayer != 1 )
    {
      for( auto & cell : lines )
      {
        cell.lines.resize( nbLineLayer, "" );
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
  tableOutput << buildCell( cell.alignment,
                            cell.lines[idxLine],
                            cell.maxDataLength );
  if( cell.cellType != CellType::Merge )
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
    for( size_t idxSubLine = 0; idxSubLine < nbLinesRow[idxLine]; idxSubLine++ )
    {
      tableOutput << GEOS_FMT( "{:<{}}", m_verticalLine, tableLayout.getBorderMargin());
      for( auto const & cell : line )
      {
        if( cell.cellType != CellType::Hidden )
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
