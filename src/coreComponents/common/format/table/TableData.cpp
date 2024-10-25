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
 * @file TableData.cpp
 */

#include "TableData.hpp"

namespace geos
{

void TableData::addRow( std::vector< string > const & row )
{
  if( m_rows.size() != 0 && row.size() != m_rows[m_rows.size() - 1].size() )
  {
    string msg = "Remarks : some cells may be missing";
    if( std::find( m_errorsMsg.begin(), m_errorsMsg.end(), msg ) == m_errorsMsg.end())
    {
      m_errorsMsg.push_back( msg );
    }
  }
  m_rows.push_back( row );
}

void TableData::clear()
{
  m_rows.clear();
  m_errorsMsg.clear();
}

std::vector< std::vector< string > > const & TableData::getTableDataRows() const
{
  return m_rows;
}

std::vector< string > const & TableData::getErrorMsgs() const
{
  return m_errorsMsg;
}

void TableData2D::collectTableValues( arraySlice1d< real64 const > rowAxisValues,
                                      arraySlice1d< real64 const > columnAxisValues,
                                      arrayView1d< real64 const > values )
{
  integer const nX = rowAxisValues.size();
  integer const nY = columnAxisValues.size();

  for( integer i = 0; i < nX; i++ )
  {
    for( integer y = 0; y < nY; y++ )
    {
      addCell( rowAxisValues[i], columnAxisValues[y], values[ y*nX + i ] );
    }
  }
}

TableData2D::TableDataHolder TableData2D::convertTable2D( arrayView1d< real64 const > const values,
                                                          units::Unit const valueUnit,
                                                          ArrayOfArraysView< real64 const > const coordinates,
                                                          string_view rowAxisDescription,
                                                          string_view columnAxisDescription )
{
  string const rowFmt = GEOS_FMT( "{} = {{}}", rowAxisDescription );
  string const columnFmt = GEOS_FMT( "{} = {{}}", columnAxisDescription );

  collectTableValues( coordinates[0], coordinates[1], values );
  return buildTableData( string( units::getDescription( valueUnit )),
                         rowFmt,
                         columnFmt );
}

TableData2D::TableDataHolder TableData2D::buildTableData( string_view targetUnit,
                                                          string_view rowFmt,
                                                          string_view columnFmt ) const
{
  TableData2D::TableDataHolder tableData1D;
  std::vector< size_t > rowsLength;

  tableData1D.headerNames.push_back( string( targetUnit ) );

  for( auto const & columnValue : m_columnValues )
  {
    tableData1D.headerNames.push_back( GEOS_FMT( columnFmt, columnValue ) );
  }

  // insert row value and row cell values
  for( auto const & [rowValue, rowMap] : m_data )
  {
    std::vector< string > currentRowValues;
    currentRowValues.reserve( rowMap.size() );
    currentRowValues.push_back( GEOS_FMT( rowFmt, rowValue ) );

    std::set< real64 >::const_iterator columnIt = m_columnValues.begin();
    for( auto const & [columnValue, cellValue] : rowMap )
    {
      // if a column value(s) is/are missing, insert empty entry(ies)
      while( columnValue > *( columnIt++ ) && columnIt != m_columnValues.end() )
      {
        currentRowValues.push_back( "" );
      }
      currentRowValues.push_back( GEOS_FMT( "{}", cellValue ) );
    }

    tableData1D.tableData.addRow( std::move( currentRowValues ) );
    rowsLength.push_back( currentRowValues.size() );
  }

  return tableData1D;
}
}
