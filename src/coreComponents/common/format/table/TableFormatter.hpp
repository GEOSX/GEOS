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
 * @file TableFormatter.hpp
 */

#ifndef GEOS_COMMON_FORMAT_TABLE_TABLEFORMATTER_HPP
#define GEOS_COMMON_FORMAT_TABLE_TABLEFORMATTER_HPP

#include "TableData.hpp"
#include "TableLayout.hpp"
#include "TableTypes.hpp"

namespace geos
{

/**
 * @brief abstract class for formatting table data
 */
class TableFormatter
{

public:
  using RowsCellInput = std::vector< std::vector< TableData::CellData > >;
  using RowsCellLayout = std::vector< std::vector< TableLayout::CellLayout > >;


protected:

  /// Layout for a table
  TableLayout m_tableLayout;

  TableFormatter() = default;

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableFormatter( TableLayout const & tableLayout );

  /**
   * @brief Destroy the Table Formatter object
   */
  virtual ~TableFormatter() = default;
};

/**
 * @brief class for CSV formatting
 */
class TableCSVFormatter : public TableFormatter
{
public:

  /**
   * @brief Construct a new Table Formatter
   */
  TableCSVFormatter():
    TableFormatter( TableLayout() )
  {}

  /**
   * @brief Construct a new Table Formatter from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableCSVFormatter( TableLayout const & tableLayout );

  /**
   * @brief Destroy the TableCSVFormatter object
   */
  virtual ~TableCSVFormatter() = default;

  /**
   * @return The string with all tableColumnData names.
   */
  string headerToString() const;

  /**
   * @brief Convert the table data to a CSV string..
   * @return The CSV string representation of the table data.
   */
  string dataToString( TableData const & tableData ) const;

  /**
   * @brief Convert a data source to a CSV string.
   * @tparam DATASOURCE The source to convert
   * @param tableData The data source to convert
   * @return The CSV string representation of a data source.
   */
  template< typename DATASOURCE >
  string toString( DATASOURCE const & tableData ) const;

};

/**
 * @brief Convert the TableData to a CSV string.
 * @param tableData The TableData to convert.
 * @return The CSV string representation of the TableData.
 */
template<>
string TableCSVFormatter::toString< TableData >( TableData const & tableData ) const;


/**
 * @brief class for log formatting
 */
class TableTextFormatter : public TableFormatter
{
public:

  /**
   * @brief Construct a new TableFormatter
   */
  TableTextFormatter():
    TableFormatter( TableLayout() )
  {}

  /**
   * @brief Construct a new TableFormatter from a tableLayout
   * @param tableLayout Contain all tableColumnData names and optionnaly the table title
   */
  TableTextFormatter( TableLayout const & tableLayout );


  /**
   * @brief Destroy the Table Text Formatter object
   */
  virtual ~TableTextFormatter() = default;

  /**
   * @return A TableLayout string representation,
   * The TableTextFormatter receives hasn't receive any data, so only the header part is returned.
   */
  string toString() const;

  /**
   * @brief Convert a data source to a table string.
   * @param tableData The data source to convert.
   * @return The table string representation of the TableData.
   */
  template< typename DATASOURCE >
  string toString( DATASOURCE const & tableData ) const;

private:

  /// symbol for separator construction
  static constexpr char m_verticalLine = '|';
  ///  for the extremity of a row
  static constexpr char m_horizontalLine = '-';


/**
 * @brief Initializes the table layout with the given table data and prepares necessary layouts for headers and data cells.
 * @param tableLayout A reference to the `TableLayout` object.
 * @param tableData A constant reference to the `TableData` object, which contains the actual data for the table.
 * @param cellsHeaderLayout A reference to a `RowsCellLayout` where the header cells will be populated.
 * @param cellsDataLayout A reference to a `RowsCellLayout` where the data cells will be populated.
 * @param sectionSeparatingLine A string that will be used as the separator line between sections in the table.
 * @param topSeparator A string that will be used as the separator at the top of the table.
 */
  void initalizeTableLayout( TableLayout & tableLayout,
                             TableData const & tableData,
                             RowsCellLayout & cellsDataLayout,
                             RowsCellLayout & cellsHeaderLayout,
                             string & sectionSeparatingLine,
                             string & topSeparator ) const;
/**
 * @brief Outputs the formatted table to the provided output stream.
 * @param tableLayout The layout of the table
 * @param tableOutput A reference to an `std::ostringstream` where the formatted table will be written.
 * @param cellsHeader The layout of the header rows
 * @param cellsData The layout of the data rows
 * @param sectionSeparatingLine The string to be used as a separator between different sections of the table.
 * @param topSeparator The string to be used as the separator at the top of the table.
 */
  void outputTable( TableLayout & tableLayout,
                    std::ostringstream & tableOutput,
                    RowsCellLayout const & cellsHeader,
                    RowsCellLayout const & cellsData,
                    string_view sectionSeparatingLine,
                    string_view topSeparator ) const;

  /**
   * @brief Sets parent-child relationships between columns and sub-columns.
   * @param columns A reference to a vector of `TableLayout::Column` objects.
   */
  void setLinks( std::vector< TableLayout::Column > & columns ) const;

  /**
   * @brief Adjusts the header layout by ensuring all header layers have consistent row sizes and formats.
   * @param tableLayout The layout of the table, containing information about columns, headers, and their layers.
   * @param cellsHeaderLayout A reference to the collection of header cells that will be updated with the gridified layout.
   */
  void gridifyHeaders( TableLayout & tableLayout,
                       RowsCellLayout & cellsDataLayout ) const;

/**
 * @brief Populates the data cells layout based on input data values.
 * @param tableLayout The layout of the table,
 * @param cellsDataLayout A reference to the layout for the data cells that will be populated.
 * @param inputDataValues A 2D vector containing the actual input data values.
 */
  void populateDataCellsLayout( TableLayout & tableLayout,
                                RowsCellLayout & cellsDataLayout,
                                RowsCellInput & inputDataValues ) const;

  /**
   * @brief Finds and sets the longest string for each column in the table.
   * @param tableLayout The layout of the table,
   * @param cellsDataLayout A reference to the collection of data cells. The function updates the maximum string
   *        length for each cell based on the longest string found in the column.
   */
  void updateColumnMaxLength( TableLayout & tableLayout,
                              RowsCellLayout & cellsDataLayout ) const;

  /**
   * @brief Computes and constructs the separator lines for the table.
   * @param tableLayout The layout of the table,
   * @param cellsDataLayout A reference to the collection of data cells that will be affected by column resizing.
   * @param sectionSeparatingLine A string reference where the separator line between sections will be stored.
   * @param topSeparator A string reference where the top separator line will be stored.
   */
  void calculateTableSeparators( TableLayout & tableLayout,
                                 RowsCellLayout & cellsDataLayout,
                                 string & sectionSeparatingLine,
                                 string & topSeparator ) const;

  /**
   * @brief Increases the size of columns to accommodate extra characters.
   * @param tableLayout The layout of the table
   * @param cellsDataLayout A reference to the collection of data cells, which will be updated with the new column sizes.
   * @param nbColumns The total number of columns in the table.
   * @param extraCharacters The total number of extra characters to be distributed across the columns.
   */
  void adjustColumnWidths( TableLayout & tableLayout,
                            RowsCellLayout & cellsHeaderLayout,
                            size_t const nbColumns,
                            size_t const extraCharacters ) const;

  /**
   * @brief Output the title row in the table
   * @param tableLayout The layout of the table
   * @param tableOutput The output stream
   * @param topSeparator The top separator string
   */
  void outputTitleRow( TableLayout & tableLayout,
                       std::ostringstream & tableOutput,
                       string_view topSeparator ) const;

  /**
   * @brief Formats a table cell and appends it to the table output.
   * @param tableLayout The layout of the table
   * @param tableOutput The output stream
   * @param cell The cell to format
   * @param idxLine The current line index used to access the specific content for the cell.
   */
  void formatCell( TableLayout & tableLayout,
                   std::ostringstream & tableOutput,
                   TableLayout::CellLayout const & cell,
                   size_t idxLine ) const;

  /**
   * @brief Outputs the formatted table lines to the output stream.
   * @param tableLayout The layout of the table
   * @param cellsLayout A collection of rows, each containing a layout of cells to be processed and formatted.
   * @param tableOutput The output stream
   * @param nbLinesRow A vector containing the number of sub-lines for each row.
   * @param sectionType The type of the section being processed (Header, Value, etc.).
   * @param sectionSeparatingLine A separator line to be added.
   */
  void outputLines( TableLayout & tableLayout,
                    RowsCellLayout const & cellsLayout,
                    std::ostringstream & tableOutput,
                    std::vector< size_t > const & nbLinesRow,
                    CellType sectionType,
                    string_view sectionSeparatingLine ) const;
};

/**
 * @brief Convert a TableData to a table string.
 * @param tableData The TableData to convert.
 * @return The table string representation of the TableData.
 */
template<>
string TableTextFormatter::toString< TableData >( TableData const & tableData ) const;
}

#endif /* GEOS_COMMON_FORMAT_TABLE_TABLEFORMATTER_HPP */
