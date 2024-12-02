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
   * @brief Prepare all the columns  with the appropriate values and formatting separator
   * @param columns  The vector containg all columns . Each tableColumnData contains its own
   *        parameters (such as name, alignment, etc.).
   * @param tableData Vector containing all rows filled with values
   * @param sectionSeparatingLine Separator string used between sections of the table
   * @param topSeparator The table top separator
   */
  void initalizeTableLayout( TableLayout & tableLayout,
                             TableData const & tableData,
                             RowsCellLayout & cellsDataLayout,
                             RowsCellLayout & cellsHeaderLayout,
                             string & sectionSeparatingLine,
                             string & topSeparator ) const;
/**
 * @brief Displays the complete table
 * @param tableOutput The output stream
 * @param columns  Vector containg all columns
 * @param sectionSeparatingLine Separator string used between sections of the table
 * @param topSeparator The table top separator
 */
  void outputTable( TableLayout & tableLayout,
                    std::ostringstream & tableOutput,
                    RowsCellLayout const & cellsHeader,
                    RowsCellLayout const & cellsData,
                    string_view sectionSeparatingLine,
                    string_view topSeparator ) const;

  void setLinks( std::vector< TableLayout::Column > & columns ) const;

  /**
   * @brief
   */
  void gridifyHeaders( TableLayout & tableLayout,
                       RowsCellLayout & cellsDataLayout ) const;

  /**
   * @brief Populate all the tableColumnData values with values extracted from TableData.
   * @param columns  Vector of columns  to populate.
   * @param tableData Vector containing all rows filled with values
   * @param isSubColumn Boolean indicating if the current tableColumnData is a subcolumn
   */
  void populateDataCellsLayout( TableLayout & tableLayout,
                                RowsCellLayout & cellsDataLayout,
                                RowsCellInput & inputDataValues ) const;

  /**
   * @brief For each tableColumnData find and set the column's longest string
   * @param tableColumnData The tableColumnData to process.
   * @note Compares the longest string from the header with the longest string from the column values.
   * If the column contains subcolumns, it recursively applies the same logic to them
   */
  void findAndSetLongestColumnString( TableLayout & tableLayout,
                                      RowsCellLayout & cellsDataLayout ) const;

  /**
   * @brief Compute the max table line length, taking into account the length of : title, columns  header/values
   * Increase the size of the columns if necessary, then builds the table's separating lines
   * @param columns Vector of tableColumnData containing containing the largest string for each tableColumnData
   * @param sectionSeparatingLine Separator string used between sections of the table
   * @param topSeparator The table top separator
   */
  void computeAndBuildTableSeparator( TableLayout & tableLayout,
                                      RowsCellLayout & cellsDataLayout,
                                      string & sectionSeparatingLine,
                                      string & topSeparator ) const;

  /**
   * @brief Increase each tableColumnData size if the title is larger than all the columns
   * @param columns Vector containing all table columns
   * @param extraCharacters ExtraCharacters to be distributed between each columns
   */
  void increaseColumnsSize( TableLayout & tableLayout,
                            RowsCellLayout & cellsHeaderLayout,
                            size_t const nbColumns,
                            size_t const extraCharacters ) const;

  /**
   * @brief Output the title row in the table
   * @param topSeparator The top separator string
   */
  void outputTitleRow( TableLayout & tableLayout,
                       std::ostringstream & tableOutput,
                       string_view topSeparator ) const;

  void formatCell( TableLayout & tableLayout,
                   std::ostringstream & tableOutput,
                   TableLayout::CellLayout const & cell,
                   size_t idxLine ) const;

  /**
   * @brief Output the values rows in the table
   * @param columns  Vector containing all table columns
   * @param tableOutput The output stream
   * @param sectionSeparatingLine Separator string used between sections of the table
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
