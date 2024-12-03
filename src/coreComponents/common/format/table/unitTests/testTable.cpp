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

// Source includes
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "dataRepository/Group.hpp"
// TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

using namespace geos;

TEST( testTable, tableEmptyRow )
{
  //table with empty row
  TableLayout const tableLayout( {"Well\nelement no.\nPV weighted\nbar",
                                  "CordX",
                                  "CoordZ",
                                  "Prev\nelement",
                                  "Next\nelement"},
                                 "InternalWellGenerator well_injector1"
                                 );

  TableData tableData;
  tableData.addRow( "value1", "[30.21543]", "3.0", 54, 0 );
  tableData.addRow( "", "", "", "", "" );
  tableData.addRow( "Duis fringilla, ligula sed porta fringilla, ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante", "[30.21543]", "30.45465142",
                    787442, 10 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString(
               tableData ),
             "\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                                                InternalWellGenerator well_injector1                                                                                 |\n"
             "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                                                                                                        Well  |       CordX  |       CoordZ  |     Prev  |     Next  |\n"
             "|                                                                                                                                 element no.  |              |               |  element  |  element  |\n"
             "|                                                                                                                                 PV weighted  |              |               |           |           |\n"
             "|                                                                                                                                         bar  |              |               |           |           |\n"
             "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                                                                                                      value1  |  [30.21543]  |          3.0  |       54  |        0  |\n"
             "|                                                                                                                                              |              |               |           |           |\n"
             "|  Duis fringilla, ligula sed porta fringilla, ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante  |  [30.21543]  |  30.45465142  |   787442  |       10  |\n"
             "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, tableClassic )
{
  TableLayout const tableLayout( {"Duis fringilla, ligula sed porta fringilla,\nligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante",
                                  "CordX",
                                  "CoordZ",
                                  "Prev\nelement",
                                  "Next\nelement"}, "InternalWellGenerator well_injector1" );

  TableData tableData;
  tableData.addRow( "value1", "[30.21543]", "3.0", 54, 0 );
  tableData.addRow( "", "", "", "", "" );
  tableData.addRow( "value23", "[30.21543]", "30.45465142", 787442, 10 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                          InternalWellGenerator well_injector1                                                           |\n"
             "-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                     Duis fringilla, ligula sed porta fringilla,  |       CordX  |       CoordZ  |     Prev  |     Next  |\n"
             "|  ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante  |              |               |  element  |  element  |\n"
             "-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                                                          value1  |  [30.21543]  |          3.0  |       54  |        0  |\n"
             "|                                                                                                  |              |               |           |           |\n"
             "|                                                                                         value23  |  [30.21543]  |  30.45465142  |   787442  |       10  |\n"
             "-----------------------------------------------------------------------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, tableColumnParamClassic )
{
  TableLayout const tableLayout( {
    TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
    TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::left},
    TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::left},
    TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
    TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::right},
    TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::right}
  }, "InternalWellGenerator well_injector1" );

  TableData tableData;
  tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n-------------------------------------------------------------------------------------------\n"
             "|                          InternalWellGenerator well_injector1                           |\n"
             "-------------------------------------------------------------------------------------------\n"
             "|  Cras egestas  |  CoordX  |  C                    |  CoordZ     |     Prev  |     Next  |\n"
             "|                |          |                       |             |  element  |  element  |\n"
             "-------------------------------------------------------------------------------------------\n"
             "|     value1     |          |  3.0                  |  3.0129877  |        2  |        1  |\n"
             "|      val1      |  v       |  [3.045,42.02,89.25]  |  3          |       10  |        3  |\n"
             "-------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, tableHiddenColumn )
{
  TableLayout const tableLayout( {
    TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
    TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
    TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::center},
    TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
    TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left, false},
    TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::center, false},
  }, "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );

  TableData tableData;
  tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n----------------------------------------------------------------------------------------------------------------\n"
             "|  Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis   |\n"
             "----------------------------------------------------------------------------------------------------------------\n"
             "|       Cras egestas       |            CoordX  |                C                |  CoordZ                    |\n"
             "----------------------------------------------------------------------------------------------------------------\n"
             "|          value1          |                    |               3.0               |  3.0129877                 |\n"
             "|           val1           |                 v  |       [3.045,42.02,89.25]       |  3                         |\n"
             "----------------------------------------------------------------------------------------------------------------\n\n" );
}

TEST( testTable, tableUniqueColumn )
{
  TableLayout const tableLayout( {
    TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
  }, "Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis" );

  TableData tableData;
  tableData.addRow( "value1" );
  tableData.addRow( "val1" );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n---------------------------------------------------------------------------------------------------------------\n"
             "|  Cras egestas ipsu a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis   |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|                                                Cras egestas                                                 |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|                                                   value1                                                    |\n"
             "|                                                    val1                                                     |\n"
             "---------------------------------------------------------------------------------------------------------------\n\n" );
}

TEST( testTable, tableEmptyTitle )
{
  TableLayout const tableLayout( {
    TableLayout::ColumnParam{{"Cras egestas"}, TableLayout::Alignment::center},
    TableLayout::ColumnParam{{"CoordX"}, TableLayout::Alignment::right},
    TableLayout::ColumnParam{{"C"}, TableLayout::Alignment::center},
    TableLayout::ColumnParam{{"CoordZ"}, TableLayout::Alignment::left},
    TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::left},
    TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::center},
  }
                                 );

  TableData tableData;
  tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n-------------------------------------------------------------------------------------------\n"
             "|  Cras egestas  |  CoordX  |           C           |  CoordZ     |  Prev     |   Next    |\n"
             "|                |          |                       |             |  element  |  element  |\n"
             "-------------------------------------------------------------------------------------------\n"
             "|     value1     |          |          3.0          |  3.0129877  |  2        |     1     |\n"
             "|      val1      |       v  |  [3.045,42.02,89.25]  |  3          |  10       |     3     |\n"
             "-------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, table2DTable )
{
  //collect
  TableData2D tableData;

  for( real64 p = 10000; p<20000; p+=5000 )
  {
    for( real64 t = 400; t>=270; t+=-50.0 )
    {
      real64 value = t/p;
      tableData.addCell( t, p, value );
    }
  }

  //convert
  string const rowFmt = GEOS_FMT( "{} = {{}}", "Temperature" );
  string const columnFmt = GEOS_FMT( "{} = {{}}", "Pression" );
  TableData2D::TableDataHolder tableconverted = tableData.buildTableData( "Viscosity kg*s",
                                                                          rowFmt,
                                                                          columnFmt );

  //format
  TableLayout const tableLayout( tableconverted.headerNames );

  //log
  TableTextFormatter const tableLog( tableLayout );
  EXPECT_EQ( tableLog.toString( tableconverted.tableData ),
             "\n---------------------------------------------------------------------\n"
             "|     Viscosity kg*s  |  Pression = 10000  |      Pression = 15000  |\n"
             "---------------------------------------------------------------------\n"
             "|  Temperature = 300  |              0.03  |                  0.02  |\n"
             "|  Temperature = 350  |             0.035  |  0.023333333333333334  |\n"
             "|  Temperature = 400  |              0.04  |   0.02666666666666667  |\n"
             "---------------------------------------------------------------------\n\n"
             );
}


TEST( testTable, table2DColumnMismatch )
{
  //test 2D table column  mismatch
  {
    // collect
    TableData2D tableData;

    tableData.addCell( 300, 10000, 0.03 );
    tableData.addCell( 300, 15000, 0.02 );
    tableData.addCell( 350, 10000, 0.035 );
    tableData.addCell( 350, 10000, 0.035 );
    tableData.addCell( 400, 10000, 0.04 );
    tableData.addCell( 400, 15000, 0.02666666666666667 );

    //convert
    string const rowFmt = GEOS_FMT( "{} = {{}}", "Temperature" );
    string const columnFmt = GEOS_FMT( "{} = {{}}", "Pression" );
    TableData2D::TableDataHolder tableConverted = tableData.buildTableData( "Viscosity kg*s",
                                                                            rowFmt,
                                                                            columnFmt );

    //format
    TableLayout const tableLayout( tableConverted.headerNames );

    //log
    TableTextFormatter const tableLog( tableLayout );
    EXPECT_EQ( tableLog.toString( tableConverted.tableData ),
               "\n--------------------------------------------------------------------\n"
               "|  Remarks : some cells may be missing                             |\n"
               "--------------------------------------------------------------------\n"
               "|     Viscosity kg*s  |  Pression = 10000  |     Pression = 15000  |\n"
               "--------------------------------------------------------------------\n"
               "|  Temperature = 300  |              0.03  |                 0.02  |\n"
               "|  Temperature = 350  |             0.035  |                       |\n"
               "|  Temperature = 400  |              0.04  |  0.02666666666666667  |\n"
               "--------------------------------------------------------------------\n\n"
               );
  }
}

TEST( testTable, layoutTable )
{
  string filename = "fluid1_phaseModel1_PhillipsBrineDensity_table";

  string log = GEOS_FMT( "The {} PVT table exceeding 500 rows.\nTo visualize the tables, go to the generated csv \n", filename );
  TableLayout const tableLayoutInfos( {TableLayout::ColumnParam{{log}, TableLayout::Alignment::left}}, filename );

  TableTextFormatter const tableLog( tableLayoutInfos );
  EXPECT_EQ( tableLog.layoutToString(),
             "\n-------------------------------------------------------------------------------------\n"
             "|                   fluid1_phaseModel1_PhillipsBrineDensity_table                   |\n"
             "-------------------------------------------------------------------------------------\n"
             "|  The fluid1_phaseModel1_PhillipsBrineDensity_table PVT table exceeding 500 rows.  |\n"
             "|  To visualize the tables, go to the generated csv                                 |\n"
             "-------------------------------------------------------------------------------------\n"
             );
}

TEST( testTable, tableSetMargin )
{
  ////////////
  //////// If setMargin used elsewhere make it public and remove comments for this test
  ////////////
  //test with tiny margin
  // {
  //   TableLayout tableLayout( {
  //     TableLayout::ColumnParam{{"Colonne 1"}, TableLayout::Alignment::center},
  //     TableLayout::ColumnParam{{"Colonne 2"}, TableLayout::Alignment::center},
  //     TableLayout::ColumnParam{{"Colonne 3"}, TableLayout::Alignment::right},
  //     TableLayout::ColumnParam{{"Colonne 4"}, TableLayout::Alignment::right},
  //     TableLayout::ColumnParam{{"Prev\nelement"}, TableLayout::Alignment::right},
  //     TableLayout::ColumnParam{{"Next\nelement"}, TableLayout::Alignment::right},
  //   }, "InternalWellGenerator well_injector1" );

  //   //tableLayout.setMargin( TableLayout::MarginValue::tiny );

  //   TableData tableData;
  //   tableData.addRow( "value 1", "long value 1", "3.0034", 3.0129877, "" , 1 );
  //   tableData.addRow(  "value 1", "long value 2", "100.45", 4.0129877, 1 , 2 );

  //   TableTextFormatter const tableText( tableLayout );
  //   EXPECT_EQ( tableText.toString( tableData ),
// "\n------------------------------------------------------------\n"
// "|           InternalWellGenerator well_injector1           |\n"
// "------------------------------------------------------------\n"
// "|Colonne 1| Colonne 2  |Colonne 3|Colonne 4|   Prev|   Next|\n"
// "|         |            |         |         |element|element|\n"
// "------------------------------------------------------------\n"
// "| value 1 |long value 1|   3.0034|3.0129877|       |      1|\n"
// "| value 1 |long value 2|   100.45|4.0129877|      1|      2|\n"
// "------------------------------------------------------------\n\n"
//              );
// }
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
