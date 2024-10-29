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
  TableLayout const tableLayout( "InternalWellGenerator well_injector1",
                                 {"Well\nelement no.\nPV weighted\nbar",
                                  "CordX",
                                  "CoordZ",
                                  "Prev\nelement",
                                  "Next\nelement"} );

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
             "|                                                                     Well                                                                     |    CordX     |    CoordZ     |   Prev    |   Next    |\n"
             "|                                                                 element no.                                                                  |              |               |  element  |  element  |\n"
             "|                                                                 PV weighted                                                                  |              |               |           |           |\n"
             "|                                                                     bar                                                                      |              |               |           |           |\n"
             "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                                                                                                      value1  |  [30.21543]  |          3.0  |       54  |        0  |\n"
             "|                                                                                                                                              |              |               |           |           |\n"
             "|  Duis fringilla, ligula sed porta fringilla, ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante  |  [30.21543]  |  30.45465142  |   787442  |       10  |\n"
             "-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, tableClassic )
{
  TableLayout const tableLayout( "InternalWellGenerator well_injector1",
                                 {"Duis fringilla, ligula sed porta fringilla,\nligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante",
                                  "CordX",
                                  "CoordZ",
                                  "Prev\nelement",
                                  "Next\nelement"} );

  TableData tableData;
  tableData.addRow( "value1", "[30.21543]", "3.0", 54, 0 );
  tableData.addRow( "", "", "", "", "" );
  tableData.addRow( "value23", "[30.21543]", "30.45465142", 787442, 10 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                          InternalWellGenerator well_injector1                                                           |\n"
             "-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                           Duis fringilla, ligula sed porta fringilla,                            |    CordX     |    CoordZ     |   Prev    |   Next    |\n"
             "|  ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante  |              |               |  element  |  element  |\n"
             "-----------------------------------------------------------------------------------------------------------------------------------------------------------\n"
             "|                                                                                          value1  |  [30.21543]  |          3.0  |       54  |        0  |\n"
             "|                                                                                                  |              |               |           |           |\n"
             "|                                                                                         value23  |  [30.21543]  |  30.45465142  |   787442  |       10  |\n"
             "-----------------------------------------------------------------------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, tableColumnParamClassic ) //TODO
{
  TableLayout tableLayout( {TableLayout::Column{{"Cras egestas"}, TableLayout::Alignment::center},
                            TableLayout::Column{{"CoordX"}, TableLayout::Alignment::left},
                            TableLayout::Column{{"C"}, TableLayout::Alignment::left},
                            TableLayout::Column{{"CoordZ"}, TableLayout::Alignment::left},
                            TableLayout::Column{{"Prev\nelement"}, TableLayout::Alignment::right},
                            TableLayout::Column{{"Next\nelement"}, TableLayout::Alignment::right}} );
  tableLayout.setValuesAlignment( TableLayout::Alignment::right );

  TableData tableData;
  tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  std::cout <<  tableText.toString( tableData ) <<std::endl;
  EXPECT_EQ( tableText.toString( tableData ),
             "\n-------------------------------------------------------------------------------------------\n"
             "|  Cras egestas  |  CoordX  |  C                    |  CoordZ     |     Prev  |     Next  |\n"
             "|                |          |                       |             |  element  |  element  |\n"
             "-------------------------------------------------------------------------------------------\n"
             "|        value1  |          |                  3.0  |  3.0129877  |        2  |        1  |\n"
             "|          val1  |       v  |  [3.045,42.02,89.25]  |          3  |       10  |        3  |\n"
             "-------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, tableHiddenColumn ) // TODO
{
  string const title = "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis";
  TableLayout tableLayout( title,
                           {TableLayout::Column{{"Cras egestas"}, TableLayout::Alignment::center},
                            TableLayout::Column{{"CoordX"}, TableLayout::Alignment::right},
                            TableLayout::Column{{"C"}, TableLayout::Alignment::center},
                            TableLayout::Column{{"CoordZ"}, TableLayout::Alignment::left},
                            TableLayout::Column{{"Prev\nelement"}, TableLayout::Alignment::left, false},
                            TableLayout::Column{{"Next\nelement"}, TableLayout::Alignment::center, false}} );
  tableLayout.setValuesAlignment( TableLayout::Alignment::left );

  TableData tableData;
  tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n---------------------------------------------------------------------------------------------------------------\n"
             "|  Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis  |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|       Cras egestas        |             CoordX  |                C                 |  CoordZ                |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|  value1                   |                     |  3.0                             |  3.0129877             |\n"
             "|  val1                     |  v                  |  [3.045,42.02,89.25]             |  3                     |\n"
             "---------------------------------------------------------------------------------------------------------------\n\n" );
}

TEST( testTable, tableUniqueColumn )
{
  string const title = "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis";
  TableLayout const tableLayout( title, {TableLayout::Column{{"Cras egestas"}, TableLayout::Alignment::center}} );

  TableData tableData;
  tableData.addRow( "value1" );
  tableData.addRow( "val1" );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n---------------------------------------------------------------------------------------------------------------\n"
             "|  Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis  |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|                                                Cras egestas                                                 |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|                                                                                                     value1  |\n"
             "|                                                                                                       val1  |\n"
             "---------------------------------------------------------------------------------------------------------------\n\n" );
}

TEST( testTable, tableEmptyTitle )
{
  TableLayout const tableLayout( {TableLayout::Column{{"Cras egestas"}, TableLayout::Alignment::center},
                                  TableLayout::Column{{"CoordX"}, TableLayout::Alignment::right},
                                  TableLayout::Column{{"C"}, TableLayout::Alignment::center},
                                  TableLayout::Column{{"CoordZ"}, TableLayout::Alignment::left},
                                  TableLayout::Column{{"Prev\nelement"}, TableLayout::Alignment::left},
                                  TableLayout::Column{{"Next\nelement"}, TableLayout::Alignment::center}} );

  TableData tableData;
  tableData.addRow( "value1", " ", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n-------------------------------------------------------------------------------------------\n"
             "|  Cras egestas  |  CoordX  |           C           |  CoordZ     |  Prev     |   Next    |\n"
             "|                |          |                       |             |  element  |  element  |\n"
             "-------------------------------------------------------------------------------------------\n"
             "|        value1  |          |                  3.0  |  3.0129877  |        2  |        1  |\n"
             "|          val1  |       v  |  [3.045,42.02,89.25]  |          3  |       10  |        3  |\n"
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
  TableLayout const tableLayout( "", tableconverted.headerNames );

  //log
  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableconverted.tableData ),
             "\n---------------------------------------------------------------------\n"
             "|   Viscosity kg*s    |  Pression = 10000  |    Pression = 15000    |\n"
             "---------------------------------------------------------------------\n"
             "|  Temperature = 300  |              0.03  |                  0.02  |\n"
             "|  Temperature = 350  |             0.035  |  0.023333333333333334  |\n"
             "|  Temperature = 400  |              0.04  |   0.02666666666666667  |\n"
             "---------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, layoutTable )
{
  string filename = "fluid1_phaseModel1_PhillipsBrineDensity_table";
  string log = GEOS_FMT( "The {} PVT table exceeding 500 rows.\nTo visualize the tables, go to the generated csv \n", filename );
  TableLayout const tableLayoutInfos( filename, {TableLayout::Column{{log}, TableLayout::Alignment::left}} );

  TableTextFormatter const tableText( tableLayoutInfos );
  EXPECT_EQ( tableText.toString(),
             "\n-------------------------------------------------------------------------------------\n"
             "|                   fluid1_phaseModel1_PhillipsBrineDensity_table                   |\n"
             "-------------------------------------------------------------------------------------\n"
             "|  The fluid1_phaseModel1_PhillipsBrineDensity_table PVT table exceeding 500 rows.  |\n"
             "|  To visualize the tables, go to the generated csv                                 |\n"
             "-------------------------------------------------------------------------------------\n"
             );
}

TEST( testTable, subColumns )
{
  {
    TableLayout const tableLayout( {" ",
                                    "Column1",
                                    TableLayout::Column{"Nodes ", TableLayout::Alignment::right, true, {"Locales", "Ghost", "Active"}},
                                    "Column3",
                                    TableLayout::Column{"Column4 ", TableLayout::Alignment::right, true, {"Locales", "Ghost"}},
                                    "Column5"} );

    TableData tableData;
    tableData.addRow( "min", "125", "375,0001", " YES", 2354654, 562, 43.0, 43.0, 562, 5 );
    tableData.addRow( "max", "360", "390,1", " YES", 383213213, 712, 48.0, 47.0, 72, 2 );

    TableTextFormatter tableText( tableLayout );

    EXPECT_EQ( tableText.toString( tableData ),
               "\n--------------------------------------------------------------------------------------------------------\n"
               "|       |  Column1  |                            Nodes   |  Column3  |           Column4   |  Column5  |\n"
               "--------------------------------------------------------------------------------------------------------\n"
               "|       |           |   Locales  |  Ghost  |     Active  |           |  Locales  |  Ghost  |           |\n"
               "--------------------------------------------------------------------------------------------------------\n"
               "|  min  |      125  |  375,0001  |    YES  |    2354654  |      562  |       43  |     43  |      562  |\n"
               "|  max  |      360  |     390,1  |    YES  |  383213213  |      712  |       48  |     47  |       72  |\n"
               "--------------------------------------------------------------------------------------------------------\n\n"
               );
  }
}

TEST( testTable, variadicTest )
{
  {
    TableLayout const layoutTest( "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis nascetur ridiculus mus",
                                  { "Rank",
                                    TableLayout::Column{"Nodes", TableLayout::Alignment::center, true, {"local", "ghost"}},
                                    "Edge",
                                    TableLayout::Column{"Faces", TableLayout::Alignment::center, true, {"local", "ghost"}},
                                    TableLayout::Column{"Elems", TableLayout::Alignment::center, true, {"local", "ghost"}} } );

    TableData tableData;
    tableData.addRow( "min(local/total)", 1, 2, 3, 4, 5, 6, 7 );
    tableData.addRow( "min(local/total)", 1, 2, 3, 4, 5, 6, 7 );
    TableTextFormatter log( layoutTest );
    EXPECT_EQ( log.toString( tableData ),
               "\n--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|  Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis nascetur ridiculus mus  |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|             Rank              |           Nodes           |      Edge      |           Faces           |           Elems           |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|                               |    local    |    ghost    |                |    local    |    ghost    |    local    |    ghost    |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|             min(local/total)  |          1  |          2  |             3  |          4  |          5  |          6  |          7  |\n"
               "|             min(local/total)  |          1  |          2  |             3  |          4  |          5  |          6  |          7  |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n\n"
               );
  }
}
TEST( testTable, testLineWrap )
{
  TableLayout tableLayout( {"Cras egestas", "CoordX", "C", "CoordZ", "Prev\nelement", "Next\nelement"} );
  tableLayout.setTitle( "title" ).setMargin( TableLayout::MarginValue::tiny ).disableLineWrap();

  TableData tableData;
  tableData.addRow( "1", "2", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "1", "2", "3.0", 3.0129877, 2.0f, 1 );

  TableTextFormatter const tableText( tableLayout );

  EXPECT_EQ( tableText.toString( tableData ),
             "---------------------------------------------------\n"
             "|                      title                      |\n"
             "---------------------------------------------------\n"
             "|Cras egestas|CoordX| C | CoordZ  | Prev  | Next  |\n"
             "|            |      |   |         |element|element|\n"
             "---------------------------------------------------\n"
             "|           1|     2|3.0|3.0129877|      2|      1|\n"
             "|           1|     2|3.0|3.0129877|      2|      1|\n"
             "---------------------------------------------------\n"
             );
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
