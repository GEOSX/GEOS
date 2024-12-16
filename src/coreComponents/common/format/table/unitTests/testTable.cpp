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
  tableData.addRow( "", " ", "", "", "" );
  tableData.addRow( "Duis fringilla, ligula sed porta fringilla, ligula wisi commodo felis,ut adipiscing felis dui in enim. Suspendisse malesuada ultrices ante", "[30.21543]", "30.45465142",
                    787442, 10 );
  TableTextFormatter const tableText( tableLayout );
  tableText.toString( tableData );
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

TEST( testTable, tableColumnParamClassic )
{
  TableLayout const tableLayout( {
    TableLayout::Column()
      .setName( "Cras egestas" ),
    TableLayout::Column()
      .setName( "CoordX" ),
    TableLayout::Column()
      .setName( "C" ),
    TableLayout::Column()
      .setName( "CoordZ" ),
    TableLayout::Column()
      .setName( "Prev\nelement" ),
    TableLayout::Column()
      .setName( "Next\nelement" )} );

  TableData tableData;
  tableData.addRow( "value1", "gaz\nwater", "3.0\n42.0", "3.0129877\n0.0123456", "2\n3", "1\n4" );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  std::cout << tableText.toString( tableData )<< std::endl;
  EXPECT_EQ( tableText.toString( tableData ),
             "\n-------------------------------------------------------------------------------------------\n"
             "|  Cras egestas  |  CoordX  |  C                    |  CoordZ     |     Prev  |     Next  |\n"
             "|                |          |                       |             |  element  |  element  |\n"
             "-------------------------------------------------------------------------------------------\n"
             "|        value1  |     gaz  |                  3.0  |  3.0129877  |        2  |        1  |\n"
             "|                |   water  |                 42.0  |  0.0123456  |        3  |        4  |\n"
             "|          val1  |       v  |  [3.045,42.02,89.25]  |          3  |       10  |        3  |\n"
             "-------------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, tableHiddenColumn )
{
  string const title = "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis";
  TableLayout tableLayout( title,
  {
    TableLayout::Column()
      .setName( "Cras egestas" )
      .setValuesAlignment( TableLayout::Alignment::left )
      .setHeaderAlignment( TableLayout::Alignment::center ),
    TableLayout::Column()
      .setName( "CoordX" )
      .setValuesAlignment( TableLayout::Alignment::left )
      .setHeaderAlignment( TableLayout::Alignment::right )
      .setVisibility( CellType::Hidden ),
    TableLayout::Column()
      .setName( "C" )
      .setValuesAlignment( TableLayout::Alignment::left )
      .setHeaderAlignment( TableLayout::Alignment::center ),
    TableLayout::Column()
      .setName( "CoordZ" )
      .setValuesAlignment( TableLayout::Alignment::left )
      .setHeaderAlignment( TableLayout::Alignment::left ),
    TableLayout::Column()
      .setName( "Prev\nelement" )
      .setVisibility( CellType::Hidden ),
    TableLayout::Column()
      .setName( "Next\nelement" )
      .setVisibility( CellType::Hidden )
  } );

  TableData tableData;
  tableData.addRow( "value1", "", "3.0", 3.0129877, 2.0f, 1 );
  tableData.addRow( "val1", "v", "[3.045,42.02,89.25]", 3.0, 10.0f, 3 );

  TableTextFormatter const tableText( tableLayout );
  EXPECT_EQ( tableText.toString( tableData ),
             "\n---------------------------------------------------------------------------------------------------------------\n"
             "|  Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis  |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|           Cras egestas           |                    C                    |  CoordZ                        |\n"
             "|                                  |                                         |                                |\n"
             "---------------------------------------------------------------------------------------------------------------\n"
             "|  value1                          |  3.0                                    |  3.0129877                     |\n"
             "|  val1                            |  [3.045,42.02,89.25]                    |  3                             |\n"
             "---------------------------------------------------------------------------------------------------------------\n\n" );
}

TEST( testTable, tableUniqueColumn )
{
  string const title = "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis";
  TableLayout const tableLayout( title, {"Cras egestas"} );

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
  TableLayout const tableLayout( {
    TableLayout::Column()
      .setName( "Cras egestas" )
      .setHeaderAlignment( TableLayout::Alignment::center ),
    TableLayout::Column()
      .setName( "CoordX" )
      .setHeaderAlignment( TableLayout::Alignment::right ),
    "C",
    TableLayout::Column()
      .setName( "CoordZ" )
      .setHeaderAlignment( TableLayout::Alignment::left ),
    TableLayout::Column()
      .setName( "Prev\nelement" )
      .setHeaderAlignment( TableLayout::Alignment::left ),
    TableLayout::Column()
      .setName( "Next\nelement" )
      .setHeaderAlignment( TableLayout::Alignment::center )} );

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
  string log = GEOS_FMT( "The {} PVT table exceeding 500 rows.\nTo visualize the tables, go to the generated csv", filename );
  TableLayout const tableLayoutInfos( filename,
  {
    TableLayout::Column()
      .setName( log )
      .setHeaderAlignment( TableLayout::Alignment::left )} );
  TableTextFormatter const tableText( tableLayoutInfos );
  EXPECT_EQ( tableText.toString(),
             "\n-------------------------------------------------------------------------------------\n"
             "|                   fluid1_phaseModel1_PhillipsBrineDensity_table                   |\n"
             "-------------------------------------------------------------------------------------\n"
             "|  The fluid1_phaseModel1_PhillipsBrineDensity_table PVT table exceeding 500 rows.  |\n"
             "|  To visualize the tables, go to the generated csv                                 |\n"
             "-------------------------------------------------------------------------------------\n\n"
             );
}

TEST( testTable, subColumns )
{
  {
    TableLayout const tableLayout( {
      " ",
      "Column1",
      TableLayout::Column()
        .setName( "Nodes" )
        .addSubColumns( {"LocalesNodes", "GhostNodes", "ActiveNodes" } ),
      "Column3",
      TableLayout::Column()
        .setName( "Column4" )
        .addSubColumns( { TableLayout::Column()
                            .setName( "Locales" ).addSubColumns( {"SubLocales1", "SubLocales2"} ),
                          TableLayout::Column()
                            .setName( "Ghost" ).addSubColumns( {"SubGhost1", "SubGhost2"} ),
                          TableLayout::Column()
                            .setName( "Active" ).addSubColumns( {"SubActive1", "SubActive2"} )
                        } ),
      "Column5"
    } );

    TableData tableData;
    tableData.addRow( "3547", "1289", "7534", "6901", "4832", "9281", "1154", "5360", "2739", "9004", "1497", "6", "7" );
    tableData.addRow( "5142", "8290", "364", "2310", "7011", "1427", "2574", "9043", "5305", "608", "980", "6", "7" );
    tableData.addRow( "3174", "8259", "6092", "1783", "7435", "2891", "914", "178", "4635", "5839", "8124", "6", "7" );
    tableData.addRow( "6193", "7481", "1305", "9037", "4306.1", "6157", "1849", "2753", "910", "2369", "9992", "6", "7" );
    tableData.addRow( "8012", "5729.2112", "6975", "3201.213", "9448", "1820", "4125", "182.12", "7453", "5069", "3912", "6", "7" );
    tableData.addRow( "4381", "6728", "5204", "8663", "2035", "7804", "6310", "9621", "4158", "789", "2537", "6", "7" );

    TableTextFormatter tableText( tableLayout );
    EXPECT_EQ( tableText.toString(
                 tableData ),
               "\n------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
               "|        |   Column1   |                     Nodes                     |  Column3  |                                         Column4                                         |  Column5  |\n"
               "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
               "|        |             |  LocalesNodes  |  GhostNodes  |  ActiveNodes  |           |            Locales            |           Ghost           |           Active            |           |\n"
               "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
               "|        |             |                |              |               |           |  SubLocales1  |  SubLocales2  |  SubGhost1  |  SubGhost2  |  SubActive1  |  SubActive2  |           |\n"
               "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
               "|  3547  |       1289  |          7534  |        6901  |         4832  |     9281  |         1154  |         5360  |       2739  |       9004  |        1497  |           6  |        7  |\n"
               "|  5142  |       8290  |           364  |        2310  |         7011  |     1427  |         2574  |         9043  |       5305  |        608  |         980  |           6  |        7  |\n"
               "|  3174  |       8259  |          6092  |        1783  |         7435  |     2891  |          914  |          178  |       4635  |       5839  |        8124  |           6  |        7  |\n"
               "|  6193  |       7481  |          1305  |        9037  |       4306.1  |     6157  |         1849  |         2753  |        910  |       2369  |        9992  |           6  |        7  |\n"
               "|  8012  |  5729.2112  |          6975  |    3201.213  |         9448  |     1820  |         4125  |       182.12  |       7453  |       5069  |        3912  |           6  |        7  |\n"
               "|  4381  |       6728  |          5204  |        8663  |         2035  |     7804  |         6310  |         9621  |       4158  |        789  |        2537  |           6  |        7  |\n"
               "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n"
               );
  }
}

TEST( testTable, variadicTest )
{
  {
    TableLayout const layoutTest( "Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis nascetur ridiculus mus",
    {
      "Rank",
      TableLayout::Column()
        .setName( "Nodes" )
        .addSubColumns( {"Locales", "Ghost" } ),
      "Edge",
      TableLayout::Column()
        .setName( "Faces" )
        .addSubColumns( {"Locales", "Ghost" } ),
      TableLayout::Column()
        .setName( "Elems" )
        .addSubColumns( {"Locales", "Ghost"} ),
    } );
    std::cout << std::endl;
    TableData tableData;
    tableData.addRow( "min(local/total)", 1, 2, 3, 4, 5, 6, 7 );
    tableData.addRow( "min(local/total)", 1, 2, 3, 4, 5, 6, 7 );
    TableTextFormatter log( layoutTest );
    EXPECT_EQ( log.toString( tableData ),
               "\n--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|  Cras egestas ipsum a nisl. Vivamus variu dolor utsisicdis parturient montes, nascetur ridiculus mus. Duis nascetur ridiculus mus  |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|           Rank            |           Nodes            |     Edge      |           Faces            |            Elems             |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|                           |   Locales    |    Ghost    |               |   Locales    |    Ghost    |   Locales    |     Ghost     |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n"
               "|         min(local/total)  |           1  |          2  |            3  |           4  |          5  |           6  |            7  |\n"
               "|         min(local/total)  |           1  |          2  |            3  |           4  |          5  |           6  |            7  |\n"
               "--------------------------------------------------------------------------------------------------------------------------------------\n\n"
               );
  }
}

TEST( testTable, testLineBreak )
{
  TableLayout tableLayout( {"Cras egestas", "CoordX", "C", "CoordZ", "Prev\nelement", "Next\nelement"} );
  tableLayout.setTitle( "title" ).setMargin( TableLayout::MarginValue::tiny ).disableLineBreak();

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

TEST( testTable, testCellMerging )
{
  TableLayout const tableLayout( {
    TableLayout::Column()
      .setName( "Cras egestas" ),
    TableLayout::Column()
      .setName( "CoordX" ),
    "C",
    TableLayout::Column()
      .setName( "CoordZ" ),
    TableLayout::Column()
      .setName( "Prev\nelement" ),
    TableLayout::Column()
      .setName( "Next\nelement" )} );

  TableData tableData;
  tableData.addRow( "ProductA", 1234, 40, "ProductName", 5678, 60 );
  tableData.addRow( "ProductA", 54, 4564575, "long size value", 5454554512, 60 );
  tableData.addSeparator();
  tableData.addRow( "ProductA", 54, 4564575, "long size value", 5454554512, 60 );
  tableData.addRow( 3.14f, 2.718f, CellType::MergeNext, 1.618f, 0.577f, CellType::MergeNext );
  tableData.addRow( "P1\nP2\nP3", "2002\n2003\n2004", CellType::MergeNext, "1212121245452145454545", 4004, CellType::MergeNext );
  tableData.addRow( "Long product size", 54, 4564575, "long size value", 5454554512, 60 );
  tableData.addRow( "ProductAfdggfd", 5445, 4565, "PrName", 5454512, 64650 );
  tableData.addRow( 3.14f, 2.718f, CellType::MergeNext, 1.618f, 0.577f, CellType::MergeNext );
  tableData.addSeparator();
  tableData.addRow( CellType::MergeNext, CellType::MergeNext, CellType::MergeNext, CellType::MergeNext, CellType::MergeNext, "Item2" );
  tableData.addSeparator();
  tableData.addRow( 1500, 2500, CellType::MergeNext, CellType::MergeNext, CellType::MergeNext, CellType::MergeNext );
  tableData.addSeparator();
  tableData.addRow( 1.23f, 4.56f, CellType::MergeNext, 7.89f, 0.12f, 40 );
  tableData.addRow( "Long product size", 54, 4564575, "long size value", 5454554512, 60 );
  tableData.addRow( "ProductA", 54, 4564575, "long size value", 5454554512, 60 );
  tableData.addSeparator();
  tableData.addRow( "P1", "2002", CellType::MergeNext,  CellType::MergeNext, CellType::MergeNext, "1212121245452145454545" );
  tableData.addSeparator();
  tableData.addRow( "Alpha", 1001, 8, "Beta\nwater", "2002\n1.0", CellType::MergeNext );

  TableTextFormatter const tableText( tableLayout );
  std::cout <<tableText.toString( tableData ) << std::endl;
  EXPECT_EQ( tableText.toString( tableData ),
             "\n----------------------------------------------------------------------------\n"
             "|  Cras egestas  |  CoordX  |  C   |    CoordZ     |   Prev    |   Next    |\n"
             "|                |          |      |               |  element  |  element  |\n"
             "----------------------------------------------------------------------------\n"
             "|      ProductA  |    1234  |  40  |  ProductName  |     5678  |       60  |\n"
             "----------------------------------------------------------------------------\n"
             "|            P1  |    2002  |  1212121245452145454545  |     4004  |           |\n"
             "|            P2  |    2003  |                      |           |           |\n"
             "|            P3  |    2004  |                      |           |           |\n"
             "|          3.14  |   2.718  |               1.618  |    0.577  |           |\n"
             "----------------------------------------------------------------------------\n"
             "|                                                                   Item2  |\n"
             "----------------------------------------------------------------------------\n"
             "|          1500  |    2500  |                                              |\n"
             "----------------------------------------------------------------------------\n"
             "|          1.23  |    4.56  |                7.89  |     0.12  |       40  |\n"
             "----------------------------------------------------------------------------\n"
             "|         Alpha  |    1001  |   8  |         Beta  |     2002  |           |\n"
             "|                |          |      |        water  |      1.0  |           |\n"
             "----------------------------------------------------------------------------\n\n"
             );
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
