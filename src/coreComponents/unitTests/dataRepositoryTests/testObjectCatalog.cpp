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

#define OBJECTCATALOGVERBOSE 2

// Source includes
#include "dataRepository/ObjectCatalog.hpp"
#include "dataRepository/DataContext.hpp"
#include "common/logger/Logger.hpp"
#include "mainInterface/initialization.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace LvArray;
using namespace geos;

//START_SPHINX_BASE
class Base
{
public:
  Base( int & junk, double const & junk2 )
  {
    GEOS_LOG( "calling Base constructor with arguments ("<<junk<<" "<<junk2<<")" );
  }

  virtual ~Base()
  {
    GEOS_LOG( "calling Base destructor" );
  }

  using CatalogInterface = dataRepository::CatalogInterface< Base, int &, double const & >;
  static CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string getCatalogName() = 0;
};
//STOP_SPHINX

//START_SPHINX_DERIVED1
class Derived1 : public Base
{
public:
  Derived1( int & junk, double const & junk2 ):
    Base( junk, junk2 )
  {
    GEOS_LOG( "calling Derived1 constructor with arguments ("<<junk<<" "<<junk2<<")" );
  }

  ~Derived1()
  {
    GEOS_LOG( "calling Derived1 destructor" );
  }
  static string catalogName() { return "derived1"; }
  string getCatalogName() { return catalogName(); }
};
REGISTER_CATALOG_ENTRY( Base, Derived1, int &, double const & )
//STOP_SPHINX

//START_SPHINX_DERIVED2
class Derived2 : public Base
{
public:
  Derived2( int & junk, double const & junk2 ):
    Base( junk, junk2 )
  {
    GEOS_LOG( "calling Derived2 constructor with arguments ("<<junk<<" "<<junk2<<")" );
  }

  ~Derived2()
  {
    GEOS_LOG( "calling Derived2 destructor" );
  }
  static string catalogName() { return "derived2"; }
  string getCatalogName() { return catalogName(); }

};
REGISTER_CATALOG_ENTRY( Base, Derived2, int &, double const & )
//STOP_SPHINX

//START_SPHINX_TEST
TEST( testObjectCatalog, testRegistration )
{
  GEOS_LOG( "EXECUTING MAIN" );
  int junk = 1;
  double junk2 = 3.14;
  dataRepository::DataFileContext const context = dataRepository::DataFileContext( "Base Test Class", __FILE__, __LINE__ );

  // allocate a new Derived1 object
  std::unique_ptr< Base >
  derived1 = Base::CatalogInterface::factory( "derived1", context,
                                              junk, junk2 );

  // allocate a new Derived2 object
  std::unique_ptr< Base >
  derived2 = Base::CatalogInterface::factory( "derived2", context,
                                              junk, junk2 );

  EXPECT_STREQ( derived1->getCatalogName().c_str(),
                Derived1::catalogName().c_str() );

  EXPECT_STREQ( derived2->getCatalogName().c_str(),
                Derived2::catalogName().c_str() );
  GEOS_LOG( "EXITING MAIN" );
}
//STOP_SPHINX

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
