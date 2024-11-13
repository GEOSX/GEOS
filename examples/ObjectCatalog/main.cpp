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



#include <iostream>

#include "common/logger/Logger.hpp"
#include "Base.hpp"
#include "Derived1.hpp"

int main( int argc, char *argv[] )
{
  GEOS_LOG( "EXECUTING MAIN" );
  int junk = 1;
  double junk2 = 3.14;
  double junk3 = 2*3.14;
  Parameter param;


  GEOS_LOG( "Attempting to create a Derived1 object" );
  std::unique_ptr<Base> derived1 = Base::CatalogInterface::Factory( "derived1", junk, junk2, param );
  GEOS_LOG( "Attempting to create a Derived2 object" );
  std::unique_ptr<Base> derived2 = Base::CatalogInterface::Factory( "derived2", junk, junk3, param );

  Base::CatalogInterface::catalog_cast<Derived1>( *(derived2.get()));
  GEOS_LOG( "EXITING MAIN" );
}
