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

#include "codingUtilities/Utilities.hpp"

#include <gtest/gtest.h>

#include <map>

using namespace geos;

TEST( Utilities, MapExtraction )
{
  std::map< string, int > const m{
    { "k0", 0 },
    { "k1", 1 },
    { "k2", 2 }
  };

  EXPECT_EQ( mapKeys( m ), std::vector< string >( { "k0", "k1", "k2" } ) );
  EXPECT_EQ( mapKeys< std::set >( m ), std::set< string >( { "k0", "k1", "k2" } ) );
  EXPECT_EQ( mapValues( m ), std::vector< int >( { 0, 1, 2 } ) );
  EXPECT_EQ( mapValues< std::set >( m ), std::set< int >( { 0, 1, 2 } ) );
}
