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

#include "gtest/gtest.h"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/TriaxialDriver.hpp"

#include "dataRepository/xmlWrapper.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

using namespace geos;
using namespace ::geos::constitutive;

template< typename POLICY >
void testTriaxialDriver()
{
  TriaxialDriver driver;
}


#ifdef GEOS_USE_DEVICE
TEST( TriaxialTests, testTriaxialDevice )
{
  testCamClayDriver< geos::parallelDevicePolicy< > >();
}
#endif
TEST( TriaxialTests, testTriaxialHost )
{
  testCamClayDriver< serialPolicy >();
}
