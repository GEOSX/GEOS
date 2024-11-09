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
#include "mainInterface/initialization.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/MobilityKernel.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::singlePhaseBaseKernels;

// Sphinx start after test mobility

TEST( SinglePhaseBaseKernels, mobility )
{
  int constexpr NTEST = 3;

  real64 const dens[NTEST]        = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[NTEST] = { 1e-5, 1e-10, 0.0    };
  real64 const visc[NTEST]        = { 5.0, 2.0, 1.0    };
  real64 const dVisc_dPres[NTEST] = { 1e-7, 0.0, 0.0    };

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );

    real64 mob;
    real64 dMob_dPres;

    MobilityKernel::compute( dens[i], dDens_dPres[i], visc[i], dVisc_dPres[i], mob, dMob_dPres );

    // compute etalon
    real64 const mob_et = dens[i] / visc[i];
    real64 const dMob_dPres_et = mob_et * (dDens_dPres[i] / dens[i] - dVisc_dPres[i] / visc[i]);

    EXPECT_DOUBLE_EQ( mob, mob_et );
    EXPECT_DOUBLE_EQ( dMob_dPres, dMob_dPres_et );
  }
}

// Sphinx end before test mobility

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
