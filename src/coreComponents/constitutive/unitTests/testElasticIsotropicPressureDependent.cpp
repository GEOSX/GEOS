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
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/ElasticIsotropicPressureDependent.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace ::geos::constitutive;

TEST( ElasticIsotropicPressureDependentTests, testThermalExpansionAndTemperature )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ElasticIsotropicPressureDependent constitutiveModel( "model", &rootGroup );

  // Create kernel updates
  ElasticIsotropicPressureDependentUpdates updates = constitutiveModel.createKernelUpdates();

  // Test that the default TEC derivative w.r.t. temperature and the default reference temperature are nil.
  EXPECT_DOUBLE_EQ( updates.m_dThermalExpansionCoefficient_dTemperature, 0 );
  EXPECT_DOUBLE_EQ( updates.m_referenceTemperature, 0 );
}
