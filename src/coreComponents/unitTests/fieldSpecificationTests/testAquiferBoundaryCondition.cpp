/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;


TEST( FieldSpecification, Aquifer )
{

  FieldSpecificationManager & fieldSpecificationManager = FieldSpecificationManager::getInstance();

  AquiferBoundaryCondition & aquiferBC = dynamicCast< AquiferBoundaryCondition & >( *fieldSpecificationManager.createChild( "Aquifer", "aquiferBoundaryCondition" ) );

  // set up the aquifer as in the simulation matched against IX

  auto & aquiferPorosity = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferPorosityString() );
  aquiferPorosity = 2e-1;

  auto & aquiferPermeability = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferPermeabilityString() );
  aquiferPermeability = 3e-13;

  auto & aquiferInitialPressure = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferInitialPressureString() );
  aquiferInitialPressure = 9e6;

  auto & aquiferWaterDensity = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferWaterDensityString() );
  aquiferWaterDensity = 962.81;

  auto & aquiferWaterViscosity = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferWaterViscosityString() );
  aquiferWaterViscosity = 0.00089;

  auto & aquiferTotalCompressibility = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferTotalCompressibilityString() );
  aquiferTotalCompressibility = 1e-10;

  auto & aquiferElevation = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferElevationString() );
  aquiferElevation = 0.0;

  auto & aquiferThickness = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferThicknessString() );
  aquiferThickness = 18.0;

  auto & aquiferInnerRadius = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferInnerRadiusString() );
  aquiferInnerRadius = 2000.0;

  auto & aquiferAngle = aquiferBC.getReference< real64 >( AquiferBoundaryCondition::viewKeyStruct::aquiferAngleString() );
  aquiferAngle = 20.0;

  aquiferBC.postInputInitializationRecursive();

  AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = aquiferBC.createKernelWrapper();

  real64 const timeAtBeginningOfStep = 0.0;
  real64 const dt = 8640.0;
  real64 const pres = 2.7212e+07;
  real64 const dPres = 0.0;
  real64 const gravCoef = -49.05;
  real64 const areaFraction = 1.0;
  real64 dAquiferVolFlux_dPres = 0.0;

  real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                          dt,
                                                          pres,
                                                          dPres,
                                                          gravCoef,
                                                          areaFraction,
                                                          dAquiferVolFlux_dPres );

  // observed flux value in the simulation matched against IX
  real64 const refAquiferVolFlux = -0.2043541482797776;

  ASSERT_NEAR( refAquiferVolFlux, aquiferVolFlux, 1e-10 );
}


int main( int argc, char * * argv )
{
  GeosxState state( basicSetup( argc, argv ) );

  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
