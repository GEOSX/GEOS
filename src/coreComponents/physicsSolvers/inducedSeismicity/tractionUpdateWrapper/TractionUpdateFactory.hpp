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

/**
 * @file tractionUpdateFactory.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_TRACTIONUPDATEFACTORY_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_TRACTIONUPDATEFACTORY_HPP

#include "SpringSliderTractionUpdate.hpp"
#include "OneWayCoupledTractionUpdate.hpp"

namespace geos
{

namespace inducedSeismicity
{

std::unique_ptr< FaultTractionUpdateBase > tractionUpdateFactory ( TractionUpdateType const & tractionUpdateType,
                                                                   string const & contactSolverName,
                                                                   string const & flowSolverName,
                                                                   QuasiDynamicEQ * qdSolver )
{

  ContactSolverBase * contactSolver = nullptr;
  FlowSolverBase * flowSolver = nullptr;
  switch( tractionUpdateType )
  {
    case TractionUpdateType::SpringSlider:
      return std::make_unique< SpringSliderTractionUpdate >( qdSolver );
      break;
    case TractionUpdateType::OneWayCoupled:
      contactSolver = &qdSolver->getParent().getGroup< ContactSolverBase >( contactSolverName );
      flowSolver = &qdSolver->getParent().getGroup< FlowSolverBase >( flowSolverName );
      return std::make_unique< OneWayCoupledTractionUpdate >( contactSolver, flowSolver );
    default:
      GEOS_ERROR( GEOS_FMT( "Traction update type {} not recognized", tractionUpdateType ) );
      return nullptr;
      break;
  }
}

ENUM_STRINGS( TractionUpdateType,
              "SpringSlider",
              "OneWayCoupled" );

} // namespace inducedSeismicity

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_TRACTIONUPDATEFACTORY_HPP */
