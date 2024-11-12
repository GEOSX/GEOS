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
 * @file CompositionalMultiphaseBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace flow
{

using array2dLayoutPhase = array2d< real64, compflow::LAYOUT_PHASE >;
using array3dLayoutPhase = array3d< real64, constitutive::relperm::LAYOUT_MOB >;
using array3dLayoutPhase_dC = array3d< real64, compflow::LAYOUT_PHASE_DC >;
using array4dLayoutPhase_dC = array4d< real64, constitutive::relperm::LAYOUT_MOB_DC >;
using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array3dLayoutComp_dC = array3d< real64, compflow::LAYOUT_COMP_DC >;
using array3dLayoutPhaseComp = array3d< real64, compflow::LAYOUT_PHASE_COMP >;

DECLARE_FIELD( globalCompDensity,
               "globalCompDensity",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Global component density" );

DECLARE_FIELD( globalCompDensity_n,
               "globalCompDensity_n",
               array2dLayoutComp,
               0,
               NOPLOT,
               NO_WRITE,
               "Global component density updates at the previous converged time step" );

DECLARE_FIELD( globalCompDensity_k,
               "globalCompDensity_k",
               array2dLayoutComp,
               0,
               NOPLOT,
               NO_WRITE,
               "Global component density updates at the previous sequential iteration" );

DECLARE_FIELD( globalCompFraction,
               "globalCompFraction",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Global component fraction" );

DECLARE_FIELD( faceGlobalCompFraction,
               "faceGlobalCompFraction",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Face global component fraction" );

DECLARE_FIELD( dGlobalCompFraction_dGlobalCompDensity,
               "dGlobalCompFraction_dGlobalCompDensity",
               array3dLayoutComp_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of global component fraction with respect to component density" );

DECLARE_FIELD( phaseVolumeFraction,
               "phaseVolumeFraction",
               array2dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase volume fraction" );

DECLARE_FIELD( dPhaseVolumeFraction,
               "dPhaseVolumeFraction",
               array3dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase volume fraction with respect to pressure, temperature, global component density" );

DECLARE_FIELD( phaseMobility,
               "phaseMobility",
               array3dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase mobility" );

DECLARE_FIELD( dPhaseMobility,
               "dPhaseMobility",
               array4dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase volume fraction with respect to pressure, temperature, global component density" );

// this is needed for time step selector
DECLARE_FIELD( phaseVolumeFraction_n,
               "phaseVolumeFraction_n",
               array2dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase volume fraction at the previous converged time step" );

DECLARE_FIELD( phaseMobility_n,
               "phaseMobility_n",
               array3dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase mobility at the previous converged time step" );

DECLARE_FIELD( phaseOutflux,
               "phaseOutflux",
               array2dLayoutPhase,
               0,
               NOPLOT,
               NO_WRITE,
               "Phase outflux" );

DECLARE_FIELD( componentOutflux,
               "componentOutflux",
               array2dLayoutComp,
               0,
               NOPLOT,
               NO_WRITE,
               "Component outflux" );

DECLARE_FIELD( phaseCFLNumber,
               "phaseCFLNumber",
               array1d< real64 >,
               0,
               LEVEL_0,
               NO_WRITE,
               "Phase CFL number" );

DECLARE_FIELD( componentCFLNumber,
               "componentCFLNumber",
               array1d< real64 >,
               0,
               LEVEL_0,
               NO_WRITE,
               "Component CFL number" );

DECLARE_FIELD( globalCompDensityScalingFactor,
               "globalCompDensityScalingFactor",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Scaling factors for global component densities" );

DECLARE_FIELD( compAmount,
               "compAmount",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Component amount" );

DECLARE_FIELD( compAmount_n,
               "compAmount_n",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Component amount at the previous converged time step" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEFIELDS_HPP_
