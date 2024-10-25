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
 * @file CompositionalMultiphaseWellFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace well
{

using array2dLayoutFluid_dC = array2d< real64, compflow::LAYOUT_FLUID_DC >;
using array2dLayoutPhase = array2d< real64, compflow::LAYOUT_PHASE >;
using array3dLayoutPhase_dC = array3d< real64, compflow::LAYOUT_PHASE_DC >;
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
               WRITE_AND_READ,
               "Global component density at the previous converged time step" );

DECLARE_FIELD( mixtureConnectionRate,
               "wellElementMixtureConnectionRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Mixture connection rate" );

DECLARE_FIELD( mixtureConnectionRate_n,
               "wellElementMixtureConnectionRate_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Mixture connection rate at the previous converged time step" );

DECLARE_FIELD( globalCompFraction,
               "globalCompFraction",
               array2dLayoutComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Global component fraction" );

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
               "Derivative of phase volume fraction with respect to pressure, temperature, and global component density" );

DECLARE_FIELD( phaseVolumeFraction_n,
               "phaseVolumeFraction_n",
               array2dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase volume fraction at the previous converged time step" );

DECLARE_FIELD( totalMassDensity,
               "totalMassDensity",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total mass density" );

DECLARE_FIELD( dTotalMassDensity,
               "dTotalMassDensity",
               array2dLayoutFluid_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of total mass density with respect to pressure, temperature, and global component density" );


DECLARE_FIELD( compPerforationRate,
               "compPerforationRate",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Component perforation rate" );

DECLARE_FIELD( dCompPerforationRate,
               "dCompPerforationRate",
               array4d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of component perforation rate with respect to pressure temperature and global component density" );



DECLARE_FIELD( globalCompDensityScalingFactor,
               "globalCompDensityScalingFactor",
               array1d< real64 >,
               1,
               NOPLOT,
               NO_WRITE,
               "Scaling factors for global component densities" );


}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLFIELDS_HPP_
