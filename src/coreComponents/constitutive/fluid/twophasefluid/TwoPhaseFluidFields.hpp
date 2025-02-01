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
 * @file TwoPhaseFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_TWOPHASEFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_TWOPHASEFLUIDFIELDS_HPP_

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "mesh/MeshFields.hpp"


namespace geos
{

namespace fields
{

namespace twophasefluid
{

using array3dLayoutPhase = array3d< real64, constitutive::multifluid::LAYOUT_PHASE >;
using array4dLayoutPhase_d = array4d< real64, constitutive::multifluid::LAYOUT_PHASE_DC >;

DECLARE_FIELD( phaseDensity,
               "phaseDensity",
               array3dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase density" );

DECLARE_FIELD( phaseDensity_n,
               "phaseDensity_n",
               array3dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase density at the previous converged time step" );

DECLARE_FIELD( dPhaseDensity,
               "dPhaseDensity",
               array4dLayoutPhase_d,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase density with respect to pressure" );

DECLARE_FIELD( phaseViscosity,
               "phaseViscosity",
               array3dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase viscosity" );

DECLARE_FIELD( dPhaseViscosity,
               "dPhaseViscosity",
               array4dLayoutPhase_d,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase viscosity with respect to pressure" );

} // namespace twophasefluid

} // namespace constitutive
} // namespace geos

#endif // GEOS_CONSTITUTIVE_FLUID_TWOPHASEFLUIDFIELDS_HPP_
