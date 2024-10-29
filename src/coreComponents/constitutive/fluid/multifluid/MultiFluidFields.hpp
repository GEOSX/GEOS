/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace multifluid
{

using array2dLayoutFluid = array2d< real64, constitutive::multifluid::LAYOUT_FLUID >;
using array3dLayoutFluid_dC = array3d< real64, constitutive::multifluid::LAYOUT_FLUID_DC >;
using array3dLayoutPhase = array3d< real64, constitutive::multifluid::LAYOUT_PHASE >;
using array4dLayoutPhase_dC = array4d< real64, constitutive::multifluid::LAYOUT_PHASE_DC >;
using array4dLayoutPhaseComp = array4d< real64, constitutive::multifluid::LAYOUT_PHASE_COMP >;
using array5dLayoutPhaseComp_dC = array5d< real64, constitutive::multifluid::LAYOUT_PHASE_COMP_DC >;

DECLARE_FIELD( phaseFraction,
               "phaseFraction",
               array3dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase fraction" );

DECLARE_FIELD( dPhaseFraction,
               "dPhaseFraction",
               array4dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase fraction with respect to pressure, temperature, and global component fractions" );

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
               array4dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase density with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseMassDensity,
               "phaseMassDensity",
               array3dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase mass density" );

DECLARE_FIELD( dPhaseMassDensity,
               "dPhaseMassDensity",
               array4dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase mass density with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseViscosity,
               "phaseViscosity",
               array3dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase viscosity" );

DECLARE_FIELD( dPhaseViscosity,
               "dPhaseViscosity",
               array4dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase viscosity with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseEnthalpy,
               "phaseEnthalpy",
               array3dLayoutPhase,
               0,
               NOPLOT, // default behavior overridden by thermal models
               NO_WRITE,
               "Phase enthalpy" );

DECLARE_FIELD( phaseEnthalpy_n,
               "phaseEnthalpy_n",
               array3dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase enthalpy at the previous converged time step" );

DECLARE_FIELD( dPhaseEnthalpy,
               "dPhaseEnthalpy",
               array4dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase enthalpy with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseInternalEnergy,
               "phaseInternalEnergy",
               array3dLayoutPhase,
               0,
               NOPLOT, // default behavior overridden by thermal models
               NO_WRITE,
               "Phase internal energy" );

DECLARE_FIELD( phaseInternalEnergy_n,
               "phaseInternalEnergy_n",
               array3dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase internal energy at the previous converged time step" );

DECLARE_FIELD( dPhaseInternalEnergy,
               "dPhaseInternalEnergy",
               array4dLayoutPhase_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase internal energy with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( phaseCompFraction,
               "phaseCompFraction",
               array4dLayoutPhaseComp,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase component fraction" );

DECLARE_FIELD( phaseCompFraction_n,
               "phaseCompFraction_n",
               array4dLayoutPhaseComp,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase component fraction at the previous converged time step" );

DECLARE_FIELD( dPhaseCompFraction,
               "dPhaseCompFraction",
               array5dLayoutPhaseComp_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase component fraction with respect to pressure, temperature, and global component fractions" );

DECLARE_FIELD( totalDensity,
               "totalDensity",
               array2dLayoutFluid,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Total density" );

DECLARE_FIELD( totalDensity_n,
               "totalDensity_n",
               array2dLayoutFluid,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Total density at the previous converged time step" );

DECLARE_FIELD( dTotalDensity,
               "dTotalDensity",
               array3dLayoutFluid_dC,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of total density with respect to pressure, temperature, and global component fractions" );

}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDFIELDS_HPP_
