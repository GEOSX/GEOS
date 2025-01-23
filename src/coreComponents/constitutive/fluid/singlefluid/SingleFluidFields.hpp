/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SingleFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUIDFIELDS_HPP_

#include "mesh/MeshFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidUtils.hpp"
namespace geos
{

namespace fields
{

namespace singlefluid
{

using array2dLayoutFluid = array2d< real64, constitutive::singlefluid::LAYOUT_FLUID >;
using array3dLayoutFluid_der = array3d< real64, constitutive::singlefluid::LAYOUT_FLUID_DER >;

DECLARE_FIELD( density,
               "density",
               array2dLayoutFluid,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Density" );

DECLARE_FIELD( dDensity,
               "dDensity",
               array3dLayoutFluid_der,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "dDensity" );

DECLARE_FIELD( density_n,
               "density_n",
               array2dLayoutFluid,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Density at the previous converged time step" );

DECLARE_FIELD( viscosity,
               "viscosity",
               array2dLayoutFluid,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Viscosity" );

DECLARE_FIELD( dViscosity,
               "dViscosity",
               array3dLayoutFluid_der,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "dViscosity" );

DECLARE_FIELD( internalEnergy,
               "internalEnergy",
               array2dLayoutFluid,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "InternalEnergy" );

DECLARE_FIELD( dInternalEnergy,
               "dInternalEnergy",
               array3dLayoutFluid_der,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "dInternalEnergy" );

DECLARE_FIELD( internalEnergy_n,
               "internalEnergy_n",
               array2dLayoutFluid,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Fluid internal energy at the previous converged step" );

DECLARE_FIELD( enthalpy,
               "enthalpy",
               array2dLayoutFluid,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Enthalpy" );

DECLARE_FIELD( dEnthalpy,
               "dEnthalpy",
               array3dLayoutFluid_der,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "dEnthalpy" );
}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_SINGLEFLUIDFIELDS_HPP_
