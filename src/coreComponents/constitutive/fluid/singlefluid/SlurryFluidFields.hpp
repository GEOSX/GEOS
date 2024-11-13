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
 * @file SlurryFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SLURRYFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SLURRYFLUIDFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace slurryfluid
{

DECLARE_FIELD( dDensity_dProppantConcentration,
               "dDens_dProppantConc",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of density with respect to proppant concentration" );

DECLARE_FIELD( dDensity_dComponentConcentration,
               "dDens_dCompConc",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of density with respect to component concentration" );

DECLARE_FIELD( componentDensity,
               "componentDensity",
               array3d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Component density" );

DECLARE_FIELD( dComponentDensity_dPressure,
               "dCompDens_dPres",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of component density with respect to pressure" );

DECLARE_FIELD( dComponentDensity_dComponentConcentration,
               "dCompDens_dCompConc",
               array4d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of component density with respect to component concentration" );

DECLARE_FIELD( fluidDensity,
               "FluidDensity",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Fluid density" );

DECLARE_FIELD( dFluidDensity_dPressure,
               "dFluidDens_dPres",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of fluid density with respect to pressure" );

DECLARE_FIELD( dFluidDensity_dComponentConcentration,
               "dFluidDens_dCompConc",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of fluid density with respect to component concentration" );

DECLARE_FIELD( dViscosity_dProppantConcentration,
               "dVisc_dProppantConc",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of viscosity with respect to proppant concentration" );

DECLARE_FIELD( dViscosity_dComponentConcentration,
               "dVisc_dCompConc",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of viscosity with respect to component concentration" );

DECLARE_FIELD( fluidViscosity,
               "FluidViscosity",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Fluid viscosity" );

DECLARE_FIELD( dFluidViscosity_dPressure,
               "dFluidVisc_dPres",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of fluid viscosity with respect to pressure" );

DECLARE_FIELD( dFluidViscosity_dComponentConcentration,
               "dFluidVisc_dCompConc",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of fluid viscosity with respect to component concentration" );

}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_SLURRYFLUIDFIELDS_HPP_
