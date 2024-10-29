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
 * @file ParticleFluidFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDFIELDS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace particlefluid
{

DECLARE_FIELD( settlingFactor,
               "settlingFactor",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Settling factor" );

DECLARE_FIELD( dSettlingFactor_dPressure,
               "dSettlingFactor_dPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of settling factor with respect to pressure" );

DECLARE_FIELD( dSettlingFactor_dProppantConcentration,
               "dSettlingFactor_dProppantConcentration",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of settling factor with respect to proppant concentration" );

DECLARE_FIELD( dSettlingFactor_dComponentConcentration,
               "dSettlingFactor_dComponentConcentration",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of settling factor with respect to component concentration" );

DECLARE_FIELD( collisionFactor,
               "collisionFactor",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Collision factor" );

DECLARE_FIELD( dCollisionFactor_dProppantConcentration,
               "dCollisionFactor_dProppantConcentration",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of collision factor with respect to proppant concentration" );

DECLARE_FIELD( proppantPackPermeability,
               "proppantPackPermeability",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Proppant pack permeability" );

}

}

}

#endif // GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDFIELDS_HPP_
