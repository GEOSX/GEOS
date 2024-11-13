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
 * @file PermeabilityFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_PERMEABILITYFIELDS_HPP_
#define GEOS_CONSTITUTIVE_PERMEABILITYFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace permeability
{

DECLARE_FIELD( permeability,
               "permeability",
               array3d< real64 >,
               -1,
               LEVEL_0,
               WRITE_AND_READ,
               "Rock permeability" );

DECLARE_FIELD( dPerm_dPressure,
               "dPerm_dPressure",
               array3d< real64 >,
               0,
               LEVEL_3,
               WRITE_AND_READ,
               "Derivative of rock permeability with respect to pressure" );

DECLARE_FIELD( dPerm_dAperture,
               "dPerm_dAperture",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of rock permeability with respect to aperture" );

DECLARE_FIELD( dPerm_dDispJump,
               "dPerm_dDispJump",
               array4d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of rock permeability with respect to displacement jump" );

DECLARE_FIELD( dPerm_dTraction,
               "dPerm_dTraction",
               array4d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of rock permeability with respect to the traction vector" );

DECLARE_FIELD( permeabilityMultiplier,
               "permeabilityMultiplier",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Rock permeability multiplier" );

}

}

}

#endif // GEOS_CONSTITUTIVE_PERMEABILITYFIELDS_HPP_
