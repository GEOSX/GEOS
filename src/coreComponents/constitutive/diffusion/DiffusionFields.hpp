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
 * @file DiffusionFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DIFFUSION_DIFFUSIONFIELDS_HPP_
#define GEOS_CONSTITUTIVE_DIFFUSION_DIFFUSIONFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace diffusion
{

DECLARE_FIELD( diffusivity,
               "diffusivity",
               array3d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Diffusivity" );

DECLARE_FIELD( dDiffusivity_dTemperature,
               "dDiffusivity_dTemperature",
               array3d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivatives of diffusivity with respect to temperature" );

DECLARE_FIELD( phaseDiffusivityMultiplier,
               "phaseDiffusivityMultiplier",
               array3d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase multipliers for the diffusivity coefficients" );


}

}

}

#endif // GEOS_CONSTITUTIVE_DIFFUSION_DIFFUSIONFIELDS_HPP_
