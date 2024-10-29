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
 * @file SinglePhaseWellFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLFIELDS_HPP_

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

DECLARE_FIELD( connectionRate,
               "connectionRate",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Connection rate" );

DECLARE_FIELD( connectionRate_n,
               "connectionRate_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Connection rate at the previous converged time step" );

DECLARE_FIELD( density_n,
               "density_n",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Density at the previous converged time step" );

DECLARE_FIELD( perforationRate,
               "perforationRate",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Perforation rate" );

DECLARE_FIELD( dPerforationRate_dPres,
               "dPerforationRate_dPres",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of perforation rate with respect to pressure" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLFIELDS_HPP_
