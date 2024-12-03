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
 * @file inducedSeismicityFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace inducedSeismicity
{

DECLARE_FIELD( initialProjectedNormalTraction,
               "initialProjectedNormalTraction",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Initial meanNormal stress acting on the fault" );

DECLARE_FIELD( projectedNormalTraction,
               "projectedNormalTraction",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanNormal stress acting on the fault" );

DECLARE_FIELD( projectedNormalTraction_n,
               "projectedNormalTraction_n",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "meanNormal stress acting on the fault at the previous converged time step" );

DECLARE_FIELD( initialProjectedShearTraction,
               "initialProjectedShearTraction",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Initial meanShear stress acting on the fault" );

DECLARE_FIELD( projectedShearTraction,
               "projectedShearTraction",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "meanShear stress acting on the fault" );

DECLARE_FIELD( projectedShearTraction_n,
               "projectedShearTraction_n",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "meanShear stress acting on the fault at the previous converged time step" );

DECLARE_FIELD( seismicityRate,
               "seismicityRate",
               array1d< real64 >,
               1.0,
               LEVEL_0,
               WRITE_AND_READ,
               "Seismicity rate" );

DECLARE_FIELD( logDenom,
               "logDenom",
               array1d< real64 >,
               0,
               LEVEL_2,
               WRITE_AND_READ,
               "Log of the denominator of the integral form of the seismicity rate" );
}

}

}

#endif // GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_
