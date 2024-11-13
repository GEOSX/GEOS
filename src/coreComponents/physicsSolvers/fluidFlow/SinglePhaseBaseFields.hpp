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
 * @file SinglePhaseBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace flow
{

DECLARE_FIELD( mobility,
               "mobility",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Mobility" );

DECLARE_FIELD( dMobility_dPressure,
               "dMobility_dPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of mobility with respect to pressure" );

DECLARE_FIELD( dMobility_dTemperature,
               "dMobility_dTemperature",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of mobility with respect to temperature" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEFIELDS_HPP_
