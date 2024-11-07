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
 * @file rateAndStateFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_RATEANDSTATEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_RATEANDSTATEFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace rateAndState
{

DECLARE_FIELD( slipRate,
               "slipRate",
               array1d< real64 >,
               1.0e-6,
               LEVEL_0,
               WRITE_AND_READ,
               "Slip rate" );

DECLARE_FIELD( stateVariable,
               "stateVariable",
               array1d< real64 >,
               0.6,
               LEVEL_0,
               WRITE_AND_READ,
               "Rate- and state-dependent friction state variable" );

DECLARE_FIELD( slipVelocity,
               "slipVelocity",
               array2d< real64 >,
               1.0e-6,
               LEVEL_0,
               WRITE_AND_READ,
               "Slip velocity" );

DECLARE_FIELD( stateVariable_n,
               "stateVariable_n",
               array1d< real64 >,
               0.6,
               NOPLOT,
               WRITE_AND_READ,
               "Rate- and state-dependent friction state variable at previous time step" );

DECLARE_FIELD( deltaSlip,
               "deltaSlip",
               array2d< real64 >,
               0.0,
               LEVEL_0,
               WRITE_AND_READ,
               "Slip increment" );


}

}

}

#endif // GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_INDUCEDSEISMICITYFIELDS_HPP_
