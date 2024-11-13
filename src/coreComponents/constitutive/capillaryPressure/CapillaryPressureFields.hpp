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
 * @file CapillaryPressureFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREFIELDS_HPP_
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREFIELDS_HPP_

#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace cappres
{

using array3dLayoutCapPressure = array3d< real64, constitutive::cappres::LAYOUT_CAPPRES >;
using array4dLayoutCapPressure_dS = array4d< real64, constitutive::cappres::LAYOUT_CAPPRES_DS >;

DECLARE_FIELD( phaseCapPressure,
               "phaseCapPressure",
               array3dLayoutCapPressure,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase capillary pressure" );

DECLARE_FIELD( dPhaseCapPressure_dPhaseVolFraction,
               "dPhaseCapPressure_dPhaseVolFraction",
               array4dLayoutCapPressure_dS,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Derivative of phase capillary pressure with respect to phase volume fraction" );

DECLARE_FIELD( jFuncMultiplier,
               "jFuncMultiplier",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Multiplier for the Leverett J-function" );

}

}

}

#endif // GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSUREFIELDS_HPP_
