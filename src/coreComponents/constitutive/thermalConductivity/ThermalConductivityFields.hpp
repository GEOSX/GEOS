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
 * @file ThermalConductivityFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYFIELDS_HPP_
#define GEOS_CONSTITUTIVE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYFIELDS_HPP_

#include "constitutive/relativePermeability/layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace thermalconductivity
{

DECLARE_FIELD( effectiveConductivity,
               "effectiveConductivity",
               array3d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Effective conductivity" );

DECLARE_FIELD( rockThermalConductivity,
               "rockThermalConductivity",
               array3d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Rock thermal conductivity" );

}

}

}

#endif // GEOS_CONSTITUTIVE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYFIELDS_HPP_
