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
 * @file SinglePhaseThermalConductivityFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYFIELDS_HPP_
#define GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYFIELDS_HPP_

#include "constitutive/relativePermeability/layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace thermalconductivity
{

DECLARE_FIELD( dEffectiveConductivity_dPorosity,
               "dEffectiveConductivity_dPorosity",
               array3d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of effective conductivity with respect to porosity" );

}

}

}

#endif // GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYFIELDS_HPP_
