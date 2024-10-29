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
 * @file RelativePermeabilityFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYFIELDS_HPP_
#define GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYFIELDS_HPP_

#include "constitutive/relativePermeability/layouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace relperm
{

using array2dLayoutPhase = array2d< real64, compflow::LAYOUT_PHASE >;
using array3dLayoutRelPerm = array3d< real64, constitutive::relperm::LAYOUT_RELPERM >;
using array4dLayoutRelPerm_dS = array4d< real64, constitutive::relperm::LAYOUT_RELPERM_DS >;

DECLARE_FIELD( phaseRelPerm,
               "phaseRelPerm",
               array3dLayoutRelPerm,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase relative permeability" );

DECLARE_FIELD( phaseRelPerm_n,
               "phaseRelPerm_n",
               array3dLayoutRelPerm,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase relative permeability at previous time" );

DECLARE_FIELD( dPhaseRelPerm_dPhaseVolFraction,
               "dPhaseRelPerm_dPhaseVolFraction",
               array4dLayoutRelPerm_dS,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase relative permeability with respect to phase volume fraction" );

DECLARE_FIELD( phaseTrappedVolFraction,
               "phaseTrappedVolFraction",
               array3dLayoutRelPerm,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase trapped volume fraction" );

DECLARE_FIELD( phaseMaxHistoricalVolFraction,
               "phaseMaxHistoricalVolFraction",
               array2dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase max historical phase volume fraction" );

DECLARE_FIELD( phaseMinHistoricalVolFraction,
               "phaseMinHistoricalVolFraction",
               array2dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase min historical phase volume fraction" );

}

}

}

#endif // GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYFIELDS_HPP_
