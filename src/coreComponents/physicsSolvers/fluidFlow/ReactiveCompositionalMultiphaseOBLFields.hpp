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
 * @file ReactiveCompositionalMultiphaseOBLFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace flow
{

using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array2dLayoutOBLOpVals = array2d< real64, compflow::LAYOUT_OBL_OPERATOR_VALUES >;
using array3dLayoutOBLOpDers = array3d< real64, compflow::LAYOUT_OBL_OPERATOR_DERIVATIVES >;

DECLARE_FIELD( globalCompFraction_n,
               "globalCompFraction_n",
               array2dLayoutComp,
               0,
               NOPLOT,
               NO_WRITE,
               "Global component fraction at the previous converged time step" );

DECLARE_FIELD( bcGlobalCompFraction,
               "bcGlobalCompFraction",
               array2dLayoutComp,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Boundary condition global component fraction" );

DECLARE_FIELD( referencePorosity,
               "referencePorosity",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Reference porosity" );

DECLARE_FIELD( referencePoreVolume,
               "referencePoreVolume",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Reference volume of porous space" );

DECLARE_FIELD( referenceRockVolume,
               "referenceRockVolume",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Reference volume of rock" );

DECLARE_FIELD( rockVolumetricHeatCapacity,
               "rockVolumetricHeatCapacity",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Volumetric rock heat capacity" );

DECLARE_FIELD( rockThermalConductivity,
               "rockThermalConductivity",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Rock thermal conductivity" );

DECLARE_FIELD( rockKineticRateFactor,
               "rockKineticRateFactor",
               array1d< real64 >,
               1, // set default value to 1 to avoid a fieldspec in the xml for most cases
               NOPLOT,
               NO_WRITE,
               "Rock kinetic rate factor" );

DECLARE_FIELD( OBLOperatorValues,
               "OBLOperatorValues",
               array2dLayoutOBLOpVals,
               0,
               NOPLOT,
               NO_WRITE,
               "Values of OBL operators" );


DECLARE_FIELD( OBLOperatorDerivatives,
               "OBLOperatorDerivatives",
               array3dLayoutOBLOpDers,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivatives of OBL operators with respect to all primary variables (Dofs)" );

DECLARE_FIELD( OBLOperatorValues_n,
               "OBLOperatorValues_n",
               array2dLayoutOBLOpVals,
               0,
               NOPLOT,
               NO_WRITE,
               "Values of OBL operators at the previous converged time step" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_REACTIVECOMPOSITIONALMULTIPHASEOBLFIELDS_HPP_
