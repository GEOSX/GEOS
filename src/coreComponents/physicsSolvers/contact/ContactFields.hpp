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
 * @file ContactFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_CONTACTFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_CONTACTFIELDS_HPP_

#include "mesh/MeshFields.hpp"
#include "common/format/EnumStrings.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace contact
{

/**
 * @struct FractureState
 *
 * A struct for the fracture states
 */
struct FractureState
{
  enum State : integer
  {
    Stick = 0,   ///< element is closed: no jump across the discontinuity.
    NewSlip = 1, ///< element just starts sliding: no normal jump across the discontinuity, but sliding is allowed.
    Slip = 2,    ///< element is sliding: no normal jump across the discontinuity, but sliding is allowed.
    Open = 3     ///< element is open: no constraints are imposed.
  };
};

DECLARE_FIELD( iterativePenalty,
               "iterativePenalty",
               array2d< real64 >,
               1.e5,
               LEVEL_0,
               WRITE_AND_READ,
               "Penalty coefficients used in the iterative procedure of the Augmented Lagrangian Method" );

DECLARE_FIELD( rotationMatrix,
               "rotationMatrix",
               array3d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "An array that holds the rotation matrices on the fracture" );

DECLARE_FIELD( dispJump,
               "displacementJump",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Displacement jump vector in the local reference system" );

DECLARE_FIELD( slip,
               "slip",
               array1d< real64 >,
               0,
               LEVEL_0,
               NO_WRITE,
               "Slip." );

DECLARE_FIELD( deltaDispJump,
               "deltaDisplacementJump",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Delta displacement jump vector" );

DECLARE_FIELD( oldDispJump,
               "oldDisplacementJump",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Displacement jump vector at the previous time-step" );

DECLARE_FIELD( traction,
               "traction",
               array2d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Fracture traction vector in the local reference system." );

DECLARE_FIELD( deltaTraction,
               "deltaTraction",
               array2d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "An array that holds the traction increments on the fracture." );

DECLARE_FIELD( dTraction_dJump,
               "dTraction_dJump",
               array3d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of the traction w.r.t. the displacement jump." );

DECLARE_FIELD( dTraction_dPressure,
               "dTraction_dPressure",
               array1d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of the traction w.r.t. to the fluid pressure." );

DECLARE_FIELD( fractureState,
               "fractureState",
               array1d< integer >,
               FractureState::Stick,
               LEVEL_0,
               WRITE_AND_READ,
               "Fracture state." );

DECLARE_FIELD( oldFractureState,
               "oldFractureState",
               array1d< integer >,
               FractureState::Stick,
               NOPLOT,
               NO_WRITE,
               "Fracture state at the previous timestep." );


ENUM_STRINGS( FractureState::State, "stick", "new_slip", "slip", "open" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_CONTACT_CONTACTFIELDS_HPP_
