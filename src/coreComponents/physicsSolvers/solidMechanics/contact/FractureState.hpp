/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FractureState.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_FRACTURESTATE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_FRACTURESTATE_HPP_

#include "common/format/EnumStrings.hpp"

namespace geos
{

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

ENUM_STRINGS( FractureState::State, "stick", "new_slip", "slip", "open" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_CONTACT_FRACTURESTATE_HPP_
