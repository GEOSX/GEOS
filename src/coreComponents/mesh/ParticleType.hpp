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
 * @file ParticleType.hpp
 */

#ifndef GEOS_MESH_PARTICLETYPE_HPP
#define GEOS_MESH_PARTICLETYPE_HPP

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

/**
 * @brief Denotes type of particle shape/interpolation scheme
 */
enum class ParticleType : integer
{
  SinglePoint,          ///< Single-point (delta dirac characteristic function)
  SinglePointBSpline,   ///< Single-point with B-spline shape functions
  CPDI,                 ///< Convected particle domain interpolation, parallelepiped domain
  CPTI,                 ///< Convected particle tetrahedral-domain interpolation, tet domain
  CPDI2,                 ///< "2nd-order" CPDI, hexahedral domain
  Count
};

/// Strings for ParticleType
ENUM_STRINGS( ParticleType,
              "SinglePoint",
              "SinglePointBSpline",
              "CPDI",
              "CPTI",
              "CPDI2" );

} // namespace geos

#endif //GEOS_MESH_PARTICLETYPE_HPP
