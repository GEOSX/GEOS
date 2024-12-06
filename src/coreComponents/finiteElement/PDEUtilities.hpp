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
 * @file PDEUtilities.hpp
 */

#ifndef GEOS_FINITEELEMENT_PDEUTILITIES_HPP_
#define GEOS_FINITEELEMENT_PDEUTILITIES_HPP_

#include "common/format/EnumStrings.hpp"

namespace geos
{

namespace PDEUtilities
{

/**
 * @brief Differential operator type.
 */
enum class DifferentialOperator : integer
{
  Divergence,
  Gradient,
  Identity,
  SymmetricGradient
};

/**
 * @brief Function space type.
 */
enum class FunctionSpace : integer
{
  invalid,
  P0,
  H1,
  H1vector
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( FunctionSpace,
              "invalid",
              "P0",
              "H1",
              "H1vector" );

} // namespace PDEUtilities

} // namespace geos

#endif //GEOS_FINITEELEMENT_PDEUTILITIES_HPP_
