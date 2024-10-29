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
 * @file traits.hpp
 */
#ifndef GEOS_LINEARALGEBRA_UTILITIES_TRAITS_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_TRAITS_HPP_

#include "codingUtilities/SFINAE_Macros.hpp"

namespace geos
{

namespace traits
{

/**
 * @brief Traits helper struct templated on the vector type.
 * @tparam VECTOR the vector type
 *
 * Suitable for detecting methods that take vectors as input.
 */
template< typename VECTOR >
struct VectorBasedTraits
{
  /// Alias for vector type
  using Vector = VECTOR;

  /// Trait for detecting applyTranspose() method on a class
  HAS_MEMBER_FUNCTION_NO_RTYPE( applyTranspose, std::declval< Vector >(), std::declval< Vector >() );
};

}

} // namespace geos

#endif //GEOS_LINEARALGEBRA_UTILITIES_TRAITS_HPP_
