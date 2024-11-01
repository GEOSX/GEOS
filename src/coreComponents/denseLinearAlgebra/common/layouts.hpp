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
 * @file common.hpp
 */
#ifndef GEOS_DENSELINEARALGEBRA_COMMON_LAYOUTS_HPP_
#define GEOS_DENSELINEARALGEBRA_COMMON_LAYOUTS_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

/**
 * Permutation info for row major and column major layouts
 */
struct MatrixLayout
{
  /// typedef LAI row major permutation consistent with RAJA order
  using ROW_MAJOR_PERM = RAJA::PERM_IJ;

  /// typedef LAI col major permutation consistent with RAJA order
  using COL_MAJOR_PERM = RAJA::PERM_JI;

  /// row major data unit stride dim
  constexpr static int const ROW_MAJOR = LvArray::typeManipulation::getStrideOneDimension( ROW_MAJOR_PERM {} );

  /// column major unit stride dim
  constexpr static int const COL_MAJOR = LvArray::typeManipulation::getStrideOneDimension( COL_MAJOR_PERM {} );
};

}

#endif //GEOS_DENSELINEARALGEBRA_COMMON_COMMON_HPP_
