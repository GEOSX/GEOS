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
 * @file ReverseCutHillMcKeeOrdering.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_REVERSECUTHILLMCKEEORDERING_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_REVERSECUTHILLMCKEEORDERING_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace reverseCutHillMcKeeOrdering
{

/**
 * @brief This function actually does the RCM ordering of a symmetric csr matrix (entire)
          The original implementation can be found in src/parcsr_ls/par_ilu.c in hypre
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[out] perm the permutation array
 */
void
computePermutation( localIndex const * const offsets,
                    globalIndex const * const columns,
                    localIndex const rankOffset,
                    arrayView1d< localIndex > const perm );


} // namespace reverseCutHillMcKeeOrdering

} // namespace geos

#endif // GEOS_LINEARALGEBRA_UTILITIES_REVERSECUTHILLMCKEEORDERING_HPP_
