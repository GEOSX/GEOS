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
 * @file CompositionalMultiphaseUtilities.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEUTILITIES_H_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEUTILITIES_H_

#include "common/DataTypes.hpp"

namespace geos
{

namespace compositionalMultiphaseUtilities
{

/**
 * @brief In each block, shift the elements from 0 to numRowsToShift-1 one position ahead and replaces the first element
 * with the sum of all elements from 0 to numRowsToShift-1 of the input block one-dimensional array of values.
 *
 * @tparam VEC type of one-dimensional array of values
 * @param numRowsToShift number of rows to shift in each block
 * @param numRowsInBlock block size
 * @param numBlocks number of blocks
 * @param v block one-dimensional array of values
 *
 * @detail numRowsToShift is different from numRowsInBlock if there is an equation that we do *not* want to
 * modify, like the energy flux for thermal simulations. In that specific case, numRowsToShift
 * is set to numRowsInBlock-1 to make sure that we don't modify the last row of each block (the energy flux)
 *
 * Let's consider the following blocks arising when computing the localFlux in the thermal compositional solver (for two-component flow)
 *
 *
 *    (mass_comp_1)_1       \
 *    (mass_comp_2)_1       |=> block 1
 *    (energy)_1            /
 *
 *    (mass_comp_1)_2       \
 *    (mass_comp_2)_2       |=> block 2
 *    (energy)_2            /
 *
 * In this case:
 *      - numRowsToShift = 2
 *      - numRowsInBlock = 3
 *      - numBlocks = 2
 */
template< typename VEC >
GEOS_HOST_DEVICE
void shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( integer const numRowsToShift,
                                                                integer const numRowsInBlock,
                                                                integer const numBlocks,
                                                                VEC && v )
{
  for( integer i = 0; i < numBlocks; ++i )
  {
    integer const ind = i * numRowsInBlock + numRowsToShift - 1;
    real64 tmp = v[ind];
    for( int j = ind - 1; j >= i * numRowsInBlock; --j )
    {
      v[j+1] = v[j];
      tmp += v[j];
    }
    v[i*numRowsInBlock] = tmp;
  }
}

/**
 * @brief Shifts all elements one position ahead and replaces the first element
 * with the sum of all elements in the input one-dimensional array of values.
 *
 * @tparam VEC type of one-dimensional array of values
 * @param numRows vector size
 * @param v one-dimensional array of values
 */
template< typename VEC >
GEOS_HOST_DEVICE
void shiftElementsAheadByOneAndReplaceFirstElementWithSum( integer const numRows,
                                                           VEC && v )
{
  shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numRows, numRows, 1, v );
}

/**
 * @brief In each block, shift the elements from 0 to numRowsToShift-1, shifts all rows one position ahead
 * and replaces the first row with the sum of all rows from 0 to numRowsToShift-1 of the input block row
 * two-dimensional array of values.
 *
 * @tparam MATRIX type of two-dimensional array of values
 * @tparam VEC type of one-dimensional array of values
 * @param numRowsToShift number of rows to shift in each block
 * @param numRowsInBlock block number of rows
 * @param numColsInBlock block number of columns
 * @param numBlocks number of row blocks
 * @param mat block row two-dimensional array of values
 * @param work one-dimensional working array of values of size numColsInBlock
 *
 * @detail numRowsToShift is different from numRowsInBlock if there is an equation that we do *not* want to
 * modify, like the energy flux for thermal simulations. In that specific case, numRowsToShift
 * is set to numRowsInBlock-1 to make sure that we don't modify the last row of each block (the energy flux)
 *
 * Let's consider the following blocks arising when computing the localFluxJacobian in the thermal compositional solver (for two-component
 * flow)
 *
 *                       p    \rho_1   \rho_2   T
 *    (mass_comp_1)_1    .       .       .      .       \
 *    (mass_comp_2)_1    .       .       .      .       |=> block 1
 *    (energy)_1         .       .       .      .       /
 *
 *    (mass_comp_1)_2    .       .       .      .       \
 *    (mass_comp_2)_2    .       .       .      .       |=> block 2
 *    (energy)_2         .       .       .      .       /
 *
 * In this case:
 *      - numRowsToShift = 2
 *      - numRowsInBlock = 3
 *      - numColsInBlock = 4
 *      - numBlocks = 2
 *
 */
template< typename MATRIX, typename VEC >
GEOS_HOST_DEVICE
void shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( integer const numRowsToShift,
                                                              integer const numRowsInBlock,
                                                              integer const numColsInBlock,
                                                              integer const numBlocks,
                                                              MATRIX && mat,
                                                              VEC && work )
{
  for( integer k = 0; k < numBlocks; ++k )
  {
    integer const ind = k * numRowsInBlock + numRowsToShift - 1;
    for( integer j = 0; j < numColsInBlock; ++j )
    {
      work[j] = mat[ind][j];
    }
    for( integer i = ind - 1; i >= k * numRowsInBlock; --i )
    {
      for( integer j = 0; j < numColsInBlock; ++j )
      {
        mat[i+1][j] = mat[i][j];
        work[j] += mat[i][j];
      }
    }
    for( integer j = 0; j < numColsInBlock; ++j )
    {
      mat[k*numRowsInBlock][j] = work[j];
    }
  }
}

/**
 * @brief Shifts all rows one position ahead and replaces the first row
 * with the sum of all rows in the input block row two-dimensional array
 * of values.
 *
 * @tparam MATRIX type of two-dimensional array of values
 * @tparam VEC type of one-dimensional array of values
 * @param numRowsInBlock number of rows
 * @param numColsInBlock number of columns
 * @param mat block row two-dimensional array of values
 * @param work one-dimensional working array of values of size numColsInBlock
 */
template< typename MATRIX, typename VEC >
GEOS_HOST_DEVICE
void shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( integer const numRowsInBlock,
                                                         integer const numColsInBlock,
                                                         MATRIX && mat,
                                                         VEC && work )
{
  shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numRowsInBlock, numRowsInBlock, numColsInBlock, 1, mat, work );
}

} // namespace compositionalMultiphaseUtilities

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEUTILITIES_HPP_
