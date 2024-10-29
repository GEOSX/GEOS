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

#include "ReverseCutHillMcKeeOrdering.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace reverseCutHillMcKeeOrdering
{

/**
 * @brief This function finds the unvisited node with the minimum degree
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] degree the degree array
 * @param[in] marker the marker array for unvisited node
 * @param[out] rootRef reference to the root (= node with the minimum degree)
 */
static void
findNodeWithMinDegree( localIndex const numRows,
                       arrayView1d< localIndex > const degree,
                       arrayView1d< localIndex > const marker,
                       localIndex & rootRef )
{
  RAJA::ReduceMinLoc< serialReduce, localIndex > minDegree( LvArray::NumericLimits< localIndex >::max, -1 );

  forAll< serialPolicy >( numRows, [=] ( localIndex const i )
  {
    if( marker[i] < 0 )
    {
      minDegree.minloc( degree[i], i );
    }
  } );
  rootRef = minDegree.getLoc();
}

/**
 * @brief This function build level structure start from root
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[in] root pointer to the root
 * @param[out] marker the marker array for unvisited node
 * @param[out] level_i points to the start/end of position on level_j, similar to CSR Matrix
 * @param[out] level_j store node number on each level
 * @param[out] numLevelsRef return the number of level on this level structure
 */
static void
buildLevel( localIndex const * const offsets,
            globalIndex const * const columns,
            localIndex const numRows,
            localIndex const rankOffset,
            localIndex const root,
            arrayView1d< localIndex > const marker,
            arrayView1d< localIndex > const level_i,
            arrayView1d< localIndex > const level_j,
            localIndex & numLevelsRef )
{
  localIndex l1 = 0;
  localIndex l2 = 0;
  localIndex lCurrent = 0;
  localIndex r1 = 0;
  localIndex r2 = 0;
  localIndex row_i = 0;
  localIndex row_j = 0;
  localIndex numLevels = 0;

  // set first level first
  level_i[0] = 0;
  level_j[0] = root;
  marker[root] = 0;
  numLevels = 1;
  l1 = 0;
  l2 = 1;
  lCurrent = l2;

  // explore nbhds of all nodes in current level
  while( l2 > l1 )
  {
    level_i[numLevels++] = l2;

    // loop through last level
    for( localIndex i = l1; i < l2; i++ )
    {
      // the node to explore
      row_i = level_j[i];
      r1 = offsets[row_i];
      r2 = offsets[row_i + 1];
      for( localIndex j = r1; j < r2; j++ )
      {
        row_j = columns[j] - rankOffset;
        if( row_j >= 0 && row_j < numRows )
        {
          if( marker[row_j] < 0 )
          {
            // unmarked row
            marker[row_j] = 0;
            level_j[lCurrent++] = row_j;
          }
        }
      }
    }
    l1 = l2;
    l2 = lCurrent;
  }
  // after this we always have a "ghost" last level
  numLevels--;

  // reset marker
  for( localIndex i = 0; i < l2; i++ )
  {
    marker[level_j[i]] = -1;
  }

  numLevelsRef = numLevels;
}


/**
 * @brief This function find a pseudo-peripheral node start from root
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[out] rootRef on return will be a end of the pseudo-peripheral
 * @param[out] marker the marker array for unvisited node
 * @param[out] level_i level structure i
 * @param[out] level_j level structure j
 */
static void
findPseudoPeripheralNode( localIndex const * const offsets,
                          globalIndex const * const columns,
                          localIndex const numRows,
                          localIndex const rankOffset,
                          localIndex & rootRef,
                          arrayView1d< localIndex > const marker,
                          arrayView1d< localIndex > const level_i,
                          arrayView1d< localIndex > const level_j )
{
  localIndex r1 = 0;
  localIndex r2 = 0;
  localIndex row = 0;
  localIndex minDegree = 0;
  localIndex levDegree = 0;
  localIndex numLevels = 0;
  localIndex newNumLevels = 0;
  localIndex root = rootRef;

  level_i.zero();
  level_j.zero();

  // build initial level structure from root
  buildLevel( offsets, columns, numRows, rankOffset, root,
              marker, level_i, level_j, newNumLevels );

  numLevels  = newNumLevels - 1;
  while( numLevels < newNumLevels )
  {
    numLevels = newNumLevels;
    r1 =  level_i[numLevels - 1];
    r2 =  level_i[numLevels];
    minDegree = numRows;
    for( localIndex i = r1; i < r2; i++ )
    {
      // select the last level, pick min-degree node
      row = level_j[i];
      levDegree = offsets[row + 1] - offsets[row];
      if( minDegree > levDegree )
      {
        minDegree = levDegree;
        root = row;
      }
    }
    buildLevel( offsets, columns, numRows, rankOffset, root,
                marker, level_i, level_j, newNumLevels );
  }
  rootRef = root;
}

/**
 * @brief This qsort is very specialized, not worth to put into utilities
 * Sort a part of array perm based on degree value (ascend)
 * That is, if degree[perm[i]] < degree[perm[j]], we should have i < j
 * @param[in] perm the perm array
 * @param[in] start start in perm
 * @param[in] end end in perm
 * @param[out] degree degree array
 */
static void
qsort( arrayView1d< localIndex > const perm,
       localIndex const start,
       localIndex const end,
       arrayView1d< localIndex > const degree )
{
  if( start >= end )
  {
    return;
  }

  std::swap( perm[start], perm[(start + end) / 2] );
  localIndex mid = start;

  // loop to split
  for( localIndex i = start + 1; i <= end; i++ )
  {
    if( degree[perm[i]] < degree[perm[start]] )
    {
      std::swap( perm[++mid], perm[i] );
    }
  }
  std::swap( perm[start], perm[mid] );
  qsort( perm, mid + 1, end, degree );
  qsort( perm, start, mid - 1, degree );
}

/**
 * @brief Last step in RCM, reverse it
 * @param[out] perm perm array
 * @param[in] start start position
 * @param[in] end end position
 */
static void
reversePermutation( arrayView1d< localIndex > const perm,
                    localIndex const start,
                    localIndex const end )
{
  localIndex i = 0;
  localIndex j = 0;
  localIndex const mid = (start + end + 1) / 2;

  for( i = start, j = end; i < mid; i++, j-- )
  {
    std::swap( perm[i], perm[j] );
  }
}


/**
 * @brief This function generate numbering for a connected component
 * @param[in] offsets row offsets in the matrix to reorder
 * @param[in] columns column indices in the matrix to reorder
 * @param[in] numRows number of rows in the matrix to reorder
 * @param[in] rankOffset offset of this rank (assuming numComps = 1)
 * @param[in] root pointer to the root
 * @param[out] marker the marker array for unvisited node
 * @param[out] perm permutation array
 * @param[out] currentNumRef number of nodes already have a perm value
 */
static void
computeNewOrdering( localIndex const * const offsets,
                    globalIndex const * const columns,
                    localIndex const numRows,
                    localIndex const rankOffset,
                    localIndex const root,
                    arrayView1d< localIndex > const marker,
                    arrayView1d< localIndex > const perm,
                    localIndex & currentNumRef )
{
  localIndex i = 0;
  localIndex j = 0;
  localIndex l1 = 0;
  localIndex l2 = 0;
  localIndex r1 = 0;
  localIndex r2 = 0;
  localIndex row_i = 0;
  localIndex row_j = 0;
  localIndex rowStart = 0;
  localIndex rowEnd = 0;
  localIndex currentNum = currentNumRef;

  marker[root]       = 0;
  l1                 = currentNum;
  perm[currentNum++] = root;
  l2                 = currentNum;

  while( l2 > l1 )
  {
    // loop through all nodes is current level
    for( i = l1; i < l2; i++ )
    {
      row_i = perm[i];
      r1 = offsets[row_i];
      r2 = offsets[row_i + 1];
      rowStart = currentNum;
      for( j = r1; j < r2; j++ )
      {
        row_j = columns[j] - rankOffset;
        if( row_j >= 0 && row_j < numRows )
        {
          if( marker[row_j] < 0 )
          {
            // save the degree in marker and add it to perm */
            marker[row_j] = offsets[row_j + 1] - offsets[row_j];
            perm[currentNum++] = row_j;
          }
        }
      }
      rowEnd = currentNum;
      qsort( perm, rowStart, rowEnd - 1, marker );
    }
    l1 = l2;
    l2 = currentNum;
  }

  // reverse
  reversePermutation( perm, currentNumRef, currentNum - 1 );
  currentNumRef = currentNum;
}

void
computePermutation( localIndex const * const offsets,
                    globalIndex const * const columns,
                    localIndex const rankOffset,
                    arrayView1d< localIndex > const perm )
{
  localIndex root = 0;
  localIndex currentNum = 0;
  localIndex const numRows = perm.size();

  // get the degree for each node
  array1d< localIndex > degree( numRows );
  array1d< localIndex > marker( numRows );
  // at most numRows levels
  array1d< localIndex > level_i( numRows+1 );
  array1d< localIndex > level_j( numRows );
  forAll< parallelHostPolicy >( numRows, [&] ( localIndex const i )
  {
    degree[i] = offsets[i + 1] - offsets[i];
    marker[i] = -1;
  } );

  // start RCM loop
  currentNum = 0;
  while( currentNum < numRows )
  {
    // Find the unvisited node with the minimum degree
    findNodeWithMinDegree( numRows, degree, marker, root );

    // This is a new connected component
    findPseudoPeripheralNode( offsets, columns, numRows, rankOffset, root, marker, level_i, level_j );

    // Numbering of this component
    computeNewOrdering( offsets, columns, numRows, rankOffset, root, marker, perm, currentNum );

  }
}

} // end namespace reverseCutHillMcKeeOrdering

} // end namespace geos
