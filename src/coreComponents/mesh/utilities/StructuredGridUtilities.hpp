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

#ifndef GEOS_MESH_UTILITIES_STRUCTUREDGRIDUTILITIES_HPP
#define GEOS_MESH_UTILITIES_STRUCTUREDGRIDUTILITIES_HPP

/**
 * @file StructuredGridUtilities.hpp
 */

namespace structuredGrid
{

/**
 * @brief Given n, compute n^dim, where dim is the spatial dimension.
 * @param[in] n the input integer whose power is computed here
 * @return the power of n
 */
template< int dim >
int dimpower( int n );

/// @cond DO_NOT_DOCUMENT

template<> inline int dimpower< 1 >( int n ) { return n; }
template<> inline int dimpower< 2 >( int n ) { return n*n; }
template<> inline int dimpower< 3 >( int n ) { return n*n*n; }

/// @endcond

/**
 * @brief Given a lexicographical index, map back to the original
 * i,j,k indices of the point. The first variation here assumes
 * a uniform number of points nnx in all coordinate directions.
 * @tparam dim the number of spatial dimensions
 * @param[in] index the lexicographical index
 * @param[in] nnx the number of points in all coordinate directions
 * @param[out] indices the original (i,j,k) indices of the point
 */
template< int dim >
void map_index( const int index,
                const int nnx,
                std::vector< int > & indices );

/// @cond DO_NOT_DOCUMENT

template<>
inline
void map_index< 1 >( const int index,
                     const int nnx,
                     std::vector< int > & indices )
{
  GEOS_ASSERT_GT( nnx, index );
  GEOS_DEBUG_VAR( nnx );
  indices[0] = index;
}

template<>
inline
void map_index< 2 >( const int index,
                     const int nnx,
                     std::vector< int > & indices )
{
  GEOS_ASSERT_GT( nnx*nnx, index );
  indices[0] = index % nnx;
  indices[1] = index / nnx;
}

template<>
inline
void map_index< 3 >( const int index,
                     const int nnx,
                     std::vector< int > & indices )
{
  GEOS_ASSERT_GT( nnx*nnx*nnx, index );
  indices[0] = index % nnx;
  indices[1] = (index / nnx) % nnx;
  indices[2] = index / (nnx*nnx);
}

/// @endcond

} /* namespace structuredGrid */

#endif /* GEOS_MESH_UTILITIES_STRUCTUREDGRIDUTILITIES_HPP */
