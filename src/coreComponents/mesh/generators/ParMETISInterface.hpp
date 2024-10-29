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
 * @file ParMETISInterface.hpp
 */

#ifndef GEOS_MESH_GENERATORS_PARMETISINTERFACE_HPP_
#define GEOS_MESH_GENERATORS_PARMETISINTERFACE_HPP_

#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{

#if defined(GEOS_USE_HIP) // still need int32 hypre for the current hip-capable build
/// Typedef to allow us to specify required parmetis integer type.
using pmet_idx_t = int32_t;
#else
/// Typedef to allow us to specify required parmetis integer type.
using pmet_idx_t = int64_t;
#endif

namespace parmetis
{

/**
 * @brief Convert a element-node mesh map into a dual (element-element) graph
 * @param elemToNodes the input mesh represented by its elem-node map
 * @param elemDist the parallel distribution of elements: element index offset on each rank
 * @param comm the MPI communicator of processes to partition over
 * @param minCommonNodes minimum number of shared nodes to create an graph edge
 * @return a graph with an edge for every pair of elements that share at least @p minCommonNodes nodes;
 *         target element indices are global with respect to offsets in @p elemDist.
 * @note elemDist must be a comm-wide exclusive scan of elemToNodes.size();
 *       the user may compute it once and reuse in a subsequent call to partition().
 */
ArrayOfArrays< pmet_idx_t, pmet_idx_t >
meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
            arrayView1d< pmet_idx_t const > const & elemDist,
            MPI_Comm comm,
            int const minCommonNodes );

/**
 * @brief Partition a mesh according to its dual graph.
 * @param graph the input graph (edges of locally owned nodes)
 * @param vertDist the parallel distribution of vertices: vertex index offset on each rank
 * @param numParts target number of partitions
 * @param comm the MPI communicator of processes to partition over
 * @param numRefinements number of partition refinement iterations
 * @return an array of target partitions for each element in local mesh
 */
array1d< pmet_idx_t >
partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
           arrayView1d< pmet_idx_t const > const & vertDist,
           pmet_idx_t const numParts,
           MPI_Comm comm,
           int const numRefinements );

} // namespace parmetis
} // namespace geos

#endif //GEOS_MESH_GENERATORS_PARMETISINTERFACE_HPP_
