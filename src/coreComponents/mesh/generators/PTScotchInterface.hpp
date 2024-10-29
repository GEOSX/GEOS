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
 * @file PTScotchInterface.hpp
 */

#ifndef GEOS_MESH_GENERATORS_PTSCOTCHINTERFACE_HPP_
#define GEOS_MESH_GENERATORS_PTSCOTCHINTERFACE_HPP_

#include "common/DataTypes.hpp"

#include "common/MpiWrapper.hpp"

namespace geos::ptscotch
{

/**
 * @brief Partition a mesh according to its dual graph.
 * @param graph the input graph (edges of locally owned nodes)
 * @param numParts target number of partitions
 * @param comm the MPI communicator of processes to partition over
 * @return an array of target partitions for each element in local mesh
 */
array1d< int64_t >
partition( ArrayOfArraysView< int64_t const, int64_t > const & graph,
           int64_t const numParts,
           MPI_Comm comm );

} // namespace geos::ptscotch

#endif //GEOS_MESH_GENERATORS_PTSCOTCHINTERFACE_HPP_
