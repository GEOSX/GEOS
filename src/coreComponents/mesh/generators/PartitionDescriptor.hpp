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
 * @file PartitionDescriptor.hpp
 */

#ifndef GEOS_MESH_PARTITIONDESCRIPTOR_H_
#define GEOS_MESH_PARTITIONDESCRIPTOR_H_

#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include <set>


namespace geos
{

/**
 * @class PartitionDescriptor
 * @brief Simple utility to retrieve partition information in case of Metis or Spatial partition.
 */
class PartitionDescriptor
{
public:
  /**
   * @brief indicate if the partition is described using a Metis Neighbor list.
   * @return A boolean indicating if the partition is described usins a Metis neighbor list.
   */
  bool hasMetisNeighborList() const { return m_hasMetisNeighborList; }

  /**
   * @brief Sets the boolean that indicates if the partition is described using a Metis Neighbor list.
   * @param hasMetisNeighborList A boolean indicating if the partition is described usins a Metis neighbor list.
   */
  void setHasMetisNeighborList( bool hasMetisNeighborList ) { m_hasMetisNeighborList = hasMetisNeighborList; }

  /**
   * @brief Gets a reference to the list of metis neighbor list.
   * @return A reference to the Metis neighbor list.
   */
  std::set< int > const & getMetisNeighborList() const { return m_metisNeighborList; }

  /**
   * @brief Sets the list of metis neighbor list.
   * @param metisNeighborList A reference to the Metis neighbor list.
   */
  void setMetisNeighborList( std::vector< int > const & metisNeighborList )
  {
    m_metisNeighborList.clear();
    m_metisNeighborList.insert( metisNeighborList.cbegin(), metisNeighborList.cend() );
  }

  /**
   * @brief indicate if the partition is described using a spatial partition.
   * @return A boolean indicating if the parition is described using a spatial partition.
   */
  bool hasSpatialPartition() const { return !m_hasMetisNeighborList; }

  /**
   * @brief Sets the boolean that indicates if the partition is described using a Metis Neighbor list.
   * @param hasSpatialPartition a boolean indicating if the parition is described using a spatial partition.
   */
  void setHasSpatialPartition( bool hasSpatialPartition ) { m_hasMetisNeighborList = !hasSpatialPartition; }

  /**
   * @brief Returns a reference to the spatialPartition
   * @return The spatial partiton.
   */
  SpatialPartition const & getSpatialPartition() const { return m_spatialPartition; }

  /**
   * @brief Sets the spatialPartition
   * @param spatialPartition The spatial partiton.
   */
  void setSpatialPartition( SpatialPartition const & spatialPartition ) { m_spatialPartition = spatialPartition; }

private:

  bool m_hasMetisNeighborList;         //< Indicate if we use metis neighbor list or spatial partition to describe the partition
  std::set< int > m_metisNeighborList;   //< The list of neighbors computed wwith metis
  SpatialPartition m_spatialPartition; //< The spatial partition
};

}
#endif /* GEOS_MESH_PARTITIONDESCRIPTOR_H_ */
