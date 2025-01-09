/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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

#include "PartitionDescriptorABC.hpp"
#include "common/DataTypes.hpp"

#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include <set>
#include <vector>

namespace geos
{

/**
 * @class PartitionDescriptor
 * @brief Simple utility to retrieve partition information in case of Metis or Spatial partition.
 */
class PartitionDescriptor : public PartitionDescriptorABC
{
public:
  PartitionDescriptor();
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
  std::set< int > const & getMetisNeighborList() const override { return m_metisNeighborList; }

  /**
   * @brief Sets the list of metis neighbor list.
   * @param metisNeighborList A reference to the Metis neighbor list.
   */
  void setMetisNeighborList( std::set< int > const & metisNeighborList )
  {
    m_metisNeighborList.clear();
    m_metisNeighborList.insert( metisNeighborList.cbegin(), metisNeighborList.cend() );
  }

  void setPartitions( array1d< int > const & partition );

  void setPartitions( unsigned int x,
                      unsigned int y,
                      unsigned int z );

  [[nodiscard]] array1d< int > getPartitions() const override
  {
    return m_partitions;
  }

  [[nodiscard]] array1d< int > getPeriodic() const override
  {
    return m_periodic;
  }

  [[nodiscard]] array1d< int > getCoords() const override
  {
    return m_coords;
  }

  [[nodiscard]] std::array< real64, 9 > getGrid() const override;

  [[nodiscard]] std::array< real64, 3 > getBlockSize() const override;

  [[nodiscard]] std::array< real64, 6 > getBoundingBox() const override;

  void setPeriodic( int i, int value )
  {
    m_periodic[i] = value;
  }

  void setPeriodic( array1d< int > periodic ) {
    m_periodic = periodic;
  }

  [[nodiscard]] bool isCoordInPartition( const real64 & coord, const int dir ) const;

  void setSizes( real64 const ( &min )[3],
                 real64 const ( &max )[3] );

  // /**
  //  * @brief indicate if the partition is described using a spatial partition.
  //  * @return A boolean indicating if the parition is described using a spatial partition.
  //  */
  // bool hasSpatialPartition() const { return !m_hasMetisNeighborList; }

  // /**
  //  * @brief Sets the boolean that indicates if the partition is described using a Metis Neighbor list.
  //  * @param hasSpatialPartition a boolean indicating if the parition is described using a spatial partition.
  //  */
  // void setHasSpatialPartition( bool hasSpatialPartition ) { m_hasMetisNeighborList = !hasSpatialPartition; }

  // /**
  //  * @brief Returns a reference to the spatialPartition
  //  * @return The spatial partiton.
  //  */
  // SpatialPartition const & getSpatialPartition() const { return m_spatialPartition; }

  // /**
  //  * @brief Sets the spatialPartition
  //  * @param spatialPartition The spatial partiton.
  //  */
  // void setSpatialPartition( SpatialPartition const & spatialPartition ) { m_spatialPartition = spatialPartition; }

private:


  /**
   * @brief Defines a distance/buffer below which we are considered in the contact zone ghosts.
   * @param bufferSize The distance.
   */
  void setContactGhostRange( const real64 bufferSize );

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   * @note Rough copy/paste of DomainPartition::AddNeighbors
   */
//  void addNeighbors( const unsigned int idim,
//                     MPI_Comm & cartcomm,
//                     int * ncoords );

  /// Size of the group associated with the MPI communicator
  int m_size;

  /// MPI rank of the current partition
  int m_rank;

  /// Ghost position (min).
  real64 m_contactGhostMin[3];
  /// Ghost position (max).
  real64 m_contactGhostMax[3];

  /// Minimum extent of partition dimensions (excluding ghost objects)
  real64 m_min[3];
  /// Maximum extent of partition dimensions (excluding ghost objects)
  real64 m_max[3];

  /// Locations of partition boundaries
  array1d< real64 > m_partitionLocations[3];

  /// Length of partition dimensions (excluding ghost objects).
  real64 m_blockSize[3];

  /// ijk partition indexes
  array1d< int > m_coords;

  /// Total length of problem dimensions (excluding ghost objects).
  real64 m_gridSize[3];
  /// Minimum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMin[3];
  /// Maximum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMax[3];

  /// Array of neighbor communicators.
//  std::vector< NeighborCommunicator > m_neighbors;

  /// The (x, y , y) MPI split (in case we need it)
  array1d< int > m_partitions;

  /// Boolean like array of length 3 (space dimensions). 1 means periodic.
  array1d< int > m_periodic;

  bool m_hasMetisNeighborList;         //< Indicate if we use metis neighbor list or spatial partition to describe the partition
  std::set< int > m_metisNeighborList;   //< The list of neighbors computed wwith metis
  // SpatialPartition m_spatialPartition; //< The spatial partition
};

}
#endif /* GEOS_MESH_PARTITIONDESCRIPTOR_H_ */
