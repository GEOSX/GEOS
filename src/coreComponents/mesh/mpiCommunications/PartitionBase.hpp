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

#ifndef GEOS_MESH_MPICOMMUNICATIONS_PARTITIONBASE_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_PARTITIONBASE_HPP_

#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

/**
 * @brief Base class for partitioning.
 */
class PartitionBase
{
public:

  /**
   * @brief Virtual empty destructor for C++ inheritance reasons
   */
  virtual ~PartitionBase();

  /**
   * @brief Checks if the point located inside the current partition in the given direction dir.
   * @param coord The point coordinates.
   * @param dir The considered direction.
   * @return The predicate result.
   */
  virtual bool isCoordInPartition( const real64 & coord, const int dir ) const = 0;

  /**
   * @brief Defines the dimensions of the grid.
   * @param min Global minimum spatial dimensions.
   * @param max Global maximum spatial dimensions.
   */
  virtual void setSizes( real64 const ( &min )[ 3 ],
                         real64 const ( &max )[ 3 ] ) = 0;

  /**
   * @brief Defines the number of partitions along the three (x, y, z) axis.
   * @param xPartitions Number of partitions along x.
   * @param yPartitions Number of partitions along y.
   * @param zPartitions Number of partitions along z.
   */
  virtual void setPartitions( unsigned int xPartitions,
                              unsigned int yPartitions,
                              unsigned int zPartitions ) = 0;

  virtual void setNeighborList( std::vector< int > const & neighbors ) {}

  /**
   * @brief Returns the number of colors.
   * @return The number of associated colors.
   */
  int getNumColors() const { return m_numColors; }

  int getColor() const { return m_color; }

  int getSize() const { return m_size; }

  int getRank() const { return m_rank; }

  void setNumColors( int const numColors ) { m_numColors = numColors; }

  void setColor( int const color ) { m_color = color; }

  void setSize( int const size ) { m_size = size; }

  real64 * getLocalMin() { return m_min; }
  real64 * getLocalMax() { return m_max; }
  real64 * getGlobalMin() { return m_gridMin; }
  real64 * getGlobalMax() { return m_gridMax; }

protected:
  /**
   * @brief default constructor.
   */
  PartitionBase();

  virtual void setColorValue() = 0;

protected:

  /**
   * @brief Array of neighbor communicators.
   */
  std::vector< NeighborCommunicator > m_neighbors;

  int const m_rank; ///< Rank of the current process

  /// Size of the group associated with the MPI communicator
  int m_size;

  /// the number of colors in the global partitioning
  int m_numColors;

  /// color of this partition
  int m_color;


  /// Minimum extent of partition dimensions (excluding ghost objects)
  real64 m_min[3];
  /// Maximum extent of partition dimensions (excluding ghost objects)
  real64 m_max[3];

  /// Total length of problem dimensions (excluding ghost objects).
  real64 m_gridSize[3];
  /// Minimum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMin[3];
  /// Maximum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMax[3];

};

}

#endif /* GEOS_MESH_MPICOMMUNICATIONS_PARTITIONBASE_HPP_ */
