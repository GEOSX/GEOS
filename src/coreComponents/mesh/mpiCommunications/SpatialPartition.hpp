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

#ifndef GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_

#include "dataRepository/Group.hpp"
#include "PartitionBase.hpp"
#include "mesh/DomainPartition.hpp"

#include <array>
#include <map>

constexpr int nsdof = 3;
namespace geos
{

// inline bool isEqual( const real64& val1, const real64& val2, const real64& tolfac=0.0 )
// {
//   realT tol = 0.0;
//   if( tolfac > 1.0e-15 )
//     tol = fabs(tolfac) * (fabs(val1)+fabs(val2))*0.5;
//   return val1<=(val2+tol) && val1>=(val2-tol);
// }

// CC: Taken from old geos
// Planar Sorter
// Sorts pairs of local and global indexes by the positions of their corresponding node points in a plane.
class PlanarSorter {

public:
	PlanarSorter(const arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD >& refPos, int dim) :
               dimension(dim), 
               refPositions(refPos) {};

	// sort operator for pairs containing local indexes (sort based on 1st element in pair)
	bool operator()(const std::pair<localIndex, localIndex>& lhs,
	                const std::pair<localIndex, localIndex>& rhs) 
  {
		bool rv = false;
		int a = 0;
		int b = 2;
		if (dimension == 0)
			a = 1;
		if (dimension == 2)
			b = 1;


		auto const & lhsVect = refPositions[lhs.first];
		auto const & rhsVect = refPositions[rhs.first];

		// const arraySlice1d<real64 const>& lhsVect = refPositions[lhs.first];
		// const arraySlice1d<real64 const>& rhsVect = refPositions[rhs.first];

		if (lhsVect[a] < rhsVect[a]) {
			rv = true;
		} else if (isEqual(lhsVect[a], rhsVect[a])
				&& (lhsVect[b] < rhsVect[b])) {
			rv = true;
		};

		return rv;
	};

	// sort operator for local indexes
	bool operator()(const localIndex& lhs, 
                  const localIndex& rhs)
  {
		bool rv = false;
		int a = 0;
		int b = 2;
		if (dimension == 0)
			a = 1;
		if (dimension == 2)
			b = 1;

		// const arraySlice1d<real64 const>& lhsVect = refPositions[lhs];
		// const arraySlice1d<real64 const>& rhsVect = refPositions[rhs];

		auto const & lhsVect = refPositions[lhs];
		auto const & rhsVect = refPositions[rhs];

		if (lhsVect[a] < rhsVect[a]) {
			rv = true;
		} else if (isEqual(lhsVect[a], rhsVect[a])
				&& (lhsVect[b] < rhsVect[b])) {
			rv = true;
		};

		return rv;
	};

private:
	int dimension;
	const arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD >& refPositions;
};

/**
 * @brief Concrete (cartesian?) partitioning.
 */
class SpatialPartition : public PartitionBase
{
public:
  SpatialPartition( string const & name,
                    Group * const parent );

  ~SpatialPartition() override;

    struct viewKeyStruct
  {
    static constexpr char const * periodicString() { return "periodic"; }
    static constexpr char const * minString() { return "min"; }
    static constexpr char const * maxString() { return "max"; }
    static constexpr char const * partitionLocationsString() { return "partitionLocations"; }
    static constexpr char const * blockSizeString() { return "blockSize"; }
    static constexpr char const * gridSizeString() { return "gridSize"; }
    static constexpr char const * gridMinString() { return "gridMin"; }
    static constexpr char const * gridMaxString() { return "gridMax"; }
    static constexpr char const * contactGhostMinString() { return "contactGhostMin"; }
    static constexpr char const * contactGhostMaxString() { return "contactGhostMax"; }
  } partitionViewKeys;

  static string catalogName() { return "SpatialPartition"; }
  
  virtual string getCatalogName() const override { return catalogName(); }

  void postInputInitialization() override; 

  bool isCoordInPartition( const real64 & coord, const int dir ) const override;

  bool isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                      const real64 & boundaryRadius ) const;

  void updateSizes( arrayView1d< real64 > const domainL,
                    real64 const dt );

//  void setSizes( real64 const ( &min )[ 3 ],
//                 real64 const ( &max )[ 3 ] ) override;

  void initializeNeighbors();

  // real64 * getLocalMin()
  array1d< real64 > const & getLocalMin()
  {
    return m_min;
  }

  // real64 * getLocalMax()
  array1d< real64 > const & getLocalMax()
  {
    return m_max;
  }

  // real64 * getGlobalMin()
  array1d< real64 > const & getGlobalMin()
  {
    return m_gridMin;
  }

  // real64 * getGlobalMax()
  array1d< real64 > const & getGlobalMax()
  {
    return m_gridMax;
  }

  void setCoords( array1d< int > coords ) {
    m_coords = coords;
  }

  /**
   * @brief Get the ijk coordinates of the partition in the domain.
   * @return An array containing number of partition in X, Y and Z directions.
   */
  array1d< int > const & getCoords() const
  {
    return m_coords;
  }
  
  void setPartitions( unsigned int xPartitions,
                      unsigned int yPartitions,
                      unsigned int zPartitions ) override;

  /**
   * @brief Get the number of domains in each dimension for a regular partition with InternalMesh.
   * @return An array containing number of partition in X, Y and Z directions.
   */
  array1d< int > const & getPartitions() const
  {
    return m_partitions;
  }

  void setPeriodic( array1d< int > periodic ) {
    m_periodic = periodic;
  }

  array1d< int > const & getPeriodic() const {
    return m_periodic;
  }

  int getColor() override;

  void repartitionMasterParticles( ParticleSubRegion & subRegion,
                                   MPI_iCommData & commData );

  void getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
                                                   MPI_iCommData & commData,
                                                   const real64 & boundaryRadius );

  //CC: overrides global indices on periodic faces so they are matched when finding neighboring nodes
  void setPeriodicDomainBoundaryObjects( MeshBody & grid,
                                         NodeManager & nodeManager,
                                         EdgeManager & edgeManager,
                                         FaceManager & faceManager );

  /**
   * @brief Send coordinates to neighbors as part of repartition.
   * @param[in] particleCoordinatesSendingToNeighbors Single list of coordinates sent to all neighbors
   * @param[in] commData Solver's MPI communicator
   * @param[in] particleCoordinatesReceivedFromNeighbors List of lists of coordinates received from each neighbor
   */
  void sendCoordinateListToNeighbors( arrayView1d< R1Tensor > const & particleCoordinatesSendingToNeighbors,
                                      MPI_iCommData & commData,
                                      std::vector< array1d< R1Tensor > > & particleCoordinatesReceivedFromNeighbors
                                      );

  template< typename indexType >
  void sendListOfIndicesToNeighbors( std::vector< array1d< indexType > > & listSendingToEachNeighbor,
                                     MPI_iCommData & commData,
                                     std::vector< array1d< indexType > > & listReceivedFromEachNeighbor );

  void sendParticlesToNeighbor( ParticleSubRegionBase & subRegion,
                                std::vector< int > const & newParticleStartingIndices,
                                std::vector< int > const & numberOfIncomingParticles,
                                MPI_iCommData & commData,
                                std::vector< array1d< localIndex > > const & particleLocalIndicesToSendToEachNeighbor );

  /**
   * @brief Get the metis neighbors indices, const version. @see DomainPartition#m_metisNeighborList
   * @return Container of global indices.
   */
  std::set< int > const & getMetisNeighborList() const
  {
    return m_metisNeighborList;
  }

  /**
   * @brief Sets the list of metis neighbor list.
   * @param metisNeighborList A reference to the Metis neighbor list.
   */
  void setMetisNeighborList( std::set< int > const & metisNeighborList )
  {
    m_metisNeighborList = metisNeighborList;
  }

  void setGrid( std::array< real64, 9 > const & grid )
  {
    m_gridSize[0] = grid[0];
    m_gridSize[1] = grid[1];
    m_gridSize[2] = grid[2];

    m_gridMin[0] = grid[3];
    m_gridMin[1] = grid[4];
    m_gridMin[2] = grid[5];

    m_gridMax[0] = grid[6];
    m_gridMax[1] = grid[7];
    m_gridMax[2] = grid[8];
  }

  void setBlockSize( std::array< real64, 3 > const & blockSize )
  {
    m_blockSize[0] = blockSize[0];
    m_blockSize[1] = blockSize[1];
    m_blockSize[2] = blockSize[2];
  }

  void setBoundingBox( std::array< real64, 6 > const & bb )
  {
    m_min[0] = bb[0];
    m_min[1] = bb[1];
    m_min[2] = bb[2];

    m_max[0] = bb[3];
    m_max[1] = bb[4];
    m_max[2] = bb[5];
  }

private:

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * @note Rough copy/paste of DomainPartition::AddNeighbors
   */
  void addNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );

  /**
   * @brief Defines a distance/buffer below which we are considered in the contact zone ghosts.
   * @param bufferSize The distance.
   */
  void setContactGhostRange( const real64 bufferSize );

  /// Minimum extent of partition dimensions (excluding ghost objects)
  array1d< real64 > m_min;

  /// Maximum extent of partition dimensions (excluding ghost objects)
  array1d< real64 > m_max;

  /// Locations of partition boundaries
  array1d< array1d< real64 > > m_partitionLocations;

  /// Length of partition dimensions (excluding ghost objects).
  array1d< real64 > m_blockSize;

  /// Minimum extent of problem dimensions (excluding ghost objects).
  array1d< real64 > m_gridMin;

  /// Maximum extent of problem dimensions (excluding ghost objects).
  array1d< real64 > m_gridMax;
  
  /// Total length of problem dimensions (excluding ghost objects).
  array1d< real64 > m_gridSize;
  
  /// ijk partition indexes
  array1d< int > m_coords;
  
  /// number of partitions
  array1d< int > m_partitions;

  /// Flag for periodicity in each direction
  array1d< int > m_periodic;

  /**
   * @brief Ghost position (min).
   */
  array1d< real64 > m_contactGhostMin;

  /**
   * @brief Ghost position (max).
   */
   array1d< real64 > m_contactGhostMax;

  /**
   * @brief Contains the global indices of the metis neighbors in case `metis` is used. Empty otherwise.
   */
  std::set< int > m_metisNeighborList;
};

}
#endif /* GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_ */
