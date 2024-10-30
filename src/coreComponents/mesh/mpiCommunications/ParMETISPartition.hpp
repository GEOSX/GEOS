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

#ifndef GEOS_MESH_MPICOMMUNICATIONS_PARMETISPARTITION_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_PARMETISPARTITION_HPP_


#include "PartitionBase.hpp"
#include "mesh/DomainPartition.hpp"


#include <map>

constexpr int nsdof = 3;
namespace geos
{

/**
 * @brief ParMETIS partitioning information.
 */
class ParMETISPartition : public PartitionBase
{
public:
  ParMETISPartition();

  ~ParMETISPartition() override;

  bool isCoordInPartition( const real64 & coord, const int dir ) const override;

  bool isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                      const real64 & boundaryRadius ) const;

  void updateSizes( arrayView1d< real64 > const domainL,
                    real64 const dt );

  void setSizes( real64 const ( &min )[ 3 ],
                 real64 const ( &max )[ 3 ] ) override;



  void setPartitions( unsigned int xPartitions,
                      unsigned int yPartitions,
                      unsigned int zPartitions ) override;

  int getColor() override;

  void repartitionMasterParticles( ParticleSubRegion & subRegion,
                                   MPI_iCommData & commData );

  void getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
                                                   MPI_iCommData & commData,
                                                   const real64 & boundaryRadius );

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
  virtual void setNeighborList( std::vector< int > const & metisNeighborList ) override
  {
    m_metisNeighborList.clear();
    m_metisNeighborList.insert( metisNeighborList.cbegin(), metisNeighborList.cend() );
  }

  /**
   * @brief Get the number of domains in each dimension for a regular partition with InternalMesh.
   * @return An array containing number of partition in X, Y and Z directions.
   */
  array1d< int > const & getPartitions() const
  {
    return m_Partitions;
  }

  /**
   * @brief Boolean like array of length 3 (space dimensions).
   *
   * 1 means periodic.
   */
  array1d< int > m_Periodic;
  /// ijk partition indexes
  array1d< int > m_coords;

private:

  /**
   * @brief Defines a distance/buffer below which we are considered in the contact zone ghosts.
   * @param bufferSize The distance.
   */
  void setContactGhostRange( const real64 bufferSize );



  /**
   * @brief Ghost position (min).
   */
  real64 m_contactGhostMin[3];

  /**
   * @brief Ghost position (max).
   */
  real64 m_contactGhostMax[3];

  /// number of partitions
  array1d< int > m_Partitions;

  /**
   * @brief Contains the global indices of the metis neighbors in case `metis` is used. Empty otherwise.
   */
  std::set< int > m_metisNeighborList;

};

}
#endif /* GEOS_MESH_MPICOMMUNICATIONS_PARMETISPARTITION_HPP_ */
