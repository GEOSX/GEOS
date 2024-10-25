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

#include "SpatialPartition.hpp"
#include "codingUtilities/Utilities.hpp"
#include "LvArray/src/genericTensorOps.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"

#include <cmath>

namespace geos
{

namespace
{

// Modulo
// returns a positive value regardless of the sign of numerator
real64 Mod( real64 num, real64 denom )
{
  if( fabs( denom )<fabs( num )*1.0e-14 )
  {
    return num;
  }

  return num - denom * std::floor( num/denom );
}

// MapValueToRange
// returns a periodic value in the range [min, max)
real64 MapValueToRange( real64 value, real64 min, real64 max )
{
  return Mod( value-min, max-min )+min;
}

}

SpatialPartition::SpatialPartition():
  PartitionBase(),
  m_Periodic( nsdof ),
  m_coords( nsdof ),
  m_min{ 0.0 },
  m_max{ 0.0 },
  m_blockSize{ 1.0 },
  m_gridSize{ 0.0 },
  m_gridMin{ 0.0 },
  m_gridMax{ 0.0 },
  m_Partitions()
{
  m_size = 0;
  m_rank = 0;
  m_numColors = 8,
  setPartitions( 1, 1, 1 );
}

SpatialPartition::~SpatialPartition()
{}

void SpatialPartition::setPartitions( unsigned int xPartitions,
                                      unsigned int yPartitions,
                                      unsigned int zPartitions )
{
  m_Partitions.resize( 3 );
  m_Partitions( 0 ) = xPartitions;
  m_Partitions( 1 ) = yPartitions;
  m_Partitions( 2 ) = zPartitions;
  m_size = 1;
  for( int i = 0; i < nsdof; i++ )
  {
    m_size *= m_Partitions( i );
  }
  setContactGhostRange( 0.0 );
}

int SpatialPartition::getColor()
{
  int color = 0;

  if( isOdd( m_coords[0] ) )
  {
    color += 1;
  }

  if( isOdd( m_coords[1] ) )
  {
    color += 2;
  }

  if( isOdd( m_coords[2] ) )
  {
    color += 4;
  }

  m_numColors = 8;

  return color;
}

void SpatialPartition::addNeighbors( const unsigned int idim,
                                     MPI_Comm & cartcomm,
                                     int * ncoords )
{

  if( idim == nsdof )
  {
    bool me = true;
    for( int i = 0; i < nsdof; i++ )
    {
      if( ncoords[i] != this->m_coords( i ))
      {
        me = false;
        break;
      }
    }
    if( !me )
    {
      int const rank = MpiWrapper::cartRank( cartcomm, ncoords );
      m_neighbors.push_back( NeighborCommunicator( rank ) );
    }
  }
  else
  {
    const int dim = this->m_Partitions( LvArray::integerConversion< localIndex >( idim ) );
    const bool periodic = this->m_Periodic( LvArray::integerConversion< localIndex >( idim ) );
    for( int i = -1; i < 2; i++ )
    {
      ncoords[idim] = this->m_coords( LvArray::integerConversion< localIndex >( idim ) ) + i;
      bool ok = true;
      if( periodic )
      {
        if( ncoords[idim] < 0 )
          ncoords[idim] = dim - 1;
        else if( ncoords[idim] >= dim )
          ncoords[idim] = 0;
      }
      else
      {
        ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
      }
      if( ok )
      {
        addNeighbors( idim + 1, cartcomm, ncoords );
      }
    }
  }
}

void SpatialPartition::updateSizes( arrayView1d< real64 > const domainL,
                                    real64 const dt )
{
  for( int i=0; i<3; i++ )
  {
    real64 ratio = 1.0 + domainL[i] * dt;
    m_min[i] *= ratio;
    m_max[i] *= ratio;
    //m_PartitionLocations[i] *= ratio; ?
    m_blockSize[i] *= ratio;
    m_gridSize[i] *= ratio;
    m_gridMin[i] *= ratio;
    m_gridMax[i] *= ratio;
    m_contactGhostMin[i] *= ratio;
    m_contactGhostMax[i] *= ratio;
  }
}

void SpatialPartition::setSizes( real64 const ( &min )[ 3 ],
                                 real64 const ( &max )[ 3 ] )
{

  {
    //get size of problem and decomposition
    m_size = MpiWrapper::commSize( MPI_COMM_GEOS );

    //check to make sure our dimensions agree
    {
      int check = 1;
      for( int i = 0; i < nsdof; i++ )
      {
        check *= this->m_Partitions( i );
      }
      GEOS_ERROR_IF_NE( check, m_size );
    }

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOS, nsdof, m_Partitions.data(), m_Periodic.data(), reorder, &cartcomm );
    }
    m_rank = MpiWrapper::commRank( cartcomm );
    MpiWrapper::cartCoords( cartcomm, m_rank, nsdof, m_coords.data());

    //add neighbors
    {
      int ncoords[nsdof];
      m_neighbors.clear();
      addNeighbors( 0, cartcomm, ncoords );
    }

    MpiWrapper::commFree( cartcomm );
  }

  // global values
  LvArray::tensorOps::copy< 3 >( m_gridMin, min );
  LvArray::tensorOps::copy< 3 >( m_gridMax, max );
  LvArray::tensorOps::copy< 3 >( m_gridSize, max );
  LvArray::tensorOps::subtract< 3 >( m_gridSize, min );

  // block values
  LvArray::tensorOps::copy< 3 >( m_blockSize, m_gridSize );

  LvArray::tensorOps::copy< 3 >( m_min, min );
  for( int i = 0; i < nsdof; ++i )
  {
    const int nloc = m_Partitions( i ) - 1;
    const localIndex nlocl = static_cast< localIndex >(nloc);
    if( m_PartitionLocations[i].empty() )
    {
      // the default "even" spacing
      m_blockSize[ i ] /= m_Partitions( i );
      m_min[ i ] += m_coords( i ) * m_blockSize[ i ];
      m_max[ i ] = min[ i ] + (m_coords( i ) + 1) * m_blockSize[ i ];

      m_PartitionLocations[i].resize( nlocl );
      for( localIndex j = 0; j < m_PartitionLocations[ i ].size(); ++j )
      {
        m_PartitionLocations[ i ][ j ] = (j+1) * m_blockSize[ i ];
      }
    }
    else if( nlocl == m_PartitionLocations[i].size() )
    {
      const int parIndex = m_coords[i];
      if( parIndex == 0 )
      {
        m_min[i] = min[i];
        m_max[i] = m_PartitionLocations[i][parIndex];
      }
      else if( parIndex == nloc )
      {
        m_min[i] = m_PartitionLocations[i][parIndex-1];
        m_max[i] = max[i];
      }
      else
      {
        m_min[i] = m_PartitionLocations[i][parIndex-1];
        m_max[i] = m_PartitionLocations[i][parIndex];
      }
    }
    else
    {
      GEOS_ERROR( "SpatialPartition::setSizes(): number of partition locations does not equal number of partitions - 1\n" );
    }
  }
}

bool SpatialPartition::isCoordInPartition( const real64 & coord, const int dir ) const
{
  bool rval = true;
  const int i = dir;
  if( m_Periodic( i ))
  {
    if( m_Partitions( i ) != 1 )
    {
      real64 localCenter = MapValueToRange( coord, m_gridMin[ i ], m_gridMax[ i ] );
      rval = rval && localCenter >= m_min[ i ] && localCenter < m_max[ i ];
    }

  }
  else
  {
    rval = rval && (m_Partitions[ i ] == 1 || (coord >= m_min[ i ] && coord < m_max[ i ]));
  }

  return rval;
}

bool SpatialPartition::isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                                      const real64 & boundaryRadius ) const
// test a point relative to a boundary box. If non-zero buffer specified, expand the box.
{
  for( int i = 0; i < nsdof; i++ )
  {
    // Is particle already in bounds of partition?
    if( !(m_Partitions( i )==1 || ( elemCenter[i] >= (m_min[i] - boundaryRadius) && elemCenter[i] <= (m_max[i] + boundaryRadius) ) ) )
    {
      // Particle not in bounds, check if direction has a periodic boundary
      if( m_Periodic( i ) && (m_coords[i] == 0 || m_coords[i] == m_Partitions[i] - 1) )
      {
        // Partition minimum boundary is periodic
        if( m_coords[i] == 0 && ( (elemCenter[i] - m_gridSize[i]) < (m_min[i] - boundaryRadius) ) )
        {
          return false;
        }
        // Partition maximum boundary is periodic
        if( m_coords[i] == m_Partitions[i] - 1 && ( (elemCenter[i] + m_gridSize[i]) > (m_max[i] + boundaryRadius) ) )
        {
          return false;
        }
      }
      else
      {
        // No periodic boundary then particle is not in partition
        return false;
      }
    }
  }
  return true;
}

void SpatialPartition::setContactGhostRange( const real64 bufferSize )
{
  LvArray::tensorOps::copy< 3 >( m_contactGhostMin, m_min );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMin, -bufferSize );

  LvArray::tensorOps::copy< 3 >( m_contactGhostMax, m_max );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMax, bufferSize );
}

void SpatialPartition::repartitionMasterParticles( ParticleSubRegion & subRegion,
                                                   MPI_iCommData & commData )
{

  /*
   * Search for any particles owned by this partition, which are no longer in the
   * partition domain.  Send a copy of these particles to their new partition, but
   * keep them as ghosts on the current partition.
   *
   * This assumes that particles have already been partitioned consistent with the background
   * grid, but does not assume any specific partition topology, only that the neighbor list
   * is complete, and each partition can evaluate its own isCoordinateInPartition() function.
   *
   * After this function, each particle should be in its correct partition, and the
   * Rank of particles that were moved from the current partition will be correct,
   * but the master-ghost map in the neighbor list still needs to updated.  Particles that
   * were ghosts of an object that has been repartitioned will need to have their ghost
   * rank updated, and may need to be deleted/added elsewhere.
   *
   */

  // (1) Identify any particles that are master on the current partition, but whose center lies
  // outside of the partition domain.  Rank() for particles is defined such that it always
  // equals the rank of the owning process. Thus a particle is master iff Rank==partition.rank
  //
  // Temporarily set the ghost rank of any particle to be moved to "-1".  If the particle is
  // requested by another partition, its ghost rank will be updated.  Any particle that still
  // has a Rank=-1 at the end of this function is lost and needs to be deleted.  This
  // should only happen if it has left the global domain (hopefully at an outflow b.c.).

  arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
  arrayView1d< localIndex > const particleRank = subRegion.getParticleRank();
  array1d< R1Tensor > outOfDomainParticleCoordinates;
  std::vector< localIndex > outOfDomainParticleLocalIndices;
  unsigned int nn = m_neighbors.size();   // Number of partition neighbors.

  forAll< serialPolicy >( subRegion.size(), [&, particleCenter, particleRank] GEOS_HOST ( localIndex const pp )
    {
      bool inPartition = true;
      R1Tensor p_x;
      for( int i=0; i<3; i++ )
      {
        p_x[i] = particleCenter[pp][i];
        inPartition = inPartition && isCoordInPartition( p_x[i], i );
      }
      if( particleRank[pp]==this->m_rank && !inPartition )
      {
        outOfDomainParticleCoordinates.emplace_back( p_x ); // Store the coordinate of the out-of-domain particle
        outOfDomainParticleLocalIndices.push_back( pp );   // Store the local index "pp" for the current coordinate.
        particleRank[pp] = -1;                             // Temporarily set particleRank of out-of-domain particle to -1 until it is
                                                           // requested by someone.
      }
    } );


  // (2) Pack the list of particle center coordinates to each neighbor, and send/receive the list to neighbors.

  std::vector< array1d< R1Tensor > > particleCoordinatesReceivedFromNeighbors( nn );

  sendCoordinateListToNeighbors( outOfDomainParticleCoordinates.toView(),       // input: Single list of coordinates sent to all neighbors
                                 commData,                                      // input: Solver MPI communicator
                                 particleCoordinatesReceivedFromNeighbors );    // output: List of lists of coordinates received from each
                                                                                // neighbor.



  // (3) check the received lists for particles that are in the domain of the
  //     current partition.  make a list of the locations in the coordinate list
  //     of the particles that are to be owned by the current partition.

  std::vector< array1d< localIndex > > particleListIndicesRequestingFromNeighbors( nn );
  for( size_t n=0; n<nn; n++ )
  {
    // Loop through the unpacked list and make a list of the index of any point in partition interior domain
    for( int pp=0; pp<particleCoordinatesReceivedFromNeighbors[n].toView().size(); pp++ )
    {
      bool inPartition = true;
      for( int j=0; j<3; j++ )
      {
        inPartition = inPartition && isCoordInPartition( particleCoordinatesReceivedFromNeighbors[n].toView()[pp][j], j );
      }
      if( inPartition )
      {
        // Request particle to be transferred, and take ownership
        particleListIndicesRequestingFromNeighbors[n].emplace_back( pp );
      }
    }
  }


  // (4) Pack and send/receive list of requested indices to each neighbor.  These are the indices
  //     in the list of coordinates, not the LocalIndices on the sending processor. Unpack it
  //     and store the request list.

  std::vector< array1d< localIndex > > particleListIndicesRequestedFromNeighbors( nn );

  sendListOfIndicesToNeighbors< localIndex >( particleListIndicesRequestingFromNeighbors,
                                              commData,
                                              particleListIndicesRequestedFromNeighbors );


  // (5) Update the ghost rank of the out-of-domain particles to be equal to the rank
  //     of the partition requesting to own the particle.

  std::vector< array1d< localIndex > > particleLocalIndicesRequestedFromNeighbors( nn );
  {
    unsigned int numberOfRequestedParticles = 0;
    std::vector< int > outOfDomainParticleRequests( outOfDomainParticleLocalIndices.size(), 0 );

    for( size_t n=0; n<nn; n++ )
    {
      int ni = particleListIndicesRequestedFromNeighbors[n].size();
      numberOfRequestedParticles += ni;

      // The corresponding local index for each item in the request list is stored here:
      particleLocalIndicesRequestedFromNeighbors[n].resize( ni );

      forAll< serialPolicy >( ni, [&, particleRank] GEOS_HOST ( localIndex const k )
        {
          int const i = particleListIndicesRequestedFromNeighbors[n][k];
          outOfDomainParticleRequests[i] += 1;
          localIndex pp = outOfDomainParticleLocalIndices[i];

          particleLocalIndicesRequestedFromNeighbors[n][k] = pp;
          // Set ghost rank of the particle equal to neighbor rank.
          particleRank[pp] = m_neighbors[n].neighborRank();
        } );
    }

    // Check that there is exactly one processor requesting each out-of-domain particle.
    if( numberOfRequestedParticles != outOfDomainParticleLocalIndices.size())
    {
      std::cout << "Rank " << m_rank << " has requests for " << numberOfRequestedParticles << " out of " << outOfDomainParticleLocalIndices.size() << " out-of-domain particles" << std::endl;
    }
    for( size_t i=0; i<outOfDomainParticleRequests.size(); i++ )
    {
      if( outOfDomainParticleRequests[i] != 1 )
      {
        std::cout << "Rank " << m_rank << " particle as " << outOfDomainParticleRequests[i] << " != 1 requests!" << std::endl;
      }
    }
  }


  // (5.1) Resize particle subRegion to accommodate incoming particles.
  //       Keep track of the starting indices and number of particles coming from each neighbor.

  int oldSize = subRegion.size();
  int newSize = subRegion.size();
  std::vector< int > newParticleStartingIndices( nn );
  std::vector< int > numberOfIncomingParticles( nn );
  for( size_t n=0; n<nn; n++ )
  {
    numberOfIncomingParticles[n] = particleListIndicesRequestingFromNeighbors[n].size();
    newParticleStartingIndices[n] = newSize;
    newSize += numberOfIncomingParticles[n];
  }
  if( newSize > oldSize )
  {
    subRegion.resize( newSize ); // TODO: Does this handle constitutive fields owned by the subRegion?
  }


  // (6) Pack a buffer for the particles to be sent to each neighbor, and send/receive

  //int sizeBeforeParticleSend = subRegion.size(); // subregion size changes after this, so we need this here to use to size the deletion
  // loop
  sendParticlesToNeighbor( subRegion,
                           newParticleStartingIndices,
                           numberOfIncomingParticles,
                           commData,
                           particleLocalIndicesRequestedFromNeighbors );


  // (7) Delete any out-of-domain particles that were not requested by a neighbor.  These particles
  //     will still have Rank=-1. This should only happen if the particle has left the global domain.
  //     which will hopefully only occur at outflow boundary conditions.  If it happens for a particle in
  //     the global domain, print a warning.

  arrayView2d< real64 > const particleCenterAfter = subRegion.getParticleCenter();
  arrayView1d< int > const particleRankAfter = subRegion.getParticleRank();
  std::set< localIndex > indicesToErase;
  forAll< serialPolicy >( subRegion.size(), [&, particleRankAfter, particleCenterAfter] GEOS_HOST ( localIndex const p )
    {
      if( particleRankAfter[p] == -1 )
      {
        GEOS_LOG_RANK( "Deleting orphan out-of-domain particle during repartition at p_x = " << particleCenterAfter[p] );
        indicesToErase.insert( p );
      }
      else if( particleRankAfter[p] != m_rank )
      {
        indicesToErase.insert( p );
      }
    } );
  subRegion.erase( indicesToErase );

  // Resize particle region owning this subregion
  ParticleRegion & region = dynamicCast< ParticleRegion & >( subRegion.getParent().getParent() );
  int newRegionSize = region.getNumberOfParticles();
  region.resize( newRegionSize );
}


void SpatialPartition::getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
                                                                   MPI_iCommData & commData,
                                                                   const real64 & boundaryRadius )
{

  /*
   * Make a list of the coordinates and global IDs of all non-ghost objects on the current
   * partition.  These should all be interior to the partition domain (excluding the ghost
   * region).  Send this list to each neighbor.  Have each neighbor check the list for
   * objects within its bounding box (including ghost radius) that do not already exist
   * on the processor.  Send a request list for those objects.  Mark all ghost objects as
   * potentially abandoned.
   *
   */

  // MPM-specific code where we assume there are 2 mesh bodies and only one of them has particles
  dataRepository::Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody1 = meshBodies.getGroup< MeshBody >( 0 );
  MeshBody & meshBody2 = meshBodies.getGroup< MeshBody >( 1 );
  MeshBody & particles = meshBody1.hasParticles() ? meshBody1 : meshBody2;
  ParticleManager & particleManager = particles.getBaseDiscretization().getParticleManager();

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // (1) Identify any particles that are master on the current partition, but whose center lies
    // outside of the partition domain.  Rank() for particles is defined such that it always
    // equals the rank of the owning process. Thus a particle is master iff Rank==partition.rank
    //
    // Temporarily set the ghost rank of any particle to be moved to "-1".  If the particle is
    // requested by another partition, its ghost rank will be updated.  Any particle that still
    // has a Rank=-1 at the end of this function is lost and needs to be deleted.  This
    // should only happen if it has left the global domain (hopefully at an outflow b.c.).

    arrayView2d< real64 > const particleCenter = subRegion.getParticleCenter();
    arrayView1d< localIndex > const particleRank = subRegion.getParticleRank();
    arrayView1d< globalIndex > const particleGlobalID = subRegion.getParticleID();
    array1d< R1Tensor > inDomainMasterParticleCoordinates;   // Theoretically the same as particle position evaluated at
                                                             // subRegion.nonGhostIndices()?
    std::vector< globalIndex > inDomainMasterParticleGlobalIndices;
    unsigned int nn = m_neighbors.size();   // Number of partition neighbors.

    forAll< serialPolicy >( subRegion.size(), [&, particleCenter, particleRank, particleGlobalID] GEOS_HOST ( localIndex const p )
      {
        bool inPartition = true;
        R1Tensor p_x;
        for( int i=0; i<3; i++ )
        {
          p_x[i] = particleCenter[p][i];
          inPartition = inPartition && isCoordInPartition( p_x[i], i );
        }
        if( particleRank[p]==this->m_rank && inPartition )
        {
          inDomainMasterParticleCoordinates.emplace_back( p_x );  // Store the coordinate of the out-of-domain particle
          inDomainMasterParticleGlobalIndices.push_back( particleGlobalID[p] );     // Store the local index "pp" for the current
                                                                                    // coordinate.
        }
      } );


    // (2) Pack the list of particle center coordinates to each neighbor, and send/receive the list to neighbors.

    std::vector< array1d< R1Tensor > > particleCoordinatesReceivedFromNeighbors( nn );

    sendCoordinateListToNeighbors( inDomainMasterParticleCoordinates.toView(),          // input: Single list of coordinates sent to all
                                                                                        // neighbors
                                   commData,                                      // input: Solver MPI communicator
                                   particleCoordinatesReceivedFromNeighbors    // output: List of lists of coordinates received from each
                                                                               // neighbor.
                                   );


    // (3) Pack the list of particle global indices to each neighbor, and send the list to neighbors.

    std::vector< array1d< globalIndex > > particleGlobalIndicesSendingToNeighbors( nn );
    std::vector< array1d< globalIndex > > particleGlobalIndicesReceivedFromNeighbors( nn );

    for( size_t n=0; n<nn; ++n )
    {
      particleGlobalIndicesSendingToNeighbors[n].resize( inDomainMasterParticleGlobalIndices.size() );
      for( size_t i=0; i<inDomainMasterParticleGlobalIndices.size(); ++i )
      {
        particleGlobalIndicesSendingToNeighbors[n][i] = inDomainMasterParticleGlobalIndices[i];
      }
    }
    sendListOfIndicesToNeighbors< globalIndex >( particleGlobalIndicesSendingToNeighbors,
                                                 commData,
                                                 particleGlobalIndicesReceivedFromNeighbors );


    // (4) Look through the received coordinates and make a list of the particles that are within
    //     the bounding box (including ghost region) and for which the globalID does not already
    //     exist on the current partition.  Make a request list of the index (order in the list)
    //     of the particle.  This will be sent from the master as a new ghost on the current
    //     partition.

    std::vector< array1d< globalIndex > > particleGlobalIndicesRequestingFromNeighbors( nn );

    for( size_t n=0; n<nn; ++n )
    {
      for( localIndex i=0; i<particleCoordinatesReceivedFromNeighbors[n].size(); ++i )
      {
        if( isCoordInPartitionBoundingBox( particleCoordinatesReceivedFromNeighbors[n][i], boundaryRadius ) )
        {
          globalIndex gI = particleGlobalIndicesReceivedFromNeighbors[n][i];

          // This particle should be a ghost on the current processor. See if it already exists here.
          bool alreadyHere = false;
          forAll< serialPolicy >( subRegion.size(), [&, particleGlobalID, particleRank] GEOS_HOST ( localIndex const p )
            {
              if( gI == particleGlobalID[p] )
              {
                // The particle already exists as a ghost on this partition, so we should update its rank.
                particleRank[p] = m_neighbors[n].neighborRank();
                alreadyHere = true;
              }
            } );
          if( !alreadyHere )
          {
            // The global index is not represented on this partition, so we should add the particle.
            particleGlobalIndicesRequestingFromNeighbors[n].emplace_back( gI );
          }
        }
      }
    }


    // (5) Pack and send request list of Global Indices to each neighbor, receive and unpack
    //     this into a list requested from each neighbor.

    std::vector< array1d< globalIndex > > particleGlobalIndicesRequestedFromNeighbors( nn );

    sendListOfIndicesToNeighbors< globalIndex >( particleGlobalIndicesRequestingFromNeighbors,
                                                 commData,
                                                 particleGlobalIndicesRequestedFromNeighbors );


    // (6) Temporarily set the ghost rank of all ghosts to "-1".  After ghosts are unpacked from the
    //     masters, the ghost rank will be overwritten.  At the end of this function, any ghosts that
    //     still have ghostRank=-1 are orphans and need to be deleted.

    int partitionRank = this->m_rank;
    forAll< parallelHostPolicy >( subRegion.size(), [=] GEOS_HOST ( localIndex const p )   // TODO: Worth moving to device?
      {
        if( particleRank[p] != partitionRank )
        {
          particleRank[p] = -1;
        }
      } );


    // (6.1) Resize particle subRegion to accommodate incoming particles.
    //       Keep track of the starting indices and number of particles coming from each neighbor.

    int oldSize = subRegion.size();
    int newSize = subRegion.size();
    std::vector< int > newParticleStartingIndices( nn );
    std::vector< int > numberOfIncomingParticles( nn );
    for( size_t n=0; n<nn; n++ )
    {
      numberOfIncomingParticles[n] = particleGlobalIndicesRequestingFromNeighbors[n].size();
      newParticleStartingIndices[n] = newSize;
      newSize += numberOfIncomingParticles[n];
    }
    if( newSize > oldSize )
    {
      subRegion.resize( newSize );   // TODO: Does this handle constitutive fields owned by the subregion's parent region?
    }


    // (7) Pack/Send/Receive/Unpack particles to be sent to each neighbor.

    {
      std::vector< array1d< localIndex > > particleLocalIndicesRequestedFromNeighbors( nn );

      for( size_t n=0; n<nn; n++ )
      {
        // make a list of the local indices corresponding to the global indices in the request list for the current neighbor.
        particleLocalIndicesRequestedFromNeighbors[n].resize( particleGlobalIndicesRequestedFromNeighbors[n].size() );
        for( localIndex i=0; i<particleLocalIndicesRequestedFromNeighbors[n].size(); ++i )
        {
          particleLocalIndicesRequestedFromNeighbors[n][i] = subRegion.globalToLocalMap( particleGlobalIndicesRequestedFromNeighbors[n][i] );
        }
      }
      // Send/receive the particles, which will add missing ghosts on the partition.
      sendParticlesToNeighbor( subRegion,
                               newParticleStartingIndices,
                               numberOfIncomingParticles,
                               commData,
                               particleLocalIndicesRequestedFromNeighbors );
    }


    // (8) Delete any particles that have ghostRank=-1.  These will be ghosts from
    //     a previous step for which the master is no longer in the ghost domain,
    // std::set< localIndex > indicesToErase;
    // arrayView1d< localIndex > const particleRankNew = subRegion.getParticleRank();
    // forAll< serialPolicy >( subRegion.size(), [=, &indicesToErase] GEOS_HOST ( localIndex const p ) // parallelize with atomics
    // {
    //   if( particleRankNew[p] == -1 )
    //   {
    //     indicesToErase.insert(p);
    //   }
    // } );
    // subRegion.erase(indicesToErase);

  } );
}

/**
 * @brief Send coordinates to neighbors as part of repartition.
 * @param[in] particleCoordinatesSendingToNeighbors Single list of coordinates sent to all neighbors
 * @param[in] commData Solver's MPI communicator
 * @param[in] particleCoordinatesReceivedFromNeighbors List of lists of coordinates received from each neighbor
 */
void SpatialPartition::sendCoordinateListToNeighbors( arrayView1d< R1Tensor > const & particleCoordinatesSendingToNeighbors,
                                                      MPI_iCommData & commData,
                                                      std::vector< array1d< R1Tensor > > & particleCoordinatesReceivedFromNeighbors
                                                      )
{
  // Number of neighboring partitions
  unsigned int nn = m_neighbors.size();

  // The send buffer is identical for each neighbor, and contains the packed coordinate list.
  unsigned int sizeToBePacked = 0;                                                        // size of the outgoing data with packing=false
                                                                                          // (we need to run through it first without
                                                                                          // packing so we can size the buffer)
  unsigned int sizeOfPacked = 0;                                                          // size of the outgoing data with packing=true
  buffer_unit_type * junk;                                                                // junk buffer, stores nothing since we're just
                                                                                          // getting the buffer size on the first pass
  sizeToBePacked += bufferOps::Pack< false >( junk,
                                              particleCoordinatesSendingToNeighbors );    // get the size of the coordinate list
  buffer_type sendBuffer( sizeToBePacked );                                               // the actual sized buffer that we pack into
  buffer_unit_type * sendBufferPtr = sendBuffer.data();                                   // get a pointer to the buffer
  sizeOfPacked += bufferOps::Pack< true >( sendBufferPtr,
                                           particleCoordinatesSendingToNeighbors );       // pack the coordinate list into the buffer
  GEOS_ERROR_IF_NE( sizeToBePacked, sizeOfPacked );                                       // make sure the packer is self-consistent

  // Declare the receive buffers
  std::vector< unsigned int > sizeOfReceived( nn );
  std::vector< buffer_type > receiveBuffer( nn );

  // send the coordinate list to each neighbor.  Using an asynchronous send,
  // the mpi request will be different for each send, but the buffer is the same
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    // Send/receive the size of the packed buffer
    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      m_neighbors[n].mpiISendReceive( &sizeOfPacked,
                                      1,
                                      sendRequest[n],
                                      &(sizeOfReceived[n]),
                                      1,
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Send/receive the buffer containing the list of coordinates
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      receiveBuffer[n].resize( sizeOfReceived[n] );
      m_neighbors[n].mpiISendReceive( sendBuffer.data(),   // TODO: This can't be sendBufferPtr, why not? I guess cuz sendBufferPtr gets
                                                           // incremented (moved) during packing.
                                      sizeOfPacked,
                                      sendRequest[n],
                                      receiveBuffer[n].data(),
                                      sizeOfReceived[n],
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Unpack the received coordinate list from each neighbor
  for( size_t n=0; n<nn; n++ )
  {
    // Unpack the buffer to an array of coordinates.
    const buffer_unit_type * receiveBufferPtr = receiveBuffer[n].data();  // needed for const cast
    bufferOps::Unpack( receiveBufferPtr, particleCoordinatesReceivedFromNeighbors[n] );
  }
}

template< typename indexType >
void SpatialPartition::sendListOfIndicesToNeighbors( std::vector< array1d< indexType > > & listSendingToEachNeighbor,
                                                     MPI_iCommData & commData,
                                                     std::vector< array1d< indexType > > & listReceivedFromEachNeighbor )
{
  // Number of neighboring partitions
  unsigned int nn = m_neighbors.size();

  // Pack the outgoing lists of local indices
  std::vector< unsigned int > sizeOfPacked( nn );
  std::vector< buffer_type > sendBuffer( nn );
  for( size_t n=0; n<nn; n++ )
  {
    unsigned int sizeToBePacked = 0;                                                  // size of the outgoing data with packing=false (we
                                                                                      // need to run through it first without packing so we
                                                                                      // can size the buffer)
    sizeOfPacked[n] = 0;                                                              // size of the outgoing data with packing=true
    buffer_unit_type * junk;                                                          // junk buffer, stores nothing since we're just
                                                                                      // getting the buffer size on the first pass
    sizeToBePacked += bufferOps::Pack< false >( junk,
                                                listSendingToEachNeighbor[n] );       // get the size of the list of local indices
    sendBuffer[n].resize( sizeToBePacked );                                           // the actual sized buffer that we pack into
    buffer_unit_type * sendBufferPtr = sendBuffer[n].data();                          // get a pointer to the buffer
    sizeOfPacked[n] += bufferOps::Pack< true >( sendBufferPtr,
                                                listSendingToEachNeighbor[n] );       // pack the list of local indices into the buffer
    GEOS_ERROR_IF_NE( sizeToBePacked, sizeOfPacked[n] );                             // make sure the packer is self-consistent
  }

  // Declare the receive buffers
  std::vector< unsigned int > sizeOfReceived( nn ); // TODO: decide if these number-of-neighbor-sized arrays should be array1d, std::vector
                                                    // or std::array
  std::vector< buffer_type > receiveBuffer( nn );

  // send the list of local indices to each neighbor using an asynchronous send
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    // Send/receive the size of the packed buffer
    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      m_neighbors[n].mpiISendReceive( &(sizeOfPacked[n]),
                                      1,
                                      sendRequest[n],
                                      &(sizeOfReceived[n]),
                                      1,
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Send/receive the buffer containing the list of local indices
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      receiveBuffer[n].resize( sizeOfReceived[n] );
      m_neighbors[n].mpiISendReceive( sendBuffer[n].data(),   // TODO: This can't be sendBufferPtr, why not? I guess cuz sendBufferPtr gets
                                                              // incremented (moved) during packing.
                                      sizeOfPacked[n],
                                      sendRequest[n],
                                      receiveBuffer[n].data(),
                                      sizeOfReceived[n],
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Unpack the received list of local indices from each neighbor
  for( size_t n=0; n<nn; n++ )
  {
    // Unpack the buffer to an array of coordinates.
    const buffer_unit_type * receiveBufferPtr = receiveBuffer[n].data();  // needed for const cast
    bufferOps::Unpack( receiveBufferPtr, listReceivedFromEachNeighbor[n] );
  }
}

void SpatialPartition::sendParticlesToNeighbor( ParticleSubRegionBase & subRegion,
                                                std::vector< int > const & newParticleStartingIndices,
                                                std::vector< int > const & numberOfIncomingParticles,
                                                MPI_iCommData & commData,
                                                std::vector< array1d< localIndex > > const & particleLocalIndicesToSendToEachNeighbor )
{
  unsigned int nn = m_neighbors.size();

  // Pack the send buffer for the particles being sent to each neighbor
  std::vector< buffer_type > sendBuffer( nn );
  std::vector< unsigned int > sizeOfPacked( nn );

  for( size_t n=0; n<nn; n++ )
  {
    sizeOfPacked[n] = subRegion.particlePack( sendBuffer[n], particleLocalIndicesToSendToEachNeighbor[n].toView(), false );
    sendBuffer[n].resize( sizeOfPacked[n] );
    unsigned int sizeCheck = subRegion.particlePack( sendBuffer[n], particleLocalIndicesToSendToEachNeighbor[n].toView(), true );
    GEOS_ERROR_IF_NE( sizeCheck, sizeOfPacked[n] );
  }

  // Declare the receive buffers
  std::vector< unsigned int > sizeOfReceived( nn ); // TODO: decide if these number-of-neighbor-sized arrays should be array1d, std::vector
                                                    // or std::array
  std::vector< buffer_type > receiveBuffer( nn );

  // send/receive the size of the packed particle data to each neighbor using an asynchronous send
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      m_neighbors[n].mpiISendReceive( &(sizeOfPacked[n]),
                                      1,
                                      sendRequest[n],
                                      &(sizeOfReceived[n]),
                                      1,
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data() );
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data() );
  }

  // Send/receive the buffer containing the list of local indices
  {
    array1d< MPI_Request > sendRequest( nn );
    array1d< MPI_Status >  sendStatus( nn );
    array1d< MPI_Request > receiveRequest( nn );
    array1d< MPI_Status >  receiveStatus( nn );

    for( size_t n=0; n<nn; n++ )
    {
      // Initialize to null
      sendRequest[n] = MPI_REQUEST_NULL;
      receiveRequest[n] = MPI_REQUEST_NULL;

      // Perform comms
      receiveBuffer[n].resize( sizeOfReceived[n] );
      m_neighbors[n].mpiISendReceive( sendBuffer[n].data(),
                                      sizeOfPacked[n],
                                      sendRequest[n],
                                      receiveBuffer[n].data(),
                                      sizeOfReceived[n],
                                      receiveRequest[n],
                                      commData.commID(),
                                      MPI_COMM_GEOS );
    }
    MPI_Waitall( nn, sendRequest.data(), sendStatus.data());
    MPI_Waitall( nn, receiveRequest.data(), receiveStatus.data());
  }

  // Unpack the received particle data.
  for( size_t n=0; n<nn; n++ )
  {
    // Unpack the buffer to an array of coordinates.
    subRegion.particleUnpack( receiveBuffer[n], newParticleStartingIndices[n], numberOfIncomingParticles[n] );
  }

}

}
