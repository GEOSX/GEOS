/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "PartitionDescriptor.hpp"

#include "common/MpiWrapper.hpp"

#include "LvArray/src/tensorOps.hpp"

namespace geos
{

PartitionDescriptor::PartitionDescriptor():
  m_min{ 0.0 },
  m_max{ 0.0 },
  m_blockSize{ 1.0 },
  m_coords( 3 ),
  m_gridSize{ 0.0 },
  m_gridMin{ 0.0 },
  m_gridMax{ 0.0 },
  m_partitions(),
  m_periodic( 3 )
{
  m_size = 0;
  m_rank = 0;
  setPartitions( 1, 1, 1 );
}

std::array< real64, 9 > PartitionDescriptor::getGrid() const
{
  return { m_gridSize[0], m_gridSize[1], m_gridSize[2], m_gridMin[0], m_gridMin[1], m_gridMin[2], m_gridMax[0], m_gridMax[1], m_gridMax[2] };
}

std::array< real64, 3 > PartitionDescriptor::getBlockSize() const
{
  return { m_blockSize[0], m_blockSize[1], m_blockSize[2] };
}

std::array< real64, 6 > PartitionDescriptor::getBoundingBox() const
{
  return { m_min[0], m_min[1], m_min[2], m_max[0], m_max[1], m_max[2] };
}

void PartitionDescriptor::setContactGhostRange( const real64 bufferSize )
{
  LvArray::tensorOps::copy< 3 >( m_contactGhostMin, m_min );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMin, -bufferSize );

  LvArray::tensorOps::copy< 3 >( m_contactGhostMax, m_max );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMax, bufferSize );
}


void PartitionDescriptor::setPartitions( array1d< int > const & partition )
{
  m_partitions.resize( 3 );
  m_partitions( 0 ) = partition[0];
  m_partitions( 1 ) = partition[1];
  m_partitions( 2 ) = partition[2];
  m_size = 1;
  for( int i = 0; i < 3; i++ )
  {
    m_size *= m_partitions( i );
  }
  setContactGhostRange( 0.0 );
}


void PartitionDescriptor::setPartitions( unsigned int x,
                                         unsigned int y,
                                         unsigned int z )
{
  array1d< int > p( 3 );
  p[0] = x;
  p[1] = y;
  p[2] = z;
  setPartitions( p );
}


namespace
{

// Modulo
// returns a positive value regardless of the sign of numerator
real64 Mod( real64 num, real64 denom )
{
  if( fabs( denom ) < fabs( num ) * 1.0e-14 )
  {
    return num;
  }

  return num - denom * std::floor( num / denom );
}

// MapValueToRange
// returns a periodic value in the range [min, max)
real64 MapValueToRange( real64 value,
                        real64 min,
                        real64 max )
{
  return Mod( value - min, max - min ) + min;
}

}


bool PartitionDescriptor::isCoordInPartition( const real64 & coord,
                                              const int dir ) const
{
  bool rval = true;
  const int i = dir;
  if( m_periodic( i ) )
  {
    if( m_partitions( i ) != 1 )
    {
      real64 localCenter = MapValueToRange( coord, m_gridMin[i], m_gridMax[i] );
      rval = rval && localCenter >= m_min[i] && localCenter < m_max[i];
    }
  }
  else
  {
    rval = rval && ( m_partitions[i] == 1 || ( coord >= m_min[i] && coord < m_max[i] ) );
  }

  return rval;
}

//void PartitionDescriptor::addNeighbors( const unsigned int idim,
//                                        MPI_Comm & cartcomm,
//                                        int * ncoords )
//{
//
//  if( idim == 3 )
//  {
//    bool me = true;
//    for( int i = 0; i < 3; i++ )
//    {
//      if( ncoords[i] != this->m_coords( i ))
//      {
//        me = false;
//        break;
//      }
//    }
//    if( !me )
//    {
//      int const rank = MpiWrapper::cartRank( cartcomm, ncoords );
//      m_neighbors.push_back( NeighborCommunicator( rank ) );
//    }
//  }
//  else
//  {
//    const int dim = this->m_partitions( LvArray::integerConversion< localIndex >( idim ) );
//    const bool periodic = this->m_periodic( LvArray::integerConversion< localIndex >( idim ) );
//    for( int i = -1; i < 2; i++ )
//    {
//      ncoords[idim] = this->m_coords( LvArray::integerConversion< localIndex >( idim ) ) + i;
//      bool ok = true;
//      if( periodic )
//      {
//        if( ncoords[idim] < 0 )
//          ncoords[idim] = dim - 1;
//        else if( ncoords[idim] >= dim )
//          ncoords[idim] = 0;
//      }
//      else
//      {
//        ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
//      }
//      if( ok )
//      {
//        addNeighbors( idim + 1, cartcomm, ncoords );
//      }
//    }
//  }
//}

void PartitionDescriptor::setSizes( real64 const ( &min )[3],
                                    real64 const ( &max )[3] )
{

  {
    //get size of problem and decomposition
    m_size = MpiWrapper::commSize( MPI_COMM_GEOSX );

    //check to make sure our dimensions agree
    {
      int check = 1;
      for( int i = 0; i < 3; i++ )
      {
        check *= this->m_partitions( i );
      }
      GEOS_ERROR_IF_NE( check, m_size );
    }

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOSX, 3, m_partitions.data(), m_periodic.data(), reorder, &cartcomm );
    }
    m_rank = MpiWrapper::commRank( cartcomm );
    MpiWrapper::cartCoords( cartcomm, m_rank, 3, m_coords.data());

    //add neighbors
//    {
//      int ncoords[3];
//      m_neighbors.clear();
//      addNeighbors( 0, cartcomm, ncoords );
//    }

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
  for( int i = 0; i < 3; ++i )
  {
    const int nloc = m_partitions( i ) - 1;
    const localIndex nlocl = static_cast< localIndex >(nloc);
    if( m_partitionLocations[i].empty() )
    {
      // the default "even" spacing
      m_blockSize[ i ] /= m_partitions( i );
      m_min[ i ] += m_coords( i ) * m_blockSize[ i ];
      m_max[ i ] = min[ i ] + (m_coords( i ) + 1) * m_blockSize[ i ];

      m_partitionLocations[i].resize( nlocl );
      for( localIndex j = 0; j < m_partitionLocations[ i ].size(); ++j )
      {
        m_partitionLocations[ i ][ j ] = (j+1) * m_blockSize[ i ];
      }
    }
    else if( nlocl == m_partitionLocations[i].size() )
    {
      const int parIndex = m_coords[i];
      if( parIndex == 0 )
      {
        m_min[i] = min[i];
        m_max[i] = m_partitionLocations[i][parIndex];
      }
      else if( parIndex == nloc )
      {
        m_min[i] = m_partitionLocations[i][parIndex-1];
        m_max[i] = max[i];
      }
      else
      {
        m_min[i] = m_partitionLocations[i][parIndex-1];
        m_max[i] = m_partitionLocations[i][parIndex];
      }
    }
    else
    {
      GEOS_ERROR( "PartitionDescriptor::setSizes(): number of partition locations does not equal number of partitions - 1\n" );
    }
  }
}

} // end of namespace geos
