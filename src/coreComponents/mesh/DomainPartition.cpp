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
 * @file DomainPartition.cpp
 */

#include "DomainPartition.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"



namespace geos
{
using namespace dataRepository;

DomainPartition::DomainPartition( string const & name,
                                  Group * const parent ):
  Group( name, parent )
{
  this->registerWrapper( "Neighbors", &m_neighbors ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( false );

  // this->registerWrapper< SpatialPartition, PartitionBase >( keys::partitionManager ).
  //   setRestartFlags( RestartFlags::WRITE ). //RestartFlags::NO_WRITE ).
  //   setSizedFromParent( false );

  registerGroup< SpatialPartition >( groupKeys.partitionManager );
  registerGroup( groupKeys.meshBodies );
  registerGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager );
}


DomainPartition::~DomainPartition()
{}

void DomainPartition::initializationOrder( string_array & order )
{
  SortedArray< string > usedNames;
  {
    order.emplace_back( string( groupKeysStruct::constitutiveManagerString() ) );
    usedNames.insert( groupKeysStruct::constitutiveManagerString() );
  }

  {
    order.emplace_back( string( groupKeysStruct::meshBodiesString() ) );
    usedNames.insert( groupKeysStruct::meshBodiesString() );
  }


  for( auto const & subGroup : this->getSubGroups() )
  {
    if( usedNames.count( subGroup.first ) == 0 )
    {
      order.emplace_back( subGroup.first );
    }
  }
}

void DomainPartition::setupBaseLevelMeshGlobalInfo()
{
  GEOS_MARK_FUNCTION;

#if defined(GEOS_USE_MPI)
  // PartitionBase & partition1 = getReference< PartitionBase >( groupKeys.partitionManager ); //keys::partitionManager );
  // SpatialPartition & partition = dynamic_cast< SpatialPartition & >(partition1);
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >( getGroup( groupKeys.partitionManager ) );

  const std::set< int > metisNeighborList = partition.getMetisNeighborList();
  if( metisNeighborList.empty() )
  {
    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOSX, 3, partition.getPartitions().data(), partition.getPeriodic().data(), reorder, &cartcomm );
      GEOS_ERROR_IF( cartcomm == MPI_COMM_NULL, "Fail to run MPI_Cart_create and establish communications" );
    }
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int nsdof = 3;

    MpiWrapper::cartCoords( cartcomm, rank, nsdof, partition.getCoords().data() );

    int ncoords[3];
    addNeighbors( 0, cartcomm, ncoords );

    MpiWrapper::commFree( cartcomm );
  }
  else
  {
    for( integer const neighborRank : metisNeighborList )
    {
      m_neighbors.emplace_back( neighborRank );
    }
  }

  // Create an array of the first neighbors.
  array1d< int > firstNeighborRanks;
  for( NeighborCommunicator const & neighbor : m_neighbors )
  {
    firstNeighborRanks.emplace_back( neighbor.neighborRank() );
  }

  int neighborsTag = 54;

  // Send this list of neighbors to all neighbors.
  std::vector< MPI_Request > requests( m_neighbors.size(), MPI_REQUEST_NULL );

  for( std::size_t i = 0; i < m_neighbors.size(); ++i )
  {
    MpiWrapper::iSend( firstNeighborRanks.toView(), m_neighbors[ i ].neighborRank(), neighborsTag, MPI_COMM_GEOSX, &requests[ i ] );
  }

  // This set will contain the second (neighbor of) neighbors ranks.
  std::set< int > secondNeighborRanks;

  array1d< int > neighborOfNeighborRanks;
  for( std::size_t i = 0; i < m_neighbors.size(); ++i )
  {
    MpiWrapper::recv( neighborOfNeighborRanks, m_neighbors[ i ].neighborRank(), neighborsTag, MPI_COMM_GEOSX, MPI_STATUS_IGNORE );

    // Insert the neighbors of the current neighbor into the set of second neighbors.
    secondNeighborRanks.insert( neighborOfNeighborRanks.begin(), neighborOfNeighborRanks.end() );
  }

  // Remove yourself and all the first neighbors from the second neighbors.
  secondNeighborRanks.erase( MpiWrapper::commRank() );
  for( NeighborCommunicator const & neighbor : m_neighbors )
  {
    secondNeighborRanks.erase( neighbor.neighborRank() );
  }

  for( integer const neighborRank : secondNeighborRanks )
  {
    m_neighbors.emplace_back( neighborRank );
  }

  MpiWrapper::waitAll( requests.size(), requests.data(), MPI_STATUSES_IGNORE );

#endif

  forMeshBodies( [&]( MeshBody & meshBody )
  {
    if( !meshBody.hasParticles() ) // Currently, particle-based mesh bodies do not construct their
                                   // own domain decomposition. MPM borrows that of the grid.
    {
      MeshLevel & meshLevel = meshBody.getBaseDiscretization();

      NodeManager & nodeManager = meshLevel.getNodeManager();
      FaceManager & faceManager = meshLevel.getFaceManager();
      EdgeManager & edgeManager = meshLevel.getEdgeManager();

      nodeManager.setMaxGlobalIndex();
      for( NeighborCommunicator const & neighbor : m_neighbors )
      {
        neighbor.addNeighborGroupToMesh( meshLevel );
      }

      CommunicationTools::getInstance().assignGlobalIndices( faceManager,
                                                             nodeManager,
                                                             m_neighbors );

      CommunicationTools::getInstance().assignGlobalIndices( edgeManager,
                                                             nodeManager,
                                                             m_neighbors );

      // CC: Add check if there even are any periodic boundaries?
      // TODO: If GEOS_USE_MPI flag is set off this will throw compile error since partition is not included
      partition.setPeriodicDomainBoundaryObjects( meshBody,
                                                  nodeManager,
                                                  edgeManager,
                                                  faceManager );


      CommunicationTools::getInstance().findMatchedPartitionBoundaryObjects( faceManager,
                                                                             m_neighbors );

      CommunicationTools::getInstance().findMatchedPartitionBoundaryObjects( edgeManager,
                                                                             m_neighbors );

      // w.r.t. edges and faces, finding the matching nodes between partitions is a bit trickier.
      // Because for contact mechanics and fractures, some nodes can be collocated.
      // And the fracture elements will point to those nodes.
      // While they are not the _same_ nodes (which is the criterion for edges and faces),
      // we still want those collocated nodes to be exchanged between the ranks.
      // This is why we gather some additional information: what are those collocated nodes
      // and also what are the nodes that we require but are not present on the current rank!
      std::set< std::set< globalIndex > > collocatedNodesBuckets;
      std::set< globalIndex > requestedNodes;
      meshLevel.getElemManager().forElementSubRegions< FaceElementSubRegion >(
        [&, g2l = &nodeManager.globalToLocalMap()]( FaceElementSubRegion const & subRegion )
      {
        ArrayOfArraysView< array1d< globalIndex > const > const buckets = subRegion.get2dElemToCollocatedNodesBuckets();
        for( localIndex e2d = 0; e2d < buckets.size(); ++e2d )
        {
          for( integer ni = 0; ni < buckets.sizeOfArray( e2d ); ++ni )
          {
            array1d< globalIndex > const & bucket = buckets( e2d, ni );
            std::set< globalIndex > tmp( bucket.begin(), bucket.end() );
            collocatedNodesBuckets.insert( tmp );

            for( globalIndex const gni: bucket )
            {
              auto const it = g2l->find( gni );
              if( it == g2l->cend() )
              {
                requestedNodes.insert( gni );
              }
            }
          }
        }
      } );

      CommunicationTools::getInstance().findMatchedPartitionBoundaryNodes( nodeManager,
                                                                           m_neighbors,
                                                                           collocatedNodesBuckets,
                                                                           requestedNodes );
    }
  } );
}


void DomainPartition::setupCommunications( bool use_nonblocking )
{
  forMeshBodies( [&]( MeshBody & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel & meshLevel )
    {
      if( !meshBody.hasParticles() ) // Currently, particle-based mesh bodies do not construct their
                                     // own domain decomposition. MPM borrows that of the grid.
      {
        if( meshLevel.getName() == MeshBody::groupStructKeys::baseDiscretizationString() )
        {
          NodeManager & nodeManager = meshLevel.getNodeManager();
          FaceManager & faceManager = meshLevel.getFaceManager();

          CommunicationTools::getInstance().setupGhosts( meshLevel, m_neighbors, use_nonblocking );
          faceManager.sortAllFaceNodes( nodeManager, meshLevel.getElemManager() );
          faceManager.computeGeometry( nodeManager );
        }
        else if( !meshLevel.isShallowCopyOf( meshBody.getMeshLevels().getGroup< MeshLevel >( 0 )) )
        {
          for( NeighborCommunicator const & neighbor : m_neighbors )
          {
            neighbor.addNeighborGroupToMesh( meshLevel );
          }
          NodeManager & nodeManager = meshLevel.getNodeManager();
          FaceManager & faceManager = meshLevel.getFaceManager();

          CommunicationTools::getInstance().findMatchedPartitionBoundaryObjects( faceManager, m_neighbors );
          CommunicationTools::getInstance().findMatchedPartitionBoundaryObjects( nodeManager, m_neighbors );
          CommunicationTools::getInstance().setupGhosts( meshLevel, m_neighbors, use_nonblocking );
        }
        else
        {
          GEOS_LOG_LEVEL_RANK_0( 3, "No communication setup is needed since it is a shallow copy of the base discretization." );
        }
      }
    } );
  } );
}

void DomainPartition::addNeighbors( const unsigned int idim,
                                    MPI_Comm & cartcomm,
                                    int * ncoords )
{  
  // PartitionBase & partition1 = getReference< PartitionBase >( keys::partitionManager );
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >( getGroup( groupKeys.partitionManager ) ); // partition1);

  if( idim == nsdof )
  {
    bool me = true;
    for( int i = 0; i < nsdof; i++ )
    {
      if( ncoords[i] != partition.getCoords()( i ))
      {
        me = false;
        break;
      }
    }
    int const neighborRank = MpiWrapper::cartRank( cartcomm, ncoords );
    if( !me && !std::any_of( m_neighbors.begin(), m_neighbors.end(), [=]( NeighborCommunicator const & nn ) { return nn.neighborRank( ) == neighborRank; } ) )
    {
      m_neighbors.emplace_back( NeighborCommunicator( neighborRank ) );
    }
  }
  else
  {
    const int dim = partition.getPartitions()( LvArray::integerConversion< localIndex >( idim ));
    const bool periodic = partition.getPeriodic()( LvArray::integerConversion< localIndex >( idim ));
    for( int i = -1; i < 2; i++ )
    {
      ncoords[idim] = partition.getCoords()( LvArray::integerConversion< localIndex >( idim )) + i;
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

} /* namespace geos */
