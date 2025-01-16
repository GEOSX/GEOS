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
 * @file DomainPartition.cpp
 */

#include "DomainPartition.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"

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

  this->registerWrapper< SpatialPartition, PartitionBase >( keys::partitionManager ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( false );

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
  PartitionBase & partition1 = getReference< PartitionBase >( keys::partitionManager );
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(partition1);

  const std::set< int > metisNeighborList = partition.getMetisNeighborList();
  if( metisNeighborList.empty() )
  {

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOS, 3, partition.getPartitions().data(), partition.m_Periodic.data(), reorder, &cartcomm );
      GEOS_ERROR_IF( cartcomm == MPI_COMM_NULL, "Fail to run MPI_Cart_create and establish communications" );
    }
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
    int nsdof = 3;

    MpiWrapper::cartCoords( cartcomm, rank, nsdof, partition.m_coords.data() );

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
    MpiWrapper::iSend( firstNeighborRanks.toView(), m_neighbors[ i ].neighborRank(), neighborsTag, MPI_COMM_GEOS, &requests[ i ] );
  }

  // This set will contain the second (neighbor of) neighbors ranks.
  std::set< int > secondNeighborRanks;

  array1d< int > neighborOfNeighborRanks;
  for( std::size_t i = 0; i < m_neighbors.size(); ++i )
  {
    MpiWrapper::recv( neighborOfNeighborRanks, m_neighbors[ i ].neighborRank(), neighborsTag, MPI_COMM_GEOS, MPI_STATUS_IGNORE );

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
  PartitionBase & partition1 = getReference< PartitionBase >( keys::partitionManager );
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(partition1);

  if( idim == nsdof )
  {
    bool me = true;
    for( int i = 0; i < nsdof; i++ )
    {
      if( ncoords[i] != partition.m_coords( i ))
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
    const bool periodic = partition.m_Periodic( LvArray::integerConversion< localIndex >( idim ));
    for( int i = -1; i < 2; i++ )
    {
      ncoords[idim] = partition.m_coords( LvArray::integerConversion< localIndex >( idim )) + i;
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

void DomainPartition::outputPartitionInformation() const
{
  using stringutilities::addCommaSeparators;

  struct RankMeshStats
  {
    // Table rows will follow this enum ordering
    enum StatIndex
    {
      Node = 0, Edge, Face, Elem, Count
    };
    std::array< globalIndex, 4 > localCount = {};
    std::array< globalIndex, 4 > ghostCount = {};
    std::array< double, 4 > ratio = {};
  };

  auto fillStats = []( RankMeshStats & stat,
                       RankMeshStats::StatIndex statIndex,
                       ObjectManagerBase const & objectManager )
  {
    stat.localCount[ statIndex ] += objectManager.getNumberOfLocalIndices();
    stat.ghostCount[ statIndex ] += objectManager.getNumberOfGhosts();
  };

  auto computeRatios = []( RankMeshStats & stat )
  {
    for( size_t i = 0; i < RankMeshStats::StatIndex::Count; ++i )
    {
      stat.ratio[i] = stat.localCount[i] + stat.ghostCount[i] == 0 ? 0 :
                      (double)stat.localCount[i] / (double)(stat.localCount[i] + stat.ghostCount[i]);
    }
  };

  auto addLocalGhostRow = []( TableData & tableData, RankMeshStats const & stat, string_view heading )
  {
    tableData.addRow( heading,
                      addCommaSeparators( stat.localCount[0] ), addCommaSeparators( stat.ghostCount[0] ),
                      addCommaSeparators( stat.localCount[0] + stat.ghostCount[0] ),
                      addCommaSeparators( stat.localCount[1] ), addCommaSeparators( stat.ghostCount[1] ),
                      addCommaSeparators( stat.localCount[1] + stat.ghostCount[1] ),
                      addCommaSeparators( stat.localCount[2] ), addCommaSeparators( stat.ghostCount[2] ),
                      addCommaSeparators( stat.localCount[2] + stat.ghostCount[2] ),
                      addCommaSeparators( stat.localCount[3] ), addCommaSeparators( stat.ghostCount[3] ),
                      addCommaSeparators( stat.localCount[3] + stat.ghostCount[3] ) );
  };

  auto addSummaryRow = []( TableData & tableData, std::array< double, 4 > stats, string_view heading )
  {
    tableData.addRow( heading,
                      CellType::MergeNext, CellType::MergeNext, stats[0],
                      CellType::MergeNext, CellType::MergeNext, stats[1],
                      CellType::MergeNext, CellType::MergeNext, stats[2],
                      CellType::MergeNext, CellType::MergeNext, stats[3] );
  };

  GEOS_LOG_RANK_0( "MPI Partitioning information:" );

  forMeshBodies( [&]( MeshBody const & meshBody )
  {
    meshBody.getMeshLevels().forSubGroupsIndex< MeshLevel >( [&]( int const level, MeshLevel const & meshLevel )
    {
      if( level!=0 )
      {
        // formatting is done on rank 0
        std::vector< RankMeshStats > allRankStats;
        allRankStats.resize( MpiWrapper::commSize() );

        { // Compute stats of the current rank, then gather it on rank 0
          RankMeshStats rankStats{};
          fillStats( rankStats, RankMeshStats::Node, meshLevel.getNodeManager() );
          fillStats( rankStats, RankMeshStats::Edge, meshLevel.getEdgeManager() );
          fillStats( rankStats, RankMeshStats::Face, meshLevel.getFaceManager() );

          meshLevel.getElemManager().forElementSubRegions< CellElementSubRegion >(
            [&]( CellElementSubRegion const & subRegion )
          {
            fillStats( rankStats, RankMeshStats::Elem, subRegion );
          } );

          computeRatios( rankStats );

          MpiWrapper::gather( rankStats, allRankStats, 0 );
        }


        if( MpiWrapper::commRank() == 0 )
        {
          TableLayout const layout( "Mesh partitioning over ranks",
                                    {TableLayout::Column()
                                       .setName( "Ranks" ),
                                     TableLayout::Column()
                                       .setName( "Nodes" )
                                       .addSubColumns( {  "Local", "Ghost", "Total" } ),
                                     TableLayout::Column()
                                       .setName( "Edges" )
                                       .addSubColumns( {  "Local", "Ghost", "Total" } ),
                                     TableLayout::Column()
                                       .setName( "Faces" )
                                       .addSubColumns( {  "Local", "Ghost", "Total" } ),
                                     TableLayout::Column()
                                       .setName( "Elems" )
                                       .addSubColumns( {  "Local", "Ghost", "Total" } )} );
          TableData tableData;

          for( int rankId = 0; rankId < MpiWrapper::commSize(); ++rankId )
          {
            if( rankId == 1 )
              tableData.addSeparator();

            addLocalGhostRow( tableData, allRankStats[rankId], std::to_string( rankId ) );
          }

          RankMeshStats sumStats{};
          RankMeshStats minStats{};
          RankMeshStats maxStats{};

          for( size_t statId = 0; statId < RankMeshStats::Count; ++statId )
          {
            minStats.localCount[statId] = std::numeric_limits< globalIndex >::max();
            minStats.ghostCount[statId] = std::numeric_limits< globalIndex >::max();
            minStats.ratio[statId] = std::numeric_limits< double >::max();

            maxStats.localCount[statId] = std::numeric_limits< globalIndex >::min();
            maxStats.ghostCount[statId] = std::numeric_limits< globalIndex >::min();
            maxStats.ratio[statId] = std::numeric_limits< double >::min();
          }

          for( int rankId = 0; rankId < MpiWrapper::commSize(); ++rankId )
          {
            for( size_t statId = 0; statId < RankMeshStats::Count; ++statId )
            {
              sumStats.localCount[statId] += allRankStats[rankId].localCount[statId];
              sumStats.ghostCount[statId] += allRankStats[rankId].ghostCount[statId];

              minStats.localCount[statId] = std::min( minStats.localCount[statId], allRankStats[rankId].localCount[statId] );
              minStats.ghostCount[statId] = std::min( minStats.ghostCount[statId], allRankStats[rankId].ghostCount[statId] );
              minStats.ratio[statId] = std::min( minStats.ratio[statId], allRankStats[rankId].ratio[statId] );

              maxStats.localCount[statId] = std::max( maxStats.localCount[statId], allRankStats[rankId].localCount[statId] );
              maxStats.ghostCount[statId] = std::max( maxStats.ghostCount[statId], allRankStats[rankId].ghostCount[statId] );
              maxStats.ratio[statId] = std::max( maxStats.ratio[statId], allRankStats[rankId].ratio[statId] );
            }
          }

          tableData.addSeparator();
          addLocalGhostRow( tableData, sumStats, "sum" );
          addLocalGhostRow( tableData, minStats, "min" );
          addLocalGhostRow( tableData, maxStats, "max" );

          std::array< double, 4 > localTotalMinRatio;
          std::array< double, 4 > localTotalMaxRatio;

          for( size_t statId = 0; statId < RankMeshStats::Count; ++statId )
          {
            localTotalMinRatio[statId] = std::numeric_limits< globalIndex >::max();
            localTotalMaxRatio[statId] = std::numeric_limits< globalIndex >::min();
          }

          for( size_t statId = 0; statId < RankMeshStats::Count; ++statId )
          {
            localTotalMinRatio[statId] = std::min( localTotalMinRatio[statId], minStats.ratio[statId] );
            localTotalMaxRatio[statId] = std::max( localTotalMinRatio[statId], maxStats.ratio[statId] );
          }
          tableData.addSeparator();
          addSummaryRow( tableData, localTotalMinRatio, "min(local/total)" );
          addSummaryRow( tableData, localTotalMaxRatio, "max(local/total)" );

          TableTextFormatter logPartition( layout );
          GEOS_LOG_RANK_0( logPartition.toString( tableData ));
        }

      }
    } );
  } );

}

} /* namespace geos */
