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

#include "NewGhosting.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_RowMatrixTransposer.h>


#include <NamedType/named_type.hpp>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDataSetSurfaceFilter.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include <algorithm>
#include <utility>

namespace geos
{

template< typename OUTPUT, typename INPUT >
inline LVARRAY_HOST_DEVICE
OUTPUT intConv( INPUT input )
{
  return LvArray::integerConversion< OUTPUT >( input );
}

}  // end of namespace


namespace geos::ghosting
{

using NodeLocIdx = fluent::NamedType< localIndex, struct NodeLocIdxTag, fluent::Comparable, fluent::Printable >;
using NodeGlbIdx = fluent::NamedType< globalIndex, struct NodeGlbIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeLocIdx = fluent::NamedType< localIndex, struct EdgeLocIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeGlbIdx = fluent::NamedType< globalIndex, struct EdgeGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::Subtractable, fluent::PreIncrementable >;
using FaceLocIdx = fluent::NamedType< localIndex, struct FaceLocIdxTag, fluent::Comparable, fluent::Printable >;
using FaceGlbIdx = fluent::NamedType< globalIndex, struct FaceGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::Subtractable, fluent::PreIncrementable >;
using CellLocIdx = fluent::NamedType< localIndex, struct CellLocIdxTag, fluent::Comparable, fluent::Printable >;
using CellGlbIdx = fluent::NamedType< globalIndex, struct CellGlbIdxTag, fluent::Comparable, fluent::Printable >;

using MpiRank = fluent::NamedType< int, struct MpiRankTag, fluent::Comparable, fluent::Printable, fluent::Addable >;

EdgeGlbIdx operator "" _egi( unsigned long long int i )
{
  return EdgeGlbIdx{ EdgeGlbIdx::UnderlyingType( i ) };
}

FaceGlbIdx operator "" _fgi( unsigned long long int i )
{
  return FaceGlbIdx{ FaceGlbIdx::UnderlyingType( i ) };
}

MpiRank operator "" _mpi( unsigned long long int i )
{
  return MpiRank{ MpiRank::UnderlyingType( i ) };
}

void to_json( json & j,
              const MpiRank & v )
{
  j = v.get();
}

void from_json( const json & j,
                MpiRank & v )
{
  v = MpiRank{ j.get< MpiRank::UnderlyingType >() };  // TODO use a `traits` instead
}

void from_json( const json & j,
                NodeGlbIdx & v )
{
  v = NodeGlbIdx{ j.get< NodeGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const NodeGlbIdx & v )
{
  j = v.get();
}

void to_json( json & j,
              const EdgeGlbIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                EdgeGlbIdx & v )
{
  v = EdgeGlbIdx{ j.get< EdgeGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const FaceGlbIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                FaceGlbIdx & v )
{
  v = FaceGlbIdx{ j.get< FaceGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const CellGlbIdx & v )
{
  j = v.get();
}

using Edge = std::tuple< NodeGlbIdx, NodeGlbIdx >;
using Face = std::vector< NodeGlbIdx >;

struct Exchange
{
  std::set< NodeGlbIdx > nodes;
  std::set< Edge > edges;
  std::set< Face > faces;
};

void to_json( json & j,
              const Exchange & v )
{
  j = json{ { "nodes", v.nodes },
            { "edges", v.edges },
            { "faces", v.faces } };
}

void from_json( const json & j,
                Exchange & v )
{
  v.nodes = j.at( "nodes" ).get< std::set< NodeGlbIdx > >();
  v.edges = j.at( "edges" ).get< std::set< Edge > >();
  v.faces = j.at( "faces" ).get< std::set< Face > >();
}

/**
 * @brief Extract the cells at the boundary of the mesh.
 * @param mesh The vtk mesh.
 * @return The vtk cell ids.
 */
std::set< vtkIdType > extractBoundaryCells( vtkSmartPointer< vtkDataSet > mesh )
{
  // TODO Better handle the boundary information, forgetting about the 3d cells and simply handling the outside shell.
  auto f = vtkDataSetSurfaceFilter::New();
  f->PassThroughCellIdsOn();
  f->PassThroughPointIdsOff();
  f->FastModeOff();

  string const originalCellsKey = "ORIGINAL_CELLS";
  f->SetOriginalCellIdsName( originalCellsKey.c_str() );
  auto boundaryMesh = vtkPolyData::New();
  f->UnstructuredGridExecute( mesh, boundaryMesh );
  vtkIdTypeArray const * originalCells = vtkIdTypeArray::FastDownCast( boundaryMesh->GetCellData()->GetArray( originalCellsKey.c_str() ) );

  std::set< vtkIdType > boundaryCellIdxs;
  for( auto i = 0; i < originalCells->GetNumberOfTuples(); ++i )
  {
    boundaryCellIdxs.insert( originalCells->GetValue( i ) );
  }

  return boundaryCellIdxs;
}

/**
 * @brief Order the nodes of the faces in a way that can be reproduced across the MPI ranks.
 * @param nodes The list of nodes as provided by the mesh.
 * @return A face with the nodes in the appropriate order
 * @details The nodes will be ordered in the following way.
 * First, we look for the lowest node index. It will become the first node.
 * Then we must pick the second node. We have to choices: just before or just after the first.
 * (we do not want to shuffle the nodes completely, we need to keep track of the order).
 * To do this, we select the nodes with the lowest index as well.
 * This also defines a direction in which we'll pick the other nodes.
 * For example, the face <tt>[2, 3, 1, 5, 9, 8]</tt> will become <tt>[1, 3, 2, 8, 9, 5]</tt>
 * because we'll start with @c 1 and then select the @c 3 over the @c 5.
 * Which defines the direction @c 2, @c 8, @c 9, @c 5.
 * @note This is the same pattern that we apply for edges.
 * Except that edges having only two nodes, it's not necessary to implement a dedicated function
 * and <tt>std::minmax</tt> is enough.
 */
Face reorderFaceNodes( std::vector< NodeGlbIdx > const & nodes )
{
  std::size_t const n = nodes.size();

  // Handles negative values of `i`.
  auto const modulo = [n]( integer const & i ) -> std::size_t
  {
    integer mod = i % n;
    if( mod < 0 )
    {
      mod += n;
    }
    return mod;
  };

  Face f;
  f.reserve( n );

  auto const it = std::min_element( nodes.cbegin(), nodes.cend() );
  std::size_t const minIdx = std::distance( nodes.cbegin(), it );
  int const increment = nodes[modulo( minIdx - 1 )] < nodes[modulo( minIdx + 1 )] ? -1 : 1;  // TODO based on increment, I can say if the face is flipped or not.
  integer i = minIdx;
  for( std::size_t count = 0; count < n; ++count, i = i + increment )
  {
    f.emplace_back( nodes.at( modulo( i ) ) );
  }

  return f;
}

Exchange buildSpecificData( vtkSmartPointer< vtkDataSet > mesh,
                            std::set< vtkIdType > const & cellIds )
{
  vtkIdTypeArray * gids = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );
  Span< vtkIdType > const globalPtIds( (vtkIdType *) gids->GetPointer( 0 ), gids->GetNumberOfTuples() );

  std::vector< Edge > tmpEdges;
  std::vector< Face > tmpFaces;

  // Pre-allocation of the temporary vectors.
  {
    std::size_t numEdges = 0;
    std::size_t numFaces = 0;
    for( vtkIdType const & c: cellIds )
    {
      vtkCell * cell = mesh->GetCell( c );
      numEdges += cell->GetNumberOfEdges();
      numFaces += cell->GetNumberOfFaces();
    }
    tmpEdges.reserve( numEdges );
    tmpFaces.reserve( numFaces );
  }

  // Filling the temporary.
  for( auto c = 0; c < mesh->GetNumberOfCells(); ++c )
  {
    vtkCell * cell = mesh->GetCell( c );

    for( auto e = 0; e < cell->GetNumberOfEdges(); ++e )
    {
      vtkCell * edge = cell->GetEdge( e );
      vtkIdType const ln0 = edge->GetPointId( 0 );
      vtkIdType const ln1 = edge->GetPointId( 1 );
      vtkIdType const gn0 = globalPtIds[ln0];
      vtkIdType const gn1 = globalPtIds[ln1];
      tmpEdges.emplace_back( std::minmax( { NodeGlbIdx{ gn0 }, NodeGlbIdx{ gn1 } } ) );
    }

    for( auto f = 0; f < cell->GetNumberOfFaces(); ++f )
    {
      vtkCell * face = cell->GetFace( f );
      vtkIdList * pids = face->GetPointIds();
      std::vector< NodeGlbIdx > nodes( pids->GetNumberOfIds() );
      for( std::size_t i = 0; i < nodes.size(); ++i )
      {
        vtkIdType const lni = face->GetPointId( i );
        vtkIdType const gni = globalPtIds[lni];
        nodes[i] = NodeGlbIdx{ gni };
      }
      tmpFaces.emplace_back( reorderFaceNodes( nodes ) );
    }
  }

  std::set< NodeGlbIdx > nodes;
  std::transform( std::cbegin( globalPtIds ), std::cend( globalPtIds ),
                  std::inserter( nodes, std::end( nodes ) ), []( vtkIdType const & id )
                  { return NodeGlbIdx{ id }; } );

  // Removing the duplicates by copying into a `std::set`.
  std::set< Edge > edges{ tmpEdges.cbegin(), tmpEdges.cend() };  // SortedArray requires the input to be sorted already.
  std::set< Face > faces{ tmpFaces.cbegin(), tmpFaces.cend() };

  return { std::move( nodes ), std::move( edges ), std::move( faces ) };
}

array1d< std::uint8_t > convertExchange( Exchange const & exchange )
{
  std::vector< std::uint8_t > const tmp = json::to_cbor( json( exchange ) );
  array1d< std::uint8_t > result;
  result.reserve( tmp.size() );
  for( std::uint8_t const & t: tmp )
  {
    result.emplace_back( t );
  }
  return result;
}

Exchange convertExchange( array1d< std::uint8_t > const & input )
{
  std::vector< std::uint8_t > const tmp( std::cbegin( input ), std::cend( input ) );
  return json::from_cbor( tmp, false ).template get< Exchange >();
}

Exchange buildFullData( vtkSmartPointer< vtkDataSet > mesh )
{
  std::set< vtkIdType > cellIds;
  for( vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i )
  {
    cellIds.insert( cellIds.end(), i );
  }

  return buildSpecificData( mesh, cellIds );
}


Exchange buildExchangeData( vtkSmartPointer< vtkDataSet > mesh )
{
  return buildSpecificData( mesh, extractBoundaryCells( mesh ) );
}

struct Buckets
{
  std::map< std::set< MpiRank >, std::set< NodeGlbIdx > > nodes;
  std::map< std::set< MpiRank >, std::set< Edge > > edges;
  std::map< std::set< MpiRank >, std::set< Face > > faces;
};

//void to_json( json & j,
//              const Buckets & v )
//{
//  j = json{ { "nodes", v.nodes }, { "edges", v.edges }, { "faces", v.faces }  };
//}

/**
 * @brief
 * @tparam T Typically NodeGloIdx or a container of NodeGlbIdx (for edges and faces)
 * @param exchanged
 * @param curRank
 * @param neighbors
 * @return
 */
template< typename T >
std::map< std::set< MpiRank >, std::set< T > >
buildIntersectionBuckets( std::map< MpiRank, std::set< T > const & > const & exchanged,
                          MpiRank curRank,
                          std::set< MpiRank > const & neighbors )
{
  std::map< T, std::set< MpiRank > > counts;  // TODO Use better intersection algorithms?
  // We "register" all the edges of the current rank: they are the only one we're interested in.
  for( T const & node: exchanged.at( curRank ) )
  {
    counts.emplace_hint( counts.end(), node, std::set< MpiRank >{ curRank } );
  }

  // We now loop on the neighbor edges.
  // If a neighbor has an edge in common with the current rank, they we store it.
  for( MpiRank const & neighborRank: neighbors )  // This does not include the current rank.
  {
    for( T const & node: exchanged.at( neighborRank ) )
    {
      auto it = counts.find( node );
      if( it != counts.cend() )  // TODO Extract `counts.cend()` out of the loop.
      {
        it->second.insert( neighborRank );
      }
    }
  }

  std::map< std::set< MpiRank >, std::set< T > > nodeBuckets;
  for( auto const & [node, ranks]: counts )
  {
    if( ranks.find( curRank ) != ranks.cend() )
    {
      nodeBuckets[ranks].insert( node );
    }
  }
  return nodeBuckets;
}


/**
 * @brief Compute the intersection between for the ranks based on the information @p exchanged.
 * @param exchanged Geometrical information provided by the ranks (including the current rank).
 * @param curRank The current MPI rank.
 * @param neighbors excluding current rank
 * @return The intersection buckets.
 */
Buckets buildIntersectionBuckets( std::map< MpiRank, Exchange > const & exchanged,
                                  MpiRank curRank,
                                  std::set< MpiRank > const & neighbors )
{
  std::map< MpiRank, std::set< NodeGlbIdx > const & > nodeInfo;
  for( auto const & [rank, exchange]: exchanged )
  {
    nodeInfo.emplace( rank, exchange.nodes );
  }
  std::map< std::set< MpiRank >, std::set< NodeGlbIdx > > const nodeBuckets = buildIntersectionBuckets( nodeInfo, curRank, neighbors );

  std::map< MpiRank, std::set< Edge > const & > edgeInfo;
  for( auto const & [rank, exchange]: exchanged )
  {
    edgeInfo.emplace( rank, exchange.edges );
  }
  std::map< std::set< MpiRank >, std::set< Edge > > const edgeBuckets = buildIntersectionBuckets( edgeInfo, curRank, neighbors );

  // For faces, the algorithm can be a tad simpler because faces can be shared by at most 2 ranks.
  std::map< std::set< MpiRank >, std::set< Face > > faceBuckets;
  std::set< Face > curFaces = exchanged.at( curRank ).faces;
  for( MpiRank const & neighborRank: neighbors )  // This does not include the current rank.
  {
    for( Face const & face: exchanged.at( neighborRank ).faces )
    {
      auto it = curFaces.find( face );
      if( it != curFaces.cend() )
      {
        faceBuckets[{ curRank, neighborRank }].insert( *it );
        curFaces.erase( it );
      }
    }
  }
  faceBuckets[{ curRank }] = curFaces;

//  // Checking if neighbors is too wide...  // TODO do we care here?
//  std::set< MpiRank > usefulNeighbors;
//  for( auto const & [ranks, nodes]: nodeBuckets )
//  {
//    if( not nodes.empty() )
//    {
//      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
//    }
//  }
//  for( auto const & [ranks, edges]: edgeBuckets )
//  {
//    if( not edges.empty() )
//    {
//      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
//    }
//  }
//  std::vector< MpiRank > uselessNeighbors;
//  std::set_difference( neighbors.cbegin(), neighbors.cend(), usefulNeighbors.cbegin(), usefulNeighbors.cend(), std::back_inserter( uselessNeighbors ) );
//  // TODO... Remove the neighbors?

  return { nodeBuckets, edgeBuckets, faceBuckets };
}


// TODO Duplicated
std::map< MpiRank, Exchange > exchange( int commId,
                                        std::vector< NeighborCommunicator > & neighbors,
                                        Exchange const & data )
{
  MPI_iCommData commData( commId );
  integer const numNeighbors = LvArray::integerConversion< integer >( neighbors.size() );
  commData.resize( numNeighbors );


  array1d< std::uint8_t > const cv = convertExchange( data );

  for( integer i = 0; i < numNeighbors; ++i )
  {
    neighbors[i].mpiISendReceiveSizes( cv,
                                       commData.mpiSendBufferSizeRequest( i ),
                                       commData.mpiRecvBufferSizeRequest( i ),
                                       commId,
                                       MPI_COMM_GEOSX );
  }

  MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
  MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus() );

  array1d< array1d< std::uint8_t > > rawExchanged( neighbors.size() );

  for( integer i = 0; i < numNeighbors; ++i )
  {
    neighbors[i].mpiISendReceiveData( cv,
                                      commData.mpiSendBufferRequest( i ),
                                      rawExchanged[i],
                                      commData.mpiRecvBufferRequest( i ),
                                      commId,
                                      MPI_COMM_GEOSX );
  }
  MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );
  MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus() );

  std::map< MpiRank, Exchange > output;
  for( auto i = 0; i < numNeighbors; ++i )
  {
    output[MpiRank{ neighbors[i].neighborRank() }] = convertExchange( rawExchanged[i] );
  }
  return output;
}


/**
 * @brief Estimate an upper bound of the serialized size of the size @p buckets.
 * @param buckets
 * @return Size in bytes.
 * @details @c MPI_Scan requires a fixed buffer size.
 * An a-priori estimation of the maximum size of the buffer leads to crazy amounts of memory
 * because of the combinatorial pattern of the rank intersections.
 * Therefore, we use an initial MPI communication to get a better estimation.
 * @node We don't care about the @c nodes because their global numbering is already done.
 */
std::size_t buildMaxBufferSize( Buckets const & buckets )  // TODO give the size buckets instead?
{
  auto const f = []( auto const & bucket ) -> std::size_t
  {
    std::size_t size = std::size( bucket );
    for( auto const & [ranks, _]: bucket )
    {
      size += std::size( ranks ) + 1 + 1;  // One `+1` for the size of the ranks set, the other one for the offset
    }
    return size;
  };

  return MpiWrapper::sum( f( buckets.edges ) + f( buckets.faces ) );  // Add the ScannedOffsets::{edge, face}Restart
}

struct BucketSizes
{
  using mapping = std::map< std::set< MpiRank >, localIndex >;
  mapping edges;
  mapping faces;
};

void to_json( json & j,
              const BucketSizes & v )
{
  j = json{ { "edges", v.edges },
            { "faces", v.faces } };
}

void from_json( const json & j,
                BucketSizes & v )
{
  v.edges = j.at( "edges" ).get< BucketSizes::mapping >();
  v.faces = j.at( "faces" ).get< BucketSizes::mapping >();
}

struct BucketOffsets
{
  std::map< std::set< MpiRank >, EdgeGlbIdx > edges;
  std::map< std::set< MpiRank >, FaceGlbIdx > faces;
};

void to_json( json & j,
              const BucketOffsets & v )
{
  j = json{ { "edges", v.edges },
            { "faces", v.faces } };
}

void from_json( const json & j,
                BucketOffsets & v )
{
  v.edges = j.at( "edges" ).get< std::map< std::set< MpiRank >, EdgeGlbIdx > >();
  v.faces = j.at( "faces" ).get< std::map< std::set< MpiRank >, FaceGlbIdx > >();
}

/**
 * @brief Compute the sizes of the intersection buckets.
 * @param buckets The intersection buckets.
 * @return The buckets of sizes.
 */
BucketSizes getBucketSize( Buckets const & buckets )
{
  BucketSizes output;

  for( auto const & [ranks, edges]: buckets.edges )
  {
    output.edges.emplace_hint( output.edges.end(), ranks, std::size( edges ) );
  }

  for( auto const & [ranks, faces]: buckets.faces )
  {
    output.faces.emplace_hint( output.faces.end(), ranks, std::size( faces ) );
  }

  return output;
}

/**
 * @brief
 * @tparam GLB_IDX The global index type of the geometrical quantity considered (typically @c EdgeGlbIdx or @c FaceGlbIdx).
 * @param sizes The sizes of the intersection bucket (from the current MPI rank).
 * @param offsets The offsets for each intersection bucket (from the MPI_Scan process, i.e. the previous ranks).
 * @param curRank Current MPI rank.
 * @return The offset buckets, updated with the sizes of the current MPI rank.
 */
template< typename GLB_IDX >
std::map< std::set< MpiRank >, GLB_IDX > updateBucketOffsets( std::map< std::set< MpiRank >, localIndex > const & sizes,
                                                              std::map< std::set< MpiRank >, GLB_IDX > const & offsets,
                                                              MpiRank curRank )
{
  std::map< std::set< MpiRank >, GLB_IDX > reducedOffsets;

  // We only keep the offsets that are still relevant to the current and higher ranks.
  // Note that `reducedOffsets` will be used by the _current_ rank too.
  // Therefore, we'll send information that may not be useful to higher ranks.
  // This is surely acceptable because the information is tiny.
  for( auto const & [ranks, offset]: offsets )
  {
    MpiRank const maxConcerned = *std::max_element( ranks.cbegin(), ranks.cend() );
    // TODO invert
//    if( maxConcerned >= curRank )
//    {
//      reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, offset );
//    }
    if( maxConcerned < curRank )
    {
      continue;
    }
    reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, offset );
  }

  // Add the offsets associated to the new size buckets from the current rank.
  GLB_IDX nextOffset{ 0 };
  for( auto const & [ranks, size]: sizes )
  {
    auto const it = reducedOffsets.find( ranks );
    if( it == reducedOffsets.end() )
    {
      reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, nextOffset );
      nextOffset += GLB_IDX{ size };
    }
    else
    {
      nextOffset = it->second + GLB_IDX{ size };  // Define the new offset from the last
    }
  }

  // Add an extra entry based for the following rank
  reducedOffsets.emplace_hint( reducedOffsets.end(), std::set< MpiRank >{ curRank + 1_mpi }, nextOffset );

  return reducedOffsets;
}

std::vector< std::uint8_t > serialize( BucketSizes const & sizes )
{
  return json::to_cbor( json( sizes ) );
}

std::vector< std::uint8_t > serialize( BucketOffsets const & offsets )
{
  return json::to_cbor( json( offsets ) );
}

/**
 * @brief
 * @tparam V Container of std::uint8_t
 * @param data
 * @return
 */
template< class V >
BucketOffsets deserialize( V const & data )
{
  return json::from_cbor( data, false ).template get< BucketOffsets >();
}

/**
 * @brief Custom MPI reduction function. Merges the sizes of the bucket from current rank into numbering offsets.
 * @param[in] in Contains the reduced result from the previous ranks. Needs to be unpacked into a @c BucketOffsets instance.
 * @param[inout] inout As an @e input, contains a pointer to the @c BucketSizes instance from the current rank.
 * As an @e output, will contained the (packed) instance of @c BucketOffsets after the @c reduction.
 * The @e output instance will be used by both current and next ranks.
 * @param[in] len Number of data.
 * @param[in] dataType Type of data. Mean to be @c MPI_BYTE for the current case.
 */
void f( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dataType )
{
  GEOS_ASSERT_EQ( *dataType, MPI_BYTE );

  MpiRank const curRank{ MpiWrapper::commRank() };

  // offsets provided by the previous rank(s)
  BucketOffsets const offsets = deserialize( Span< std::uint8_t >( (std::uint8_t *) in, *len ) );

  // Sizes provided by the current rank, under the form of a pointer to the data.
  // No need to serialize, since we're on the same rank.
  std::uintptr_t addr;
  std::memcpy( &addr, inout, sizeof( std::uintptr_t ) );
  BucketSizes const * sizes = reinterpret_cast<BucketSizes const *>(addr);

  BucketOffsets updatedOffsets;
  updatedOffsets.edges = updateBucketOffsets< EdgeGlbIdx >( sizes->edges, offsets.edges, curRank );
  updatedOffsets.faces = updateBucketOffsets< FaceGlbIdx >( sizes->faces, offsets.faces, curRank );

  // Serialize the updated offsets, so they get sent to the next rank.
  std::vector< std::uint8_t > const serialized = serialize( updatedOffsets );
  std::memcpy( inout, serialized.data(), serialized.size() );
}

struct MaxGlbIdcs
{
  NodeGlbIdx nodes;
  EdgeGlbIdx edges;
  FaceGlbIdx faces;
  CellGlbIdx cells;
};

void to_json( json & j,
              const MaxGlbIdcs & mo )
{
  j = json{ { "nodes", mo.nodes },
            { "edges", mo.edges },
            { "faces", mo.faces },
            { "cells", mo.cells } };
}

void g( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dataType )
{
  GEOS_ASSERT_EQ( *len, 1 );

  MaxGlbIdcs const * i = reinterpret_cast<MaxGlbIdcs const *>(in);
  MaxGlbIdcs * io = reinterpret_cast<MaxGlbIdcs *>(inout);

  io->nodes = std::max( i->nodes, io->nodes );
  io->edges = std::max( i->edges, io->edges );
  io->faces = std::max( i->faces, io->faces );
  io->cells = std::max( i->cells, io->cells );
}

MaxGlbIdcs gatherOffset( vtkSmartPointer< vtkDataSet > mesh,
                         EdgeGlbIdx const & maxEdgeId,
                         FaceGlbIdx const & maxFaceId )
{
  MaxGlbIdcs offsets{ NodeGlbIdx{ 0 }, maxEdgeId, maxFaceId, CellGlbIdx{ 0 } };

  auto const extract = []( vtkDataArray * globalIds ) -> vtkIdType
  {
    vtkIdTypeArray * gids = vtkIdTypeArray::FastDownCast( globalIds );
    Span< vtkIdType > const s( (vtkIdType *) gids->GetPointer( 0 ), gids->GetNumberOfTuples() );
    return *std::max_element( s.begin(), s.end() );
  };

  offsets.nodes = NodeGlbIdx{ extract( mesh->GetPointData()->GetGlobalIds() ) };
  offsets.cells = CellGlbIdx{ extract( mesh->GetCellData()->GetGlobalIds() ) };

  // Otherwise, use `MPI_Type_create_struct`.
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, EdgeGlbIdx::UnderlyingType > );
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, FaceGlbIdx::UnderlyingType > );
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, CellGlbIdx::UnderlyingType > );

//  MPI_Datatype const underlying = internal::getMpiType< NodeGlbIdx::UnderlyingType >();
  MPI_Datatype t;
//  MPI_Type_contiguous( sizeof( MatrixOffsets ) / sizeof( underlying ), underlying, &t );
  MPI_Type_contiguous( 4, MPI_LONG_LONG_INT, &t );
  MPI_Type_commit( &t );

  MPI_Op op;
  MPI_Op_create( g, true, &op );

  MaxGlbIdcs result( offsets );

  MPI_Allreduce( &offsets, &result, 1, t, op, MPI_COMM_WORLD );

  return result;
}

struct MeshGraph  // TODO add the local <-> global mappings here?
{
  std::map< CellGlbIdx, std::set< FaceGlbIdx > > c2f;  // TODO What about the metadata (e.g. flip the face)
  std::map< FaceGlbIdx, std::set< EdgeGlbIdx > > f2e;
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n; // TODO use Edge here?
  std::set< NodeGlbIdx > nodes;
  std::set< FaceGlbIdx > otherFaces;  // Faces that are there but not owned.
  std::set< EdgeGlbIdx > otherEdges;  // Edges that are there but not owned.
  std::set< NodeGlbIdx > otherNodes;  // Nodes that are there but not owned.
  // TODO add the nodes here?
  // TODO add all types of connections here? How?
};

void to_json( json & j,
              const MeshGraph & v )  // For display
{
  j = json{ { "c2f", v.f2e },
            { "f2e", v.f2e },
            { "e2n", v.e2n } };
}


/**
 * @brief Builds the graph information for the owned elements only.
 * @param mesh
 * @param buckets
 * @param offsets
 * @param curRank
 * @return
 */
MeshGraph buildMeshGraph( vtkSmartPointer< vtkDataSet > mesh,  // TODO give a sub-mesh?
                          Buckets const & buckets,
                          BucketOffsets const & offsets,
                          MpiRank curRank )
{
  MeshGraph result;

  auto const isCurrentRankOwning = [&curRank]( std::set< MpiRank > const & ranks ) -> bool
  {
    return curRank == *std::min_element( std::cbegin( ranks ), std::cend( ranks ) );
  };

  for( auto const & [ranks, ns]: buckets.nodes )
  {
    std::set< NodeGlbIdx > & nodes = isCurrentRankOwning( ranks ) ? result.nodes : result.otherNodes;
    nodes.insert( std::cbegin( ns ), std::cend( ns ) );
  }

  // The `e2n` is a mapping for all the geometrical entities, not only the one owned like `result.e2n`.
  // TODO check that it is really useful.
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n;
  for( auto const & [ranks, edges]: buckets.edges )
  {
    bool const isOwning = isCurrentRankOwning( ranks );
    auto & m = isOwning ? result.e2n : e2n;
    EdgeGlbIdx i = offsets.edges.at( ranks );  // TODO hack
    for( Edge const & edge: edges )
    {
      m[i] = edge;
      if( not isOwning )
      {
        result.otherEdges.insert( i );  // TODO use the keys of e2n instead?
      }
      ++i;
    }
  }
  e2n.insert( std::cbegin( result.e2n ), std::cend( result.e2n ) );

  // Simple inversion
  std::map< std::tuple< NodeGlbIdx, NodeGlbIdx >, EdgeGlbIdx > n2e;
  for( auto const & [e, n]: e2n )
  {
    n2e[n] = e;  // TODO what about ownership?
  }

  for( auto const & [ranks, faces]: buckets.faces )
  {
    bool const isOwning = isCurrentRankOwning( ranks );

    if( isOwning )
    {
      FaceGlbIdx i = offsets.faces.at( ranks );
      for( Face face: faces )  // Intentional copy for the future `emplace_back`.
      {
        face.emplace_back( face.front() );  // Trick to build the edges.
        for( std::size_t ii = 0; ii < face.size() - 1; ++ii )
        {
          NodeGlbIdx const & n0 = face[ii], & n1 = face[ii + 1];
          std::pair< NodeGlbIdx, NodeGlbIdx > const p0 = std::make_pair( n0, n1 );
          std::pair< NodeGlbIdx, NodeGlbIdx > const p1 = std::minmax( n0, n1 );
          result.f2e[i].insert( n2e.at( p1 ) );
          bool const flipped = p0 != p1;  // TODO store somewhere.
        }
        ++i;
      }
    }
    else
    {
      FaceGlbIdx const size = FaceGlbIdx{ intConv< FaceGlbIdx::UnderlyingType >( std::size( faces ) ) };
      for( FaceGlbIdx ii = offsets.faces.at( ranks ); ii < size; ++ii )
      {
        result.otherFaces.insert( ii );  // TODO insert iota
      }
    }
  }

  std::map< std::vector< NodeGlbIdx >, FaceGlbIdx > n2f;
  for( auto const & [ranks, faces]: buckets.faces )
  {
    FaceGlbIdx i = offsets.faces.at( ranks );
    for( Face const & face: faces )  // TODO hack
    {
      n2f[face] = i;
      ++i;
    }
  }

  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );
  vtkIdTypeArray const * globalCellIds = vtkIdTypeArray::FastDownCast( mesh->GetCellData()->GetGlobalIds() );  // TODO do the mapping beforehand
  for( vtkIdType c = 0; c < mesh->GetNumberOfCells(); ++c )
  {
    vtkCell * cell = mesh->GetCell( c );
    CellGlbIdx const gci{ globalCellIds->GetValue( c ) };
    // TODO copy paste?
    for( auto f = 0; f < cell->GetNumberOfFaces(); ++f )
    {
      vtkCell * face = cell->GetFace( f );
      vtkIdList * pids = face->GetPointIds();
      std::vector< NodeGlbIdx > faceNodes( pids->GetNumberOfIds() );
      for( std::size_t i = 0; i < faceNodes.size(); ++i )
      {
        vtkIdType const lni = face->GetPointId( i );
        vtkIdType const gni = globalPtIds->GetValue( lni );
        faceNodes[i] = NodeGlbIdx{ gni };
      }
      std::vector< NodeGlbIdx > const reorderedFaceNodes = reorderFaceNodes( faceNodes );
      result.c2f[gci].insert( n2f.at( reorderedFaceNodes ) );
      // TODO... bool const flipped = ... compare nodes and reorderedNodes. Or ask `reorderFaceNodes` to tell
    }
  }

  return result;
}


void assembleAdjacencyMatrix( MeshGraph const & graph,
                              MaxGlbIdcs const & gis,
                              MpiRank curRank )
{
  std::size_t const edgeOffset = gis.nodes.get() + 1;
  std::size_t const faceOffset = edgeOffset + gis.edges.get() + 1;
  std::size_t const cellOffset = faceOffset + gis.faces.get() + 1;

  std::size_t const n = cellOffset + gis.cells.get() + 1;  // Total number of entries in the graph.

  std::size_t const numOwnedNodes = std::size( graph.nodes );
  std::size_t const numOwnedEdges = std::size( graph.e2n );
  std::size_t const numOwnedFaces = std::size( graph.f2e );
  std::size_t const numOwnedCells = std::size( graph.c2f );
  std::size_t const numOwned = numOwnedNodes + numOwnedEdges + numOwnedFaces + numOwnedCells;

  std::size_t const numOtherNodes = std::size( graph.otherNodes );
  std::size_t const numOtherEdges = std::size( graph.otherEdges );
  std::size_t const numOtherFaces = std::size( graph.otherFaces );
  std::size_t const numOther = numOtherNodes + numOtherEdges + numOtherFaces;

  std::vector< int > ownedGlbIdcs, numEntriesPerRow;  // TODO I couldn't use a vector of `std::size_t`
  std::vector< std::vector< int > > indices;
  ownedGlbIdcs.reserve( numOwned );
  numEntriesPerRow.reserve( numOwned );
  indices.reserve( numOwned );

  std::vector< int > otherGlbIdcs;  // TODO I couldn't use a vector of `std::size_t`
  otherGlbIdcs.reserve( numOther );
  for( NodeGlbIdx const & ngi: graph.otherNodes )
  {
    otherGlbIdcs.emplace_back( ngi.get() );
  }
  for( EdgeGlbIdx const & egi: graph.otherEdges )
  {
    otherGlbIdcs.emplace_back( egi.get() + edgeOffset );
  }
  for( FaceGlbIdx const & fgi: graph.otherFaces )
  {
    otherGlbIdcs.emplace_back( fgi.get() + faceOffset );
  }
  GEOS_ASSERT_EQ( numOther, std::size( otherGlbIdcs ) );

  for( NodeGlbIdx const & ngi: graph.nodes )
  {
    auto const i = ngi.get();
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( 1 );
    std::vector< int > const tmp( 1, ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }
  for( auto const & [egi, nodes]: graph.e2n )
  {
    auto const i = egi.get() + edgeOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::tuple_size_v< decltype( nodes ) > + 1 );  // `+1` comes from the identity
    std::vector< int > const tmp{ int( std::get< 0 >( nodes ).get() ), int( std::get< 1 >( nodes ).get() ), ownedGlbIdcs.back() };
    indices.emplace_back( tmp );
  }
  for( auto const & [fgi, edges]: graph.f2e )
  {
    auto const i = fgi.get() + faceOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::size( edges ) + 1 );  // `+1` comes from the identity
    std::vector< int > tmp;
    tmp.reserve( numEntriesPerRow.back() );
    for( EdgeGlbIdx const & egi: edges )
    {
      tmp.emplace_back( egi.get() + edgeOffset );
    }
    tmp.emplace_back( ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }
  for( auto const & [cgi, faces]: graph.c2f )
  {
    auto const i = cgi.get() + cellOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::size( faces ) + 1 );  // `+1` comes from the identity
    std::vector< int > tmp;
    tmp.reserve( numEntriesPerRow.back() );
    for( FaceGlbIdx const & fgi: faces )
    {
      tmp.emplace_back( fgi.get() + faceOffset );
    }
    tmp.emplace_back( ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }

  GEOS_ASSERT_EQ( numOwned, std::size( ownedGlbIdcs ) );
  GEOS_ASSERT_EQ( numOwned, std::size( numEntriesPerRow ) );
  GEOS_ASSERT_EQ( numOwned, std::size( indices ) );
  for( std::size_t i = 0; i < numOwned; ++i )
  {
    GEOS_ASSERT_EQ( indices[i].size(), std::size_t( numEntriesPerRow[i] ) );
  }

//  GEOS_LOG_RANK( "n = " << n );
//  GEOS_LOG_RANK( "ownedGlbIdcs = " << json( ownedGlbIdcs ) );

  Epetra_MpiComm const & comm = Epetra_MpiComm( MPI_COMM_GEOSX );
  Epetra_Map const rowMap( n, numOwned, ownedGlbIdcs.data(), 0, comm );

  Epetra_CrsMatrix adj( Epetra_DataAccess::Copy, rowMap, numEntriesPerRow.data(), true );

  for( std::size_t i = 0; i < numOwned; ++i )
  {
    std::vector< int > const & rowIndices = indices[i];
    std::vector< double > const rowValues( std::size( rowIndices ), 1. );
    GEOS_ASSERT_EQ( std::size( rowIndices ), std::size_t( numEntriesPerRow[i] ) );
    adj.InsertGlobalValues( ownedGlbIdcs[i], std::size( rowIndices ), rowValues.data(), rowIndices.data() );
  }

  adj.FillComplete();
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/adj.mat", adj );

  // Now let's build the domain indicator matrix.
  // It's rectangular, one dimension being the number of MPI ranks, the other the number of nodes in the mesh graph.
  Epetra_Map const mpiMap( MpiWrapper::commSize(), 0, comm );  // Rows
  Epetra_CrsMatrix indicator( Epetra_DataAccess::Copy, mpiMap, numOwned + numOther, true );

  std::vector< double > const ones( n, 1. );
  indicator.InsertGlobalValues( curRank.get(), numOwned, ones.data(), ownedGlbIdcs.data() );
  indicator.InsertGlobalValues( curRank.get(), numOther, ones.data(), otherGlbIdcs.data() );

  Epetra_Map const graphNodeMap( int( n ), 0, comm );  // Columns
  indicator.FillComplete( graphNodeMap, mpiMap );

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "indicator.NumGlobalCols() = " << indicator.NumGlobalCols() );
    GEOS_LOG_RANK( "indicator.NumGlobalRows() = " << indicator.NumGlobalRows() );
  }
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/indicator.mat", indicator );

  int NN( MpiWrapper::commSize() );
  GEOS_LOG_RANK("A");
  auto multiply = [&]()-> Epetra_CrsMatrix
  {
    // Upward (n -> e -> f -> c)

    Epetra_CrsMatrix result0( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, false, indicator, true, result0, false );
    result0.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-0.mat", result0 );

    Epetra_CrsMatrix result1( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, false, result0, false, result1, false );
    result1.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-1.mat", result1 );

    Epetra_CrsMatrix result2( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, false, result1, false, result2, false );
    result2.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-2.mat", result2 );

    // Unneeded step
    Epetra_CrsMatrix result3( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, false, result2, false, result3, false );
    result3.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-3.mat", result3 );

    // Downward (c -> f -> e -> n)

    Epetra_CrsMatrix result4( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, true, result3, false, result4, false );
    result4.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-4.mat", result4 );

    Epetra_CrsMatrix result5( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, true, result4, false, result5, false );
    result5.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-5.mat", result5 );

    Epetra_CrsMatrix result6( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, true, result5, false, result6, false );
    result6.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-6.mat", result6 );

    // Unneeded step
    Epetra_CrsMatrix result7( Epetra_DataAccess::Copy, rowMap, NN, false );
    EpetraExt::MatrixMatrix::Multiply( adj, true, result6, false, result7, false );
    result7.FillComplete( mpiMap, graphNodeMap );
    EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-7.mat", result7 );
    return result7;
  };

  Epetra_CrsMatrix ghosted( multiply() );
  ghosted.FillComplete( mpiMap, graphNodeMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ghosted.mat", ghosted );

  // Transposing in order to get easy access to the row...
  Epetra_RowMatrixTransposer transposer( &ghosted );
  Epetra_CrsMatrix * ptGhosted = nullptr;  // The transposer returns a pointer we must handle ourselves.
  transposer.CreateTranspose( true, ptGhosted );
  GEOS_LOG_RANK( "ptGhosted->NumGlobalCols() = " << ptGhosted->NumGlobalCols() );
  GEOS_LOG_RANK( "ptGhosted->NumGlobalRows() = " << ptGhosted->NumGlobalRows() );

  // My test says that `tGhosted` is filled here.
  int extracted = 0;
  std::vector< double > extractedValues( n );
  std::vector< int > extractedIndices( n );
  ptGhosted->ExtractGlobalRowCopy( curRank.get(), int( n ), extracted, extractedValues.data(), extractedIndices.data() );
  extractedValues.resize( extracted );
  extractedIndices.resize( extracted );
  GEOS_LOG_RANK( "extracted = " << extracted );
  {
    std::vector< int > interest;
    for( int i = 0; i < extracted; ++i )
    {
      int const & index = extractedIndices[i];
      double const & val = extractedValues[i];
      if( val > 0 and index < int( edgeOffset ) )
      {
//        interest.push_back( index - cellOffset );
        interest.push_back( index );
      }
    }
    GEOS_LOG_RANK( "ghost interest = " << json( interest ) );
  }

  GEOS_LOG_RANK("B");

  delete ptGhosted;
}


void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< MpiRank > const & neighbors )
{
  // Now we exchange the data with our neighbors.
  MpiRank const curRank{ MpiWrapper::commRank() };

  std::vector< NeighborCommunicator > ncs;
  ncs.reserve( neighbors.size() );
  for( MpiRank const & rank: neighbors )
  {
    ncs.emplace_back( rank.get() );
  }

  CommID const commId = CommunicationTools::getInstance().getCommID();
  std::map< MpiRank, Exchange > exchanged = exchange( int( commId ), ncs, buildExchangeData( mesh ) );
  exchanged[MpiRank{ MpiWrapper::commRank() }] = buildFullData( mesh );

//  std::cout << "exchanged on rank " << curRank << " -> " << json( exchanged[MpiRank{ MpiWrapper::commRank() }] ) << std::endl;

  Buckets const buckets = buildIntersectionBuckets( exchanged, curRank, neighbors );

  std::size_t const maxBufferSize = 100 * buildMaxBufferSize( buckets );

  std::vector< std::uint8_t > sendBuffer( maxBufferSize );
  std::vector< std::uint8_t > recvBuffer( maxBufferSize );

  BucketSizes const sizes = getBucketSize( buckets );
  BucketOffsets offsets;
  if( curRank == 0_mpi )
  {
    // The `MPI_Scan` process will not call the reduction operator for rank 0.
    // So we need to reduce ourselves for ourselves.
    offsets.edges = updateBucketOffsets< EdgeGlbIdx >( sizes.edges, { { { 0_mpi, }, 0_egi } }, curRank );
    offsets.faces = updateBucketOffsets< FaceGlbIdx >( sizes.faces, { { { 0_mpi, }, 0_fgi } }, curRank );
    // Still we need to send this reduction to the following rank, by copying to it to the send buffer.
    std::vector< std::uint8_t > const bytes = serialize( offsets );
    std::memcpy( sendBuffer.data(), bytes.data(), bytes.size() );
  }
  else
  {
    // For the other ranks, the reduction operator will be called during the `Mpi_Scan` process.
    // So unlike for rank 0, we do not have to do it ourselves.
    // In order to provide the `sizes` to the reduction operator, since `sizes` will only be used on the current ranl,
    // we'll provide the information as a pointer to the instance.
    // The reduction operator will then compute the new offsets and send them to the following rank.
    std::uintptr_t const addr = reinterpret_cast<std::uintptr_t>(&sizes);
    std::memcpy( sendBuffer.data(), &addr, sizeof( std::uintptr_t ) );
  }

  MPI_Op op;
  MPI_Op_create( f, false, &op );

  MPI_Scan( sendBuffer.data(), recvBuffer.data(), maxBufferSize, MPI_BYTE, op, MPI_COMM_WORLD );

  if( curRank != 0_mpi )
  {
    offsets = deserialize( recvBuffer );
  }

  std::cout << "offsets on rank " << curRank << " -> " << json( offsets ) << std::endl;

  MpiRank const nextRank = curRank + 1_mpi;
  MaxGlbIdcs const matrixOffsets = gatherOffset( mesh, offsets.edges.at( { nextRank } ) - 1_egi, offsets.faces.at( { nextRank } ) - 1_fgi );
  std::cout << "matrixOffsets on rank " << curRank << " -> " << json( matrixOffsets ) << std::endl;

  MeshGraph const graph = buildMeshGraph( mesh, buckets, offsets, curRank );  // TODO change into buildOwnedMeshGraph?
//  if( curRank == MpiRank{ 1 } )
//  {
//    std::cout << "My graph is " << json( graph ) << std::endl;
//  }

  assembleAdjacencyMatrix( graph, matrixOffsets, curRank );
}

void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< int > const & neighbors )
{
  std::set< MpiRank > neighbors_;
  for( int const & rank: neighbors )
  {
    neighbors_.insert( MpiRank{ rank } );
  }

  return doTheNewGhosting( mesh, neighbors_ );
}

}  // end of namespace geos::ghosting
