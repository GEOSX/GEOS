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
 * @file VTKMeshGeneratorTools.hpp
 */

#include "VTKMeshGeneratorTools.hpp"

#include <vtkAppendFilter.h>
#include <vtkDIYGhostUtilities.h>
#include <vtkDIYUtilities.h>

// NOTE: do NOT include anything from GEOS here.
// See full explanation in VTKMeshGeneratorTools.hpp.

namespace geos::vtk
{

vtkSmartPointer< vtkUnstructuredGrid >
redistribute( vtkPartitionedDataSet & localParts,
              MPI_Comm mpiComm )
{
  // The code below is modified from vtkDIYKdTreeUtilities::Exchange():
  // https://gitlab.kitware.com/vtk/vtk/-/blob/7037a148605bf9628710d8b729c22f27dd0ede93/Filters/ParallelDIY2/vtkDIYKdTreeUtilities.cxx#L289
  // We cannot call that function directly because vtkDIYKdTreeUtilities.hpp is a private header in VTK.
  // We also make simplifying assumptions about the nature of input (e.g. one partition per target rank).

  diy::mpi::communicator comm( mpiComm );
  assert( static_cast< int >( localParts.GetNumberOfPartitions() ) == comm.size() );

  using BlockType = std::vector< vtkSmartPointer< vtkUnstructuredGrid > >;

  diy::Master master( comm, 1, -1,
                      [] { return static_cast< void * >( new BlockType() ); },
                      []( void * b ) { delete static_cast< BlockType * >( b ); } );

  diy::ContiguousAssigner const assigner( comm.size(), comm.size() );
  diy::RegularDecomposer< diy::DiscreteBounds > decomposer( 1, diy::interval( 0, comm.size() - 1 ), comm.size() );
  decomposer.decompose( comm.rank(), assigner, master );
  assert( master.size() == 1 );

  int const myRank = comm.rank();
  diy::all_to_all( master, assigner, [myRank, &localParts]( BlockType * block, diy::ReduceProxy const & reduceProxy )
  {
    if( reduceProxy.in_link().size() == 0 )
    {
      // enqueue blocks to send.
      block->reserve( localParts.GetNumberOfPartitions() );
      for( unsigned int partId = 0; partId < localParts.GetNumberOfPartitions(); ++partId )
      {
        if( auto part = vtkUnstructuredGrid::SafeDownCast( localParts.GetPartition( partId ) ) )
        {
          int const targetRank = static_cast< int >( partId );
          if( targetRank == myRank )
          {
            // short-circuit messages to self.
            block->push_back( part );
          }
          else
          {
            reduceProxy.enqueue< vtkDataSet * >( reduceProxy.out_link().target( targetRank ), part );
          }
        }
      }
    }
    else
    {
      for( int i = 0; i < reduceProxy.in_link().size(); ++i )
      {
        int const gid = reduceProxy.in_link().target( i ).gid;
        while( reduceProxy.incoming( gid ) )
        {
          vtkDataSet * ptr = nullptr;
          reduceProxy.dequeue< vtkDataSet * >( reduceProxy.in_link().target( i ), ptr );

          vtkSmartPointer< vtkUnstructuredGrid > sptr;
          sptr.TakeReference( vtkUnstructuredGrid::SafeDownCast( ptr ) );
          block->push_back( sptr );
        }
      }
    }
  } );

  // At this point of the process, it is legitimate to have ranks with no cells for the cases with fractures.
  // But this leaves us with a technical problem since `vtkAppendFilter`
  // (that will be used to merge the different pieces of the meshes)
  // discards the empty the data sets it merges.
  // The definition of "empty" in its context is having no points nor cells...
  //
  // However, some other information which was defined in the discarded data sets gets lost too!
  // In particular, the cell, points and field data were _defined_, but _legitimately_ _empty_.
  // After the `vtkAppendFilter` processing, the definition of those fields is no more available.
  //
  // This leaves us with a specific case when importing fields from vtk,
  // since we need to take extra care of the empty data sets, while we should not have to do this.
  // To circumvent this issue, we gather the point, cell and field data by hand,
  // before registering them by hand again into the final `vtkUnstructuredGrid`.

  // This little structure stores the information we'll need to register back into
  // the cell, points and field data into the final vtkUnstructuredGrid.
  // We should not need it outside of this function, but if we needed, make it a little more solid.
  struct FieldMetaInfo
  {
    enum Location
    {
      CELL,
      POINT,
      FIELD
    };
    std::string name;
    int numComponents;
    int dataType;
    Location location;

    bool operator<( FieldMetaInfo const & other ) const
    {
      return std::tie( name, numComponents, dataType, location ) < std::tie( other.name, other.numComponents, other.dataType, other.location );
    }
  };
  // First step is to gather all the field information
  std::set< FieldMetaInfo > fieldMetaInfo;
  for( unsigned int i = 0; i < master.size(); ++i )
  {
    for( vtkUnstructuredGrid * ug: *master.block< BlockType >( i ) )
    {
      if( !ug )
      {
        break;
      }
      for( int c = 0; c < ug->GetCellData()->GetNumberOfArrays(); ++c )
      {
        auto array = ug->GetCellData()->GetArray( c );
        fieldMetaInfo.insert( { array->GetName(), array->GetNumberOfComponents(), array->GetDataType(), FieldMetaInfo::Location::CELL } );
      }
      for( int c = 0; c < ug->GetPointData()->GetNumberOfArrays(); ++c )
      {
        auto array = ug->GetPointData()->GetArray( c );
        fieldMetaInfo.insert( { array->GetName(), array->GetNumberOfComponents(), array->GetDataType(), FieldMetaInfo::Location::POINT } );
      }
      for( int c = 0; c < ug->GetFieldData()->GetNumberOfArrays(); ++c )
      {
        auto array = ug->GetFieldData()->GetArray( c );
        fieldMetaInfo.insert( { array->GetName(), array->GetNumberOfComponents(), array->GetDataType(), FieldMetaInfo::Location::FIELD } );
      }
    }
  }

  vtkNew< vtkAppendFilter > appender;
  appender->MergePointsOn();
  for( unsigned int i = 0; i < master.size(); ++i )
  {
    for( vtkUnstructuredGrid * ug: *master.block< BlockType >( i ) )
    {
      appender->AddInputDataObject( ug );
    }
  }
  appender->Update();

  vtkUnstructuredGrid * result = vtkUnstructuredGrid::SafeDownCast( appender->GetOutputDataObject( 0 ) );
  // Now we register back the field info.
  if( result->GetNumberOfCells() == 0 )
  {
    for( FieldMetaInfo const & info: fieldMetaInfo )
    {
      vtkAbstractArray * array = vtkAbstractArray::CreateArray( info.dataType );
      array->SetNumberOfComponents( info.numComponents );
      array->SetNumberOfTuples( 0 );
      array->SetName( info.name.c_str() );
      if( info.location == 0 )
      {
        result->GetCellData()->AddArray( array );
      }
      if( info.location == 1 )
      {
        result->GetPointData()->AddArray( array );
      }
      if( info.location == 2 )
      {
        result->GetFieldData()->AddArray( array );
      }
    }
  }

  return result;
}

std::vector< vtkBoundingBox >
exchangeBoundingBoxes( vtkDataSet & dataSet, MPI_Comm mpiComm )
{
  // The code below is modified from vtkDIYGhostUtilities::ExchangeBoundingBoxes():
  // https://gitlab.kitware.com/vtk/vtk/-/blob/1f0e4b2d0be7cd328795131642b5bf7984f681c1/Parallel/DIY/vtkDIYGhostUtilities.txx#L300
  // It makes some simplifications (e.g. just one input dataset per rank).

  using BlockType = std::map< int, vtkBoundingBox >;

  diy::mpi::communicator comm( mpiComm );
  diy::Master master( comm, 1, -1,
                      [] { return static_cast< void * >( new BlockType() ); },
                      []( void * b ) { delete static_cast< BlockType * >( b ); } );

  diy::ContiguousAssigner const assigner( comm.size(), comm.size() );
  diy::RegularDecomposer< diy::DiscreteBounds > decomposer( 1, diy::interval( 0, comm.size() - 1 ), comm.size() );
  decomposer.decompose( comm.rank(), assigner, master );
  assert( master.size() == 1 );

  diy::all_to_all( master, assigner, [&dataSet]( BlockType * block, diy::ReduceProxy const & rp )
  {
    int myBlockId = rp.gid();
    if( rp.round() == 0 )
    {
      vtkBoundingBox bb( dataSet.GetBounds() );
      for( int i = 0; i < rp.out_link().size(); ++i )
      {
        diy::BlockID const blockId = rp.out_link().target( i );
        if( blockId.gid != myBlockId )
        {
          rp.enqueue( blockId, bb.GetMinPoint(), 3 );
          rp.enqueue( blockId, bb.GetMaxPoint(), 3 );
        }
      }
    }
    else
    {
      double minPoint[3], maxPoint[3];
      for( int i = 0; i < static_cast< int >( rp.in_link().size() ); ++i )
      {
        diy::BlockID const blockId = rp.in_link().target( i );
        if( blockId.gid != myBlockId )
        {
          rp.dequeue( blockId, minPoint, 3 );
          rp.dequeue( blockId, maxPoint, 3 );

          block->emplace( blockId.gid, vtkBoundingBox( minPoint[0], maxPoint[0], minPoint[1],
                                                       maxPoint[1], minPoint[2], maxPoint[2] ) );
        }
      }
    }
  } );

  BlockType & boxMap = *master.block< BlockType >( 0 );
  boxMap.emplace( comm.rank(), vtkBoundingBox( dataSet.GetBounds() ) );
  assert( static_cast< int >( boxMap.size() ) == comm.size() );

  std::vector< vtkBoundingBox > boxes;
  boxes.reserve( boxMap.size() );
  for( auto const & rankBox : boxMap )
  {
    boxes.push_back( rankBox.second );
  }
  return boxes;
}

} // namespace geos::vtk
