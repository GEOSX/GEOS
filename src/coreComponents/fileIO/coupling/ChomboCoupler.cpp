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

#include "ChomboCoupler.hpp"
#include "hdf5_interface/coupler.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGeneratorFields.hpp"

#include <cstdint>
#include <tuple>
#include <cstdio>

namespace geos
{

ChomboCoupler::ChomboCoupler( MPI_Comm const comm, const string & outputPath, const string & inputPath, MeshLevel & mesh ):
  m_comm( comm ),
  m_outputPath( outputPath ),
  m_inputPath( inputPath ),
  m_face_offset( -1 ),
  m_n_faces_written( -1 ),
  m_node_offset( -1 ),
  m_n_nodes_written( -1 ),
  m_mesh( mesh ),
  m_counter( 0 )
{
  m_mesh.getFaceManager().registerWrapper< array1d< real64 > >( "ChomboPressure" );
}

void ChomboCoupler::write( double dt )
{
  ++m_counter;
  FaceManager const & faces = m_mesh.getFaceManager();
  ElementRegionManager const & elemRegionManager = m_mesh.getElemManager();


  ArrayOfArraysView< localIndex const > const & face_connectivity = faces.nodeList().toViewConst();
  FaceManager::ElemMapType const & toElementRelation = faces.toElementRelation();
  arrayView2d< localIndex const > const & faceToElementRegionIndex = toElementRelation.m_toElementRegion.toViewConst();

  localIndex const n_faces = face_connectivity.size();

  /* Copy the face connectivity into a contiguous array. */
  std::int64_t * connectivity_array = new std::int64_t[4 * n_faces];
  for( localIndex i = 0; i < n_faces; ++i )
  {
    for( localIndex j = 0; j < 4; ++j )
    {
      connectivity_array[4 * i + j] = face_connectivity( i, j );
    }
  }

  arrayView1d< integer const > const & ruptureState =
    faces.getField< fields::surfaceGeneration::ruptureState >();
  arrayView1d< integer const > const & ghostRank = faces.ghostRank();

  localIndex voidRegionIndex = -1;
  elemRegionManager.forElementRegionsComplete( [&]( localIndex const elemRegionIndex,
                                                    ElementRegionBase const & elemRegion )
  {
    if( elemRegion.getName() == "void" )
    {
      voidRegionIndex = elemRegionIndex;
    }
  } );

  bool * faceMask = new bool[n_faces];
  for( localIndex i = 0; i < n_faces; ++i )
  {
    bool isVoid = (faceToElementRegionIndex[i][0] == voidRegionIndex) ||
                  (faceToElementRegionIndex[i][1] == voidRegionIndex);
    faceMask[i] = (ruptureState[i] > 1) && (ghostRank[i] < 0) && (!isVoid);
  }



  /* Build the face FieldMap. */
  FieldMap_in face_fields;
  real64 const * pressure_ptr = faces.getReference< real64_array >( "ChomboPressure" ).data();
  face_fields["Pressure"] = std::make_tuple( H5T_NATIVE_DOUBLE, 1, pressure_ptr );

  /* Build the node FieldMap. */
  copyNodalData();

  FieldMap_in node_fields;
  node_fields["position"] = std::make_tuple( H5T_NATIVE_DOUBLE, 3, m_referencePositionCopy.data() );
  node_fields["displacement"] = std::make_tuple( H5T_NATIVE_DOUBLE, 3, m_displacementCopy.data() );
  node_fields["velocity"] = std::make_tuple( H5T_NATIVE_DOUBLE, 3, m_velocityCopy.data() );

  writeBoundaryFile( m_comm, m_outputPath.data(), dt, faceMask,
                     m_face_offset, m_n_faces_written, n_faces, connectivity_array, face_fields,
                     m_node_offset, m_n_nodes_written, m_referencePositionCopy.size( 0 ), node_fields );
  GEOS_LOG_RANK_0( "Wrote file: " << m_outputPath );
  delete[] connectivity_array;
  delete[] faceMask;
}

void ChomboCoupler::read( bool usePressures )
{
  GEOS_LOG_RANK_0( "Waiting for file existence: " << m_inputPath );
  waitForFileExistence( m_comm, m_inputPath.data() );

  if( usePressures )
  {
    FaceManager & faces = m_mesh.getFaceManager();
    NodeManager & nodes = m_mesh.getNodeManager();

    const localIndex n_faces = faces.size();
    const localIndex n_nodes = nodes.size();

    /* Build the face FieldMap. */
    FieldMap_out face_fields;
    real64 * pressure_ptr = faces.getReference< real64_array >( "ChomboPressure" ).data();
    face_fields["Pressure"] = std::make_tuple( H5T_NATIVE_DOUBLE, 1, pressure_ptr );

    FieldMap_out node_fields;
    node_fields["position"] = std::make_tuple( H5T_NATIVE_DOUBLE, 3, m_referencePositionCopy.data() );

    readBoundaryFile( m_comm, m_inputPath.data(),
                      m_face_offset, m_n_faces_written, n_faces, face_fields,
                      m_node_offset, m_n_nodes_written, n_nodes, node_fields );

    arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & reference_pos = nodes.referencePosition();
    for( localIndex i = 0; i < n_nodes; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        reference_pos( i, j ) = m_referencePositionCopy( i, j );
      }
    }
  }

  int rank;
  MPI_Comm_rank( m_comm, &rank );
  if( rank == 0 )
  {
    std::remove( m_inputPath.data() );
  }
}

void ChomboCoupler::copyNodalData()
{
  NodeManager const & nodes = m_mesh.getNodeManager();
  localIndex const numNodes = nodes.size();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePos = nodes.referencePosition();
  fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const & displacement =
    nodes.getField< fields::solidMechanics::totalDisplacement >();
  fields::solidMechanics::arrayViewConst2dLayoutVelocity const & velocity =
    nodes.getField< fields::solidMechanics::velocity >();

  GEOS_ERROR_IF_NE( referencePos.size( 0 ), numNodes );
  GEOS_ERROR_IF_NE( referencePos.size( 1 ), 3 );
  GEOS_ERROR_IF_NE( displacement.size( 0 ), numNodes );
  GEOS_ERROR_IF_NE( displacement.size( 1 ), 3 );
  GEOS_ERROR_IF_NE( velocity.size( 0 ), numNodes );
  GEOS_ERROR_IF_NE( velocity.size( 1 ), 3 );

  m_referencePositionCopy.resizeWithoutInitializationOrDestruction( numNodes, 3 );
  m_displacementCopy.resizeWithoutInitializationOrDestruction( numNodes, 3 );
  m_velocityCopy.resizeWithoutInitializationOrDestruction( numNodes, 3 );

  for( localIndex i = 0; i < numNodes; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      m_referencePositionCopy( i, j ) = referencePos( i, j );
      m_displacementCopy( i, j ) = displacement( i, j );
      m_velocityCopy( i, j ) = velocity( i, j );
    }
  }
}


} // namespace geos
