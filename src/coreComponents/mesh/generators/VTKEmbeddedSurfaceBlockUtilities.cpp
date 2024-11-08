/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "VTKEmbeddedSurfaceBlockUtilities.hpp"
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include "CellBlockManager.hpp"

namespace geos::vtk
{

void importEmbeddedFractureNetwork( string const & embeddedSurfaceBlockName,
                                    vtkSmartPointer< vtkDataSet > embeddedSurfaceMesh,
                                    vtkSmartPointer< vtkDataSet > mesh,
                                    CellBlockManager & cellBlockManager )
{

  // Ouassim change
  vtkIdType numEdfmFracs = embeddedSurfaceMesh->GetNumberOfCells();
  //1.Get EDFM vertices
  vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast( embeddedSurfaceMesh );
  vtkPoints * const nodes = grid->GetPoints();
  vtkIdType numNodes = nodes->GetNumberOfPoints();

  ArrayOfArrays< localIndex > elem2dToNodes( numEdfmFracs );
  for( vtkIdType i = 0; i < numEdfmFracs; ++i )
  {
    auto val = grid->GetCell( i )->GetPointIds();
    elem2dToNodes.emplaceBack( i, val->GetId( 0 ));
    elem2dToNodes.emplaceBack( i, val->GetId( 1 ));
    elem2dToNodes.emplaceBack( i, val->GetId( 2 ));
    elem2dToNodes.emplaceBack( i, val->GetId( 3 ));
  }


  ArrayOfArrays< real64 > embeddedSurfaceElemNodes( LvArray::integerConversion< localIndex >( numNodes ) );
  for( vtkIdType i = 0; i < numNodes; ++i )
  {
    double point[3];
    nodes->GetPoint( i, point );
    embeddedSurfaceElemNodes.emplaceBack( i, point[0] );
    embeddedSurfaceElemNodes.emplaceBack( i, point[1] );
    embeddedSurfaceElemNodes.emplaceBack( i, point[2] );
  }

  //2. Get EDFM normal vectors
  ArrayOfArrays< real64 > embeddedSurfaceNormalVectors( numEdfmFracs );
  const char * const normal_str = "normal_vectors";
  vtkDataArray * normals =grid->GetCellData()->GetArray( normal_str );
  for( vtkIdType i = 0; i < numEdfmFracs; ++i )
  {
    auto val = normals->GetTuple3( i );
    embeddedSurfaceNormalVectors.emplaceBack( i, val[0] );
    embeddedSurfaceNormalVectors.emplaceBack( i, val[1] );
    embeddedSurfaceNormalVectors.emplaceBack( i, val[2] );
  }

  //3. Get EDFM  tangential length vectors (horizontal)
  ArrayOfArrays< real64 > embeddedSurfaceTangentialLengthVectors( numEdfmFracs );
  const char * const tl_str = "tangential_length_vectors";
  vtkDataArray * tangential_len =grid->GetCellData()->GetArray( tl_str );
  for( vtkIdType i = 0; i < numEdfmFracs; ++i )
  {
    auto val = tangential_len->GetTuple3( i );
    embeddedSurfaceTangentialLengthVectors.emplaceBack( i, val[0] );
    embeddedSurfaceTangentialLengthVectors.emplaceBack( i, val[1] );
    embeddedSurfaceTangentialLengthVectors.emplaceBack( i, val[2] );
  }

  //4. Get EDFM  tangential width vectors (vertical)
  ArrayOfArrays< real64 > embeddedSurfaceTangentialWidthVectors( numEdfmFracs );
  const char * const tw_str = "tangential_width_vectors";
  vtkDataArray * tangential_width =grid->GetCellData()->GetArray( tw_str );
  for( vtkIdType i = 0; i < numEdfmFracs; ++i )
  {
    auto val = tangential_width->GetTuple3( i );
    embeddedSurfaceTangentialWidthVectors.emplaceBack( i, val[0] );
    embeddedSurfaceTangentialWidthVectors.emplaceBack( i, val[1] );
    embeddedSurfaceTangentialWidthVectors.emplaceBack( i, val[2] );
  }

  //4. Get EDFM aperture
  array1d< real64 > embeddedSurfaceAperture( numEdfmFracs );
  const char * const aperture_str = "aperture";
  vtkDataArray * apertures =grid->GetCellData()->GetArray( aperture_str );
  for( vtkIdType i = 0; i < numEdfmFracs; ++i )
  {
    auto val =apertures->GetTuple1( i );
    embeddedSurfaceAperture[i] = val;
  }

  //5. Get EDFM permeability
  array1d< real64 > embeddedSurfacePermeability( numEdfmFracs );
  const char * const perm_str = "permeability";
  vtkDataArray * perms =grid->GetCellData()->GetArray( perm_str );
  for( vtkIdType i = 0; i < numEdfmFracs; ++i )
  {
    auto val =perms->GetTuple1( i );
    embeddedSurfaceAperture[i] = val;
  }

  //6. Get edfm to matrix cell mapping
  const char * const fid_to_mid_mapping  = "fracture_to_parent_matrix_cell_mapping";
  vtkDataArray * fid_mid =grid->GetCellData()->GetArray( fid_to_mid_mapping );

  ArrayOfArrays< localIndex > toBlockIndex( numEdfmFracs );
  ArrayOfArrays< localIndex > toCellIndex( numEdfmFracs );
  for( vtkIdType i = 0; i < numEdfmFracs; ++i )
  {
    auto fidmid =fid_mid->GetTuple2( i );
    toBlockIndex.emplaceBack( fidmid[0], 0 );// cell block is set to 0 for now
    toCellIndex.emplaceBack( fidmid[0], fidmid[1] );
  }
  ToCellRelation< ArrayOfArrays< localIndex > > elem2dTo3d( std::move( toBlockIndex ), std::move( toCellIndex ) );

  EmbeddedSurfaceBlock & embeddedSurfBlock = cellBlockManager.registerEmbeddedSurfaceBlock( embeddedSurfaceBlockName );

  embeddedSurfBlock.setNumEmbeddedSurfElem( numEdfmFracs );
  embeddedSurfBlock.setEmbeddedSurfElemNodes( std::move( embeddedSurfaceElemNodes ));
  embeddedSurfBlock.setEmbeddedSurfElemNormalVectors( std::move( embeddedSurfaceNormalVectors ));
  embeddedSurfBlock.setEmbeddedSurfElemTangentialLengthVectors( std::move( embeddedSurfaceTangentialLengthVectors ));
  embeddedSurfBlock.setEmbeddedSurfElemTangentialWidthVectors( std::move( embeddedSurfaceTangentialWidthVectors ));
  embeddedSurfBlock.setEmbeddedSurfElemAperture( std::move( embeddedSurfaceAperture ));
  embeddedSurfBlock.setEmbeddedSurfElemPermeability( std::move( embeddedSurfacePermeability ));
  embeddedSurfBlock.setEmbeddedSurfElemTo3dElem( std::move( elem2dTo3d ));
  embeddedSurfBlock.setEmbeddedSurfElemToNodes( std::move( elem2dToNodes ));
}
}
