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

/*
 * @file VTKWellGenerator.cpp
 *
 */

#include "VTKWellGenerator.hpp"

#include "mesh/generators/VTKUtilities.hpp"

#include "LvArray/src/tensorOps.hpp"

#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPolyLine.h>
#include <vtkCellArrayIterator.h>

namespace geos
{
using namespace dataRepository;

VTKWellGenerator::VTKWellGenerator( string const & name, Group * const parent ):
  WellGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Path to the well file" );

}

void VTKWellGenerator::fillPolylineDataStructure( )
{
  GEOS_MARK_FUNCTION;

  vtkSmartPointer< vtkMultiProcessController > controller = vtk::getController();
  vtkMultiProcessController::SetGlobalController( controller );

  GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': reading well from {}", catalogName(), getName(), m_filePath ) );
  {
    GEOS_LOG_LEVEL_RANK_0( 2, "  reading the dataset..." );
    vtk::AllMeshes allMeshes = vtk::loadAllMeshes( m_filePath, "main", array1d< string >());
    vtkSmartPointer< vtkDataSet > loadedMesh = allMeshes.getMainMesh();
    controller->Broadcast( loadedMesh, 0 );

    vtkSmartPointer< vtkPolyData > polyData = vtkPolyData::SafeDownCast( loadedMesh );

    // load points
    vtkPoints * points = polyData->GetPoints();
    m_polyNodeCoords.resize( points->GetNumberOfPoints(), 3 );
    globalIndex ipoint = 0;
    for( vtkIdType c = 0; c < points->GetNumberOfPoints(); ++c, ++ipoint )
    {
      real64 point[3];
      points->GetPoint( c, point );
      LvArray::tensorOps::copy< 3 >( m_polyNodeCoords[ipoint], point );
    }

    GEOS_ERROR_IF( polyData->GetLines()->GetNumberOfCells() == 0, GEOS_FMT( "{}: Error! Your VTK file {} doesn't contain any well",
                                                                            this->getName(), m_filePath ));

    GEOS_LOG_RANK_0_IF( polyData->GetLines()->GetNumberOfCells() > 1, GEOS_FMT( "{}: Warning! Your VTK file {} contains multiple wells. Only the first one will be read",
                                                                                this->getName(), m_filePath ));

    // load edges
    // polyData->GetLines()->InitTraversal();
    // vtkNew< vtkIdList > idList;
    // polyData->GetLines()->GetNextCell( idList );

    // vtkNew< vtkIdList > idList2;
    // polyData->GetLines()->GetCell( 0, idList2 );

    //newer version of vtk prefer local thread-safe iterator
    //TODO deal with multiple lines and add a for loop to handle them instead of accessing line 0
    vtkNew< vtkIdList > cellPts;
    auto iter = ::vtk::TakeSmartPointer( polyData->GetLines()->NewIterator());
    iter->GetCellAtId( 0, cellPts );

    const globalIndex nbSegments = cellPts->GetNumberOfIds() - 1;
    m_segmentToPolyNodeMap.resizeDimension< 0 >( nbSegments );
    m_segmentToPolyNodeMap.resizeDimension< 1 >( m_numNodesPerElem );

    globalIndex iseg = 0;
    for( vtkIdType pointId = 0; pointId < nbSegments; ++pointId )
    {
      m_segmentToPolyNodeMap[iseg][0] = cellPts->GetId( pointId );
      m_segmentToPolyNodeMap[iseg][1] = cellPts->GetId( pointId + 1 );
      ++iseg;
    }
  }
}

REGISTER_CATALOG_ENTRY( WellGeneratorBase, VTKWellGenerator, string const &, Group * const )
}
