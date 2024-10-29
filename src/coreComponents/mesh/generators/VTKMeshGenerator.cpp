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
 * @file VTKMeshGenerator.cpp
 */

#include "VTKMeshGenerator.hpp"

#include "VTKFaceBlockUtilities.hpp"
#include "VTKMeshGeneratorTools.hpp"
#include "CellBlockManager.hpp"

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"

#include <vtkXMLUnstructuredGridWriter.h>

namespace geos
{
using namespace dataRepository;


VTKMeshGenerator::VTKMeshGenerator( string const & name,
                                    Group * const parent )
  : ExternalMeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::regionAttributeString(), &m_attributeName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "attribute" ).
    setDescription( "Name of the VTK cell attribute to use as region marker" );

  registerWrapper( viewKeyStruct::nodesetNamesString(), &m_nodesetNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the VTK nodesets to import" );

  registerWrapper( viewKeyStruct::mainBlockNameString(), &m_mainBlockName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( "main" ).
    setDescription( "For multi-block files, name of the 3d mesh block." );

  registerWrapper( viewKeyStruct::faceBlockNamesString(), &m_faceBlockNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "For multi-block files, names of the face mesh block." );

  registerWrapper( viewKeyStruct::partitionRefinementString(), &m_partitionRefinement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Number of partitioning refinement iterations (defaults to 1, recommended value)."
                    "A value of 0 disables graph partitioning and keeps simple kd-tree partitions (not recommended). "
                    "Values higher than 1 may lead to slightly improved partitioning, but yield diminishing returns." );

  registerWrapper( viewKeyStruct::partitionMethodString(), &m_partitionMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Method (library) used to partition the mesh" );

  registerWrapper( viewKeyStruct::useGlobalIdsString(), &m_useGlobalIds ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Controls the use of global IDs in the input file for cells and points."
                    " If set to 0 (default value), the GlobalId arrays in the input mesh are used if available, and generated otherwise."
                    " If set to a negative value, the GlobalId arrays in the input mesh are not used, and generated global Ids are automatically generated."
                    " If set to a positive value, the GlobalId arrays in the input mesh are used and required, and the simulation aborts if they are not available" );
}

void VTKMeshGenerator::fillCellBlockManager( CellBlockManager & cellBlockManager, array1d< int > const & GEOS_UNUSED_PARAM( partition ) )
{
  // TODO refactor void MeshGeneratorBase::generateMesh( DomainPartition & domain )
  GEOS_MARK_FUNCTION;

  MPI_Comm const comm = MPI_COMM_GEOSX;
  vtkSmartPointer< vtkMultiProcessController > controller = vtk::getController();
  vtkMultiProcessController::SetGlobalController( controller );

  GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': reading mesh from {}", catalogName(), getName(), m_filePath ) );
  {
    GEOS_LOG_LEVEL_RANK_0( 2, "  reading the dataset..." );
    vtk::AllMeshes allMeshes = vtk::loadAllMeshes( m_filePath, m_mainBlockName, m_faceBlockNames );
    GEOS_LOG_LEVEL_RANK_0( 2, "  redistributing mesh..." );
    vtk::AllMeshes redistributedMeshes =
      vtk::redistributeMeshes( getLogLevel(), allMeshes.getMainMesh(), allMeshes.getFaceBlocks(), comm, m_partitionMethod, m_partitionRefinement, m_useGlobalIds );
    m_vtkMesh = redistributedMeshes.getMainMesh();
    m_faceBlockMeshes = redistributedMeshes.getFaceBlocks();
    GEOS_LOG_LEVEL_RANK_0( 2, "  finding neighbor ranks..." );
    std::vector< vtkBoundingBox > boxes = vtk::exchangeBoundingBoxes( *m_vtkMesh, comm );
    std::vector< int > const neighbors = vtk::findNeighborRanks( std::move( boxes ) );
    m_partition.setMetisNeighborList( std::move( neighbors ) );
    GEOS_LOG_LEVEL_RANK_0( 2, "  done!" );
  }
  GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': generating GEOSX mesh data structure", catalogName(), getName() ) );


  GEOS_LOG_LEVEL_RANK_0( 2, "  preprocessing..." );
  m_cellMap = vtk::buildCellMap( *m_vtkMesh, m_attributeName );

  GEOS_LOG_LEVEL_RANK_0( 2, "  writing nodes..." );
  cellBlockManager.setGlobalLength( writeNodes( getLogLevel(), *m_vtkMesh, m_nodesetNames, cellBlockManager, this->m_translate, this->m_scale ) );

  GEOS_LOG_LEVEL_RANK_0( 2, "  writing cells..." );
  writeCells( getLogLevel(), *m_vtkMesh, m_cellMap, cellBlockManager );

  GEOS_LOG_LEVEL_RANK_0( 2, "  writing surfaces..." );
  writeSurfaces( getLogLevel(), *m_vtkMesh, m_cellMap, cellBlockManager );

  GEOS_LOG_LEVEL_RANK_0( 2, "  building connectivity maps..." );
  cellBlockManager.buildMaps();

  for( auto const & [name, mesh]: m_faceBlockMeshes )
  {
    vtk::importFractureNetwork( name, mesh, m_vtkMesh, cellBlockManager );
  }

  GEOS_LOG_LEVEL_RANK_0( 2, "  done!" );
  vtk::printMeshStatistics( *m_vtkMesh, m_cellMap, comm );
}

void VTKMeshGenerator::importVolumicFieldOnArray( string const & cellBlockName,
                                                  string const & meshFieldName,
                                                  bool isMaterialField,
                                                  dataRepository::WrapperBase & wrapper ) const
{
  for( auto const & typeRegions: m_cellMap )
  {
    // Restrict data import to 3D cells
    if( getElementDim( typeRegions.first ) == 3 )
    {
      for( auto const & regionCells: typeRegions.second )
      {
        string const currentCellBlockName = vtk::buildCellBlockName( typeRegions.first, regionCells.first );
        // We don't know how the user mapped cell blocks to regions, so we must check all of them
        if( cellBlockName != currentCellBlockName )
        {
          continue;
        }

        vtkDataArray * vtkArray = vtk::findArrayForImport( *m_vtkMesh, meshFieldName );
        if( isMaterialField )
        {
          return vtk::importMaterialField( regionCells.second, vtkArray, wrapper );
        }
        else
        {
          return vtk::importRegularField( regionCells.second, vtkArray, wrapper );
        }
      }
    }
  }

  GEOS_ERROR( "Could not import field \"" << meshFieldName << "\" from cell block \"" << cellBlockName << "\"." );
}


void VTKMeshGenerator::importSurfacicFieldOnArray( string const & faceBlockName,
                                                   string const & meshFieldName,
                                                   dataRepository::WrapperBase & wrapper ) const
{
  // Note that there is no additional work w.r.t. the cells on which we want to import the fields,
  // because the face blocks are heterogeneous.
  // We always take the whole data, we do not select cell type by cell type.
  vtkSmartPointer< vtkDataSet > faceMesh = m_faceBlockMeshes.at( faceBlockName );

  // I've noticed that there may be some issues when reading empty arrays (empty, not nulls).
  // It looks like we may be reading above the limits of the array; ghosting is surely at stake here.
  if( faceMesh->GetNumberOfCells() == 0 )
  {
    return;
  }

  if( vtk::hasArray( *faceMesh, meshFieldName ) )
  {
    vtkDataArray * vtkArray = vtk::findArrayForImport( *faceMesh, meshFieldName );
    return vtk::importRegularField( vtkArray, wrapper );
  }

  GEOS_ERROR( "Could not import field \"" << meshFieldName << "\" from face block \"" << faceBlockName << "\"." );
}


void VTKMeshGenerator::importFieldOnArray( Block block,
                                           string const & blockName,
                                           string const & meshFieldName,
                                           bool isMaterialField,
                                           dataRepository::WrapperBase & wrapper ) const
{
  GEOS_ASSERT_MSG( m_vtkMesh, "Must call generateMesh() before importFields()" );

  switch( block )
  {
    case MeshGeneratorBase::Block::VOLUMIC:
      return importVolumicFieldOnArray( blockName, meshFieldName, isMaterialField, wrapper );
    case MeshGeneratorBase::Block::SURFACIC:
    case MeshGeneratorBase::Block::LINEIC:
      return importSurfacicFieldOnArray( blockName, meshFieldName, wrapper );
  }
}

void VTKMeshGenerator::freeResources()
{
  m_vtkMesh = nullptr;
  m_cellMap.clear();
  m_faceBlockMeshes.clear();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )

} // namespace geos
