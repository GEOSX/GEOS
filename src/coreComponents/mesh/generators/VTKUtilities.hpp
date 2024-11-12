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

/**
 * @file VTKUtilities.hpp
 */

#ifndef GEOS_MESH_GENERATORS_VTKUTILITIES_HPP
#define GEOS_MESH_GENERATORS_VTKUTILITIES_HPP

#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"
#include "mesh/generators/CellBlockManager.hpp"

#include <vtkDataSet.h>
#include <vtkMultiProcessController.h>
#include <vtkSmartPointer.h>

#include <numeric>
#include <unordered_set>

namespace geos
{
namespace vtk
{

/**
 * @brief Choice of advanced mesh partitioner
 */
enum class PartitionMethod : integer
{
  parmetis, ///< Use ParMETIS library
  ptscotch, ///< Use PTScotch library
};

/// Strings for VTKMeshGenerator::PartitionMethod enumeration
ENUM_STRINGS( PartitionMethod,
              "parmetis",
              "ptscotch" );

/**
 * @brief Type of map used to store cell lists.
 *
 * This should be an unordered_map, but some outdated standard libraries on some systems
 * do not provide std::hash specialization for enums. This is not performance critical though.
 */
using CellMapType = std::map< ElementType, std::unordered_map< int, std::vector< vtkIdType > > >;

/**
 * @brief Return a VTK controller for multiprocessing.
 * @return Return a VTK controller for multiprocessing.
 */
vtkSmartPointer< vtkMultiProcessController > getController();

/**
 * @brief Gathers all the vtk meshes together.
 */
class AllMeshes
{
public:
  AllMeshes() = default;

  /**
   * @brief Builds the compound from values.
   * @param main The main 3d mesh (the matrix).
   * @param faceBlocks The fractures meshes.
   */
  AllMeshes( vtkSmartPointer< vtkDataSet > const & main,
             std::map< string, vtkSmartPointer< vtkDataSet > > const & faceBlocks )
    : m_main( main ),
    m_faceBlocks( faceBlocks )
  { }

  /**
   * @return the main 3d mesh for the simulation.
   */
  vtkSmartPointer< vtkDataSet > getMainMesh()
  {
    return m_main;
  }

  /**
   * @return a mapping linking the name of each face block to its mesh.
   */
  std::map< string, vtkSmartPointer< vtkDataSet > > & getFaceBlocks()
  {
    return m_faceBlocks;
  }

  /**
   * @brief Defines the main 3d mesh for the simulation.
   * @param main The new 3d mesh.
   */
  void setMainMesh( vtkSmartPointer< vtkDataSet > main )
  {
    m_main = main;
  }

  /**
   * @brief Defines the face blocks/fractures.
   * @param faceBlocks A map which connects each name of the face block to its mesh.
   */
  void setFaceBlocks( std::map< string, vtkSmartPointer< vtkDataSet > > const & faceBlocks )
  {
    m_faceBlocks = faceBlocks;
  }

private:
  /// The main 3d mesh (namely the matrix).
  vtkSmartPointer< vtkDataSet > m_main;

  /// The face meshes (namely the fractures).
  std::map< string, vtkSmartPointer< vtkDataSet > > m_faceBlocks;
};

/**
 * @brief Load the VTK file into the VTK data structure
 * @param[in] filePath the Path of the file to load
 * @param[in] mainBlockName The name of the block to import (will be considered for multi-block files only).
 * @param[in] faceBlockNames The names of the face blocks to import  (will be considered for multi-block files only).
 * @return The compound of the main mesh and the face block meshes.
 */
AllMeshes loadAllMeshes( Path const & filePath,
                         string const & mainBlockName,
                         array1d< string > const & faceBlockNames );

/**
 * @brief Compute the rank neighbor candidate list.
 * @param[in] boundingBoxes the bounding boxes used by the VTK partitioner for all ranks
 * @return the list of neighboring MPI ranks, will be updated
 */
std::vector< int >
findNeighborRanks( std::vector< vtkBoundingBox > boundingBoxes );

/**
 * @brief Generate global point/cell IDs and redistribute the mesh among MPI ranks.
 * @param[in] logLevel the log level
 * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
 * @param[in] namesToFractures the fracture meshes
 * @param[in] comm the MPI communicator
 * @param[in] method the partitionning method
 * @param[in] partitionRefinement number of graph partitioning refinement cycles
 * @param[in] useGlobalIds controls whether global id arrays from the vtk input should be used
 * @return the vtk grid redistributed
 */
AllMeshes
redistributeMeshes( integer const logLevel,
                    vtkSmartPointer< vtkDataSet > loadedMesh,
                    std::map< string, vtkSmartPointer< vtkDataSet > > & namesToFractures,
                    MPI_Comm const comm,
                    PartitionMethod const method,
                    int const partitionRefinement,
                    int const useGlobalIds );

/**
 * @brief Collect lists of VTK cell indices organized by type and attribute value.
 * @param[in] mesh the vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] attributeName name of the VTK data array containing the attribute, if any
 * @return A map from element type to a map of attribute to the associated cell ids for the current rank.
 *         The map contains entries for all types and attribute values across all MPI ranks,
 *         even if there are no cells on current rank (then the list will be empty).
 */
CellMapType buildCellMap( vtkDataSet & mesh,
                          string const & attributeName );

/**
 * @brief Print statistics for a vtk mesh
 *
 * @param[in] mesh an input mesh
 * @param cellMap a map of cell lists grouped by type
 * @param comm the MPI communicator
 */
void printMeshStatistics( vtkDataSet & mesh,
                          CellMapType const & cellMap,
                          MPI_Comm const comm );

/**
 * @brief Collect the data to be imported.
 * @param[in] mesh an input mesh
 * @param[in] sourceName a field name
 * @return A VTK data array pointer.
 */
vtkDataArray * findArrayForImport( vtkDataSet & mesh, string const & sourceName );

/**
 * @brief Check if the vtk mesh as a cell data field of name @p sourceName
 * @param[in] mesh an input mesh
 * @param[in] sourceName a field name
 * @return The boolean result.
 * @note No check is performed to see if the type of the vtk array matches any requirement.
 */
bool hasArray( vtkDataSet & mesh, string const & sourceName );

/**
 * @brief build cell block name from regionId and cellType
 * @param[in] type The type of element in the region
 * @param[in] regionId The region considered
 * @return The cell block name
 */
string buildCellBlockName( ElementType const type, int const regionId );

/**
 * @brief Imports 2d and 3d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importMaterialField( std::vector< vtkIdType > const & cellIds,
                          vtkDataArray * vtkArray,
                          dataRepository::WrapperBase & wrapper );

/**
 * @brief Imports 1d and 2d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importRegularField( std::vector< vtkIdType > const & cellIds,
                         vtkDataArray * vtkArray,
                         dataRepository::WrapperBase & wrapper );

/**
 * @brief Imports 1d and 2d arrays from @p vtkArray to @p wrapper, for all the elements/cells of the provided wrapper.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importRegularField( vtkDataArray * vtkArray,
                         dataRepository::WrapperBase & wrapper );


} // namespace vtk

/**
 * @brief Build all the vertex blocks.
 * @param[in] logLevel the log level
 * @param[in] mesh The vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] nodesetNames
 * @param[in] cellBlockManager The instance that stores the vertex blocks.
 * @param[in] translate translate the dataset
 * @param[in] scale scale the dataset
 * @return size of the dataset on x-axis
 */
real64 writeNodes( integer const logLevel,
                   vtkDataSet & mesh,
                   string_array & nodesetNames,
                   CellBlockManager & cellBlockManager,
                   const geos::R1Tensor & translate,
                   const geos::R1Tensor & scale );

/**
 * @brief Build all the cell blocks.
 * @param[in] logLevel the log level
 * @param[in] mesh The vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] cellMap Map from the surfaces index to the list of cells in this surface in this rank.
 * @param[in] cellBlockManager The instance that stores the cell blocks.
 */
void writeCells( integer const logLevel,
                 vtkDataSet & mesh,
                 const geos::vtk::CellMapType & cellMap,
                 CellBlockManager & cellBlockManager );

/**
 * @brief Build the "surface" node sets from the surface information.
 * @param[in] logLevel the log level
 * @param[in] mesh The vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in] cellMap Map from the surfaces index to the list of cells in this surface in this rank.
 * @param[out] cellBlockManager The instance that stores the node sets.
 * @note @p surfacesIdsToCellsIds will contain all the surface ids across all the MPI ranks, but only its cell ids.
 * If the current MPI rank has no cell id for a given surface, then an empty set will be created.
 */
void writeSurfaces( integer const logLevel,
                    vtkDataSet & mesh,
                    const geos::vtk::CellMapType & cellMap,
                    CellBlockManager & cellBlockManager );

} // namespace geos

#endif /* GEOS_MESH_GENERATORS_VTKUTILITIES_HPP */
