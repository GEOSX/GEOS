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
 * @file VTKMeshGenerator.hpp
 */

#ifndef GEOS_MESH_GENERATORS_VTKMESHGENERATOR_HPP
#define GEOS_MESH_GENERATORS_VTKMESHGENERATOR_HPP

#include "mesh/generators/ExternalMeshGeneratorBase.hpp"
#include "mesh/generators/VTKUtilities.hpp"

#include <vtkDataSet.h>

namespace geos
{

/**
 *  @class VTKMeshGenerator
 *  @brief The VTKMeshGenerator class provides a class implementation of VTK generated meshes.
 */
class VTKMeshGenerator : public ExternalMeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the VTKMeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  VTKMeshGenerator( const string & name,
                    Group * const parent );

/**
 * @brief Return the name of the VTKMeshGenerator in object Catalog.
 * @return string that contains the key name to VTKMeshGenerator in the Catalog
 */
  static string catalogName() { return "VTKMesh"; }

  /**
   * @brief Generate the mesh using the VTK library.
   * @param[inout] cellBlockManager the CellBlockManager that will receive the meshing information
   * @param[in] partition the number of domain in each direction (x,y,z) for InternalMesh only, not used here
   * @details This method leverages the VTK library to load the meshes.
   * The supported formats are the official VTK ones dedicated to
   * unstructured grids (.vtu, .pvtu and .vtk) and structured grids (.vts, .vti and .pvts).
   *
   * Please note that this mesh generator works only with a number of MPI processes than
   * can be decomposed into a power of 2.
   *
   * - If a .vtu, .vts, .vti or .vtk file is used, the root MPI process will load it.
   *   The mesh will be then redistribute among all the available MPI processes
   * - If a .pvtu or .pvts file is used, it means that the mesh is pre-partionned in the file system.
   *   The available MPI processes will load the pre-partionned mesh. The mesh will be then
   *   redistributed among ALL the available MPI processes.
   *
   * The properties on the mesh will be also and redistributed. The only compatible types are double and float.
   * The properties can be multi-dimensional.
   * The name of the properties has to have the right name in order to be used by GEOSX. For instance,
   * the property that stored the input porosity in GEOSX is named "referencePorosity", so the mesh has to have
   * a property names "referencePorosity".
   *
   * The regions are defined using a property called "attribute" that can be defined in the input mesh. This property
   * will be held by each volume elements. This method will created several CellBlocks, named using the combination
   * of the attribute index and the type of the element.
   * For instance, the cells of a mesh with two regions will hold the attribute "1", or "2". The CellBlocks will
   * be instantiated according to the attribute and the type of the cells. If the region "1" has wedges, tetrahedron
   * and hexahedron, three CellBlocks will be created names 1_tetrahedron, 1_wedges and 1_hexahedron.
   * The ElementRegions have to be be defined in the XML file.
   *
   * The node sets of surface are defined in the same way, using the same
   * property names "attribute" defined in the input mesh. The node sets will
   * hold a name that is just the attribute index. For instance, if a mesh has three
   * surfaces of interest, with triangles and/or quads holding an attribute value
   * of 1, 2 or 3, three node sets named "1", "2" and "3" will be instantiated by this method
   */
  virtual void fillCellBlockManager( CellBlockManager & cellBlockManager, array1d< int > const & partition ) override;

  void importFieldOnArray( Block block,
                           string const & blockName,
                           string const & meshFieldName,
                           bool isMaterialField,
                           dataRepository::WrapperBase & wrapper ) const override;

  virtual void freeResources() override;

private:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * regionAttributeString() { return "regionAttribute"; }
    constexpr static char const * mainBlockNameString() { return "mainBlockName"; }
    constexpr static char const * faceBlockNamesString() { return "faceBlocks"; }
    constexpr static char const * nodesetNamesString() { return "nodesetNames"; }
    constexpr static char const * partitionRefinementString() { return "partitionRefinement"; }
    constexpr static char const * partitionMethodString() { return "partitionMethod"; }
    constexpr static char const * useGlobalIdsString() { return "useGlobalIds"; }
  };
  /// @endcond

  void importVolumicFieldOnArray( string const & cellBlockName,
                                  string const & meshFieldName,
                                  bool isMaterialField,
                                  dataRepository::WrapperBase & wrapper ) const;

  void importSurfacicFieldOnArray( string const & faceBlockName,
                                   string const & meshFieldName,
                                   dataRepository::WrapperBase & wrapper ) const;

  /**
   * @brief The VTK mesh to be imported into GEOSX.
   * @note We keep this smart pointer as a member for use in @p importFields().
   */
  // TODO can we use unique_ptr to hold mesh?
  vtkSmartPointer< vtkDataSet > m_vtkMesh;

  /// Name of VTK dataset attribute used to mark regions
  string m_attributeName;

  /// Name of the main block to be imported (for multi-block files).
  string m_mainBlockName;

  /// Name of the face blocks to be imported (for multi-block files).
  array1d< string > m_faceBlockNames;

  /// Maps the face block name to its vtk mesh instance.
  std::map< string, vtkSmartPointer< vtkDataSet > > m_faceBlockMeshes;

  /// Names of VTK nodesets to import
  string_array m_nodesetNames;

  /// Number of graph partitioning refinement iterations
  integer m_partitionRefinement = 0;

  /// Whether global id arrays should be used, if available
  integer m_useGlobalIds = 0;

  /// Method (library) used to partition the mesh
  vtk::PartitionMethod m_partitionMethod = vtk::PartitionMethod::parmetis;

  /// Lists of VTK cell ids, organized by element type, then by region
  vtk::CellMapType m_cellMap;
};

} // namespace geos

#endif /* GEOS_MESH_GENERATORS_VTKMESHGENERATOR_HPP */
