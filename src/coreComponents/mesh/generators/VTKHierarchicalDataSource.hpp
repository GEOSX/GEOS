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
 * @file VTKHierarchicalDataRepository.hpp
 */

#ifndef GEOS_MESH_GENERATORS_VTKHIERARCHICALDATAREPOSITORY_HPP
#define GEOS_MESH_GENERATORS_VTKHIERARCHICALDATAREPOSITORY_HPP

#include "dataRepository/Group.hpp"
#include "mesh/ExternalDataRepositoryBase.hpp"

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkDataAssembly.h>
#include <vtkPartitionedDataSet.h>
#include <vtkPartitionedDataSetCollection.h>

namespace geos
{

/**
 * @class VTKHierarchicalDataRepository
 * @brief This class provides an API to access VTKPartitionedDataSetCollection through a vtkDataAssembly
 */
class VTKHierarchicalDataRepository : public ExternalDataRepositoryBase
{
public:

  /**
   * @brief Main constructor for VTKHierarchicalDataRepository base class.
   * @param[in] name of the VTKHierarchicalDataRepository object
   * @param[in] parent the parent Group pointer for the VTKHierarchicalDataRepository object
   */
  VTKHierarchicalDataRepository( string const & name, Group * const parent );

  virtual ~VTKHierarchicalDataRepository() override = default;

  /**
   * @brief Return the name of the MeshGenerator in object catalog.
   * @return string that contains the catalog name of the MeshGenerator
   */
  static string catalogName() { return "VTKHierarchicalDataRepository"; }

  /**
   * @brief Opens a vtkPartitionedDataSetCollection and gets the colletion and the associated dataAssembly
   *
   */
  void open() override;

  /**
   * @brief Performs a search in the dataAssembly to find a node of PartitionedDataSets
   *
   * @param path the path in the data assembly tree
   * @return the found dataset
   */
  vtkSmartPointer< vtkPartitionedDataSet > search( string const & path );

private:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
  };
  /// @endcond

  /// Path to the mesh file
  Path m_filePath;

  /// DataAssembly to query the dataset collection
  vtkSmartPointer< vtkDataAssembly > m_dataAssembly;

  /// Collection of datasets
  vtkSmartPointer< vtkPartitionedDataSetCollection > m_collection;
};

}

#endif
