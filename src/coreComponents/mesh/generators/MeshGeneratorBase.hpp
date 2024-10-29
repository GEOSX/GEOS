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
 * @file MeshGeneratorBase.hpp
 */

#ifndef GEOS_MESH_GENERATORS_MESHGENERATORBASE_HPP
#define GEOS_MESH_GENERATORS_MESHGENERATORBASE_HPP

#include "PartitionDescriptor.hpp"
#include "PartitionDescriptorABC.hpp"

#include "dataRepository/Group.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

// This forward declaration prevents from exposing the internals of the module,
// which are only accessed through some private functions signatures.
// In order to avoid this forward declaration, we could expose an ABC
// instead of exposing the MeshGeneratorBase implementation.
class CellBlockManager;

/**
 *  @class MeshGeneratorBase
 *  @brief The MeshGeneratorBase class provides an abstract base class implementation for different mesh types.
 *	   The MeshGeneratorBase is the Group specialization for different type of mesh handling.
 */
class MeshGeneratorBase : public dataRepository::Group
{
public:

  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the MeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  explicit MeshGeneratorBase( string const & name,
                              Group * const parent );

  /// This function is used to expand any catalogs in the data structure
  void expandObjectCatalogs() override;

  /// using alias for templated Catalog meshGenerator type
  using CatalogInterface = dataRepository::CatalogInterface< MeshGeneratorBase, string const &, Group * const >;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Accessor for the singleton Catalog object
   * @return a static reference to the Catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Generate the mesh object the input mesh object.
   * @param parent The parent group of the CellBlockManager.
   * @param[in] partition The reference to spatial partition
   */
  void generateMesh( Group & parent, array1d< int > const & partition );

  /**
   * @brief Describe which kind of block must be considered.
   */
  enum struct Block
  {
    VOLUMIC,
    SURFACIC,
    LINEIC
  };

  /**
   * @brief import field from the mesh on the array accessible via the given wrapper.
   * @param block Type of block to import from.
   * @param blockName name of the block to copy data from.
   * @param meshFieldName name of the field in the meshd
   * @param isMaterialField Indicate if we want to import material or regular fields
   * @param wrapper Wrapper to access the array
   */
  virtual void importFieldOnArray( Block block,
                                   string const & blockName,
                                   string const & meshFieldName,
                                   bool isMaterialField,
                                   dataRepository::WrapperBase & wrapper ) const = 0;

  /**
   * @brief Free internal resources associated with mesh/data import.
   *
   * This is relevant for mesh generators that load external mesh files.
   * Once this method is called, they can release any memory allocated.
   */
  virtual void freeResources() {}

  /**
   * @brief Get the name mapping between mesh volumic field names and internal GEOSX volumic field names.
   * @return The string to string mapping of field names.
   */
  std::map< string, string > const & getVolumicFieldsMapping() const { return m_volumicFields; }

  /**
   * @brief Get the name mapping between mesh surfacic field names and internal GEOSX surfacic field names.
   * @return The string to string mapping of field names.
   */
  std::map< string, string > const & getSurfacicFieldsMapping() const { return m_surfacicFields; }

  PartitionDescriptorABC const & getPartitionDescriptor() const
  {
    return m_partition;
  }

protected:
  /// Mapping from volumic field source to GEOSX field.
  std::map< string, string > m_volumicFields;

  /// Mapping from surfacic field source to GEOSX field.
  std::map< string, string > m_surfacicFields;

  /// The partition information
  PartitionDescriptor m_partition;

private:
  /**
   * @brief Fill the cellBlockManager object .
   * @param[inout] cellBlockManager the CellBlockManager that will receive the meshing information
   * @param[in] partition The (x, y , y) MPI split (in case we need it)
   */
  virtual void fillCellBlockManager( CellBlockManager & cellBlockManager, array1d< int > const & partition ) = 0;

  void attachWellInfo( CellBlockManager & cellBlockManager );
};
}

#endif /* GEOS_MESH_GENERATORS_MESHGENERATORBASE_HPP */
