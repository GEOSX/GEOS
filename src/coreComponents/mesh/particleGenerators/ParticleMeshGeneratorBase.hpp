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

#ifndef GEOS_MESH_PARTICLEGENERATORS_PARTICLEMESHGENERATORBASE_HPP
#define GEOS_MESH_PARTICLEGENERATORS_PARTICLEMESHGENERATORBASE_HPP

#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include "dataRepository/Group.hpp"
#include "dataRepository/WrapperBase.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

class ParticleMeshGeneratorBase : public dataRepository::Group
{
public:

  /**
   * @brief Main constructor for ParticleMeshGeneratorBase base class.
   * @param[in] name of the ParticleMeshGeneratorBase object
   * @param[in] parent the parent Group pointer for the ParticleMeshGeneratorBase object
   */
  explicit ParticleMeshGeneratorBase( string const & name,
                                      Group * const parent );

  /**
   * @brief Return the name of the MeshGenerator in object catalog.
   * @return string that contains the catalog name of the MeshGenerator
   */
  static string catalogName() { return "ParticleMeshGeneratorBase"; }

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

  /// using alias for templated Catalog meshGenerator type
  using CatalogInterface = dataRepository::CatalogInterface< ParticleMeshGeneratorBase, string const &, Group * const >;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

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
  void generateMesh( Group & parent, SpatialPartition & partition );

private:

  /**
   * @brief Fill the particleBlockManager object .
   * @param[inout] particleBlockManager the particleBlockManager that will receive the meshing information
   * @param[in] particleManager The reference to the particle manager
   * @param[in] partition The reference to spatial partition
   */
  virtual void fillParticleBlockManager( ParticleBlockManager & particleBlockManager, ParticleManager & particleManager, SpatialPartition const & partition ) = 0;
};
}

#endif /* GEOS_MESH_PARTICLEGENERATORS_PARTICLEMESHGENERATORBASE_HPP */
