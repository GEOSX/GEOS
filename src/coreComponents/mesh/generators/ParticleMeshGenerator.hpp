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

/**
 * @file ParticleMeshGenerator.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_PARTICLEMESHGENERATOR_HPP
#define GEOSX_MESH_GENERATORS_PARTICLEMESHGENERATOR_HPP

#include "mesh/generators/MeshGeneratorBase.hpp"

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

class ParticleManager;
class SpatialPartition;

/**
 * @class ParticleMeshGenerator
 * @brief The ParticleMeshGenerator class is a class handling import of particle data from an externel particle file.
 */
class ParticleMeshGenerator : public MeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for ParticleMeshGenerator.
   * @param[in] name of the ParticleMeshGenerator
   * @param[in] parent point to the parent Group of the ParticleMeshGenerator
   */
  ParticleMeshGenerator( const string & name, Group * const parent );

  virtual ~ParticleMeshGenerator() override = default;

  /**
   * @brief Return the name of the ParticleMeshGenerator in object Catalog.
   * @return string that contains the key name to ParticleMeshGenerator in the Catalog
   */
  static string catalogName() { return "ParticleMesh"; }

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void fillParticleBlockManager( ParticleBlockManager & particleBlockManager, ParticleManager & particleManager, SpatialPartition const & partition ) override;

  void importFieldOnArray( Block block,
                           string const & blockName,
                           string const & meshFieldName,
                           bool isMaterialField,
                           dataRepository::WrapperBase & wrapper ) const override;

protected:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * particleFilePathString() { return "particleFile"; }
    constexpr static char const * headerFilePathString() { return "headerFile"; }
    constexpr static char const * particleBlockNamesString() { return "particleBlockNames"; }
    constexpr static char const * particleTypesString() { return "particleTypes"; }
  };
  /// @endcond

  void postInputInitialization() override;

  /// Mesh number of dimension
  int m_dim;

  /// Minimum extent of particle coordinates
  real64 m_min[3];

  /// Maximum extent of particle coordinates
  real64 m_max[3];

private:

  /// Path to the particle file
  Path m_particleFilePath;

  /// Path to the header file (problem and particle metadata)
  Path m_headerFilePath;

  /// String array of particle block names associated with the particle mesh
  array1d< string > m_blockNames;

  /// String array listing the particle types present
  array1d< string > m_particleType;

public:

};

} /* namespace geos */

#endif /* GEOSX_MESH_GENERATORS_PARTICLEMESHGENERATOR_HPP */
