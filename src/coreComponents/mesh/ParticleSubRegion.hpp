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


#ifndef GEOS_MESH_PARTICLEELEMENTSUBREGION_HPP_
#define GEOS_MESH_PARTICLEELEMENTSUBREGION_HPP_

#include "mesh/generators/ParticleBlockABC.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "ParticleSubRegionBase.hpp"


namespace geos
{

class MeshLevel;

/**
 * @class ParticleSubRegion
 * Class specializing the particle subregion for a cell particle.
 * This is the class used in the physics solvers to represent a collection of mesh cell particles
 */
class ParticleSubRegion : public ParticleSubRegionBase
{
public:

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName()
  { return "ParticleSubRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return ParticleSubRegion::catalogName(); }

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleSubRegion( string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~ParticleSubRegion() override;

  ///@}

  /**
   * @name Helpers for ParticleSubRegion construction
   */
  ///@{

  /**
   * @brief Set the ghost rank of particles in this subregion
   * @param rank the mpi rank to which all particles in this subregion (in this partition) will have their ghost rank set to
   */
  void setParticleRank( int const rank );

  /**
   * @brief Fill the ParticleSubRegion by copying those of the source ParticleBlock
   * @param particleBlock the ParticleBlock which properties (connectivity info) will be copied.
   */
  void copyFromParticleBlock( ParticleBlockABC & particleBlock );

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Helper function to apply a lambda function over all constructive groups
   * @tparam LAMBDA the type of the lambda function
   * @param lambda the lambda function
   */
  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {
    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  ///@}

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ParticleSubRegionBase::viewKeyStruct
  {
    /// @return String key for the constitutive grouping
    static constexpr char const * constitutiveGroupingString() { return "ConstitutiveGrouping"; }
    /// @return String key for the constitutive map
    static constexpr char const * constitutiveMapString() { return "ConstitutiveMap"; }

    /// ViewKey for the constitutive grouping
    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString() };
    /// ViewKey for the constitutive map
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString() };
  }
  /// viewKey struct for the ParticleSubRegion class
  m_ParticleBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_ParticleBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_ParticleBlockSubRegionViewKeys; }

private:

  /// Map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

};

} /* namespace geos */

#endif /* GEOS_MESH_CELLELEMENTSUBREGION_HPP_ */
