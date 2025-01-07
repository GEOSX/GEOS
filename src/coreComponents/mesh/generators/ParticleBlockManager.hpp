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
 * @file ParticleBlockManager.hpp
 */

#ifndef GEOS_MESH_PARTICLEBLOCKMANAGER_H_
#define GEOS_MESH_PARTICLEBLOCKMANAGER_H_

#include "mesh/generators/ParticleBlock.hpp"
#include "mesh/generators/ParticleBlockManagerABC.hpp"

namespace geos
{

/**
 * @class ParticleBlockManager
 * @brief The ParticleBlockManager class provides an interface to ObjectManagerBase in order to manage ParticleBlock data.
 */
class ParticleBlockManager : public ParticleBlockManagerABC
{
public:

  /**
   * @brief Constructor for ParticleBlockManager object.
   * @param name name of this instantiation of ParticleBlockManager
   * @param parent pointer to the parent Group of this instantiation of ParticleBlockManager
   */
  ParticleBlockManager( string const & name, Group * const parent );

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  using Group::resize;

  /**
   * @brief Set the number of particles for a set of particle regions.
   * @param numParticles list of the new particle numbers
   * @param regionNames list of the particle region names
   */
  void resize( integer_array const & numParticles,
               string_array const & regionNames );

  /**
   * @brief Get particle block by name.
   * @param[in] name Name of the particle block.
   * @return Reference to the particle block instance.
   */
  ParticleBlock & getParticleBlock( string const & name )
  {
    return this->getGroup( viewKeyStruct::particleBlocks() ).getGroup< ParticleBlock >( name );
  }

  const Group & getParticleBlocks() const override;

  Group & getParticleBlocks() override;

  /**
   * @brief Registers and returns a particle block of name @p name.
   * @param name The name of the created particle block.
   * @return A reference to the new particle block. The ParticleBlockManager owns this new instance.
   */
  ParticleBlock & registerParticleBlock( string name );

  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forParticleSubRegions( LAMBDA lambda )
  {
    this->getGroup( viewKeyStruct::particleBlocks ).forSubGroups< ParticleBlock >( lambda );
  }

private:

  struct viewKeyStruct
  {
    /// Particle blocks key
    static constexpr char const * particleBlocks() { return "particleBlocks"; }
  };

  /**
   * @brief Get particle block at index @p iParticleBlock.
   * @param[in] iParticleBlock The particle block index.
   * @return Const reference to the instance.
   *
   * @note Mainly useful for iteration purposes.
   */
  ParticleBlock const & getParticleBlock( localIndex const blockIndex ) const;

  /**
   * @brief Returns the number of particles blocks
   * @return Number of particle blocks
   */
  localIndex numParticleBlocks() const;

};

}
#endif /* GEOS_MESH_PARTICLEBLOCKMANAGER_H_ */
