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

#ifndef GEOSX_MESH_PARTICLEBLOCK_HPP_
#define GEOSX_MESH_PARTICLEBLOCK_HPP_

#include "dataRepository/Group.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/generators/ParticleBlockABC.hpp"
#include "mesh/ParticleType.hpp"

namespace geos
{

/**
 * This implementation of ParticleBlockABC mainly use the cell patterns/shapes
 * to build all the particle to nodes, faces and edges mappings.
 */
class ParticleBlock : public ParticleBlockABC
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  ParticleBlock() = delete;

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleBlock( string const & name, Group * const parent );

  /**
   * @brief Copy constructor.
   * @param[in] init the source to copy
   */
  ParticleBlock( const ParticleBlock & init ) = delete;

  /**
   * @brief Default destructor.
   */
  ~ParticleBlock();

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  ///@}
  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Defines the underlying particle type (hex, tet...)
   * @param[in] particleType the particle type
   *
   * @note Allocates the values of the particle to nodes, edges, faces accordingly.
   */
  void setParticleType( ParticleType particleType );

  ParticleType getParticleType() const override
  { return m_particleType; }

  array1d< globalIndex > getParticleID() const override
  { return m_particleID; }

  /**
   * @brief Sets the global IDs of particles in this subregion.
   * @param particleID The input list of global IDs
   */
  void setParticleID( array1d< globalIndex > const particleID )
  { m_particleID = particleID; }

  array2d< real64 > getParticleCenter() const override
  { return m_particleCenter; }

  /**
   * @brief Set the list of particle center locations in this subregion.
   * @param particleCenter The input list of particle center coordinates
   */
  void setParticleCenter( array2d< real64 > const particleCenter )
  { m_particleCenter = particleCenter; }

  array2d< real64 > getParticleVelocity() const override
  { return m_particleVelocity; }

  /**
   * @brief Set the list of particle velocities in this subregion.
   * @param particleVelocity The input list of velocities
   */
  void setParticleVelocity( array2d< real64 > const particleVelocity )
  { m_particleVelocity = particleVelocity; }

  array2d< real64 > getParticleAcceleration() const override
  { return m_particleAcceleration; }

  /**
   * @brief Set the list of particle accelerations in this subregion.
   * @param particleAcceleration The input list of accelerations
   */
  void setParticleAcceleration( array2d< real64 > const particleAcceleration )
  { m_particleAcceleration = particleAcceleration; }  

  array2d< real64 > getParticleInitialMaterialDirection() const override
  { return m_particleInitialMaterialDirection; }

  /**
   * @brief Set the list of material directions in this subregion.
   * @param particleInitialMaterialDirection The input list of initial directions
   */
  void setParticleInitialMaterialDirection( array2d< real64 > const particleInitialMaterialDirection )
  { m_particleInitialMaterialDirection = particleInitialMaterialDirection; }

  array2d< real64 > getParticleMaterialDirection() const override
  { return m_particleMaterialDirection; }

  /**
   * @brief Set the list of material directions in this subregion.
   * @param particleMaterialDirection The input list of directions
   */
  void setParticleMaterialDirection( array2d< real64 > const particleMaterialDirection )
  { m_particleMaterialDirection = particleMaterialDirection; }

  array1d< int > getParticleGroup() const override
  { return m_particleGroup; }

  /**
   * @brief Set the list of particle group numbers (for contact) in this subregion.
   * @param particleGroup The input list of contact group numbers.
   */
  void setParticleGroup( array1d< int > const particleGroup )
  { m_particleGroup = particleGroup; } // TODO: Rename to particleContactGroup

  array1d< int > getParticleSurfaceFlag() const override
  { return m_particleSurfaceFlag; }

  void setParticleSurfaceFlag( array1d< int > const particleSurfaceFlag )
  { m_particleSurfaceFlag = particleSurfaceFlag; }

  array1d< real64 > getParticleDamage() const override
  { return m_particleDamage; }

  /**
   * @brief Set the list of particle damage values in this subregion.
   * @param particleDamage The input list of damage values
   */
  void setParticleDamage( array1d< real64 > const particleDamage )
  { m_particleDamage = particleDamage; }

  array1d< real64 > getParticleStrengthScale() const override
  { return m_particleStrengthScale; }

  /**
   * @brief Set the list of particle strength scale values in this subregion.
   * @param particleDamage The input list of strength scale values
   */
  void setParticleStrengthScale( array1d< real64 > const particleStrengthScale )
  { m_particleStrengthScale = particleStrengthScale; }

  array1d< real64 > getParticleVolume() const override
  { return m_particleVolume; }

  /**
   * @brief Set the list of particle volumes in this subregion.
   * @param particleVolume The input list of volumes
   */
  void setParticleVolume( array1d< real64 > const particleVolume )
  { m_particleVolume = particleVolume; }

  array3d< real64 > getParticleRVectors() const override
  { return m_particleRVectors; }

  /**
   * @brief Set the list of particle r-vectors in this subregion.
   * @param particleRVectors The input list of r-vectors
   */
  void setParticleRVectors( array3d< real64 > const particleRVectors )
  { m_particleRVectors = particleRVectors; }

  bool hasRVectors() const override
  { return m_hasRVectors; }

  array2d< real64 > getParticleInitialSurfaceNormal() const override
  { return m_particleInitialSurfaceNormal; }

  /**
   * @brief Set the list of surface normals in this subregion.
   * @param particleInitialSurfaceNormal The input list of initial normals
   */
  void setParticleInitialSurfaceNormal( array2d< real64 > const particleInitialSurfaceNormal )
  { m_particleInitialSurfaceNormal = particleInitialSurfaceNormal; }

  array2d< real64 > getParticleSurfaceNormal() const override
  { return m_particleSurfaceNormal; }

  /**
   * @brief Set the list of surface normals in this subregion.
   * @param particleSurfaceNormal The input list of normals
   */
  void setParticleSurfaceNormal( array2d< real64 > const particleSurfaceNormal )
  { m_particleSurfaceNormal = particleSurfaceNormal; }

  localIndex numParticles() const override
  { return size(); }

  /**
   * @brief Get local to global map, non-const version.
   * @return The mapping relationship as a array.
   *
   * @deprecated This accessor is meant to be used like a setter even though it's a bit like having public attribute...
   * Use a real setter instead.
   */
  arrayView1d< globalIndex > localToGlobalMap()
  { return m_localToGlobalMap; }

  array1d< globalIndex > localToGlobalMap() const override
  { return m_localToGlobalMap; }

  /**
   * @brief Resize the cell block to hold @p numParticles
   * @param numParticles The new number of particles.
   */
  void resize( dataRepository::indexType const numParticles ) override final;

  ///@}

  /**
   * @name Properties
   */
  ///@{

  ///@}

private:

  /// Contains the global index of each object.
  array1d< globalIndex > m_localToGlobalMap;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /// Member level field for the particle global ID.
  array1d< globalIndex > m_particleID;

  /// Member level field for the particle contact group.
  array1d< int > m_particleGroup;

  /// Member level field for the particle surface flag.
  array1d< int > m_particleSurfaceFlag;

  /// Member level field for the particle damage.
  array1d< real64 > m_particleDamage;

  /// Member level field for the particle strength scale.
  array1d< real64 > m_particleStrengthScale;

  /// Member level field for the particle center.
  array2d< real64 > m_particleCenter;

  /// Member level field for the particle velocity.
  array2d< real64 > m_particleVelocity;

  /// Member level field for the particle acceleration.
  array2d< real64 > m_particleAcceleration;  

  /// Member level field for the particle initial material direction.
  array2d< real64 > m_particleInitialMaterialDirection;

  /// Member level field for the particle material direction.
  array2d< real64 > m_particleMaterialDirection;

  /// Member level field for the particle volume.
  array1d< real64 > m_particleVolume;

  /// Type of particles in this subregion.
  ParticleType m_particleType;

  /// Bool flag for whether particles in this block have r-vectors defining its domain
  bool m_hasRVectors;

  /// R-vectors
  array3d< real64 > m_particleRVectors;

  /// Member level field for the particle initial surface normal.
  array2d< real64 > m_particleInitialSurfaceNormal;

  /// Member level field for the particle surface normal.
  array2d< real64 > m_particleSurfaceNormal;

  std::list< dataRepository::WrapperBase * > getExternalProperties() override
  {
    std::list< dataRepository::WrapperBase * > result;
    for( string const & externalPropertyName : m_externalPropertyNames )
    {
      result.push_back( &this->getWrapperBase( externalPropertyName ) );
    }
    return result;
  }
};

}

#endif /* GEOSX_MESH_CELLBLOCK_HPP_ */
