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
 * @file ParticleSubRegionBase.hpp
 */

#ifndef GEOS_MESH_PARTICLESUBREGIONBASE_HPP_
#define GEOS_MESH_PARTICLESUBREGIONBASE_HPP_

#include "mesh/ParticleType.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "ToParticleRelation.hpp"

namespace geos
{

class ParticleManager;
class MeshLevel;

namespace constitutive
{
class ConstitutiveBase;
}

/**
 * @class ParticleSubRegionBase
 * Abstract class for a collection of mesh particles that
 * will be derived and specialized for particles
 */
class ParticleSubRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleSubRegionBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Destructor.
   */
  ~ParticleSubRegionBase();

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get whether particle has r-vectors.
   * @return bool indicating whether particle has r-vectors.
   */
  bool hasRVectors() const
  { return m_hasRVectors; }

  /**
   * @brief Set whether particle has r-vectors.
   * @param[in] hasRVectors bool indicating whether particle has r-vectors.
   */
  void setHasRVectors( bool hasRVectors )
  { m_hasRVectors = hasRVectors; }

  /**
   * @brief Get the global ID of each particle in this subregion.
   * @return an arrayView1d of const particle global IDs
   */
  arrayView1d< globalIndex const > getParticleID() const
  { return m_particleID; }

  /**
   * @copydoc getParticleID() const
   */
  arrayView1d< globalIndex > getParticleID()
  { return m_particleID; }

  /**
   * @brief Get the contact group of each particle in this subregion.
   * @return an arrayView1d of const particle contact groups
   */
  arrayView1d< int const > getParticleGroup() const
  { return m_particleGroup; }

  /**
   * @copydoc getParticleGroup() const
   */
  arrayView1d< int > getParticleGroup()
  { return m_particleGroup; }

  /**
   * @brief Get the contact group of each particle in this subregion.
   * @return an arrayView1d of const particle surface flags
   */
  arrayView1d< int const > getParticleSurfaceFlag() const
  { return m_particleSurfaceFlag; }

  /**
   * @copydoc getParticleSurfaceFlag() const
   */
  arrayView1d< int > getParticleSurfaceFlag()
  { return m_particleSurfaceFlag; }

  /**
   * @brief Get the damage of each particle in this subregion.
   * @return an arrayView1d of const particle damage
   */
  arrayView1d< real64 const > getParticleDamage() const
  { return m_particleDamage; }

  /**
   * @copydoc getParticleDamage() const
   */
  arrayView1d< real64 > getParticleDamage()
  { return m_particleDamage; }

  /**
   * @brief Get the porosity of each particle in this subregion.
   * @return an arrayView1d of const particle porosity
   */
  arrayView1d< real64 const > getParticlePorosity() const
  { return m_particlePorosity; }

  /**
   * @copydoc getParticlePorosity() const
   */
  arrayView1d< real64 > getParticlePorosity()
  { return m_particlePorosity; }

    /**
   * @brief Get the temperature of each particle in this subregion.
   * @return an arrayView1d of const particle temperature
   */
  arrayView1d< real64 const > getParticleTemperature() const
  { return m_particleTemperature; }

  /**
   * @copydoc getParticleTemperature() const
   */
  arrayView1d< real64 > getParticleTemperature()
  { return m_particleTemperature; }

  /**
   * @brief Get the strength scale of each particle in this subregion.
   * @return an arrayView1d of const particle strength scale
   */
  arrayView1d< real64 const > getParticleStrengthScale() const
  { return m_particleStrengthScale; }

  /**
   * @copydoc getParticleStrengthScale() const
   */
  arrayView1d< real64 > getParticleStrengthScale()
  { return m_particleStrengthScale; }


  /**
   * @brief Get the ghost rank of each particle in this subregion.
   * @return an arrayView1d of const particle ghost ranks
   */
  arrayView1d< int const > getParticleRank() const
  { return m_particleRank; }

  /**
   * @copydoc getParticleRank() const
   */
  arrayView1d< int > getParticleRank()
  { return m_particleRank; }

  /**
   * @brief Get the center of each particle in this subregion.
   * @return an arrayView1d of const particle centers
   */
  arrayView2d< real64 const > getParticleCenter() const
  { return m_particleCenter; }

  /**
   * @copydoc getParticleCenter() const
   */
  arrayView2d< real64 > getParticleCenter()
  { return m_particleCenter; }

  /**
   * @brief Set the center of each particle in this subregion.
   * @param particleCenter the list of coordinates to be used as particle centers
   */
  void setParticleCenter( array2d< real64 > particleCenter )
  { m_particleCenter = particleCenter; }

  /**
   * @brief Get the velocity of each particle in this subregion.
   * @return an arrayView1d of const particle velocities
   */
  arrayView2d< real64 const > getParticleVelocity() const
  { return m_particleVelocity; }

  /**
   * @copydoc getParticleVelocity() const
   */
  arrayView2d< real64 > getParticleVelocity()
  { return m_particleVelocity; }

  /**
   * @brief Get the material direction of each particle in this subregion.
   * @return an arrayView1d of const particle material direction
   */
  arrayView2d< real64 const > getParticleMaterialDirection() const
  { return m_particleMaterialDirection; }

  /**
   * @copydoc getParticleMaterialDirection() const
   */
  arrayView2d< real64 > getParticleMaterialDirection()
  { return m_particleMaterialDirection; }

  /**
   * @brief Get the current volume of each particle in this subregion.
   * @return an arrayView1d of const current particle volumes
   */
  arrayView1d< real64 const > getParticleVolume() const
  { return m_particleVolume; }

  /**
   * @copydoc getParticleVolume() const
   */
  arrayView1d< real64 > getParticleVolume()
  { return m_particleVolume; }

  /**
   * @brief Get the r-vectors of each particle in this subregion.
   * @return an arrayView3d of const particle r-vectors
   */
  arrayView3d< real64 const > getParticleRVectors() const
  { return m_particleRVectors; }

  /**
   * @copydoc getParticleRVectors() const
   */
  arrayView3d< real64 > getParticleRVectors()
  { return m_particleRVectors; }

  /**
   * @brief Get the surface normal of each particle in this subregion.
   * @return an arrayView1d of const particle surface normal
   */
  arrayView2d< real64 const > getParticleSurfaceNormal() const
  { return m_particleSurfaceNormal; }

  /**
   * @copydoc getParticleSurfaceNormal() const
   */
  arrayView2d< real64 > getParticleSurfaceNormal()
  { return m_particleSurfaceNormal; }

  /**
   * @brief Get the surface position of each particle in this subregion.
   * @return an arrayView1d of const particle surface position
   */
  arrayView2d< real64 const > getParticleSurfacePosition() const
  { return m_particleSurfacePosition; }

  /**
   * @copydoc getParticleSurfacePosition() const
   */
  arrayView2d< real64 > getParticleSurfacePosition()
  { return m_particleSurfacePosition; }

  /**
   * @brief Get the surface traction of each particle in this subregion.
   * @return an arrayView1d of const particle surface traction
   */
  arrayView2d< real64 const > getParticleSurfaceTraction() const
  { return m_particleSurfaceTraction; }

  /**
   * @copydoc getParticleSurfaceTraction() const
   */
  arrayView2d< real64 > getParticleSurfaceTraction()
  { return m_particleSurfaceTraction; }

  /**
   * @brief Get the group in which the constitutive models of this subregion are registered.
   * @return a pointer to the const group in which the constitutive models are registered
   */
  dataRepository::Group const & getConstitutiveModels() const
  { return m_constitutiveModels; }

  /**
   * @copydoc getConstitutiveModels() const
   */
  dataRepository::Group & getConstitutiveModels()
  { return m_constitutiveModels; }

  /**
   * @brief Get a pointer to the constitutive model.
   * @tparam T The type of the constitutive model.
   * @param name The name of the constitutive model.
   * @return A pointer to the constitutive model.
   */
  template< typename T = constitutive::ConstitutiveBase >
  T const & getConstitutiveModel( string const & name ) const
  { return m_constitutiveModels.getGroup< T >( name ); }

  /**
   * @copydoc getConstitutiveModel( string const & ) const
   */
  template< typename T = constitutive::ConstitutiveBase >
  T & getConstitutiveModel( string const & name )
  { return m_constitutiveModels.getGroup< T >( name ); }

  /**
   * @brief Get the type of particle in this subregion.
   * @return the type of particle in this subregion
   */
  ParticleType getParticleType() const
  { return m_particleType; }

  /**
   * @brief Set the type of particle in this subregion.
   * @param[in] particleType the particle type
   */
  virtual void setParticleType( ParticleType const particleType )
  { m_particleType = particleType; }

  /**
   * @brief Get the number of vertices for particles in this subregion.
   * @return the type of particle in this subregion
   */
  int numberOfVerticesPerParticle() const
  { return m_numVerticesPerParticle; }

  /**
   * @brief Set the number of vertices for particles in this subregion.
   * @param[in] particleType the particle type
   */
  virtual void setNumberOfVerticesPerParticle( ParticleType const particleType )
  {
    switch( particleType )
    {
      case ParticleType::SinglePoint:
      {
        m_numVerticesPerParticle = 1;
        break;
      }
      case ParticleType::CPDI:
      {
        m_numVerticesPerParticle = 8;
        break;
      }
      default:
      {
        GEOS_ERROR( "Particle type \"" << m_particleType << "\" is not yet supported." );
        break;
      }
    }
  }

  /**
   * @brief Provide an immutable accessor to the particle neighbor list.
   * @return const reference to neighbor list
   */
  OrderedVariableToManyParticleRelation const & neighborList() const
  { return m_neighborList; }

  /**
   * @brief Get a mutable accessor to the particle neighbor list.
   * @return reference to neighbor list
   */
  OrderedVariableToManyParticleRelation & neighborList()
  { return m_neighborList; }

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the member level field for the particle global ID.
    static constexpr char const * particleIDString() { return "particleID"; }

    /// @return String key for the member level field for the particle contact group.
    static constexpr char const * particleGroupString() { return "particleGroup"; }

    /// @return String key for the member level field for the particle surface flag.
    static constexpr char const * particleSurfaceFlagString() { return "particleSurfaceFlag"; }

    /// @return String key for the member level field for the particle damage.
    static constexpr char const * particleDamageString() { return "particleDamage"; }

    /// @return String key for the member level field for the particle porosity.
    static constexpr char const * particlePorosityString() { return "particlePorosity"; }

        /// @return String key for the member level field for the particle temperature.
    static constexpr char const * particleTemperatureString() { return "particleTemperature"; }

    /// @return String key for the member level field for the particle strength scale.
    static constexpr char const * particleStrengthScaleString() { return "particleStrengthScale"; }

    /// @return String key for the member level field for the particle ghost rank.
    static constexpr char const * particleRankString() { return "particleRank"; }

    /// @return String key for the member level field for the particle center.
    static constexpr char const * particleCenterString() { return "particleCenter"; }

    /// @return String key for the member level field for the particle velocity.
    static constexpr char const * particleVelocityString() { return "particleVelocity"; }

    /// @return String key for the member level field for the particle material direction.
    static constexpr char const * particleMaterialDirectionString() { return "particleMaterialDirection"; }

    /// @return String key for the member level field for the current particle volume.
    static constexpr char const * particleVolumeString() { return "particleVolume"; }

    /// @return String key for the member level field for the particle volume.
    static constexpr char const * particleRVectorsString() { return "particleRVectors"; }
  
    /// @return String key for the member level field for the particle surface normal.
    static constexpr char const * particleSurfaceNormalString() { return "particleSurfaceNormal"; }

    /// @return String key for the member level field for the particle surface position.
    static constexpr char const * particleSurfacePositionString() { return "particleSurfacePosition"; }
  
    /// @return String key for the member level field for the particle surface traction.
    static constexpr char const * particleSurfaceTractionString() { return "particleSurfaceTraction"; }
  };

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// @return String key for the group in which the constitutive models of this subregion are registered.
    static constexpr auto constitutiveModelsString() { return "ConstitutiveModels"; }
  };

  /**
   * @brief Packs the data of particles on this subregion
   * @param buffer The buffer to pack into
   * @param localIndices the local indices of the particles to pack
   * @param doPack whether to actually pack (1) or simply check the size of the packed contents (0)
   * @return the size of the packed contents
   */
  unsigned int particlePack( buffer_type & buffer,
                             arrayView1d< localIndex > const & localIndices,
                             bool doPack ) const;

  /**
   * @brief Unpacks the data of particles on this subregion
   * @param buffer The buffer to unpack from
   * @param startingIndex The local index where unpacking should begin
   * @param numberOfIncomingParticles the number of incoming particles
   */
  void particleUnpack( buffer_type & buffer,
                       int const & startingIndex,
                       int const & numberOfIncomingParticles );

  /**
   * @brief Erases particle field data at the provided local indices
   * @param indicesToErase the local particle indices whose data should be erased
   */
  void erase( std::set< localIndex > const & indicesToErase );

  /**
   * @brief Identifies the local indices of non-ghost particles
   */
  void setActiveParticleIndices();

  /**
   * @brief Returns the local indices of all non-ghost particles
   * @return the local indices of all non-ghost particles
   */
  SortedArrayView< localIndex const > const activeParticleIndices() const
  {
    return m_activeParticleIndices.toView();
  }

  /**
   * @brief Returns the local indices of all ghost particles
   * @return the local indices of all ghost particles
   */
  SortedArrayView< localIndex const > const inactiveParticleIndices() const
  {
    return m_inactiveParticleIndices.toView();
  }

  /**
   * @brief Updates the globalToLocal and localToGlobal maps
   */
  void updateMaps();

private:
  /// Group in which the constitutive models of this subregion are registered
  dataRepository::Group m_constitutiveModels;

protected:
  /// Boolean indicating whether the particle subregion contains particles needing r-vectors defining their domain extent.
  bool m_hasRVectors;

  /// The number of vertices each particle has
  int m_numVerticesPerParticle;

  /// Member level field for particle ghost ranks
  array1d< int > m_particleRank;

  /// Member level field for the particle global ID.
  array1d< globalIndex > m_particleID;

  /// Member level field for the particle contact group.
  array1d< int > m_particleGroup;

  /// Member level field for the particle surface flag.
  array1d< int > m_particleSurfaceFlag;

  /// Member level field for the particle damage.
  array1d< real64 > m_particleDamage;

  /// Member level field for the particle porosity.
  array1d< real64 > m_particlePorosity;

  /// Member level field for the particle temperature.
  array1d< real64 > m_particleTemperature;

  /// Member level field for the particle strength scale.
  array1d< real64 > m_particleStrengthScale;

  /// Member level field for the particle center.
  array2d< real64 > m_particleCenter;

  /// Member level field for the particle velocity.
  array2d< real64 > m_particleVelocity;

  /// Member level field for the particle material direction.
  array2d< real64 > m_particleMaterialDirection;

  /// Member level field for the current particle volume.
  array1d< real64 > m_particleVolume;

  /// Type of particles in this subregion.
  ParticleType m_particleType;

  /// current half-R-vectors (center to face)
  array3d< real64 > m_particleRVectors;

  /// Member level field for the particle surface normal.
  array2d< real64 > m_particleSurfaceNormal;

  /// Member level field for the particle surface position.
  array2d< real64 > m_particleSurfacePosition;

  /// Member level field for the particle surface traction.
  array2d< real64 > m_particleSurfaceTraction;

  /// Indices of particles that are not ghosts
  SortedArray< localIndex > m_activeParticleIndices;

  /// Indices of ghost particles
  SortedArray< localIndex > m_inactiveParticleIndices;

  /// Neighbor list
  OrderedVariableToManyParticleRelation m_neighborList;

};


} /* namespace geos */

#endif /* GEOS_MESH_PARTICLESUBREGIONBASE_HPP_ */
