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

#ifndef GEOS_PARTICLEBLOCKABC_HPP
#define GEOS_PARTICLEBLOCKABC_HPP

#include "mesh/ParticleType.hpp"
#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"

#include <vector>

namespace geos
{

/**
 * Abstract base class defining the information provided by any particle block implementation.
 * Mainly particle type-related information (TODO put an example here),
 * a local to global mapping and some kind of accessors to external properties.
 *
 * The particles of ParticleBlockABC are all of the same kind.
 *
 * It's noteworthy that the ParticleBlockABC is immutable oriented.
 * The derived implementations need to have the modification/creation capabilities.
 */
class ParticleBlockABC : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  ParticleBlockABC( string const & name,
                    Group * const parent )
    :
    Group( name, parent )
  { }

  /**
   * @brief Get the type of particle in this subregion.
   * @return an enum specifying the type of particle in this subregion
   */
  virtual ParticleType getParticleType() const = 0;

  /**
   * @brief Get the list of particle global IDs in this subregion.
   * @return the list of particle global IDs in this subregion.
   */
  virtual array1d< globalIndex > getParticleID() const = 0;

  /**
   * @brief Get the list of particle group numbers (for contact) in this subregion.
   * @return the list of particle group numbers (for contact) in this subregion.
   */
  virtual array1d< int > getParticleGroup() const = 0;

  /**
   * @brief Get the list of particle surface flag values in this subregion.
   * @return the list of particle surface flag values in this subregion.
   */
  virtual array1d< int > getParticleSurfaceFlag() const = 0;

  /**
   * @brief Get the list of particle damage values in this subregion.
   * @return the list of particle damage values in this subregion.
   */
  virtual array1d< real64 > getParticleDamage() const = 0;

  /**
   * @brief Get the list of particle porosity values in this subregion.
   * @return the list of particle porosity values in this subregion.
   */
  virtual array1d< real64 > getParticlePorosity() const = 0;

  /**
   * @brief Get the list of particle temperature values in this subregion.
   * @return the list of particle temperature values in this subregion.
   */
  virtual array1d< real64 > getParticleTemperature() const = 0;

  /**
   * @brief Get the list of particle strength scale values in this subregion.
   * @return the list of particle strength scale values in this subregion.
   */
  virtual array1d< real64 > getParticleStrengthScale() const = 0;

  /**
   * @brief Get the list of particle center locations in this subregion.
   * @return the list of particle center locations in this subregion.
   */
  virtual array2d< real64 > getParticleCenter() const = 0;

  /**
   * @brief Get the list of particle velocities in this subregion.
   * @return the list of particle velocities in this subregion.
   */
  virtual array2d< real64 > getParticleVelocity() const = 0;

  /**
   * @brief Get the list of particle material directions in this subregion.
   * @return the list of particle material directions in this subregion.
   */
  virtual array2d< real64 > getParticleInitialMaterialDirection() const = 0;

  /**
   * @brief Get the list of particle material directions in this subregion.
   * @return the list of particle material directions in this subregion.
   */
  virtual array2d< real64 > getParticleMaterialDirection() const = 0;

  /**
   * @brief Get the list of particle volumes in this subregion.
   * @return the list of particle volumes in this subregion.
   */
  virtual array1d< real64 > getParticleVolume() const = 0;

  /**
   * @brief Get the list of particle r-vectors in this subregion.
   * @return the list of particle r-vectors in this subregion.
   */
  virtual array3d< real64 > getParticleRVectors() const = 0;

  /**
   * @brief Query whether this subregion has particles of a type that have r-vectors that depend on the deformation gradient
   * @return a bool indicating whether particles on this subregion have r-vectors that depend on the deformation gradient
   */
  virtual bool hasRVectors() const = 0;

  /**
   * @brief Get the list of particle initial surface normal in this subregion.
   * @return the list of particle initial surface normal in this subregion.
   */
  virtual array2d< real64 > getParticleInitialSurfaceNormal() const = 0;

  /**
   * @brief Get the list of particle surface normal in this subregion.
   * @return the list of particle surface normal in this subregion.
   */
  virtual array2d< real64 > getParticleSurfaceNormal() const = 0;

  /**
   * @brief Get the list of particle initial surface position in this subregion.
   * @return the list of particle initial surface position in this subregion.
   */
  virtual array2d< real64 > getParticleInitialSurfacePosition() const = 0;

  /**
   * @brief Get the list of particle surface position in this subregion.
   * @return the list of particle surface position in this subregion.
   */
  virtual array2d< real64 > getParticleSurfacePosition() const = 0;

  /**
   * @brief Get the list of particle initial surface traction in this subregion.
   * @return the list of particle initial surface traction in this subregion.
   */
  virtual array2d< real64 > getParticleInitialSurfaceTraction() const = 0;

  /**
   * @brief Get the list of particle surface traction in this subregion.
   * @return the list of particle surface traction in this subregion.
   */
  virtual array2d< real64 > getParticleSurfaceTraction() const = 0;

  /**
   * @brief Get the number of particles.
   * @return number of particles in the particle block
   */
  virtual localIndex numParticles() const = 0;

  /**
   * @brief Get local to global map.
   * @return The mapping relationship as an array.
   */
  virtual array1d< globalIndex > localToGlobalMap() const = 0;

  /**
   * @brief Helper function to apply a lambda function over all the external properties of the subregion
   * @tparam LAMBDA the type of the lambda function
   * @param lambda lambda function that is applied to the wrappers of external properties
   *
   * @note Unlike the other member functions of this class, this current member function is not abstract,
   * mainly because it's a template method. The abstraction is delegated to private method @p getExternalProperties.
   * @note There is some `constness` concern for this member function.
   * @see getExternalProperties()
   */
  template< typename LAMBDA >
  void forExternalProperties( LAMBDA && lambda )
  {
    for( auto * wrapperBase: this->getExternalProperties() )
    {
      lambda( *wrapperBase );
    }
  }

private:
  /**
   * @brief Returns the external properties under the form of WrapperBase pointers.
   * @return An iterable of pointers
   *
   * In order not to expose the implementation details (external properties being stored as WrapperBase),
   * this abstract member function is made private.
   * Thus the list of pointers shall not be used anyhow by end-users that must use @p forExternalProperties.
   * @note There is some `constness` concern for this member function.
   * @see forExternalProperties(LAMBDA && lambda)
   */
  virtual std::list< dataRepository::WrapperBase * > getExternalProperties() = 0;
};

}

#endif //GEOS_CELLBLOCKABC_HPP
