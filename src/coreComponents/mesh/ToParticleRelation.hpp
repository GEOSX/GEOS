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
 * @file ToParticleRelation.hpp
 */

#ifndef GEOS_MESH_TOPARTICLERELATION_HPP_
#define GEOS_MESH_TOPARTICLERELATION_HPP_

#include "InterObjectRelation.hpp"

namespace geos
{

class ParticleManager;

/**
 * @brief A relationship to a particle.
 * @tparam BASETYPE The underlying relation type to use to
 *                  store the relationsip information.
 */
template< typename BASETYPE >
class ToParticleRelation
{
public:

  /// The type of the underlying relationship storage object.
  using base_type = BASETYPE;

  /**
   * @brief Resize the underlying relationship storage.
   * @tparam DIMS The types of each dimensions resize parameter.
   * @param newdims A parameter pack of appropriate size to resize each
   *                dimension of the relationship storage.
   */
  template< typename ... DIMS >
  void resize( DIMS... newdims )
  {
    m_toParticleRegion.resize( newdims ... );
    m_toParticleSubRegion.resize( newdims ... );
    m_toParticleIndex.resize( newdims ... );
    m_numParticles.resizeDefault( m_toParticleRegion.size(), 0 );
  }

  /**
   * @brief Get the current size of the relationship storage.
   * @return The current size of the relationship storage.
   */
  localIndex size() const
  {
    return m_toParticleRegion.size();
  }

  /**
   * @brief Get the size of a specific dimension of the relationship storage.
   * @param dim The dimension to get the storage size of.
   * @return The dimension size
   */
  localIndex size( int const dim ) const
  {
    return m_toParticleRegion.size( dim );
  }

  /**
   * @brief Set the ParticleRegionManager.
   * @param input The ParticleRegionManager to set.
   */
  void setParticleManager( ParticleManager const & input )
  {
    m_particleManager = &input;
  }

  /**
   * @brief Get the ParticleRegionManager.
   * @return The current ParticleRegionManager.
   */
  ParticleManager const * getParticleManager() const
  {
    return m_particleManager;
  }

  /// The number of particles associated with the object
  array1d< localIndex > m_numParticles;
  /// The relationship between object indices and particle regions.
  BASETYPE m_toParticleRegion;
  /// The relationship between object indices and particle subregions.
  BASETYPE m_toParticleSubRegion;
  /// The relationship between object indices and particle indices.
  BASETYPE m_toParticleIndex;

  /// The current ParticleRegionManager
  ParticleManager const * m_particleManager{};
};

/// @brief A ToParticleRelation where each object is related to the same number of particles.
typedef ToParticleRelation< array2d< localIndex > > FixedToManyParticleRelation;

/// @brief A ToParticleRelation where each object is related to an arbitrary number of particles.
typedef ToParticleRelation< ArrayOfArrays< localIndex > > OrderedVariableToManyParticleRelation;

/**
 * @brief Remove a particle relation from an object in the relation.
 * @param relation The relationship mapping to remove a single particle relation from
 *                 a single object from.
 * @param firstIndex The object index to remove a particle relation from.
 * @param er The particle region to remove.
 * @param esr The particle subregion to remove.
 * @param ei The particle index to remove.
 */
void erase( OrderedVariableToManyParticleRelation & relation,
            localIndex const firstIndex,
            localIndex const er,
            localIndex const esr,
            localIndex const ei );

/**
 * @brief Insert a particle relation for an object in the relation. Checks for
 *        existing membership.
 * @param relation The relationship mapping to insert a single particle relation for
 *                 a single object into.
 * @param firstIndex The object index to insert a particle relation from.
 * @param er The particle region to insert.
 * @param esr The particle subregion to insert.
 * @param ei The particle index to insert.
 */
void insert( OrderedVariableToManyParticleRelation & relation,
             localIndex const firstIndex,
             localIndex const er,
             localIndex const esr,
             localIndex const ei );

/**
 * @brief Insert a particle relation for an object in the relation. This is slightly
 *        faster than "insert" because it does NOT check for existing membership.
 * @param relation The relationship mapping to insert a single particle relation for
 *                 a single object into.
 * @param firstIndex The object index to insert a particle relation from.
 * @param er The particle region to insert.
 * @param esr The particle subregion to insert.
 * @param ei The particle index to insert.
 */
void fastInsert( OrderedVariableToManyParticleRelation & relation,
                 localIndex const firstIndex,
                 localIndex const er,
                 localIndex const esr,
                 localIndex const ei );

/**
 * @brief Insert a particle relation for several objects in the relation.
 *        Does NOT check for existing membership.
 * @param relation The relationship mapping to insert a single particle relation for
 *                 a single object into.
 * @param firstIndex The object index to insert a particle relation from.
 * @param erArray The array of particle regions to insert.
 * @param esrArray The array of particle subregions to insert.
 * @param eiArray The array of particle indices to insert.
 */
void insertMany( OrderedVariableToManyParticleRelation & relation,
                 localIndex const firstIndex,
                 std::vector< localIndex > const & erArray,
                 std::vector< localIndex > const & esrArray,
                 std::vector< localIndex > const & eiArray );

/**
 * @brief Reserve a set number of entities for a particle to relate to
 * @param relation The relationship mapping to insert a single particle relation for
 *                 a single object into.
 * @param numToReserve Number of entities to reserve
 */
void reserveNeighbors( OrderedVariableToManyParticleRelation & relation,
                       int const numToReserve );


} /* namespace geos */

#endif /* GEOS_MESH_TOPARTICLERELATION_HPP_ */
