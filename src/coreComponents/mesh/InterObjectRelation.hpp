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
 * @file InterObjectRelation.hpp
 */

#ifndef GEOS_MESH_INTEROBJECTRELATION_HPP_
#define GEOS_MESH_INTEROBJECTRELATION_HPP_


#include <map>
#include "mesh/ObjectManagerBase.hpp"

namespace geos
{

/**
 * @tparam BASETYPE The base class to provide the implementation
 *         of the relationship mapping.
 */
template< typename BASETYPE >
class InterObjectRelation : public BASETYPE
{
public:

  /// The type of the base class
  using base_type = BASETYPE;

  /**
   * @brief A forwarding constructor.
   * @tparam ARGS The types of the arguments to forward to the BASETYPE constructor.
   * @param args A parameter pack of arguments to forward to the BASETYPE constructor.
   */
  template< typename ... ARGS >
  InterObjectRelation( ARGS && ... args ):
    BASETYPE( std::forward< ARGS >( args )... )
  {}

  /**
   * @brief copy assignement operator
   *
   * @param rhs the relation to copy
   * @return InterObjectRelation< BASETYPE >&
   */
  InterObjectRelation< BASETYPE > & operator=( InterObjectRelation< BASETYPE > const & rhs )
  {
    base_type::operator=( rhs.base() );
    m_relatedObject = rhs.m_relatedObject;
    return *this;
  }

  /**
   * @brief Get a reference to this object cast to BASETYPE const.
   * @return A reference to this object cast to BASETYPE const.
   */
  const base_type & base() const { return static_cast< const BASETYPE & >(*this); }

  /**
   * @brief Get a reference to this object cast to BASETYPE.
   * @return A reference to this object cast to BASETYPE.
   */
  base_type & base() { return dynamic_cast< BASETYPE & >(*this); }

  /**
   * @brief Set the related object.
   * @param relatedObject The related object to use for mapping.
   */
  void setRelatedObject( ObjectManagerBase const & relatedObject )
  { m_relatedObject = &relatedObject; }

  /**
   * @brief Get the related object.
   * @return The related object.
   */
  const ObjectManagerBase * relatedObject() const
  { return m_relatedObject; }

  /**
   * @brief Get the LocalToGlobal mapping from the related object.
   * @return The LocalToGlobal mapping from the related object.
   */
  arrayView1d< globalIndex const > relatedObjectLocalToGlobal() const
  { return this->m_relatedObject->localToGlobalMap(); }

  /**
   * @brief Get the GlobalToLocal mapping from the related object.
   * @return The GlobalToLocal mapping from the related object.
   */
  unordered_map< globalIndex, localIndex > const & relatedObjectGlobalToLocal() const
  { return this->m_relatedObject->globalToLocalMap(); }

private:
  ObjectManagerBase const * m_relatedObject = nullptr;
};

/**
 * @brief A relationship from single objects to many other objects, where
 *        each object is related to the same number of objects.
 **/
typedef InterObjectRelation< array2d< localIndex > > FixedOneToManyRelation;
}

#endif /* GEOS_MESH_INTEROBJECTRELATION_HPP_ */
