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

#ifndef GEOS_CELLBLOCKABC_HPP
#define GEOS_CELLBLOCKABC_HPP

#include "dataRepository/Group.hpp"
#include "mesh/ElementType.hpp"
#include "common/DataTypes.hpp"

#include <vector>

namespace geos
{

/**
 * Abstract base class defining the information provided by any cell block implementation.
 * Mainly element type related information (number of faces per element...),
 * some element to nodes, faces and edges mappings, a local to global mapping and
 * some kind of accessors to external properties.
 *
 * The cells of CellBlockABC are all of the same kind.
 *
 * It's noteworthy that the CellBlockABC is immutable oriented.
 * The derived implementations need to have the modification/creation capabilities.
 */
class CellBlockABC : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  CellBlockABC( string const & name,
                Group * const parent )
    :
    Group( name, parent )
  { }

  /**
   * @brief Get the type of element in this subregion.
   * @return the type of element in this subregion
   *
   * See class FiniteElementBase for possible element type.
   */
  virtual ElementType getElementType() const = 0;

  /**
   * @brief Get the number of nodes per element.
   * @return number of nodes per element
   */
  virtual localIndex numNodesPerElement() const = 0;

  /**
   * @brief Get the number of edges per element.
   * @return number of edges per element
   */
  virtual localIndex numEdgesPerElement() const = 0;

  /**
   * @brief Get the number of faces per element.
   * @return number of faces per element
   */
  virtual localIndex numFacesPerElement() const = 0;

  /**
   * @brief Get the number of elements.
   * @return number of elements in the cell block
   */
  virtual localIndex numElements() const = 0;

  /**
   * @brief Get the element-to-nodes map.
   * @return The mapping relationship as a 2d-array.
   */
  virtual array2d< localIndex, cells::NODE_MAP_PERMUTATION > getElemToNodes() const = 0;

  /**
   * @brief Get the element-to-edges map.
   * @return The mapping relationship as a 2d-array.
   */
  virtual array2d< localIndex > getElemToEdges() const = 0;

  /**
   * @brief Get the element-to-faces map.
   * @return The mapping relationship as a 2d-array.
   */
  virtual array2d< localIndex > getElemToFaces() const = 0;

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
  void forExternalProperties( LAMBDA && lambda ) const
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
  virtual std::list< dataRepository::WrapperBase const * > getExternalProperties() const = 0;
};

}

#endif //GEOS_CELLBLOCKABC_HPP
