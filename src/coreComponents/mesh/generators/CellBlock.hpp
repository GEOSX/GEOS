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

#ifndef GEOS_MESH_CELLBLOCK_HPP_
#define GEOS_MESH_CELLBLOCK_HPP_

#include "dataRepository/Group.hpp"
#include "mesh/generators/CellBlockABC.hpp"
#include "mesh/ElementType.hpp"

namespace geos
{

/**
 * This implementation of CellBlockABC mainly use the cell patterns/shapes
 * to build all the element to nodes, faces and edges mappings.
 */
class CellBlock : public CellBlockABC
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
  CellBlock( string const & name, Group * const parent );

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
   * @brief Defines the underlying element type (hex, tet...)
   * @param[in] elementType the element type
   *
   * @note Allocates the values of the element to nodes, edges, faces accordingly.
   */
  void setElementType( ElementType elementType );

  ElementType getElementType() const override
  { return m_elementType; }

  localIndex numNodesPerElement() const override
  { return m_numNodesPerElement; }

  localIndex numEdgesPerElement() const override
  { return m_numEdgesPerElement; }

  localIndex numFacesPerElement() const override
  { return m_numFacesPerElement; }

  localIndex numElements() const override
  { return size(); }

  /**
   * @brief Get the maximum number of nodes comprising a face.
   * @return maximum number of nodes comprising a face.
   */
  localIndex maxNodesPerFace() const
  { return m_maxNodesPerFace; }

  /**
   * @brief Puts the nodes of face @p faceNum of element @p cellIndex inside vector @p nodesInFaces.
   * @param[in] cellIndex The element index.
   * @param[in] faceNum the face number within an element
   * @param[out] nodesInFace space for result, must have enough space to fit all nodes
   * @return the number of nodes written into @p nodesInFace
   *
   * @p nodesInFaces is sorted from lower to larger node indices values.
   * @p nodesInFaces is exactly the size of the number of nodes.
   */
  localIndex getFaceNodes( localIndex const cellIndex,
                           localIndex const faceNum,
                           Span< localIndex > nodesInFace ) const;

  /**
   * @brief Get the element to nodes mapping, non-const version.
   * @return The mapping relationship as an array.
   *
   * @deprecated This accessor is meant to be used like a setter even though it's a bit like having public attribute...
   * Use a real setter instead.
   */
  arrayView2d< localIndex, cells::NODE_MAP_USD > getElemToNode()
  { return m_elementsToNodes; }

  /**
   * @brief Get the element to nodes mapping, const version.
   * @return The mapping relationship as an array view
   */
  arrayView2d< localIndex const, cells::NODE_MAP_USD > getElemToNode() const
  { return m_elementsToNodes; }

  array2d< localIndex, cells::NODE_MAP_PERMUTATION > getElemToNodes() const override
  { return m_elementsToNodes; }

  array2d< localIndex > getElemToEdges() const override
  { return m_elementsToEdges; }

  array2d< localIndex > getElemToFaces() const override
  { return m_elementsToFaces; }

  /**
   * @brief Get the element-to-faces map.
   * @return A const view of the mapping.
   */
  arrayView2d< localIndex const > getElemToFacesConstView() const
  { return m_elementsToFaces.toViewConst(); }

  /**
   * @brief Sets an entry in the element to faces mapping.
   * @param[in] cellIndex Index of the element
   * @param[in] faceNum Local index of the face of the element @p iElement (typically from 0 to 5 for an hexahedron).
   * @param[in] faceIndex The face index.
   */
  void setElementToFaces( localIndex const cellIndex,
                          localIndex const faceNum,
                          localIndex const faceIndex )
  {
    m_elementsToFaces( cellIndex, faceNum ) = faceIndex;
  }

  /**
   * @brief Sets an entry in the element to edges mapping.
   * @param[in] cellIndex Index of the element.
   * @param[in] edgeNum Local index of the edge of the element @p iElement (typically from 0 to 11 for an hexahedron).
   * @param[in] edgeIndex Index of the edge.
   *
   * In the element to edges mapping, element @p iElement has a given number of edges (typically 12 for a hexahedron).
   * Then edge @p edgeNum of this local indexing (typically 0 to 11) is meant to have global indexing of @p edgeIndex.
   */
  void setElementToEdges( localIndex const cellIndex,
                          localIndex const edgeNum,
                          localIndex const edgeIndex )
  {
    m_elementsToEdges( cellIndex, edgeNum ) = edgeIndex;
  }

  /**
   * @brief Checks if edge @p edgeIndex of element @p cellIndex has been defined.
   * @param[in] cellIndex Index of the element
   * @param[in] edgeNum Index of the edge of the element @p cellIndex (typically from 0 to 11 for an hexahedron).
   * @param[in] edgeIndex The edge index.
   * @return True if the entry is already there in the mapping. False otherwise.
   */
  bool hasElementToEdges( localIndex const cellIndex,
                          localIndex const edgeNum,
                          localIndex const edgeIndex ) const
  {
    return m_elementsToEdges( cellIndex, edgeNum ) == edgeIndex;
  }

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
   * @brief Get local to global map, const view version.
   * @return The mapping relationship as a array.
   */
  arrayView1d< globalIndex const > localToGlobalMapConstView() const
  { return m_localToGlobalMap.toViewConst(); }

  /**
   * @brief Resize the cell block to hold @p numElements
   * @param numElements The new number of elements.
   */
  void resize( dataRepository::indexType const numElements ) override final;

  /**
   * @brief Resize the cell block to hold @p numnodes
   * @param numNodes The new number of nodes.
   */
  void resizeNumNodes ( dataRepository::indexType const numNodes );
  ///@}

  /**
   * @name Properties
   */
  ///@{

  ///@}

private:

  /// Number of nodes per element in this block.
  localIndex m_numNodesPerElement = -1;

  /// Number of edges per element in this block.
  localIndex m_numEdgesPerElement = -1;

  /// Number of faces per element in this block.
  localIndex m_numFacesPerElement = -1;

  /// Maximum number of face nodes in this block.
  localIndex m_maxNodesPerFace = -1;

  /// Element-to-node relation
  array2d< localIndex, cells::NODE_MAP_PERMUTATION > m_elementsToNodes;

  /// Element-to-edges relation
  array2d< localIndex > m_elementsToEdges;

  /// Element-to-node relation
  array2d< localIndex > m_elementsToFaces;

  /// Contains the global index of each object.
  array1d< globalIndex > m_localToGlobalMap;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /// Type of element in this block.
  ElementType m_elementType;

  std::list< dataRepository::WrapperBase const * > getExternalProperties() const override
  {
    std::list< dataRepository::WrapperBase const * > result;
    for( string const & externalPropertyName : m_externalPropertyNames )
    {
      result.push_back( &this->getWrapperBase( externalPropertyName ) );
    }
    return result;
  }
};

}

#endif /* GEOS_MESH_CELLBLOCK_HPP_ */
