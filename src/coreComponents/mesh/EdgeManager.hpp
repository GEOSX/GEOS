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
 * @file EdgeManager.hpp
 */

#ifndef GEOS_MESH_EDGEMANAGER_HPP_
#define GEOS_MESH_EDGEMANAGER_HPP_

#include "mesh/ObjectManagerBase.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"
#include "LvArray/src/tensorOps.hpp"


namespace geos
{
class FaceManager;
class NodeManager;

/**
 * @class EdgeManager
 * @brief This class provides an interface to ObjectManagerBase in order to manage edge data.
 * @details An edge is bounded by two nodes. This class aims to manage all the edges of every
 * cell elements in a DomainPartition.
 */

class EdgeManager : public ObjectManagerBase
{
public:

  /// EdgeToNode map type
  using NodeMapType = FixedOneToManyRelation;

  /// EdgeToFace map type
  using FaceMapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @return the string representing the edge manager name in the catalog
   */
  static string catalogName()
  { return "EdgeManager"; }


  /**
   * @brief Getter used to access the edge manager catalog name.
   * @return the edge manager catalog name
   */
  virtual string getCatalogName() const override
  { return catalogName(); }

  ///@}

  /**
   * @brief Oversize the Face mapping by this amount for each edge (hardcoded)
   * @return extra space for each edge of the EdgetoFace map
   *
   * @note Value forwarding is due to refactoring.
   */
  static constexpr localIndex faceMapOverallocation()
  { return CellBlockManagerABC::faceMapExtraSpacePerEdge(); }

  /**
   * @name Constructors/destructors
   */
  ///@{

  /**
   * @brief main constructor for EdgeManager Objects
   * @param[in] name the name of the EdgeManager object in the repository
   * @param[in] parent the parent group of the EdgeManager object being constructed
   */
  EdgeManager( string const & name,
               Group * const parent );

  ///@}

  /**
   * @brief Resize the EdgeManager object and all it members.
   * @details the size of the EdgeMananager is the number of edges
   * in the domain.
   * @param[in] newSize the new number of edges.
   */
  virtual void resize( localIndex const newSize ) override;

  /**
   * @brief Set the node of the domain boundary object.
   * @param[in] faceManager The reference of the face manager.
   */
  void setDomainBoundaryObjects( FaceManager const & faceManager );

  /**
   * @brief Set external edges.
   * @details external edges are the edges on the faces which are external
   * @param[in] faceManager the face manager to obtain external face
   */
  void setIsExternal( FaceManager const & faceManager );

  /**
   * @brief Build sets from the node sets
   * @param[in] nodeManager The node manager that will provide the node sets.
   */
  void buildSets( NodeManager const & nodeManager );

  /**
   * @brief Copies the edges to (nodes|faces) mappings from @p cellBlockManager.
   * @param[in] cellBlockManager Provides the mappings.
   * @param[in] isBaseMeshLevel flag that indicates if we are operating on the base mesh level or on another mesh level
   */
  void setGeometricalRelations( CellBlockManagerABC const & cellBlockManager, bool isBaseMeshLevel );

  /**
   * @brief Link the current manager to other managers.
   * @param[in] nodeManager The node manager instance.
   * @param[in] faceManager The face manager instance.
   *
   * @note the @p EdgeManager do not hold any information related to the regions nor to the elements.
   * This is why the element region manager is not provided.
   */
  void setupRelatedObjectsInRelations( NodeManager const & nodeManager,
                                       FaceManager const & faceManager );

  /**
   * @brief Build faces-to-edges and nodes-to-edges relation maps.
   * @param[in] numNodes number of nodes
   * @param[in] faceToNodeMap manager face-to-nodes map
   * @param[in] faceToEdgeMap manager face-to-edges map
   */
  void buildEdges( localIndex const numNodes,
                   ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                   ArrayOfArrays< localIndex > & faceToEdgeMap );


  /**
   * @brief Build a vector containing all the global indices of each nodes of each edges
   * @param[in] nodeManager the nodeManager object.
   * @return an array of pairs of each edge's nodes globalIndices
   */
  virtual ArrayOfSets< globalIndex >
  extractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & nodeManager ) override;

  /**
   * @brief Compute the future size of a packed list.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in] packList the list of edge indices to be packed
   * @return The size of the packed list
   */
  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Pack an array of edge indices in a buffer.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in,out] buffer a buffer to pack the edge index data into
   * @param[in] packList the indices of edge to pack
   * @return The size of the packed array
   */
  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Unpack a buffer to an array of indices that has been packed.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in] buffer buffer containing the packed indices
   * @param[in,out] packList the array of local index to be unpacked
   * @param[in] overwriteUpMaps boolean to state if the Up maps should be overwritten
   * @param[in] overwriteDownMaps boolean to state if the Downs maps should be overwritten
   * @return The size of the unpacked array
   */
  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  /**
   * @brief Call fixUpDownMaps of the class ObjectManagerBase for nodes-to-edges and nodes-to-faces relation maps.
   * @param[in] clearIfUnmapped boolean to indicate if the unmaped indices should be removed or not
   */
  void fixUpDownMaps( bool const clearIfUnmapped );

  /**
   * @brief Compress all nodes-to-faces relation maps
   * so that the values of each array are contiguous with no extra capacity in between.
   */
  void compressRelationMaps();

  /**
   * @brief Clean up the edges-to-faces mapping with respect to a new list of edges.
   * @param[in] receivedEdges the new list of edges indices
   * @param[in] facesToEdges map to go from faces to edges
   */
  void depopulateUpMaps( std::set< localIndex > const & receivedEdges,
                         ArrayOfArraysView< localIndex const > const & facesToEdges );

  /**
   * @brief Check if edge \p edgeIndex contains node \p nodeIndex
   * @param[in] edgeIndex local index of the edge
   * @param[in] nodeIndex local index of the node
   * @return boolean : true if the node \p nodeIndex is part of the edge \p edgeIndex; false otherwise
   */
  bool hasNode( const localIndex edgeIndex, const localIndex nodeIndex ) const;

  /**
   * @brief Calculate the center of an edge given its index.
   * @tparam OUT_VECTOR type of output
   * @param edgeIndex index of the edge
   * @param X array view of nodes associated with the edge
   * @param edgeCenter output vector containing coordinate of the edge center
   */
  template< typename OUT_VECTOR >
  void calculateCenter( localIndex const edgeIndex,
                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                        OUT_VECTOR && edgeCenter ) const;

  /**
   * @brief Calculate the edge segment.
   * @tparam OUT_VECTOR type of output
   * @param edgeIndex index of the edge
   * @param X array view of the nodes associated with the NodeManager of current domain
   * @param edgeVector output vector containing edge segment
   */
  template< typename OUT_VECTOR >
  void calculateLength( localIndex const edgeIndex,
                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                        OUT_VECTOR && edgeVector ) const;

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{

  /**
   * @struct viewKeyStruct
   * @brief Container of keys needed to access the data of the class member
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @cond DO_NOT_DOCUMENT
    static constexpr char const * nodeListString() { return "nodeList"; }
    static constexpr char const * faceListString() { return "faceList"; }
    static constexpr char const * elementRegionListString() { return "elemRegionList"; }
    static constexpr char const * elementSubRegionListString() { return "elemSubRegionList"; }
    static constexpr char const * elementListString() { return "elemList"; }

    dataRepository::ViewKey nodesList             = { nodeListString() };
    dataRepository::ViewKey faceList              = { faceListString() };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString() };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString() };
    dataRepository::ViewKey elementList           = { elementListString() };
    /// @endcond
  }
  /// viewKeys
  viewKeys;

  ///}@

  /**
   * @name Getters for stored value.
   */
  ///@{

  /**
   * @brief Get a node list.
   * @return reference to the list of edge to node relation
   */
  NodeMapType & nodeList()       { return m_toNodesRelation; }

  /**
   * @brief Provide an accessor to an immutable node list.
   * @return reference to the list of edge to node relation
   */
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

  /**
   * @brief Get a reference to a node list
   * @param[in] edgeIndex local index of the edge
   * @param[in] nodeIndex local index of the node
   * @return a reference of the index of the edge to node relation
   */
  localIndex & nodeList( localIndex const edgeIndex, localIndex const nodeIndex )
  { return m_toNodesRelation( edgeIndex, nodeIndex ); }

  /**
   * @brief Get a node list
   * @param[in] edgeIndex local index of the edge
   * @param[in] nodeIndex local index of the node
   * @return a reference of the index of the edge to node relation
   */
  localIndex nodeList( localIndex const edgeIndex, localIndex const nodeIndex ) const
  { return m_toNodesRelation( edgeIndex, nodeIndex ); }

  /**
   * @brief Get the edge-to-face relation
   * @return A FaceMapType
   */
  FaceMapType & faceList()       { return m_toFacesRelation; }

  /**
   * @brief Get an immutable accessor to the edge-to-face relation
   * @return An immutable FaceMapType
   */
  FaceMapType const & faceList() const { return m_toFacesRelation; }

  ///@}

private:

  /// Map for the edges-to-nodes relation
  NodeMapType m_toNodesRelation;

  /// Map for the edges-to-faces relation
  FaceMapType m_toFacesRelation;

  /// Unmaped edge indices (those that are not in the nodes-to-edges relation map)
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;

  /// Unmaped edge indices (those that are not in the faces-to-edges relation map)
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  /**
   * @brief function to pack the upward and downward pointing maps.
   * @tparam DO_PACKING template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                    packed and the function returns the size of the packing that would have occurred if set to TRUE.
   * @param buffer the buffer to pack data into
   * @param packList the indices of nodes that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DO_PACKING >
  localIndex packUpDownMapsImpl( buffer_unit_type * & buffer,
                                 arrayView1d< localIndex const > const & packList ) const;

};

template< typename OUT_VECTOR >
inline void EdgeManager::calculateCenter( localIndex const edgeIndex,
                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                          OUT_VECTOR && edgeCenter ) const
{
  LvArray::tensorOps::copy< 3 >( edgeCenter, X[m_toNodesRelation[edgeIndex][0]] );
  LvArray::tensorOps::add< 3 >( edgeCenter, X[m_toNodesRelation[edgeIndex][1]] );
  LvArray::tensorOps::scale< 3 >( edgeCenter, 0.5 );
}

template< typename OUT_VECTOR >
inline void EdgeManager::calculateLength( localIndex const edgeIndex,
                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                          OUT_VECTOR && edgeSegment ) const
{
  LvArray::tensorOps::copy< 3 >( edgeSegment, X[m_toNodesRelation[edgeIndex][1]] );
  LvArray::tensorOps::subtract< 3 >( edgeSegment, X[m_toNodesRelation[edgeIndex][0]] );
}

}
#endif /* GEOS_MESH_EDGEMANAGER_HPP_ */
