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
 * @file EmbeddedSurfaceNodeManager.hpp
 */

#ifndef GEOS_MESH_EMBEDDEDSURFACENODEMANAGER_HPP_
#define GEOS_MESH_EMBEDDEDSURFACENODEMANAGER_HPP_

#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "ToElementRelation.hpp"

namespace geos
{

class EdgeManager;
class ElementRegionManager;

/**
 * @class EmbeddedSurfaceNodeManager
 * @brief The EmbeddedSurfaceNodeManager class provides an interface to ObjectManagerBase in order to manage node data.
 *
 * The EmbeddedSurfaceNodeManagerT class manages the embedded surfaces node data using the
 * ObjectDataStructureBaseT as a data manager.
 * This means that each field is stored in an array where each array entry
 * corresponds to a node of an embedded surface.
 */
class EmbeddedSurfaceNodeManager : public ObjectManagerBase
{
public:
  /// nodeToEdge map type
  using EdgeMapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /// nodeToElement map type
  using ElemMapType = OrderedVariableToManyElementRelation;

  /**
   * @brief return default size of the value array in the node-to-edge mapping
   * @return default size of value array in the node-to-edge mapping
   *
   * @note Value forwarding is due to refactoring.
   */
  static constexpr localIndex edgeMapOverallocation()
  { return CellBlockManagerABC::edgeMapExtraSpacePerNode(); }

  /**
   * @brief return default size of the value array in the node-to-element mapping
   * @return default size of value array in the node-to-element mapping
   */
  static constexpr localIndex elemMapOverallocation()
  { return 8; }


/**
 * @name Constructors/destructor
 */
  ///@{

  /**
   * @brief Main constructor for EmbeddedSurfaceNodeManager Objects.
   * @param [in] name the name of this instantiation of EmbeddedSurfaceNodeManager
   * @param [in] parent the parent group of this instantiation of EmbeddedSurfaceNodeManager
   */
  EmbeddedSurfaceNodeManager( string const & name,
                              dataRepository::Group * const parent );

  ///@}

  /**
   * @brief Resize the EmbeddedSurfaceNodeManager, and all its member vectors that relate nodes to faces, to edges, and to elements.
   * @details the size of the EmbeddedSurfaceNodeManager is the number of nodes
   * @param[in] newsize the new size of the EmbeddedSurfaceNodeManager
   */
  virtual void resize( localIndex const newsize ) override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Return the name of the node manager in the object catalog.
   * @return string that contains the EmbeddedSurfaceNodeManager catalog name
   */
  static string catalogName()
  { return "EmbeddedSurfaceNodeManager"; }

  /**
   * @brief Provide a virtual access to catalogName().
   * @return string that contains the EmbeddedSurfaceNodeManager catalog name
   */
  string getCatalogName() const override final
  { return catalogName(); }

  ///@}

  /**
   * @brief Link the EdgeManager \p edgeManager to the EmbeddedSurfaceNodeManager, and performs the node-to-edge mapping.
   * @param [in] edgeManager the edgeManager to assign this EmbeddedSurfaceNodeManager
   */
  void setEdgeMaps( EdgeManager const & edgeManager );


  /**
   * @brief Assign the ElementRegionManager \p elementRegionManager to the EmbeddedSurfaceNodeManager, and performs the node-to-element
   * mapping
   * @param [in] elementRegionManager the ElementRegionManager to assign this EmbeddedSurfaceNodeManager
   */
  void setElementMaps( ElementRegionManager const & elementRegionManager );

  /**
   * @brief Compress all EmbeddedSurfaceNodeManager member arrays so that the values of each array are contiguous with no extra capacity
   * inbetween.
   * @note The method used here on each arrays (compress) does not free any memory.
   */
  void compressRelationMaps();

  /**
   * @brief appends a node to the embeSurfaceNodeManager.
   * @param pointCoord node location
   * @param pointGhostRank ghost rank of the node
   */
  void appendNode( arraySlice1d< real64 const > const & pointCoord,
                   integer const & pointGhostRank );

  /**
   * @name Packing methods
   */
  ///@{

  /**
   * @brief Computes the pack size of the global maps elements in the @ packList.
   * @param packList The element we want packed.
   * @return The packed size.
   */
  localIndex packNewNodesGlobalMapsSize( arrayView1d< localIndex const > const & packList ) const;

  /**
   * @brief Packs the global maps for the nodes in the @ packList.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @return The packed size.
   */
  localIndex packNewNodesGlobalMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const;

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;


  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Unpacks the global maps from @p buffer and cleans them up from duplicated nodes.
   * @param buffer The buffer containing the packed data.
   * @param packList The (un)packed element.
   * @return The unpacked size.
   */
  localIndex unpackNewNodesGlobalMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList );

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  ///@}

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{

  /**
   *  @brief contains the added view access keys to be bound with class data member.
   *  @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @return String to access the reference position
    static constexpr char const * referencePositionString() { return "referencePosition"; }

    /// @return String to access the edge map
    static constexpr char const * edgeListString() { return "edgeList"; }

    /// @return String to access the element region map
    static constexpr char const * elementRegionListString() { return "elemRegionList"; }

    /// @return String to access the element subregion map
    static constexpr char const * elementSubRegionListString() { return "elemSubRegionList"; }

    /// @return String to access the element map
    static constexpr char const * elementListString() { return "elemList"; }

    /// @return String to access the parent edge globalIndex
    static constexpr char const * parentEdgeGlobalIndexString() { return "parentEdgeGlobalIndex"; }

    /// Accessor to reference position
    dataRepository::ViewKey referencePosition       = { referencePositionString() };

    /// Accessor to edge map
    dataRepository::ViewKey edgeList                = { edgeListString() };

    /// Accessor to element region map
    dataRepository::ViewKey elementRegionList       = { elementRegionListString() };

    /// Accessor to element subregion map
    dataRepository::ViewKey elementSubRegionList    = { elementSubRegionListString() };

    /// Accessor to element map
    dataRepository::ViewKey elementList             = { elementListString() };

    /// Accessor to element map
    dataRepository::ViewKey parentEdgeGlobalIndex   = { parentEdgeGlobalIndexString() };
  }
  /// viewKeys
  viewKeys;

  ///@}

  /**
   * @name Accessors for EmbeddedSurfaceNodeManager fixed data
   */
  ///@{

  /**
   * @brief Provide an immutable accessor to the nodes-to-edges relation.
   * @return const reference to  nodes-to-edges relation
   */
  EdgeMapType const & edgeList() const { return m_toEdgesRelation; }

  /**
   * @brief Get a mutable accessor to the node-to-edges relation.
   * @return reference to nodes-to-edges relation
   */
  EdgeMapType & edgeList() { return m_toEdgesRelation; }

  /**
   * @brief Get the nodes-to-elements relation.
   * @return reference to nodes-to-elements relation
   */
  ElemMapType & toElementRelation() {return m_toElements;}

  /**
   * @brief Provide a const accessor to the nodes-to-elements relation.
   * @details The returned ElemMapType gives access, for one node
   * to the element index, the element sub region, and the element region
   * in relation with a node
   * @return const reference to nodes-to-elements relation
   */
  ElemMapType const & toElementRelation() const {return m_toElements;}

  /**
   * @brief Get the mutable nodes-to-elements-regions relation.
   * @return reference to nodes-to-elements-regions relation
   */
  ArrayOfArrays< localIndex > & elementRegionList() { return m_toElements.m_toElementRegion; }

  /**
   * @brief Provide an immutable arrayView to the nodes-to-elements-regions relation.
   * @return const reference to nodes-to-elements-regions relation
   */
  ArrayOfArraysView< localIndex const > elementRegionList() const { return m_toElements.m_toElementRegion.toViewConst(); }

  /**
   * @brief Get the mutable nodes-to-elements-subregions relation.
   * @return reference to nodes-to-elements-subregions relation
   */
  ArrayOfArrays< localIndex > & elementSubRegionList() { return m_toElements.m_toElementSubRegion; }

  /**
   * @brief Provide an immutable arrayView to the nodes-to-elements-subregions relation.
   * @return const reference to nodes-to-elements-subregions relation
   */
  ArrayOfArraysView< localIndex const > elementSubRegionList() const { return m_toElements.m_toElementSubRegion.toViewConst(); }

  /**
   * @brief Get the mutable nodes-to-elements indices.
   * @return reference to nodes-to-elements indices
   */
  ArrayOfArrays< localIndex > & elementList() { return m_toElements.m_toElementIndex; }

  /**
   * @brief Provide an immutable arrayView to the nodes-to-elements indices.
   * @return const reference to nodes-to-elements indices
   */

  ArrayOfArraysView< localIndex const > elementList() const
  { return m_toElements.m_toElementIndex.toViewConst(); }

  /**
   * @brief Get the mutable reference position array. This table will contain all the node coordinates.
   * @return reference position array
   */
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & referencePosition() { return m_referencePosition; }

  /**
   * @brief Provide an immutable arrayView of the reference position. This table will contain all the node coordinates.
   * @return an immutable arrayView of the reference position.
   */
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > referencePosition() const
  { return m_referencePosition; }

  /**
   * @brief Provide an immutable arrayView of the parent Edge global index position.
   * @return an immutable arrayView of the parent edge global index.
   */
  array1d< globalIndex > & getParentEdgeGlobalIndex()
  {
    return m_parentEdgeGlobalIndex;
  }

  /**
   * @brief Provide an immutable arrayView of the parent Edge global index position.
   * @return an immutable arrayView of the parent edge global index.
   */
  arrayView1d< globalIndex const > getParentEdgeGlobalIndex() const
  {
    return m_parentEdgeGlobalIndex;
  }

  ///@}

private:

  /**
   * @brief Pack the upward and downward pointing maps into a buffer.
   * @tparam DO_PACKING template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                    packed and the function returns the size of the packing that would have occured if set to TRUE.
   * @param buffer the buffer to pack data into
   * @param packList the indices of nodes that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DO_PACKING >
  localIndex packUpDownMapsImpl( buffer_unit_type * & buffer,
                                 arrayView1d< localIndex const > const & packList ) const;

  /**
   * @brief Packing global maps.
   * @tparam DO_PACKING A template parameter to discriminate between actually packing or only computing the packing size.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @return The packed size.
   */
  template< bool DO_PACKING >
  localIndex packNewNodesGlobalMapsImpl( buffer_unit_type * & buffer,
                                         arrayView1d< localIndex const > const & packList ) const;

  /**
   * @brief Checks if a node already exists on this rank and returns its local index.
   * @param[in] nodeCoord the coordinates of the node.
   * @param[out] nodeIndexOnThisRank the local index of the node (if it exists).
   */
  void nodeExistsOnThisRank( real64 const (&nodeCoord)[3],
                             localIndex & nodeIndexOnThisRank ) const;

  /// reference position of the nodes
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_referencePosition;

  /// nodes-to-edges relation
  EdgeMapType m_toEdgesRelation;

  /// nodes-to-element relation
  ElemMapType m_toElements;

  /// parent edge global index
  array1d< globalIndex > m_parentEdgeGlobalIndex;

};
}

#endif // MESH_EMBEDDEDSURFACENODEMANAGER_HPP_
