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


#ifndef GEOS_MESH_CELLELEMENTSUBREGION_HPP_
#define GEOS_MESH_CELLELEMENTSUBREGION_HPP_

#include "mesh/generators/CellBlockABC.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "ElementSubRegionBase.hpp"


namespace geos
{

class MeshLevel;

/**
 * @class CellElementSubRegion
 * Class specializing the element subregion for a cell element.
 * This is the class used in the physics solvers to represent a collection of mesh cell elements
 */
class CellElementSubRegion : public ElementSubRegionBase
{
public:

  /// Alias for the type of the element-to-node map
  using NodeMapType = InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;
  /// Alias for the type of the element-to-edge map
  using EdgeMapType = FixedOneToManyRelation;
  /// Alias for the type of the element-to-face map
  using FaceMapType = FixedOneToManyRelation;
  /// Type of map between cell blocks and embedded elements
  using EmbSurfMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName()
  { return "CellElementSubRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return catalogName(); }

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  CellElementSubRegion( string const & name, Group * const parent );

  ///@}

  /**
   * @name Helpers for CellElementSubRegion construction
   */
  ///@{

  /**
   * @brief Fill the CellElementSubRegion by copying those of the source CellBlock
   * @param cellBlock the CellBlock which properties (connectivity info) will be copied.
   */
  void copyFromCellBlock( CellBlockABC const & cellBlock );

  ///@}

  /**
   * @brief Add fractured element to list and relative entries to the map.
   * @param cellElemIndex cell element index
   * @param embSurfIndex embedded surface element index
   */
  void addFracturedElement( localIndex const cellElemIndex,
                            localIndex const embSurfIndex );

  virtual void resizePerElementValues( localIndex const numNodesPerElement,
                                       localIndex const numEdgesPerElement,
                                       localIndex const numFacesPerElement ) override;


  /**
   * @name Overriding packing / Unpacking functions
   */
  ///@{

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  /**
   * @brief Computes the pack size of the set of fractured elements and the element-to-embeddedSurfaces map of the elements in the @
   * packList.
   * @param packList The element we want packed.
   * @param  embeddedSurfacesLocalToGlobal LocaltoGlobal map of the embedded surfaces.
   * @return The packed size.
   */
  localIndex packFracturedElementsSize( arrayView1d< localIndex const > const & packList,
                                        arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const;
  /**
   * @brief Packs the set of fractured elements and the element-to-embeddedSurfaces map of the elements in the @ packList.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @param embeddedSurfacesLocalToGlobal LocaltoGlobal map of the embedded surfaces.
   * @return The packed size.
   */
  localIndex packFracturedElements( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList,
                                    arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const;

  /**
   * @brief Unpacks the set of fractured elemetn and the element-to-embeddedSurfaces map from @p buffer.
   * @param buffer The buffer containing the packed data.
   * @param packList The (un)packed element.
   * @param embeddedSurfacesGlobalToLocal GlobalToLocal map of the embedded surfaces.
   * @return The unpacked size.
   */
  localIndex unpackFracturedElements( buffer_unit_type const * & buffer,
                                      localIndex_array & packList,
                                      unordered_map< globalIndex, localIndex > const & embeddedSurfacesGlobalToLocal );

  virtual void fixUpDownMaps( bool const clearIfUnmapped ) final override;

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Helper function to apply a lambda function over all constructive groups
   * @tparam LAMBDA the type of the lambda function
   * @param lambda the lambda function
   */
  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {

    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  ///@}

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementSubRegionBase::viewKeyStruct
  {
    /// @return String key for the constitutive point volume fraction
    static constexpr char const * constitutivePointVolumeFractionString() { return "ConstitutivePointVolumeFraction"; }
    /// @return String key for the derivatives of the shape functions with respect to the reference configuration
    static constexpr char const * dNdXString() { return "dNdX"; }
    /// @return String key for the derivative of the jacobian.
    static constexpr char const * detJString() { return "detJ"; }
    /// @return String key for the constitutive grouping
    static constexpr char const * constitutiveGroupingString() { return "ConstitutiveGrouping"; }
    /// @return String key for the constitutive map
    static constexpr char const * constitutiveMapString() { return "ConstitutiveMap"; }
    /// @return String key to embSurfMap
    static constexpr char const * toEmbSurfString() { return "ToEmbeddedSurfaces"; }
    /// @return String key to fracturedCells
    static constexpr char const * fracturedCellsString() { return "fracturedCells"; }
    /// @return String key to bubbleCells
    static constexpr char const * bubbleCellsString() { return "bubbleCells"; }
    /// @return String key to toFaceElements
    static constexpr char const * toFaceElementsString() { return "toFaceElements"; }

    /// ViewKey for the constitutive grouping
    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString() };
    /// ViewKey for the constitutive map
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString() };
  }
  /// viewKey struct for the CellElementSubRegion class
  m_CellBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_CellBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockSubRegionViewKeys; }

  /**
   * @brief Get the local indices of the nodes in a face of the element.
   * @param[in] elementIndex The local index of the target element.
   * @param[in] localFaceIndex The local index of the target face in the element (this will be [0, numFacesInElement[)
   * @param[out] nodeIndices Memory to which node indices for the face will be written, must have sufficient size.
   * @return tne number of values written into @p nodeIndices
   * @deprecated This method will be removed soon.
   */
  localIndex getFaceNodes( localIndex const elementIndex,
                           localIndex const localFaceIndex,
                           Span< localIndex > const nodeIndices ) const;

  /**
   * @brief Get the element-to-node map.
   * @return a reference to the element-to-node map
   */
  NodeMapType & nodeList() { return m_toNodesRelation; }

  /**
   * @copydoc nodeList()
   */
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

  /**
   * @brief Get the local index of the a-th node of the k-th element.
   * @param[in] k the index of the element
   * @param[in] a the index of the node in the element
   * @return a reference to the local index of the node
   */
  localIndex & nodeList( localIndex k, localIndex a ) { return m_toNodesRelation( k, a ); }

  /**
   * @copydoc nodeList( localIndex const k, localIndex a )
   */
  localIndex const & nodeList( localIndex k, localIndex a ) const { return m_toNodesRelation( k, a ); }

  /**
   * @brief Get the element-to-edge map.
   * @return a reference to element-to-edge map
   */
  FixedOneToManyRelation & edgeList() { return m_toEdgesRelation; }

  /**
   * @copydoc edgeList()
   */
  FixedOneToManyRelation const & edgeList() const { return m_toEdgesRelation; }

  /**
   * @brief Get the element-to-face map.
   * @return a reference to the element to face map
   */
  FixedOneToManyRelation & faceList() { return m_toFacesRelation; }

  /**
   * @copydoc faceList()
   */
  FixedOneToManyRelation const & faceList() const { return m_toFacesRelation; }

  /**
   * @brief Get the local index of the a-th face of the k-th element.
   * @param[in] k the index of the element
   * @param[in] a the index of the face in the element
   * @return a const reference to the local index of the face
   */
  localIndex faceList( localIndex k, localIndex a ) const { return m_toFacesRelation( k, a ); }

  /**
   * @brief @return The array of shape function derivatives.
   */
  array4d< real64 > & dNdX()
  { return m_dNdX; }

  /**
   * @brief @return The array of shape function derivatives.
   */
  arrayView4d< real64 const > dNdX() const
  { return m_dNdX; }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  array2d< real64 > & detJ()
  { return m_detJ; }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  arrayView2d< real64 const > detJ() const
  { return m_detJ; }

  /**
   * @brief @return The sorted array of local fractured elements.
   */
  SortedArray< localIndex > & fracturedElementsList()
  { return m_fracturedCells; }

  /**
   * @brief @return The sorted array view of local fractured elements.
   */
  SortedArrayView< localIndex const > const fracturedElementsList() const
  { return m_fracturedCells.toViewConst(); }

  /**
   * @brief @return The map to the embedded surfaces
   */
  EmbSurfMapType & embeddedSurfacesList() { return m_toEmbeddedSurfaces; }

  /**
   * @brief @return The map to the embedded surfaces
   */
  EmbSurfMapType const & embeddedSurfacesList() const { return m_toEmbeddedSurfaces; }

  /**
   * @brief Set the array of bubble elements.
   * @param[in] bubbleCells The array of bubble elements.
   */
  void setBubbleElementsList( arrayView1d< localIndex const > const bubbleCells )
  { m_bubbleCells = bubbleCells; }

  /**
   * @brief @return The array view of bubble elements.
   */
  arrayView1d< localIndex const > const bubbleElementsList() const
  { return m_bubbleCells.toViewConst(); }

  /**
   * @brief Set the array of neighboring face elements.
   * @param[in] faceElemsList The array of neighboring face elements.
   */
  void setFaceElementsList( arrayView2d< localIndex const > const faceElemsList )
  { m_toFaceElements = faceElemsList; }

  /**
   * @brief @return  The array of neighboring face elements.
   */
  arrayView2d< localIndex const > const faceElementsList() const
  { return m_toFaceElements.toViewConst(); }

  /**
   * @brief Compute the center of each element in the subregion.
   * @param[in] X an arrayView of (const) node positions
   */
  void calculateElementCenters( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
  {
    ElementSubRegionBase::calculateElementCenters( m_toNodesRelation, X );
  }

  void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                            FaceManager const & faceManager ) override;

private:

  /// Map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// Array of constitutive point volume fraction
  array3d< real64 > m_constitutivePointVolumeFraction;

  /// Element-to-node relation
  NodeMapType m_toNodesRelation;

  /// Element-to-edge relation
  EdgeMapType m_toEdgesRelation;

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /**
   * @brief Compute the volume of the k-th element in the subregion.
   * @param[in] k the index of the element in the subregion
   * @param[in] X an arrayView of (const) node positions
   */
  void calculateCellVolumesKernel( localIndex const k,
                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const;

  void calculateElementCenterAndVolume( localIndex const k,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const;

  /// The array of shape function derivaties.
  array4d< real64 > m_dNdX;

  /// The array of jacobian determinantes.
  array2d< real64 > m_detJ;

  /// Map of unmapped global indices in the element-to-node map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInNodelist;

  /// Map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInEdgelist;

  /// Map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInFacelist;

  /// List of the fractured elements on this rank
  SortedArray< localIndex > m_fracturedCells;

  /// Map from local Cell Elements to Embedded Surfaces
  EmbSurfMapType m_toEmbeddedSurfaces;

  /// List of the elements adjacent to faceElement (useful for bubble elements)
  array1d< localIndex > m_bubbleCells;

  /// List of face IDs adjacent to faceElement (useful for bubble elements).
  /// A 2D array is needed to store both the face ID in the physical space
  /// and the corresponding face ID in the parent element space.
  array2d< localIndex > m_toFaceElements;

  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam DO_PACKING the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DO_PACKING >
  localIndex packUpDownMapsImpl( buffer_unit_type * & buffer,
                                 arrayView1d< localIndex const > const & packList ) const;

  /**
   * @brief Links the managers to their mappings.
   * @param[in] mesh Holds the node/edge/face managers.
   *
   * Defines links the element to nodes, edges and faces to their respective node/edge/face managers.
   */
  void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;

  template< bool DO_PACKING >
  localIndex packFracturedElementsImpl( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList,
                                        arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const;
};

} /* namespace geos */

#endif /* GEOS_MESH_CELLELEMENTSUBREGION_HPP_ */
