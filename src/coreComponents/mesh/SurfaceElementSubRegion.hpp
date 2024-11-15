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
 * @file SurfaceElementSubRegion.hpp
 */

#ifndef GEOS_MESH_SURFACEELEMENTSUBREGION_HPP_
#define GEOS_MESH_SURFACEELEMENTSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"
#include "EdgeManager.hpp"
#include "CellElementSubRegion.hpp"

namespace geos
{

/**
 * @class SurfaceElementSubRegion
 *
 * The SurfaceElementSubRegion class contains the functionality to support the concept of a
 * surface element that can be either and surface element or a face element.
 */
class SurfaceElementSubRegion : public ElementSubRegionBase
{
public:

  /// Surface element to nodes map type
  using NodeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /// Surface element to edges map type
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  static string catalogName()
  { return "SurfaceElementSubRegion"; }

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  virtual string getCatalogName() const override
  {
    return catalogName();
  }

  ///@}

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name the group name
   * @param parent the parent group
   */
  SurfaceElementSubRegion( string const & name,
                           dataRepository::Group * const parent );


  /// @brief Destructor
  virtual ~SurfaceElementSubRegion() override;

  /**
   * @brief Compute the center of each element in the subregion.
   * @param[in] X an arrayView of (const) node positions
   */
  void calculateElementCenters( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
  {
    ElementSubRegionBase::calculateElementCenters( m_toNodesRelation, X );
  }

  ///@}

  /**
   * @name Relation getters
   * @brief Getter functions for the various inter-object relations
   */
  ///@{

  /**
   * @brief Get the face element to nodes map.
   * @return the face element to node map
   */
  NodeMapType const & nodeList() const
  {
    return m_toNodesRelation;
  }

  /**
   * @copydoc nodeList() const
   */
  NodeMapType & nodeList()
  {
    return m_toNodesRelation;
  }

  /**
   * @brief Get the local index of the a-th node of the k-th element.
   * @param[in] k the index of the element
   * @param[in] a the index of the node in the element
   * @return a reference to the local index of the node
   */
  localIndex & nodeList( localIndex const k, localIndex a ) { return m_toNodesRelation( k, a ); }

  /**
   * @copydoc nodeList( localIndex const k, localIndex a )
   */
  localIndex const & nodeList( localIndex const k, localIndex a ) const { return m_toNodesRelation( k, a ); }

  /**
   * @brief Get the surface element to edges map.
   * @return The face element to edge map
   */
  EdgeMapType const & edgeList() const
  {
    return m_toEdgesRelation;
  }

  /**
   * @copydoc edgeList() const
   */
  EdgeMapType & edgeList()
  {
    return m_toEdgesRelation;
  }

  using ElementSubRegionBase::numNodesPerElement;

  localIndex numNodesPerElement( localIndex const k ) const final
  { return m_toNodesRelation[k].size(); }

  ///@}


  /**
   * @name Properties Getters
   * @brief Getters to surface element properties.
   */
  ///@{

  /**
   * @brief Get face element aperture.
   * @return the aperture of the face elements
   */
  arrayView1d< real64 > getElementAperture() { return m_elementAperture; }

  /**
   * @copydoc getElementAperture()
   */
  arrayView1d< real64 const > getElementAperture() const { return m_elementAperture; }

  /**
   * @brief Get face element surface area.
   * @return the surface area of the face element
   */
  arrayView1d< real64 > getElementArea() { return m_elementArea; }

  /**
   * @copydoc getElementArea()
   */
  arrayView1d< real64 const > getElementArea() const { return m_elementArea; }

  /**
   * @brief Const accessor to the normal vectors.
   * @return a const view to the array of normal vectors.
   */
  arrayView2d< real64 const > getNormalVector() const { return m_normalVector; }

  /**
   * @brief Non const accessor to the normal vectors.
   * @return a non const view to the array of normal vectors.
   */
  arrayView2d< real64 > getNormalVector() { return m_normalVector; }

  /**
   * @brief Get normal vector of a specific surface element.
   * @param k index of the surface element
   * @return the normal vector of a specific surface element
   */
  arraySlice1d< real64 const > getNormalVector( localIndex k ) const { return m_normalVector[k]; }

  /**
   * @brief Get an array of the first tangent vector of the surface elements.
   * @return an array of the first tangent vector of the surface elements
   */
  arrayView2d< real64 const > getTangentVector1() const { return m_tangentVector1; }

  /**
   * @brief Get an array of the first tangent vector of the surface elements.
   * @return a non const view to the array of the first tangent vector of the surface elements
   */
  arrayView2d< real64 > getTangentVector1() { return m_tangentVector1; }

  /**
   * @brief Get the first tangent vector of a specific surface element.
   * @param k index of the surface element
   * @return the first tangent vector of a specific surface element
   */
  arraySlice1d< real64 const > getTangentVector1( localIndex const k ) const { return m_tangentVector1[k]; }

  /**
   * @brief Get an array of the second tangent vector of the surface elements.
   * @return aconst view to the array of the second tangent vector of the surface elements
   */
  arrayView2d< real64 const > getTangentVector2() const { return m_tangentVector2.toViewConst(); }

  /**
   * @brief Get an array of the first tangent vector of the surface elements.
   * @return a non const view to the array of the second tangent vector of the surface elements
   */
  arrayView2d< real64 > getTangentVector2() { return m_tangentVector2; }

  /**
   * @brief Get the second tangent vector of a specific surface element.
   * @param k index of the surface element
   * @return the second tangent vector of a specific surface element
   */
  arraySlice1d< real64 const > getTangentVector2( localIndex const k ) const { return m_tangentVector2[k];}


  ///@}

  /**
   * @brief Struct containing the keys to all surface element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ElementSubRegionBase::viewKeyStruct
  {
    /// @return Face element to cell regions map string.
    static constexpr char const * surfaceElementsToCellRegionsString() { return "fractureElementsToCellRegions"; }

    /// @return Face element to cell subregions map string.
    static constexpr char const * surfaceElementsToCellSubRegionsString() { return "fractureElementsToCellSubRegions"; }

    /// @return Face element to cell indices map string.
    static constexpr char const * surfaceElementsToCellIndexString() { return "fractureElementsToCellIndices"; }

    /// @return surface element to parent plane string.
    constexpr static char const * surfaceElementToParentPlaneString() { return "surfaceElementToParentPlane"; }
  };

protected:

  /// Unmapped surface elements to nodes map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;

  /// list of nodes
  NodeMapType m_toNodesRelation;

  /// list of edges
  EdgeMapType m_toEdgesRelation;

  /// Member level field for the element center
  array1d< real64 > m_elementAperture;

  /// Member level field for the element center
  array1d< real64 > m_elementArea;

  /// Normal vector to the surface element
  array2d< real64 > m_normalVector;

  /// Unit vector indicating the first tangential direction
  array2d< real64 > m_tangentVector1;

  /// Unit vector indicating the second tangential direction
  array2d< real64 > m_tangentVector2;

};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_MESH_SURFACEELEMENTSUBREGION_HPP_ */
