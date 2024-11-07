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
 * @file EmbeddedSurfaceSubRegion.hpp
 */

#ifndef GEOS_MESH_EMBEDDEDSURFACESUBREGION_HPP_
#define GEOS_MESH_EMBEDDEDSURFACESUBREGION_HPP_

#include "SurfaceElementSubRegion.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"
#include "EdgeManager.hpp"
#include "EmbeddedSurfaceNodeManager.hpp"
#include "CellElementSubRegion.hpp"
#include "mesh/generators/EmbeddedSurfaceBlockABC.hpp"
//Do we really need this include Rectangle?
#include "simpleGeometricObjects/Rectangle.hpp"

namespace geos
{

/**
 * @brief Struct defining an embedded element which has at least on node which is a ghost on this rank
 * @struct surfaceWithGhostNodes
 */
struct surfaceWithGhostNodes
{
  /// local index of the surface element
  localIndex surfaceIndex;
  /// index of the parent edge of each node
  std::vector< globalIndex > parentEdgeIndex;
  ///number of nodes of the element
  localIndex numOfNodes;

  /**
   * @brief Constructor
   */
  surfaceWithGhostNodes():
    surfaceIndex(),
    parentEdgeIndex(),
    numOfNodes( 0 ){}
  /**
   * @brief insert a new node
   * @param edgeIndex global index of the parent edge
   */
  void insert ( globalIndex const & edgeIndex );
};

/**
 * @class EmbeddedSurfaceSubRegion
 *
 * The EmbeddedSurfaceSubRegion class contains the functionality to support the concept of an embedded
 * surface element. It consists of a 2D surface that cuts a 3D matrix cell.
 */
class EmbeddedSurfaceSubRegion : public SurfaceElementSubRegion
{
public:

  /// Embedded surface element to faces map type
  using FaceMapType = FixedOneToManyRelation;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  static string catalogName()
  { return "EmbeddedSurfaceSubRegion"; }

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
  EmbeddedSurfaceSubRegion( string const & name,
                            dataRepository::Group * const parent );

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  virtual void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & facemanager ) override final;

  /**
   * @brief Function to compute the geometric quantities of a specific embedded surface element.
   * @param intersectionPoints array containing the nodes defining the embedded surface elements
   * @param k index of the embedded surface element
   */
  void calculateElementGeometricQuantities( arrayView2d< real64 const > const intersectionPoints,
                                            localIndex k );
  /**
   * @brief Function to add a new embedded surface element.
   * @param cellIndex cell element index
   * @param regionIndex cell element region index
   * @param subRegionIndex cell element subregion index
   * @param nodeManager the nodemanager group
   * @param embSurfNodeManager the embSurfNodeManager group
   * @param edgeManager the edgemanager group
   * @param cellToEdges cellElement to edges map
   * @param fracture pointer to the bounded plane which is defining the embedded surface element
   * @return boolean defining whether the embedded element was added or not
   */
  bool addNewEmbeddedSurface( localIndex const cellIndex,
                              localIndex const regionIndex,
                              localIndex const subRegionIndex,
                              NodeManager const & nodeManager,
                              EmbeddedSurfaceNodeManager & embSurfNodeManager,
                              EdgeManager const & edgeManager,
                              FixedOneToManyRelation const & cellToEdges,
                              PlanarGeometricObject const * fracture );

  /**
   * @brief Function that maps each node of an embedded surface element to its parent edge on the 3D element
   * @param elemNodesLocations geometric coordinates of each embedded surface element node (no duplicates for shared nodes)
   * @param elemToNodes embedded surface node index to the array of its node indices map
   * @param elemTo3dElem embedded element to its parent 3D cell element map
   * @param cellToEdges 3D cells to its edges map
   * @param edgeToNodes edge to its node indices map
   * @param nodesCoord geometric coordinates of 3D cells
   * @return embedded node to 3D cell edge mapping
   */
  array1d< localIndex > getEdfmNodeParentEdgeIndex( ArrayOfArraysView< real64 > const & elemNodesLocations,
                                                    ArrayOfArraysView< localIndex > const & elemToNodes,
                                                    ToCellRelation< ArrayOfArrays< localIndex > > const & elemTo3dElem,
                                                    FixedOneToManyRelation const & cellToEdges,
                                                    arrayView2d< localIndex const > const & edgeToNodes,
                                                    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord );

  /**
   * @brief Function to add a all the embedded surface elements (loaded from VTK).
   * @param regionIndex cell element region index
   * @param subRegionIndex cell element subregion index
   * @param nodeManager the nodemanager group
   * @param embSurfNodeManager the embSurfNodeManager group
   * @param edgeManager the edgemanager group
   * @param cellToEdges cellElement to edges map
   * @param embeddedSurfaceBlock embeddedSurfaceBloc holds the information for all embedded element network
   * @return boolean defining whether all the embedded elements were added or not
   */

  bool addAllEmbeddedSurfaces( localIndex const regionIndex,
                               localIndex const subRegionIndex,
                               NodeManager const & nodeManager,
                               EmbeddedSurfaceNodeManager & embSurfNodeManager,
                               EdgeManager const & edgeManager,
                               FixedOneToManyRelation const & cellToEdges,
                               EmbeddedSurfaceBlockABC const & embeddedSurfaceBlock );


  /**
   * @brief inherit ghost rank from cell elements.
   * @param cellGhostRank cell element ghost ranks
   */
  void inheritGhostRank( array1d< array1d< arrayView1d< integer const > > > const & cellGhostRank );

  /**
   * @brief Given the coordinates of a node, it computes the Heaviside function iside a cut element with respect to the fracture element.
   * @param nodeCoord coordinate of the node
   * @param k embedded surface cell index
   * @return value of the Heaviside
   */
  real64 computeHeavisideFunction( ArraySlice< real64 const, 1, nodes::REFERENCE_POSITION_USD - 1 > const nodeCoord,
                                   localIndex const k ) const;

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  ///@}

  /**
   * @brief Struct containing the keys to all embedded surface element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : SurfaceElementSubRegion::viewKeyStruct
  {
    /// @return Connectivity index string
    static constexpr char const * connectivityIndexString() { return "connectivityIndex"; }

    /// @return surfaces with ghost nodes list string
    static constexpr char const * surfaceWithGhostNodesString() { return "surfaceWithGhostNodes"; }
  }
  /// viewKey struct for the EmbeddedSurfaceSubRegion class
  viewKeys;

  virtual void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;

  /**
   * @name Properties Getters
   * @brief Getters to embedded surface elements properties.
   */
  ///@{

  /**
   * @brief Get number of jump enrichments.
   * @return a reference to the number of jump enrichments
   */
  localIndex & numOfJumpEnrichments()       {return m_numOfJumpEnrichments;}

  /**
   * @brief Get number of jump enrichments.
   * @return  a constant reference to the number of jump enrichments
   */
  localIndex const & numOfJumpEnrichments() const {return m_numOfJumpEnrichments;}

  /**
   * @brief Get the name of the bounding plate that was used to generate fracture element k.
   * @param k the index of the embedded surface element
   * @return the name of the bounded plane, the element was generated from
   */
  string const & getFractureName( localIndex k ) const { return m_parentPlaneName[k]; }

  /**
   * @brief Get the connectivity index of the  embedded surface element.
   * @return the connectivity index
   */
  arrayView1d< real64 > getConnectivityIndex()   { return m_connectivityIndex.toView();}

  /**
   * @copydoc getConnectivityIndex()
   */
  arrayView1d< real64 const > getConnectivityIndex() const { return m_connectivityIndex;}

  /**
   * @brief accessor to the m_surfaceWithGhostNodes list
   * @return the list of surfaces with at least one ghost node.
   */
  std::vector< struct surfaceWithGhostNodes > surfaceWithGhostNodes() { return m_surfaceWithGhostNodes; }

  ///@}

private:

  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DO_PACKING >
  localIndex packUpDownMapsImpl( buffer_unit_type * & buffer,
                                 arrayView1d< localIndex const > const & packList ) const;

  /// The number of jump enrichments
  localIndex m_numOfJumpEnrichments;

  /// The CI of the cells
  array1d< real64 > m_connectivityIndex;

  // Indices of geometric objects the element belongs to
  array1d< string > m_parentPlaneName;

  /// Surfaces with ghost nodes
  std::vector< struct surfaceWithGhostNodes > m_surfaceWithGhostNodes;
};


} /* namespace geos */

#endif /* GEOS_MESH_EMBEDDEDSURFACESUBREGION_HPP_ */
