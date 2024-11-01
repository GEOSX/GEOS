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
 * @file MeshLevel.hpp
 */

#ifndef GEOS_MESH_MESHLEVEL_HPP_
#define GEOS_MESH_MESHLEVEL_HPP_

#include "NodeManager.hpp"
#include "ParticleManager.hpp"
#include "EmbeddedSurfaceNodeManager.hpp"
#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "FaceManager.hpp"

namespace geos
{
class ElementRegionManager;

/**
 * @class MeshLevel
 * @brief Class facilitating the representation of a multi-level discretization of a MeshBody.
 * @details This contains the main components that compose a discretized mesh in GEOSX (nodes, faces, elements).
 *          In current practice, the code utilizes a single ``MeshLevel`` until such time as we
 *          implement a proper multi-level mesh capability.
 */
class MeshLevel : public dataRepository::Group
{
public:

  /**
   * @brief Constructor for the MeshLevel object.
   * @param[in] name the name of the MeshLevel object in the repository
   * @param[in] parent the parent group of the MeshLevel object being constructed
   */
  MeshLevel( string const & name,
             Group * const parent );


  /**
   * @brief Constructor to create a shallow MeshLevel.
   * @param[in] name the name of the MeshLevel object in the repository
   * @param[in] parent the parent group of the MeshLevel object being constructed
   * @param[in] source The MeshLevel to make a shallow copy of.
   */
  MeshLevel( string const & name,
             Group * const parent,
             MeshLevel & source );


  /**
   * @brief Constructor for the MeshLevel object.
   * @param[in] name the name of the MeshLevel object in the repository
   * @param[in] parent the parent group of the MeshLevel object being constructed
   * @param[in] source The source MeshLevel to build the new one from
   * @param[in] order The order of the MeshLevel
   */
  MeshLevel( string const & name,
             Group * const parent,
             MeshLevel const & source,
             int const order );

  virtual ~MeshLevel() override;

  /**
   * @brief Generate the sets for the objects within a MeshLevel
   */
  void generateSets();


  /**
   * @brief Collects the nodes, edges, faces, and elements that are adjacent to a given list of nodes.
   * @param[in] seedNodeList the input nodes
   * @param[out] nodeAdjacencyList the nodes adjacent to the input nodes of seedNodeList
   * @param[out] edgeAdjacencyList the edges adjacent to the input nodes of seedNodeList
   * @param[out] faceAdjacencyList the faces adjacent to the input nodes of seedNodeList
   * @param[out] elementAdjacencyList the elements adjacent to the input nodes of seedNodeList
   * @param[in] depth the depth of the search for adjacent quantities (first-order neighbors, neighbors of neighbors, etc)
   * @details All the additional information (nodes, edges, faces) connected to
   * the edges, faces, elements that touch the @p seedNodeList is also collected.
   * For instance, all the nodes, edges and faces that touch an element that relies on a node of the @p seedNodeList, will be considered.
   * Even if these "second level" geometrical elements do touch @p seedNodeList themselves.
   * The idea being obviously that any shared geometrical object needs to be fully defined.
   */
  void generateAdjacencyLists( arrayView1d< localIndex const > const & seedNodeList,
                               localIndex_array & nodeAdjacencyList,
                               localIndex_array & edgeAdjacencyList,
                               localIndex_array & faceAdjacencyList,
                               ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > & elementAdjacencyList,
                               integer const depth );


  virtual void initializePostInitialConditionsPostSubGroups() override;


  /// @cond DO_NOT_DOCUMENT

  struct viewStructKeys
  {
    dataRepository::ViewKey meshLevel                = { "meshLevel" };
//    static constexpr char const *  baseDiscretizationString() { return "baseDiscretization"; }

  } viewKeys;

  struct groupStructKeys
  {
    dataRepository::GroupKey vertexManager  = { "vertexManager" };
    dataRepository::GroupKey cellManager    = { "cellManager" };

    static constexpr char const * nodeManagerString() { return "nodeManager"; }
    static constexpr char const * edgeManagerString() { return "edgeManager"; }
    static constexpr char const * faceManagerString() { return "faceManager"; }

    // This key is defined in problem manager:
    static constexpr char const * elemManagerString() { return "ElementRegions"; }
    static constexpr char const * particleManagerString() { return "ParticleRegions"; }

    static constexpr auto embSurfNodeManagerString = "embeddedSurfacesNodeManager";
    static constexpr auto embSurfEdgeManagerString = "embeddedSurfacesEdgeManager";

    dataRepository::GroupKey nodeManager = {nodeManagerString()};
    dataRepository::GroupKey particleManager = {particleManagerString()};
    dataRepository::GroupKey edgeManager = {edgeManagerString()};
    dataRepository::GroupKey faceManager = {faceManagerString()};
    dataRepository::GroupKey elemManager = {elemManagerString()};
    dataRepository::GroupKey embSurfNodeManager = {embSurfNodeManagerString};
    dataRepository::GroupKey embSurfEdgeManager = {embSurfEdgeManagerString};
  } groupKeys;

  /// @endcond

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the node manager.
   * @return a reference to the nodeManager object
   */
  NodeManager const & getNodeManager() const
  { return *m_nodeManager; }

  /**
   * @copydoc getNodeManager() const
   */
  NodeManager & getNodeManager()
  { return *m_nodeManager; }

  /**
   * @brief Get the particle manager.
   * @return a reference to the particleManager object
   */
  ParticleManager const & getParticleManager() const
  { return *m_particleManager; }

  /**
   * @copydoc getParticleManager() const
   */
  ParticleManager & getParticleManager()
  { return *m_particleManager; }

  /**
   * @brief Get the edge manager.
   * @return a reference to the edgeManager object
   */
  EdgeManager const & getEdgeManager() const
  { return *m_edgeManager; }

  /**
   * @copydoc getEdgeManager() const
   */
  EdgeManager & getEdgeManager()
  { return *m_edgeManager; }

  /**
   * @brief Get the face manager.
   * @return a reference to the faceManager object
   */
  FaceManager const & getFaceManager() const
  { return *m_faceManager; }

  /**
   * @copydoc getFaceManager() const
   */
  FaceManager & getFaceManager()
  { return *m_faceManager; }

  /**
   * @brief Get the element region manager.
   * @return a reference to the elementRegionManager object
   */
  ElementRegionManager const & getElemManager() const
  { return *m_elementManager; }

  /**
   * @copydoc getElemManager() const
   */
  ElementRegionManager & getElemManager()
  { return *m_elementManager; }

  /**
   * @brief Get the node Manager of the embedded surfaces grid.
   * @return a pointer to the EmbeddedSurfaceNodeManager
   */
  EmbeddedSurfaceNodeManager const & getEmbSurfNodeManager() const
  { return *m_embSurfNodeManager; }

  /**
   * @copydoc getEmbSurfNodeManager() const
   */
  EmbeddedSurfaceNodeManager & getEmbSurfNodeManager()
  { return *m_embSurfNodeManager; }

  /**
   * @brief Get the edge Manager related to the embedded surfaces grid.
   * @return a pointer to the edgeManager related to the embedded surfaces grid
   */
  EdgeManager const & getEmbSurfEdgeManager() const
  { return *m_embSurfEdgeManager; }

  /**
   * @copydoc getEmbSurfEdgeManager() const
   */
  EdgeManager & getEmbSurfEdgeManager()
  { return *m_embSurfEdgeManager; }

  /**
   * @brief Getter for the modification timestamp
   * @return the timestamp of the last modification
   */
  Timestamp getModificationTimestamp() const
  { return m_modificationTimestamp; }

  /**
   * @brief Increment the modification timestamp if the mesh has been modified
   */
  void modified()
  { m_modificationTimestamp++; }

  /**
   * @return value of m_isShallowCopy.
   */
  bool isShallowCopy() const
  { return m_isShallowCopy; }

  /**
   * @brief Determines if this->MeshLevel is a shallow copy of the input.
   * @param comparisonLevel The MeshLevel to compare with.
   * @return Whether or not the comparison is true.
   */
  bool isShallowCopyOf( MeshLevel const & comparisonLevel ) const;


  /**
   * @brief If this is a shallow clone of another MeshLevel, then return the source MeshLevel.
   *
   * @return MeshLevel const&
   */
  MeshLevel const & getShallowParent() const;

  /**
   * @brief If this is a shallow clone of another MeshLevel, then return the source MeshLevel.
   *
   * @return MeshLevel const&
   */
  MeshLevel & getShallowParent();

  ///@}

private:

  /// Manager for node data
  NodeManager * const m_nodeManager;
  /// Manager for particle data
  ParticleManager * const m_particleManager;
  /// Manager for edge data
  EdgeManager * const m_edgeManager;
  /// Manager for face data
  FaceManager * const m_faceManager;
  /// Manager for element data
  ElementRegionManager * const m_elementManager;

  ///Manager for embedded surfaces nodes
  EmbeddedSurfaceNodeManager * const m_embSurfNodeManager;
  /// Manager for embedded surfaces edge data
  EdgeManager * const m_embSurfEdgeManager;

  /// Timestamp of the last modification of the mesh level
  Timestamp m_modificationTimestamp;

  bool const m_isShallowCopy = false;

  MeshLevel * const m_shallowParent;

};

} /* namespace geos */

#endif /* GEOS_MESH_MESHLEVEL_HPP_ */
