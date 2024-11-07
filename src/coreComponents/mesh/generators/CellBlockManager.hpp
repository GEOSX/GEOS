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
 * @file CellBlockManager.hpp
 */

#ifndef GEOS_MESH_CELLBLOCKMANAGER_H_
#define GEOS_MESH_CELLBLOCKMANAGER_H_

#include "mesh/generators/CellBlock.hpp"
#include "mesh/generators/FaceBlock.hpp"
#include "mesh/generators/InternalWellGenerator.hpp"
#include "mesh/generators/LineBlock.hpp"
#include "mesh/generators/LineBlockABC.hpp"
#include "mesh/generators/EmbeddedSurfaceBlock.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/PartitionDescriptor.hpp"

namespace geos
{

/**
 * @class CellBlockManager
 * @brief The CellBlockManager class provides an interface to ObjectManagerBase in order to manage CellBlock data.
 */
class CellBlockManager : public CellBlockManagerABC
{
public:

  /**
   * @brief Constructor for CellBlockManager object.
   * @param name name of this instantiation of CellBlockManager
   * @param parent pointer to the parent Group of this instantiation of CellBlockManager
   */
  CellBlockManager( string const & name, Group * const parent );

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Maximum number of nodes allowed (in memory) per each face.
   * @return The number as an integer.
   */
  static constexpr int maxNodesPerFace()
  { return 64; }

  array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const override;

  /**
   * @brief Returns a view to the vector holding the nodes coordinates
   * @return The reference
   *
   * @note This is meant to be used as a values setter.
   */
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > getNodePositions();

  ArrayOfArrays< localIndex > getNodeToEdges() const override;

  ArrayOfArrays< localIndex > getNodeToFaces() const override;

  ToCellRelation< ArrayOfArrays< localIndex > > getNodeToElements() const override;

  array2d< localIndex > getEdgeToNodes() const override;

  ArrayOfArrays< localIndex > getEdgeToFaces() const override;

  ArrayOfArrays< localIndex > getFaceToNodes() const override;

  ArrayOfArrays< localIndex > getFaceToEdges() const override;

  ToCellRelation< array2d< localIndex > > getFaceToElements() const override;

  array1d< globalIndex > getNodeLocalToGlobal() const override;

  /**
   * @brief Returns a view to the vector holding the node to global mapping.
   * @return The reference
   *
   * @note This is meant to be used as a values setter.
   */
  arrayView1d< globalIndex > getNodeLocalToGlobal();

  std::map< string, SortedArray< localIndex > > const & getNodeSets() const override;

  /**
   * @brief Returns a mutable reference to the node sets.
   * @return A reference to the mapping.
   *
   * The key of the map is the name of the set.
   * While the values are sorted arrays which sizes are meant to be managed by the client code.
   * This member function is meant to be used like a setter.
   */
  std::map< string, SortedArray< localIndex > > & getNodeSets();

  /**
   * @brief Defines the number of nodes and resizes some underlying arrays appropriately.
   * @param[in] numNodes The number of nodes on the MPI rank (that is per domain).
   * Nodes are not duplicated along subregion interfaces.
   *
   * The nodes coordinates and nodes local to global mappings get resized to @p numNodes.
   */
  void setNumNodes( localIndex numNodes );

  void generateHighOrderMaps( localIndex const order,
                              globalIndex const maxVertexGlobalID,
                              globalIndex const maxEdgeGlobalID,
                              globalIndex const maxFaceGlobalID,
                              arrayView1d< globalIndex const > const edgeLocalToGlobal,
                              arrayView1d< globalIndex const > const faceLocalToGlobal ) override;

  localIndex numNodes() const override;

  localIndex numEdges() const override;

  localIndex numFaces() const override;

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param[in] numElements list of the new element numbers
   * @param[in] regionNames list of the element region names
   */
  void resize( integer_array const & numElements,
               string_array const & regionNames );

  /**
   * @brief Trigger the computation of all the mappings.
   *
   * Call this member function to compute all the mappings.
   * Computations could be done lazily when calling getters.
   * But this is not yet implemented.
   */
  void buildMaps();

  /**
   * @brief Get cell block by name.
   * @param[in] name Name of the cell block.
   * @return Reference to the cell block instance.
   */
  CellBlock & getCellBlock( string const & name )
  {
    return this->getGroup( viewKeyStruct::cellBlocks() ).getGroup< CellBlock >( name );
  }

  const Group & getCellBlocks() const override;

  Group & getCellBlocks() override;

  std::map< integer, std::set< string > > const & getRegionAttributesCellBlocks() const override
  { return m_regionAttributesCellBlocks; }

  Group const & getFaceBlocks() const override;

  Group & getFaceBlocks() override;

  Group const & getEmbeddedSurfaceBlocks() const override;

  Group & getEmbeddedSurfaceBlocks() override;

  LineBlockABC const & getLineBlock( string name ) const override;

  /**
   * @brief Registers and returns a cell block of name @p name.
   * @param name The name of the created cell block.
   * @return A reference to the new cell block. The CellBlockManager owns this new instance.
   */
  CellBlock & registerCellBlock( string const & name );

  /**
   * @brief Registers and returns a cell block of name @p name.
   * @param cellBlockName The name of the created cell block.
   * @param regionAttribute The region attribute of the created cell block.
   * @return A reference to the new cell block. The CellBlockManager owns this new instance.
   */
  CellBlock & registerCellBlock( string const & cellBlockName, integer regionAttribute );

  /**
   * @brief Registers and returns a face block of name @p name.
   * @param name The name of the created face block.
   * @return A reference to the new face block. The CellBlockManager owns this new instance.
   */
  FaceBlock & registerFaceBlock( string const & name );

  /**
   * @brief Registers and returns a line block of name @p name.
   * @param name The name of the created line block.
   * @return A reference to the new line block. The CellBlockManager owns this new instance.
   */
  LineBlock & registerLineBlock( string const & name );
  /**
   * @brief Registers and returns an embedded surface block of name @p name.
   * @param name The name of the created embedded surface block.
   * @return A reference to the new embedded surface block. The CellBlockManager owns this new instance.
   */
  EmbeddedSurfaceBlock & registerEmbeddedSurfaceBlock( string const & name );
  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA lambda )
  {
    this->getGroup( viewKeyStruct::cellBlocks() ).forSubGroups< CellBlock >( lambda );
  }

  real64 getGlobalLength() const override { return m_globalLength; }

  /**
   * @brief Setter for the global length
   * @param globalLength the global length
   */
  void setGlobalLength( real64 globalLength ) { m_globalLength = globalLength; }

private:

  struct viewKeyStruct
  {
    /// Cell blocks key
    static constexpr char const * cellBlocks()
    { return "cellBlocks"; }

    /// Face blocks key
    static constexpr char const * faceBlocks()
    { return "faceBlocks"; }

    /// Line blocks key
    static constexpr char const * lineBlocks()
    { return "lineBlocks"; }
    /// Embedded Surface blocks key
    static constexpr char const * embeddedSurfaceBlocks()
    { return "embeddedSurfaceBlocks"; }
  };

  /**
   * @brief Returns a group containing the well blocks as @p LineBlockABC instances.
   * @return Mutable reference to the well blocks group.
   *
   * @note It should probably be better not to expose a non-const accessor here.
   */
  Group & getLineBlocks();

  /**
   * @brief Get cell block at index @p blockIndex.
   * @param[in] blockIndex The cell block index.
   * @return Const reference to the instance.
   *
   * @note Mainly useful for iteration purposes.
   */
  CellBlock const & getCellBlock( localIndex const blockIndex ) const;

  /**
   * @brief Returns the number of cells blocks
   * @return Number of cell blocks
   */
  localIndex numCellBlocks() const;

  /**
   * @brief Trigger the node to edges mapping computation.
   */
  void buildNodeToEdges();

  /**
   * @brief Trigger the face to nodes, edges and elements mappings.
   */
  void buildFaceMaps();

  template< typename BASEMAP, typename FUNC >
  void buildToCellMap( localIndex const cellIndex,
                       ToCellRelation< BASEMAP > & toCells,
                       FUNC cellToObjectGetter,
                       localIndex const overAlloc = 0 ) const;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_nodesPositions;

  ArrayOfArrays< localIndex > m_nodeToEdges;
  ArrayOfArrays< localIndex > m_edgeToFaces;
  array2d< localIndex > m_edgeToNodes;
  ArrayOfArrays< localIndex > m_faceToNodes;
  ArrayOfArrays< localIndex > m_faceToEdges;
  ToCellRelation< array2d< localIndex > > m_faceToCells;

  array1d< globalIndex > m_nodeLocalToGlobal;

  std::map< string, SortedArray< localIndex > > m_nodeSets;

  std::map< integer, std::set< string > > m_regionAttributesCellBlocks;

  real64 m_globalLength;

  localIndex m_numNodes;
  localIndex m_numFaces;
  localIndex m_numEdges;
};

}
#endif /* GEOS_MESH_CELLBLOCKMANAGER_H_ */
