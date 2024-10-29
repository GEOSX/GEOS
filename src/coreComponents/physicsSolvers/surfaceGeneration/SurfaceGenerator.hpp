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
 * @file SurfaceGenerator.hpp
 */
#ifndef GEOS_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOR_HPP_
#define GEOS_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOR_HPP_

#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

struct ModifiedObjectLists
{
  std::set< localIndex > newNodes;
  std::set< localIndex > newEdges;
  std::set< localIndex > newFaces;
  std::set< localIndex > modifiedNodes;
  std::set< localIndex > modifiedEdges;
  std::set< localIndex > modifiedFaces;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > newElements;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > modifiedElements;

  void clearNewFromModified();

  void insert( ModifiedObjectLists const & lists );
};


class SpatialPartition;

class NodeManager;
class EdgeManager;
class FaceManager;
class ElementRegionManager;
class ElementRegionBase;

/**
 * @class SurfaceGenerator
 *
 * This solver manages the mesh topology splitting methods.
 *
 */
class SurfaceGenerator : public SolverBase
{
public:
  SurfaceGenerator( const string & name,
                    Group * const parent );
  ~SurfaceGenerator() override;


  static string catalogName() { return "SurfaceGenerator"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( Group & MeshBody ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                        DomainPartition & domain ) override
  {
    solverStep( time_n, dt, cycleNumber, domain );
    return false;
  }

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  /**@}*/

  inline string const getFractureRegionName() const { return m_fractureRegionName; }

  void postInputInitialization() override final;

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override final;
  virtual void postRestartInitialization() override final;

private:

  int separationDriver( DomainPartition & domain,
                        MeshLevel & mesh,
                        std::vector< NeighborCommunicator > & neighbors,
                        int const tileColor,
                        int const numTileColors,
                        const bool prefrac,
                        const real64 time_np1 );

  /**
   * @brief Function to generate new global indices of a simple object (node, edge, face)
   * @param[in,out] object A reference to the object that needs new global indices
   * @param[in] indexList the list of local indices that need new global indices
   */
  void assignNewGlobalIndicesSerial( ObjectManagerBase & object,
                                     std::set< localIndex > const & indexList );


  /**
   * @brief Function to generate new global indices for elements
   * @param[in,out] elementManager A reference to the ElementRegionManager that needs new global indices
   * @param[in] indexList the list of local indices that need new global indices
   */
  void
  assignNewGlobalIndicesSerial( ElementRegionManager & elementManager,
                                map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & indexList );

  /**
   * @brief Function to identify which faces are ready for rupture
   * @param domain
   * @param nodeManager Field @p SIFNode is modified.
   * @param edgeManager Fields @p SIF_I, @p SIF_II, @p SIF_III are modified.
   * @param faceManager Fields @p SIFonFace, @p RuptureState are modified.
   * @param elementManager
   * @param prefrac
   */
  void identifyRupturedFaces( DomainPartition const & domain,
                              NodeManager & nodeManager,
                              EdgeManager & edgeManager,
                              FaceManager & faceManager,
                              ElementRegionManager const & elementManager,
                              bool const prefrac );

  /**
   * @brief
   * @param edgeID
   * @param trailFaceID
   * @param nodeManager
   * @param edgeManager Fields @p SIF_I, @p SIF_II, @p SIF_III are modified.
   * @param faceManager
   * @param elementManager
   * @param vecTipNorm
   * @param vecTip
   * @return
   */
  real64 calculateEdgeSif( DomainPartition const & domain,
                           localIndex const edgeID,
                           localIndex & trailFaceID,
                           NodeManager const & nodeManager,
                           EdgeManager & edgeManager,
                           FaceManager const & faceManager,
                           ElementRegionManager const & elementManager,
                           real64 ( &vecTipNorm )[3],
                           real64 ( &vecTip )[3] );

  /**
   * @brief
   * @param nodeManager Field @p SIFNode gets modified.
   * @param edgeManager
   * @param faceManager Field @p SIFonFace gets modified.
   * @param elementManager
   * @return
   */
  void calculateNodeAndFaceSif( DomainPartition const & domain,
                                NodeManager & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager & faceManager,
                                ElementRegionManager const & elementManager );

  /**
   * @brief Function to calculate f_disconnect and f_u.
   * @param edgeID
   * @param edgeLength
   * @param nodeIndices
   * @param nodeManager
   * @param edgeManager
   * @param elementManager
   * @param vecTipNorm
   * @param fNode
   * @param GdivBeta
   * @param threeNodesPinched
   * @param calculatef_u. True: calculate f_u; False: calculate f_disconnect.
   */
  int calculateElementForcesOnEdge( DomainPartition const & domain,
                                    localIndex const edgeID,
                                    real64 edgeLength,
                                    localIndex_array & nodeIndices,
                                    NodeManager const & nodeManager,
                                    EdgeManager const & edgeManager,
                                    ElementRegionManager const & elementManager,
                                    real64 ( &vecTipNorm )[3],
                                    real64 ( &fNode )[3],
                                    real64 & GdivBeta,
                                    bool threeNodesPinched,
                                    bool calculatef_u );

  /**
   * @brief
   * @param edgeID
   * @param trailFaceID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager Field @p SIFonFace is modified. @p KIC, @p SIF_I, @p SIF_II are used.
   * @param elementManager
   * @param vecTipNorm
   * @param vecTip
   * @param modifiedObjects
   * @param edgeMode
   */
  void markRuptureFaceFromEdge( localIndex const edgeID,
                                localIndex const & trailFaceID,
                                NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager & faceManager,
                                ElementRegionManager const & elementManager,
                                real64 ( &vecTipNorm )[3],
                                real64 ( &vecTip )[3],
                                ModifiedObjectLists & modifiedObjects,
                                int const edgeMode );

  /**
   * @brief
   * @param nodeIndex
   * @param nodeManager
   * @param nodeManager
   * @param edgeManager
   * @param faceManager Field @p RuptureState gets modified.
   * @param modifiedObjects
   */
  void markRuptureFaceFromNode( localIndex const nodeIndex,
                                NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager & faceManager,
                                ElementRegionManager const & elementManager,
                                ModifiedObjectLists & modifiedObjects );

  /**
   *
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   */
  void postUpdateRuptureStates( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                ElementRegionManager const & elementManager,
                                std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                                std::vector< std::set< localIndex > > & edgesToRupturedFaces );

  /**
   *
   * @param elementManager
   * @param faceManager
   * @param iFace
   */
  int checkOrphanElement( ElementRegionManager const & elementManager,
                          FaceManager const & faceManager,
                          localIndex iFace );

  /**
   *
   * @param edgeID
   * @param nodeManager
   * @param faceManager
   * @param edgeManager
   * @param prefrac
   * @return
   */
  int checkEdgeSplitability( localIndex const edgeID,
                             NodeManager const & nodeManager,
                             FaceManager const & faceManager,
                             EdgeManager const & edgeManager,
                             const bool prefrac );


//  void UpdatePathCheckingArrays();

  /**
   * @brief check and split node in mesh
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elemManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param elementManager
   * @param modifiedObjects
   * @param prefrac
   * @return
   */
  bool processNode( const localIndex nodeID,
                    real64 const time,
                    NodeManager & nodeManager,
                    EdgeManager & edgeManager,
                    FaceManager & faceManager,
                    ElementRegionManager & elemManager,
                    std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                    std::vector< std::set< localIndex > > & edgesToRupturedFaces,
                    ElementRegionManager & elementManager,
                    ModifiedObjectLists & modifiedObjects,
                    const bool prefrac );

  /**
   * @brief Find a fracture path for surface generation
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elemManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param separationPathFaces
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool findFracturePlanes( localIndex const nodeID,
                           NodeManager const & nodeManager,
                           EdgeManager const & edgeManager,
                           FaceManager const & faceManager,
                           ElementRegionManager const & elemManager,
                           std::vector< std::set< localIndex > > const & nodesToRupturedFaces,
                           std::vector< std::set< localIndex > > const & edgesToRupturedFaces,
                           std::set< localIndex > & separationPathFaces,
                           map< localIndex, int > & edgeLocations,
                           map< localIndex, int > & faceLocations,
                           map< std::pair< CellElementSubRegion const *, localIndex >, int > & elemLocations );


  /**
   * @brief given a fracture path, split the mesh according to the path.
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param modifiedObjects
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param separationPathFaces
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   */
  void performFracture( localIndex const nodeID,
                        real64 const time_np1,
                        NodeManager & nodeManager,
                        EdgeManager & edgeManager,
                        FaceManager & faceManager,
                        ElementRegionManager & elementManager,
                        ModifiedObjectLists & modifiedObjects,
                        std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                        std::vector< std::set< localIndex > > & edgesToRupturedFaces,
                        std::set< localIndex > const & separationPathFaces,
                        map< localIndex, int > const & edgeLocations,
                        map< localIndex, int > const & faceLocations,
                        map< std::pair< CellElementSubRegion const *, localIndex >, int > const & elemLocations );

  void mapConsistencyCheck( localIndex const nodeID,
                            NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            ElementRegionManager const & elementManager,
                            map< std::pair< CellElementSubRegion const *, localIndex >, int > const & elemLocations );

  /**
   * @brief function to set which side of the fracture plane all objects are on
   * @param separationPathFaces
   * @param elemManager
   * @param faceManager
   * @param nodeToElementMaps The vector is assumed not to contain any duplicate.
   * @param localFacesToEdges
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool setLocations( std::set< localIndex > const & separationPathFaces,
                     ElementRegionManager const & elemManager,
                     FaceManager const & faceManager,
                     std::vector< std::pair< CellElementSubRegion const *, localIndex > > const & nodeToElementMaps,
                     map< localIndex, std::pair< localIndex, localIndex > > const & localFacesToEdges,
                     map< localIndex, int > & edgeLocations,
                     map< localIndex, int > & faceLocations,
                     map< std::pair< CellElementSubRegion const *, localIndex >, int > & elemLocations );

  /**
   * @brief function to set which side of the fracture plane all objects are on
   * @param location
   * @param elem
   * @param separationPathFaces
   * @param elemManager
   * @param faceManager
   * @param nodesToElements The vector is assumed not to contain any duplicate.
   * @param localFacesToEdges
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool setElemLocations( int const location,
                         std::pair< CellElementSubRegion const *, localIndex > const & elem,
                         std::set< localIndex > const & separationPathFaces,
                         ElementRegionManager const & elemManager,
                         FaceManager const & faceManager,
                         std::vector< std::pair< CellElementSubRegion const *, localIndex > > const & nodesToElements,
                         map< localIndex, std::pair< localIndex, localIndex > > const & localFacesToEdges,
                         map< localIndex, int > & edgeLocations,
                         map< localIndex, int > & faceLocations,
                         map< std::pair< CellElementSubRegion const *, localIndex >, int > & elemLocations );

  /**
   *
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  real64 calculateKinkAngle( localIndex const nodeID,
                             NodeManager const & nodeManager,
                             EdgeManager const & edgeManager,
                             FaceManager const & faceManager );

  /**
   *
   * @param faceManager
   * @param edgeManager Field @p kinkAngle is modified.
   * @param nodeManager
   * @param modifiedObjects
   * @param prefrac
   */
  void calculateKinkAngles( FaceManager const & faceManager,
                            EdgeManager & edgeManager,
                            NodeManager const & nodeManager,
                            ModifiedObjectLists const & modifiedObjects,
                            bool const prefrac );

  /**
   *
   * @param ModifiedObjectLists
   */
  void synchronizeTipSets( FaceManager & faceManager,
                           EdgeManager & edgeManager,
                           NodeManager & nodeManager,
                           ModifiedObjectLists & receivedObjects );


//  void setDegreeFromCrackTip( NodeManager & nodeManager,
//                              FaceManager & faceManager );

  /**
   *
   * @param edgeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  real64 minimumToughnessOnEdge( localIndex const edgeID,
                                 NodeManager const & nodeManager,
                                 EdgeManager const & edgeManager,
                                 FaceManager const & faceManager );

  /**
   *
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  real64 minimumToughnessOnNode( localIndex const nodeID,
                                 NodeManager const & nodeManager,
                                 EdgeManager const & edgeManager,
                                 FaceManager const & faceManager );


  real64 calculateRuptureRate( SurfaceElementRegion & faceElementRegion );

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * failCriterionString() { return "failCriterion"; }
    constexpr static char const * solidMaterialNameString() { return "solidMaterialNames"; }
    constexpr static char const * fExternalString() { return "fExternal"; }
    constexpr static char const * tipNodesString() { return "tipNodes"; }
    constexpr static char const * tipEdgesString() { return "tipEdges"; }
    constexpr static char const * tipFacesString() { return "tipFaces"; }
    constexpr static char const * trailingFacesString() { return "trailingFaces"; }
    constexpr static char const * fractureRegionNameString() { return "fractureRegion"; }
    constexpr static char const * mpiCommOrderString() { return "mpiCommOrder"; }
    constexpr static char const * isPoroelasticString() {return "isPoroelastic";}

    //TODO: rock toughness should be a material parameter, and we need to make rock toughness to KIC a constitutive
    // relation.
    constexpr static char const * rockToughnessString() { return "rockToughness"; }

//    //TODO: Once the node-based SIF criterion becomes mature and robust, remove the edge-based criterion.
    constexpr static char const * nodeBasedSIFString() { return "nodeBasedSIF"; }
  };


private:

  constexpr static real64 m_nonRuptureTime = 1e9;

  /// choice of failure criterion
  integer m_failCriterion=1;

  array1d< localIndex > m_solidMaterialFullIndex;

  int m_nodeBasedSIF;

  int m_isPoroelastic;

  real64 m_rockToughness;

  // Flag for consistent communication ordering
  int m_mpiCommOrder;

  /// set of separable faces
  SortedArray< localIndex > m_separableFaceSet;

  /// copy of the original node->face mapping prior to any separation
  ArrayOfSets< localIndex > m_originalNodetoFaces;

  /// copy of the original node->edge mapping prior to any separation
  ArrayOfSets< localIndex > m_originalNodetoEdges;

  /// copy of the original face->edge mapping prior to any separation
  ArrayOfArrays< localIndex >  m_originalFaceToEdges;

  /// collection of faces that have been used for separation of each node
  array1d< SortedArray< localIndex > > m_usedFacesForNode;

  /// copy of the original face->elemRegion mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemRegion;

  /// copy of the original face->elemSubRegion mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemSubRegion;

  /// copy of the original face->elemIndex mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemIndex;

  /// name of the element region to place all new fractures
  string m_fractureRegionName;

  SortedArray< localIndex > m_tipNodes;

  SortedArray< localIndex > m_tipEdges;

  SortedArray< localIndex > m_tipFaces;

  SortedArray< localIndex > m_trailingFaces;

  SortedArray< localIndex > m_faceElemsRupturedThisSolve;

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOR_HPP_ */
