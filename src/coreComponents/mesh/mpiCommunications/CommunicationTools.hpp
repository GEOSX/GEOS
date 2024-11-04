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
 * @file CommunicationTools.hpp
 */

#ifndef GEOS_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_

#include "CommID.hpp"

#include "common/MpiWrapper.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "mesh/FieldIdentifiers.hpp"

#include <set>

namespace geos
{


class ObjectManagerBase;
class NodeManager;
class NeighborCommunicator;
class MeshLevel;
class ElementRegionManager;

class MPI_iCommData;



class CommunicationTools
{
public:
  CommunicationTools();
  ~CommunicationTools();

  static CommunicationTools & getInstance();

  /**
   * @brief Assign the global indices to the objects managed by @p manager.
   * @param manager The manager of the objects that will receive the global indices.
   * @param compositionManager The manager of the objects on which the objects of @p manager are defined. (Mainly a node manager.)
   * @param neighbors The neighbors that may share objects of @p manager.
   *
   * @note The global indices will not continuously span the [min global index, max global index] range.
   * There will be holes in the enumeration.
   * The global indices are meant to be used as label to identify objects across ranks.
   * They are not meant to iterate.
   *
   * @example For example, faces hold by @p manager (@p FaceManager) are defined on nodes managed by @p compositionManager (@p NodeManager).
   * And two faces sharing the same nodes will be considered identical and should therefore share the same global index.
   */
  void assignGlobalIndices( ObjectManagerBase & manager,
                            NodeManager const & compositionManager,
                            std::vector< NeighborCommunicator > & neighbors );

  static void assignNewGlobalIndices( ObjectManagerBase & manager,
                                      std::set< localIndex > const & indexList );

  static void assignNewGlobalIndices( ElementRegionManager & elementManager,
                                      std::map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems );

  void setupGhosts( MeshLevel & meshLevel,
                    std::vector< NeighborCommunicator > & neighbors,
                    bool use_nonblocking );

  CommID getCommID()
  { return CommID( m_freeCommIDs ); }

  void findMatchedPartitionBoundaryObjects( ObjectManagerBase & group,
                                            std::vector< NeighborCommunicator > & allNeighbors );

  void findMatchedPartitionBoundaryNodes( NodeManager & nodeManager,
                                          std::vector< NeighborCommunicator > & allNeighbors,
                                          std::set< std::set< globalIndex > > const & collocatedNodesBuckets,
                                          std::set< globalIndex > const & requestedNodes );

  void synchronizeFields( FieldIdentifiers const & fieldsToBeSync,
                          MeshLevel & mesh,
                          std::vector< NeighborCommunicator > & allNeighbors,
                          bool onDevice );

  void synchronizePackSendRecvSizes( FieldIdentifiers const & fieldsToBeSync,
                                     MeshLevel & mesh,
                                     std::vector< NeighborCommunicator > & neighbors,
                                     MPI_iCommData & icomm,
                                     bool onDevice );

  void synchronizePackSendRecv( FieldIdentifiers const & fieldsToBeSync,
                                MeshLevel & mesh,
                                std::vector< NeighborCommunicator > & allNeighbors,
                                MPI_iCommData & icomm,
                                bool onDevice );

  void asyncPack( FieldIdentifiers const & fieldsToBeSync,
                  MeshLevel & mesh,
                  std::vector< NeighborCommunicator > & neighbors,
                  MPI_iCommData & icomm,
                  bool onDevice,
                  parallelDeviceEvents & events );

  void asyncSendRecv( std::vector< NeighborCommunicator > & neighbors,
                      MPI_iCommData & icomm,
                      bool onDevice,
                      parallelDeviceEvents & events );

  void synchronizeUnpack( MeshLevel & mesh,
                          std::vector< NeighborCommunicator > & neighbors,
                          MPI_iCommData & icomm,
                          bool onDevice );

  bool asyncUnpack( MeshLevel & mesh,
                    std::vector< NeighborCommunicator > & neighbors,
                    MPI_iCommData & icomm,
                    bool onDevice,
                    parallelDeviceEvents & events,
                    MPI_Op op=MPI_REPLACE );

  void finalizeUnpack( MeshLevel & mesh,
                       std::vector< NeighborCommunicator > & neighbors,
                       MPI_iCommData & icomm,
                       bool onDevice,
                       parallelDeviceEvents & events,
                       MPI_Op op=MPI_REPLACE );

  static void checkSendRecv( ObjectManagerBase const & objectManager,
                      std::vector< NeighborCommunicator > & neighbors );

private:
  std::set< int > m_freeCommIDs;
  static CommunicationTools * m_instance;

  /**
   * @brief Exchange the boundary objects managed by the @p manager and
   * find the objects that are equivalent in order to assign them a unique global id.
   * @param manager The instance managing the (boundary) objects.
   * @param allNeighbors The neighbors involved in the check.
   * @return For each neighbor, the global indices of boundary objects (not necessarily the matching objects).
   * @note The matched objects are stored in the @p NeighborData of the @p manager, for each neighbor.
   */
  array1d< array1d< globalIndex > >
  buildNeighborPartitionBoundaryObjects( ObjectManagerBase & manager,
                                         std::vector< NeighborCommunicator > & allNeighbors );

};



} /* namespace geos */

#endif /* GEOS_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_ */
