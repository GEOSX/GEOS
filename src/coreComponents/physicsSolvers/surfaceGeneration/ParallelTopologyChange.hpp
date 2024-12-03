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
 * @file ParallelTopologyChange.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_
#define GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_

#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"

#define PARALLEL_TOPOLOGY_CHANGE_METHOD 1

namespace geos
{
class MeshLevel;
class NeighborCommunicator;
struct ModifiedObjectLists;

namespace parallelTopologyChange
{

void synchronizeTopologyChange( MeshLevel * const mesh,
                                std::vector< NeighborCommunicator > & neighbors,
                                ModifiedObjectLists & modifiedObjects,
                                ModifiedObjectLists & receivedObjects,
                                int mpiCommOrder );



struct TopologyChangeStepData
{
  void init( ElementRegionManager const & elemManager )
  {
    m_nodes.resize( 0 );
    m_edges.resize( 0 );
    m_faces.resize( 0 );
    m_elements.resize( elemManager.numRegions() );
    m_elementsView.resize( elemManager.numRegions() );
    m_elementsData.resize( elemManager.numRegions() );
    m_size = 0;

    for( localIndex er=0; er<elemManager.numRegions(); ++er )
    {
      ElementRegionBase const & elemRegion = elemManager.getRegion( er );
      m_elements[er].resize( elemRegion.numSubRegions());
      m_elementsView[er].resize( elemRegion.numSubRegions());
      m_elementsData[er].resize( elemRegion.numSubRegions());
      for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
      {
        m_elementsData[er][esr].resize( 0 );
        m_elements[er][esr].set( m_elementsData[er][esr] );
        m_elementsView[er][esr] = m_elementsData[er][esr];
      }
    }
  }


  localIndex_array m_nodes;
  localIndex_array m_edges;
  localIndex_array m_faces;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > m_elements;
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > m_elementsView;

  array1d< array1d< localIndex_array > > m_elementsData;
  buffer_type::size_type m_size;

};

struct TopologyChangeUnpackStepData : public TopologyChangeStepData
{
  void init( buffer_type const & receiveBuffer,
             ElementRegionManager const & elemManager )
  {
    m_bufferPtr = receiveBuffer.data();
    TopologyChangeStepData::init( elemManager );
  }

  buffer_unit_type const * m_bufferPtr;
};

}

}

#endif /* GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_ */
