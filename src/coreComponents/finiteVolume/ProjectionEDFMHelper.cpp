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

#include "ProjectionEDFMHelper.hpp"
#include "mesh/MeshLevel.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/SurfaceElementStencil.hpp"
#include "finiteVolume/EmbeddedSurfaceToCellStencil.hpp"

#include <limits>  // numeric_limits

namespace geos
{

ProjectionEDFMHelper::ProjectionEDFMHelper( MeshLevel const & mesh,
                                            CellElementStencilTPFA & cellStencil,
                                            EmbeddedSurfaceToCellStencil & edfmStencil,
                                            string const & embeddedSurfaceRegionName )
  : m_elementManager( mesh.getElemManager() ),
  m_nodesCoord( mesh.getNodeManager().referencePosition() ),
  m_edgeToNodes( mesh.getEdgeManager().nodeList() ),
  m_faceToRegions( mesh.getFaceManager().elementRegionList() ),
  m_faceToSubRegions( mesh.getFaceManager().elementSubRegionList() ),
  m_faceToCells( mesh.getFaceManager().elementList() ),
  m_faceToEdges( mesh.getFaceManager().edgeList().toViewConst() ),
  m_faceToNodes( mesh.getFaceManager().nodeList().toViewConst() ),
  m_cellCenters( m_elementManager.constructArrayViewAccessor< real64, 2 >( CellElementSubRegion::viewKeyStruct::elementCenterString() ) ),
  m_cellStencil( cellStencil ),
  m_edfmStencil( edfmStencil ),
  m_embeddedSurfaceRegionName( embeddedSurfaceRegionName )
{}

void ProjectionEDFMHelper::addNonNeighboringConnections() const
{
  SurfaceElementRegion const & fractureRegion = m_elementManager.getRegion< SurfaceElementRegion >( m_embeddedSurfaceRegionName );
  EmbeddedSurfaceSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< EmbeddedSurfaceSubRegion >();

  arrayView1d< integer const > const ghostRank = fractureSubRegion.ghostRank();
  OrderedVariableToManyElementRelation const & surfaceElementsToCells = fractureSubRegion.getToCellRelation();

  for( localIndex fracElement = 0; fracElement < fractureSubRegion.size(); fracElement++ )
  {
    if( ghostRank[fracElement] < 0 )
    {
      // get host cell information
      localIndex const hostCellRegionIdx  = surfaceElementsToCells.m_toElementRegion[fracElement][0];
      localIndex const hostCellSubRegionIdx = surfaceElementsToCells.m_toElementSubRegion[fracElement][0];
      localIndex const hostCellIdx = surfaceElementsToCells.m_toElementIndex[fracElement][0];
      CellDescriptor cellID( hostCellRegionIdx, hostCellSubRegionIdx, hostCellIdx );

      // get host cell faces
      CellElementRegion const & cellRegion = m_elementManager.getRegion< CellElementRegion >( hostCellRegionIdx );
      CellElementSubRegion const & cellSubRegion = cellRegion.getSubRegion< CellElementSubRegion >( hostCellSubRegionIdx );

      // pick faces for non-neighboring pEDFM connections
      std::vector< localIndex > const faces = selectFaces( cellSubRegion.faceList(), cellID, fracElement, fractureSubRegion );
      for( localIndex const faceIdx : faces )
      {
        CellDescriptor neighborCell = otherCell( faceIdx, cellID );
        real64 fractureMatrixWeights[2];
        computeFractureMatrixWeights( neighborCell, fracElement, fractureSubRegion, faceIdx, fractureMatrixWeights );

        addNonNeighboringConnection( fracElement, neighborCell, fractureMatrixWeights, fractureSubRegion );

        // zero out matrix-matrix connections that are replaced by fracture-matrix connections
        // I assume that the m-m connection index is equal to the face id
        m_cellStencil.zero( faceIdx );
      }
    }
  }
}

std::vector< localIndex > ProjectionEDFMHelper::selectFaces( FixedOneToManyRelation const & subRegionFaces,
                                                             CellDescriptor const & hostCellID,
                                                             localIndex const fracElement,
                                                             EmbeddedSurfaceSubRegion const & fractureSubRegion ) const
{
  arraySlice1d< real64 const > const n = fractureSubRegion.getNormalVector( fracElement );
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();

  real64 fracCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( centers[ fracElement ] );
  real64 cellCenterToFracCenter[ 3 ];
  real64 distToFrac = getSignedDistanceCellCenterToFracPlane( hostCellID, n, fracCenter, cellCenterToFracCenter );

  // pick faces that intersect the fracture
  std::vector< localIndex > faces;
  for( localIndex const iface : subRegionFaces[hostCellID.index] )
  {
    if( isBoundaryFace( iface ))
    {
      continue;
    }
    // face intersected by frac and non-branching
    if( intersection( fracCenter, n, iface ) && !neighborOnSameSide( iface, distToFrac, hostCellID, fractureSubRegion ) )
    {
      faces.push_back( iface );
      continue;
    }

    // if the frac is horizontal this can be 0 so we just pick one side.
    if( distToFrac < std::numeric_limits< real64 >::epsilon() )
      distToFrac = 1.0;

    if( onLargerSide( iface, distToFrac, fracCenter, n ) )
    {
      faces.push_back( iface );
      continue;
    }
  }

  return faces;
}

bool ProjectionEDFMHelper::intersection( real64 const (&fracCenter)[3],
                                         arraySlice1d< real64 const > const & fracNormal,
                                         localIndex const faceIdx ) const
{
  for( localIndex iEdge = 0; iEdge < m_faceToEdges[faceIdx].size(); iEdge++ )
  {
    localIndex edgeIdx = m_faceToEdges[faceIdx][iEdge];

    real64 signedDistanceProduct = 1.f;
    localIndex const nEdgeVertices = 2;
    for( localIndex ivertex = 0; ivertex < nEdgeVertices; ivertex++ )
    {
      real64 tmp[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_nodesCoord[m_edgeToNodes[edgeIdx][ivertex]] );
      LvArray::tensorOps::subtract< 3 >( tmp, fracCenter );
      signedDistanceProduct *= LvArray::tensorOps::AiBi< 3 >( tmp, fracNormal );
    }

    if( signedDistanceProduct <= 0 )
    {
      return true;
    }
  }
  return false;
}

bool ProjectionEDFMHelper::isBoundaryFace( localIndex const faceIdx ) const
{
  return ( m_faceToCells[faceIdx][0] < 0 || m_faceToCells[faceIdx][1] < 0 );
}

bool ProjectionEDFMHelper::onLargerSide( localIndex const faceIdx,
                                         real64 const signedDistanceCellCenterToFrac,
                                         real64 const (&fracCenter)[3],
                                         arraySlice1d< real64 const > const & fracNormal ) const
{
  /* Check If face is on the same side of the fracture as the mass center of the cell. */
  real64 const areaTolerance = 10.f * std::numeric_limits< real64 >::epsilon();
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  // compute face center
  computationalGeometry::centroid_3DPolygon( m_faceToNodes[faceIdx], m_nodesCoord,
                                             faceCenter, faceNormal, areaTolerance );
  LvArray::tensorOps::subtract< 3 >( faceCenter, fracCenter );
  return LvArray::tensorOps::AiBi< 3 >( faceCenter, fracNormal ) * signedDistanceCellCenterToFrac > 0;
}

real64 ProjectionEDFMHelper::getSignedDistanceCellCenterToFracPlane( CellDescriptor const & hostCellID,
                                                                     arraySlice1d< real64 const > const & fracNormal,
                                                                     real64 const (&fracCenter)[3],
                                                                     real64 (& cellCenterToFracOrigin)[3] ) const
{
  LvArray::tensorOps::copy< 3 >( cellCenterToFracOrigin, m_cellCenters[ hostCellID.region ]
                                 [ hostCellID.subRegion ]
                                 [ hostCellID.index ] );
  LvArray::tensorOps::subtract< 3 >( cellCenterToFracOrigin, fracCenter );

  return LvArray::tensorOps::AiBi< 3 >( cellCenterToFracOrigin, fracNormal );
}

CellDescriptor ProjectionEDFMHelper::otherCell( localIndex const faceIdx, CellDescriptor const & hostCellID ) const
{
  auto const & neighbors = m_faceToCells[faceIdx];
  localIndex const ineighbor = (neighbors[0] == hostCellID.index) ? 1 : 0;
  return CellDescriptor ( m_faceToRegions[faceIdx][ineighbor],
                          m_faceToSubRegions[faceIdx][ineighbor],
                          m_faceToCells[faceIdx][ineighbor] );
}

bool ProjectionEDFMHelper::neighborOnSameSide( localIndex const faceIdx,
                                               real64 const signedDistanceCellCenterToFrac,
                                               CellDescriptor const & hostCellID,
                                               EmbeddedSurfaceSubRegion const & fractureSubRegion ) const
{
  // find the identification of the neighbor cell
  CellDescriptor const neighborCellID = otherCell( faceIdx, hostCellID );

  // get fracture normal and origin in the neighbor cells:
  // they might be different than those in the host cell
  CellElementRegion const & neighborRegion =
    m_elementManager.getRegion< CellElementRegion >( neighborCellID.region );
  CellElementSubRegion const & neighborSubRegion =
    neighborRegion.getSubRegion< CellElementSubRegion >( neighborCellID.subRegion );
  auto const & efracSurfaces = neighborSubRegion.embeddedSurfacesList();

  // embedded fracture elements hosted in the neighbor cell
  // (there could be more than one)
  auto const & efracElements = efracSurfaces[ neighborCellID.index ];
  if( efracElements.size() == 0 )
  {
    return false;
  }
  GEOS_ERROR_IF( efracElements.size() > 1, "pEDFM with fracture intersections is not supported yet." );

  localIndex const fracElement = efracElements[0];
  arraySlice1d< real64 const > const n = fractureSubRegion.getNormalVector( fracElement );
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();


  real64 fracCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( centers[ fracElement ] );
  real64 cellCenterToFracOrigin[ 3 ];
  real64 const signedDistanceNeighborlCenterToFrac =
    getSignedDistanceCellCenterToFracPlane( neighborCellID, n, fracCenter, cellCenterToFracOrigin );

  return signedDistanceCellCenterToFrac * signedDistanceNeighborlCenterToFrac > 0;
}

void ProjectionEDFMHelper::addNonNeighboringConnection( localIndex const fracElement,
                                                        CellDescriptor const & cell,
                                                        real64 const (&weights)[2],
                                                        EmbeddedSurfaceSubRegion const & fractureSubRegion ) const
{
  array1d< localIndex > stencilCellsRegionIndex( 2 );
  array1d< localIndex > stencilCellsSubRegionIndex( 2 );
  array1d< localIndex > stencilCellsIndex( 2 );
  array1d< real64 > stencilWeights( 2 );

  // cell data
  stencilCellsRegionIndex[0] = cell.region;
  stencilCellsSubRegionIndex[0] = cell.subRegion;
  stencilCellsIndex[0] = cell.index;
  stencilWeights[0] =  weights[0];

  // fracture data
  stencilCellsRegionIndex[1] = fractureSubRegion.getParent().getParent().getIndexInParent();
  stencilCellsSubRegionIndex[1] = fractureSubRegion.getIndexInParent();
  stencilCellsIndex[1] = fracElement;
  stencilWeights[1] = weights[1];

  m_edfmStencil.add( 2,
                     stencilCellsRegionIndex.data(),
                     stencilCellsSubRegionIndex.data(),
                     stencilCellsIndex.data(),
                     stencilWeights.data(),
                     m_edfmStencil.size() );
}

void ProjectionEDFMHelper::
  computeFractureMatrixWeights( CellDescriptor const & neighborCell,
                                localIndex const fracElement,
                                EmbeddedSurfaceSubRegion const & fractureSubRegion,
                                localIndex const faceIdx,
                                real64 (& weights)[2] ) const
{

  // TODO: should I really compute the real projection, or assuming the whole face is occupied by the projection?

  // compute face center and normal
  real64 const areaTolerance = 10.f * std::numeric_limits< real64 >::epsilon();
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  real64 faceArea = computationalGeometry::centroid_3DPolygon( m_faceToNodes[faceIdx], m_nodesCoord,
                                                               faceCenter, faceNormal, areaTolerance );

  // get cell center
  real64 cellCenter[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_cellCenters[ neighborCell.region ]
                                                           [ neighborCell.subRegion ]
                                                           [ neighborCell.index ] );

  // get frac  center
  arrayView2d< real64 const > const & centers = fractureSubRegion.getElementCenter().toViewConst();
  real64 fracCenter[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( centers[ fracElement ] );

  // Compute projection point: a point  on the face that resides  on a  line
  // connecting edfm element and a neighbor cell
  real64 faceToFracVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( faceCenter );
  real64 cellToFracVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( cellCenter );
  LvArray::tensorOps::subtract< 3 >( faceToFracVector, fracCenter );
  LvArray::tensorOps::subtract< 3 >( cellToFracVector, fracCenter );

  real64 const numerator = LvArray::tensorOps::AiBi< 3 >( faceToFracVector, faceNormal );
  real64 const denom = LvArray::tensorOps::AiBi< 3 >( cellToFracVector, faceNormal );
  real64 const t = numerator / denom;

  real64 projectionPoint[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracCenter );
  LvArray::tensorOps::scaledAdd< 3 >( projectionPoint, cellToFracVector, t );

  // Compute half geometric transmissibilities
  real64 projectionPointToCellCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( cellCenter );
  LvArray::tensorOps::subtract< 3 >( projectionPointToCellCenter, projectionPoint );
  weights[0] = faceArea * LvArray::tensorOps::l2Norm< 3 >( projectionPointToCellCenter );

  real64 projectionPointToFracCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracCenter );
  LvArray::tensorOps::subtract< 3 >( projectionPointToFracCenter, projectionPoint );
  weights[1] = faceArea * LvArray::tensorOps::l2Norm< 3 >( projectionPointToFracCenter );
}

}  // end namespace geos
