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
 * @file TwoPointFluxApproximation.cpp
 *
 */
#include "TwoPointFluxApproximation.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/MpiWrapper.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/EmbeddedSurfaceToCellStencil.hpp"
#include "finiteVolume/FaceElementToCellStencil.hpp"
#include "finiteVolume/ProjectionEDFMHelper.hpp"
#include "finiteVolume/SurfaceElementStencil.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "LvArray/src/tensorOps.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;

TwoPointFluxApproximation::TwoPointFluxApproximation( string const & name,
                                                      Group * const parent )
  : FluxApproximationBase( name, parent )
{

  registerWrapper( viewKeyStruct::meanPermCoefficientString(), &m_meanPermCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "" );

  registerWrapper< CellElementStencilTPFA >( viewKeyStruct::cellStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper< SurfaceElementStencil >( viewKeyStruct::fractureStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper< EmbeddedSurfaceToCellStencil >( viewKeyStruct::edfmStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper< FaceElementToCellStencil >( viewKeyStruct::faceToCellStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::usePEDFMString(),
                   &m_useProjectionEmbeddedFractureMethod ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void TwoPointFluxApproximation::registerCellStencil( Group & stencilGroup ) const
{
  stencilGroup.registerWrapper< CellElementStencilTPFA >( viewKeyStruct::cellStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void TwoPointFluxApproximation::computeFractureStencil( MeshLevel & mesh ) const
{
  mesh.getElemManager().forElementRegions< SurfaceElementRegion >(
    [&]( SurfaceElementRegion const & region )
  {
    // For the moment, the feature is only available for DFM,
    // which is described in GEOSX with `FaceElementSubRegion`.
    if( region.subRegionType() == SurfaceElementRegion::SurfaceSubRegionType::faceElement )
    {
      string const & regionName = region.getName();
      addFractureFractureConnectionsDFM( mesh, regionName );
      addFractureMatrixConnectionsDFM( mesh, regionName );
    }
  } );
}

void TwoPointFluxApproximation::computeCellStencil( MeshLevel & mesh ) const
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  CellElementStencilTPFA & stencil = getStencil< CellElementStencilTPFA >( mesh, viewKeyStruct::cellStencilString() );

  arrayView2d< localIndex const > const & elemRegionList = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList = faceManager.elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( m_coeffName + viewKeyStruct::transMultiplierString() );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenter =
    elemManager.constructArrayViewAccessor< real64, 2 >( CellElementSubRegion::viewKeyStruct::elementCenterString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const elemGlobalIndex =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();

  // make a list of region indices to be included
  SortedArray< localIndex > regionFilter;
  arrayView1d< string const > const targetRegions = m_targetRegions.at( mesh.getParent().getParent().getName() );
  elemManager.forElementRegionsComplete< CellElementRegion >( targetRegions,
                                                              [&]( localIndex,
                                                                   localIndex const ei,
                                                                   CellElementRegion const & )
  {
    regionFilter.insert( ei );
  } );

  stencil.reserve( faceManager.size() );

  real64 const lengthTolerance = m_lengthScale * m_areaRelTol;
  real64 const areaTolerance = lengthTolerance * lengthTolerance;

  forAll< serialPolicy >( faceManager.size(), [=, &stencil]( localIndex const kf )
  {
    // Filter out boundary faces
    if( elemList[kf][0] < 0 || elemList[kf][1] < 0 || isZero( transMultiplier[kf] ) )
    {
      return;
    }

    // Filter out faces where neither cell is locally owned
    if( elemGhostRank[elemRegionList[kf][0]][elemSubRegionList[kf][0]][elemList[kf][0]] >= 0 &&
        elemGhostRank[elemRegionList[kf][1]][elemSubRegionList[kf][1]][elemList[kf][1]] >= 0 )
    {
      return;
    }

    // Filter out faces where either of two cells is outside of target regions
    if( !( regionFilter.contains( elemRegionList[kf][0] ) && regionFilter.contains( elemRegionList[kf][1] ) ) )
    {
      return;
    }

    real64 faceCenter[ 3 ], faceNormal[ 3 ], cellToFaceVec[2][ 3 ];
    real64 const faceArea = computationalGeometry::centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );

    if( faceArea < areaTolerance )
    {
      return;
    }

    stackArray1d< localIndex, 2 > regionIndex( 2 );
    stackArray1d< localIndex, 2 > subRegionIndex( 2 );
    stackArray1d< localIndex, 2 > elementIndex( 2 );
    stackArray1d< real64, 2 > stencilWeights( 2 );
    stackArray1d< real64, 2 > stencilStabilizationWeights( 2 );
    stackArray1d< globalIndex, 2 > stencilCellsGlobalIndex( 2 );

    for( localIndex ke = 0; ke < 2; ++ke )
    {
      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];
      localIndex const ei  = elemList[kf][ke];

      regionIndex[ke] = er;
      subRegionIndex[ke] = esr;
      elementIndex[ke] = ei;
      stencilCellsGlobalIndex[ke] = elemGlobalIndex[er][esr][ei];

      LvArray::tensorOps::copy< 3 >( cellToFaceVec[ke], faceCenter );
      LvArray::tensorOps::subtract< 3 >( cellToFaceVec[ke], elemCenter[er][esr][ei] );

      real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec[ke] );

      stencilWeights[ke] = faceArea / c2fDistance;
      stencilStabilizationWeights[ke] = faceArea * c2fDistance;
    }

    real64 const sumStabilizationWeight =
      ( stencilStabilizationWeights[0] + stencilStabilizationWeights[1] );

    // Ensure elements are added to stencil in order of global indices
    if( stencilCellsGlobalIndex[0] >= stencilCellsGlobalIndex[1] )
    {
      std::swap( regionIndex[0], regionIndex[1] );
      std::swap( subRegionIndex[0], subRegionIndex[1] );
      std::swap( elementIndex[0], elementIndex[1] );
      std::swap( stencilWeights[0], stencilWeights[1] );
      std::swap( stencilStabilizationWeights[0], stencilStabilizationWeights[1] );
      std::swap( cellToFaceVec[0][0], cellToFaceVec[1][0] );
      std::swap( cellToFaceVec[0][1], cellToFaceVec[1][1] );
      std::swap( cellToFaceVec[0][2], cellToFaceVec[1][2] );
    }

    stencil.add( 2,
                 regionIndex.data(),
                 subRegionIndex.data(),
                 elementIndex.data(),
                 stencilWeights.data(),
                 kf );

    stencil.addVectors( transMultiplier[kf], sumStabilizationWeight, faceNormal, cellToFaceVec );
  } );
}

void TwoPointFluxApproximation::registerFractureStencil( Group & stencilGroup ) const
{
  stencilGroup.registerWrapper< SurfaceElementStencil >( viewKeyStruct::fractureStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );

  stencilGroup.registerWrapper< EmbeddedSurfaceToCellStencil >( viewKeyStruct::edfmStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );

  stencilGroup.registerWrapper< FaceElementToCellStencil >( viewKeyStruct::faceToCellStencilString() ).
    setRestartFlags( RestartFlags::NO_WRITE );
}

void TwoPointFluxApproximation::addFractureFractureConnectionsDFM( MeshLevel & mesh,
                                                                   string const & faceElementRegionName ) const
{
  localIndex constexpr maxElems = SurfaceElementStencil::maxStencilSize;

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

  arrayView2d< real64 const > faceCenter = faceManager.faceCenter();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > X = nodeManager.referencePosition();

  SurfaceElementStencil & fractureStencil = getStencil< SurfaceElementStencil >( mesh, viewKeyStruct::fractureStencilString() );
  fractureStencil.setMeanPermCoefficient( m_meanPermCoefficient );
  fractureStencil.move( hostMemorySpace );

  SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( faceElementRegionName );
  localIndex const fractureRegionIndex = fractureRegion.getIndexInParent();

  FaceElementSubRegion & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion.faceList();

  arrayView1d< localIndex const > const & fractureConnectorsToEdges = fractureSubRegion.m_2dFaceToEdge;

  ArrayOfArraysView< localIndex const > const & fractureConnectorsToFaceElements = fractureSubRegion.m_2dFaceTo2dElems.toViewConst();

  SortedArrayView< localIndex const > const recalculateFractureConnectorEdges = fractureSubRegion.m_recalculateConnectionsFor2dFaces.toViewConst();

  // reserve memory for the connections of this fracture
  fractureStencil.reserve( fractureStencil.size() + recalculateFractureConnectorEdges.size() );

  // add new connectors/connections between face elements to the fracture stencil
  forAll< serialPolicy >( recalculateFractureConnectorEdges.size(),
                          [ recalculateFractureConnectorEdges,
                            fractureConnectorsToFaceElements,
                            fractureConnectorsToEdges,
                            &edgeManager,
                            X,
                            &faceMap,
                            faceCenter,
                            hydraulicAperture,
                            fractureRegionIndex,
                            elemGhostRank,
                            &fractureStencil]
                            ( localIndex const k )
  {
    localIndex const fci = recalculateFractureConnectorEdges[k];
    localIndex const numElems = fractureConnectorsToFaceElements.sizeOfArray( fci );
    // only do this if there are more than one element attached to the connector
    localIndex const edgeIndex = fractureConnectorsToEdges[fci];

    // For now, we do not filter out connections for which numElems == 1 in this function.
    // Instead, the filter takes place in the single-phase FluxKernels specialized for the SurfaceElementStencil
    // (see physicsSolvers/multiphysics/SinglePhaseProppantFluxKernels.cpp).
    // The reason for doing the filtering there and not here is that the ProppantTransport solver
    // needs the connections numElems == 1 to produce correct results.

    localIndex const connectorIndex = fractureStencil.size();

    GEOS_ERROR_IF( numElems > maxElems, "Max stencil size exceeded by fracture-fracture connector " << fci );

    stackArray1d< localIndex, maxElems > stencilCellsRegionIndex( numElems );
    stackArray1d< localIndex, maxElems > stencilCellsSubRegionIndex( numElems );
    stackArray1d< localIndex, maxElems > stencilCellsIndex( numElems );
    stackArray1d< real64, maxElems > stencilWeights( numElems );
    stackArray1d< R1Tensor, maxElems > stencilCellCenterToEdgeCenters( numElems );
    stackArray1d< integer, maxElems > isGhostConnectors( numElems );

    // get edge geometry
    real64 edgeCenter[3], edgeSegment[3];
    edgeManager.calculateCenter( edgeIndex, X, edgeCenter );
    edgeManager.calculateLength( edgeIndex, X, edgeSegment );
    real64 const edgeLength = LvArray::tensorOps::l2Norm< 3 >( edgeSegment );
    bool containsLocalElement = false;

    // loop over all face elements attached to the connector and add them to the stencil
    for( localIndex kfe=0; kfe<numElems; ++kfe )
    {
      localIndex const fractureElementIndex = fractureConnectorsToFaceElements[fci][kfe];

      // use straight difference between the edge center and face center for gradient length...
      // TODO: maybe do something better here??
      real64 cellCenterToEdgeCenter[ 3 ];
      LvArray::tensorOps::copy< 3 >( cellCenterToEdgeCenter, edgeCenter );
      LvArray::tensorOps::subtract< 3 >( cellCenterToEdgeCenter, faceCenter[ faceMap[fractureElementIndex][0] ] );

      // form the CellStencil entry
      stencilCellsRegionIndex[kfe] = fractureRegionIndex;
      stencilCellsSubRegionIndex[kfe] = 0;
      stencilCellsIndex[kfe] = fractureElementIndex;
      containsLocalElement = containsLocalElement || elemGhostRank[fractureRegionIndex][0][fractureElementIndex] < 0;

      // Note: this is done solely to avoid crashes when using the SolidMechanicsLagrangeContact that builds a stencil but does not register
      // the
      // hydraulicAperture
      real64 const aperture_h =  hydraulicAperture[fractureRegionIndex][0].size() == 0 ? 1.0 : hydraulicAperture[fractureRegionIndex][0][fractureElementIndex];

      stencilWeights[kfe] = aperture_h * edgeLength / LvArray::tensorOps::l2Norm< 3 >( cellCenterToEdgeCenter );

      LvArray::tensorOps::copy< 3 >( stencilCellCenterToEdgeCenters[kfe], cellCenterToEdgeCenter );
    }

    if( !containsLocalElement )
    {
      return;
    }

    // add/overwrite the stencil for index fci
    fractureStencil.add( numElems,
                         stencilCellsRegionIndex.data(),
                         stencilCellsSubRegionIndex.data(),
                         stencilCellsIndex.data(),
                         stencilWeights.data(),
                         connectorIndex );

    fractureStencil.add( numElems,
                         stencilCellCenterToEdgeCenters.data(),
                         connectorIndex );
  } );
}

void TwoPointFluxApproximation::cleanMatrixMatrixConnectionsDFM( MeshLevel & mesh,
                                                                 string const & faceElementRegionName ) const
{
  ElementRegionManager const & elemManager = mesh.getElemManager();

  CellElementStencilTPFA & cellStencil = getStencil< CellElementStencilTPFA >( mesh, viewKeyStruct::cellStencilString() );

  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( faceElementRegionName );
  FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();
  SortedArrayView< localIndex const > const newFaceElements = fractureSubRegion.m_newFaceElements.toViewConst();

  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion.faceList();

  forAll< serialPolicy >( newFaceElements.size(),
                          [newFaceElements,
                           &cellStencil,
                           &faceMap]( localIndex const k )
  {
    localIndex const kfe = newFaceElements[k];
    // Remove cell-to-cell connections from cell stencil and add in new connections.
    // This is there to shut off previously connected cells
    // that are not connected anymore due to dynamic fracturing.
    cellStencil.zero( faceMap[kfe][0] );
  } );
}


void TwoPointFluxApproximation::addFractureMatrixConnectionsDFM( MeshLevel & mesh,
                                                                 string const & faceElementRegionName ) const
{
  localIndex constexpr maxElems = SurfaceElementStencil::maxStencilSize;

  ElementRegionManager & elemManager = mesh.getElemManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  CellElementStencilTPFA & cellStencil = getStencil< CellElementStencilTPFA >( mesh, viewKeyStruct::cellStencilString() );
  FaceElementToCellStencil & faceToCellStencil = getStencil< FaceElementToCellStencil >( mesh, viewKeyStruct::faceToCellStencilString() );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenter =
    elemManager.constructArrayViewAccessor< real64, 2 >( CellElementSubRegion::viewKeyStruct::elementCenterString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

  arrayView1d< real64 const > faceArea   = faceManager.faceArea();
  arrayView2d< real64 const > faceCenter = faceManager.faceCenter();
  arrayView2d< real64 const > faceNormal = faceManager.faceNormal();

  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( m_coeffName + viewKeyStruct::transMultiplierString() );

  SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( faceElementRegionName );
  localIndex const fractureRegionIndex = fractureRegion.getIndexInParent();
  FaceElementSubRegion & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();
  OrderedVariableToManyElementRelation const & elems2dToElems3d = fractureSubRegion.getToCellRelation();

  SortedArrayView< localIndex const > const new2dElems = fractureSubRegion.m_newFaceElements.toViewConst();
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion.faceList();

  ArrayOfArraysView< localIndex const > elemRegionList = elems2dToElems3d.m_toElementRegion.toViewConst();
  ArrayOfArraysView< localIndex const > elemSubRegionList = elems2dToElems3d.m_toElementSubRegion.toViewConst();
  ArrayOfArraysView< localIndex const > elemList = elems2dToElems3d.m_toElementIndex.toViewConst();

  // reserve memory for the connections of this region
  if( cellStencil.size() != 0 )
  {
    faceToCellStencil.reserve( faceToCellStencil.size() + elems2dToElems3d.size() );
  }

  // We store the concerned region indices,
  // in order to only define the connections for the requested regions.
  auto const & regions = elemManager.getRegions();
  SortedArray< integer > regionIndices;
  for( auto const & bodyRegionNames: m_targetRegions )
  {
    for( string const & regionName: bodyRegionNames.second )
    {
      regionIndices.insert( regions.getIndex( regionName ) );
    }
  }

  forAll< serialPolicy >( new2dElems.size(),
                          [ new2dElems,
                            &elems2dToElems3d,
                            &faceToCellStencil,
                            &faceMap,
                            elemRegionList,
                            elemSubRegionList,
                            elemList,
                            elemGhostRank,
                            faceCenter,
                            elemCenter,
                            faceNormal,
                            faceArea,
                            hydraulicAperture,
                            transMultiplier,
                            regionIndices = regionIndices.toViewConst(),
                            fractureRegionIndex ] ( localIndex const k )
  {
    localIndex const kfe = new2dElems[k];
    {
      localIndex const numElems = elems2dToElems3d.m_toElementSubRegion.sizeOfArray( kfe );

      GEOS_ERROR_IF( numElems > maxElems, "Max stencil size exceeded by fracture-cell connector " << kfe );

      real64 cellToFaceVec[ 3 ], faceNormalVector[ 3 ];

      localIndex connectorIndex = faceToCellStencil.size();

      for( localIndex ke = 0; ke < numElems; ++ke )
      {
        connectorIndex += ke;

        localIndex const faceIndex = faceMap[kfe][ke];
        localIndex const er  = elemRegionList[kfe][ke];
        localIndex const esr = elemSubRegionList[kfe][ke];
        localIndex const ei  = elemList[kfe][ke];

        // remove cell-to-cell connections from cell stencil and add in new connections
        if( !regionIndices.contains( er ) )
        {
          continue;
        }

        // Filter out entries where both fracture and cell element are ghosted
        if( elemGhostRank[fractureRegionIndex][0][kfe] >= 0 && elemGhostRank[er][esr][ei] >= 0 )
        {
          continue;
        }

        LvArray::tensorOps::copy< 3 >( faceNormalVector, faceNormal[faceIndex] );

        LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter[faceIndex] );
        LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter[er][esr][ei] );

        real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

        real64 const ht = faceArea[faceIndex] / c2fDistance;
        // Note: this is done solely to avoid crashes when using the SolidMechanicsLagrangeContact that builds a stencil but does not
        // register the
        // hydraulicAperture
        real64 const aperture_h =  hydraulicAperture[fractureRegionIndex][0].size() == 0 ? 1.0 : hydraulicAperture[fractureRegionIndex][0][kfe];

        localIndex const stencilCellsRegionIndex[2]{ er, fractureRegionIndex };
        localIndex const stencilCellsSubRegionIndex[2]{ esr, 0 };
        localIndex const stencilCellsIndex[2]{ ei, kfe };
        real64 const stencilWeights[2]{ ht, 2. * faceArea[faceIndex] / aperture_h };

        faceToCellStencil.add( 2,
                               stencilCellsRegionIndex,
                               stencilCellsSubRegionIndex,
                               stencilCellsIndex,
                               stencilWeights,
                               connectorIndex );

        faceToCellStencil.addVectors( transMultiplier[faceIndex], faceNormalVector, cellToFaceVec );
      }
    }
  } );
}

void TwoPointFluxApproximation::addToFractureStencil( MeshLevel & mesh,
                                                      string const & faceElementRegionName ) const
{
  this->addFractureFractureConnectionsDFM( mesh, faceElementRegionName );
  this->cleanMatrixMatrixConnectionsDFM( mesh, faceElementRegionName );
  this->addFractureMatrixConnectionsDFM( mesh, faceElementRegionName );
}

void TwoPointFluxApproximation::addFractureMatrixConnectionsEDFM( MeshLevel & mesh,
                                                                  string const & embeddedSurfaceRegionName ) const
{
  // Add connections EmbeddedSurface to/from CellElements.
  ElementRegionManager & elemManager = mesh.getElemManager();

  EmbeddedSurfaceToCellStencil & edfmStencil = getStencil< EmbeddedSurfaceToCellStencil >( mesh, viewKeyStruct::edfmStencilString() );
  edfmStencil.move( hostMemorySpace );

  SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( embeddedSurfaceRegionName );
  localIndex const fractureRegionIndex = fractureRegion.getIndexInParent();

  EmbeddedSurfaceSubRegion & fractureSubRegion = fractureRegion.getUniqueSubRegion< EmbeddedSurfaceSubRegion >();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

  arrayView1d< real64 const > const faceArea = fractureSubRegion.getElementArea();

  OrderedVariableToManyElementRelation const & surfaceElementsToCells = fractureSubRegion.getToCellRelation();

  arrayView1d< real64 const >     const & connectivityIndex = fractureSubRegion.getConnectivityIndex();

  arrayView1d< integer const > const ghostRank = fractureSubRegion.ghostRank();

  // start from last connectorIndex from surface-To-cell connections
  localIndex connectorIndex = edfmStencil.size();
  localIndex constexpr MAX_NUM_ELEMS = EmbeddedSurfaceToCellStencil::maxStencilSize;

  // reserve memory for the connections of this fracture
  edfmStencil.reserve( edfmStencil.size() + fractureSubRegion.size() );

  // loop over the embedded surfaces and add connections to cellStencil
  for( localIndex kes = 0; kes < fractureSubRegion.size(); kes++ )
  {
    if( ghostRank[kes] < 0 )
    {
      localIndex const numElems = 2;   // there is a 1 to 1 relation

      GEOS_ERROR_IF( numElems > MAX_NUM_ELEMS, "Max stencil size exceeded by fracture-cell connector " << kes );

      stackArray1d< localIndex, MAX_NUM_ELEMS > stencilCellsRegionIndex( numElems );
      stackArray1d< localIndex, MAX_NUM_ELEMS > stencilCellsSubRegionIndex( numElems );
      stackArray1d< localIndex, MAX_NUM_ELEMS > stencilCellsIndex( numElems );
      stackArray1d< real64, MAX_NUM_ELEMS > stencilWeights( numElems );

      localIndex const er  = surfaceElementsToCells.m_toElementRegion[kes][0];
      localIndex const esr = surfaceElementsToCells.m_toElementSubRegion[kes][0];
      localIndex const ei  = surfaceElementsToCells.m_toElementIndex[kes][0];

      stencilCellsRegionIndex[0] = er;
      stencilCellsSubRegionIndex[0] = esr;
      stencilCellsIndex[0] = ei;
      stencilWeights[0] = connectivityIndex[kes];

      stencilCellsRegionIndex[1] = fractureRegionIndex;
      stencilCellsSubRegionIndex[1] = 0;
      stencilCellsIndex[1] = kes;
      stencilWeights[1] = 4. * faceArea[fractureRegionIndex] / hydraulicAperture[fractureRegionIndex][0][kes];

      edfmStencil.add( 2,
                       stencilCellsRegionIndex.data(),
                       stencilCellsSubRegionIndex.data(),
                       stencilCellsIndex.data(),
                       stencilWeights.data(),
                       connectorIndex );

      connectorIndex++;
    }
  }

}

void TwoPointFluxApproximation::addFractureFractureConnectionsEDFM( MeshLevel & mesh,
                                                                    string const & embeddedSurfaceRegionName ) const
{

  EdgeManager const & embSurfEdgeManager = mesh.getEmbSurfEdgeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();
  EmbeddedSurfaceNodeManager & nodeManager = mesh.getEmbSurfNodeManager();

  // Get the stencil
  SurfaceElementStencil & fractureStencil = getStencil< SurfaceElementStencil >( mesh, viewKeyStruct::fractureStencilString() );
  fractureStencil.move( hostMemorySpace );

  SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( embeddedSurfaceRegionName );
  localIndex const fractureRegionIndex = fractureRegion.getIndexInParent();

  EmbeddedSurfaceSubRegion & fractureSubRegion = fractureRegion.getUniqueSubRegion< EmbeddedSurfaceSubRegion >();
  arrayView2d< real64 const > const fractureElemCenter = fractureSubRegion.getElementCenter();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

  EdgeManager::FaceMapType const & edgeToEmbSurfacesMap = embSurfEdgeManager.faceList();

  localIndex constexpr maxElems = SurfaceElementStencil::maxStencilSize;

  localIndex connectorIndex = 0;

  // reserve memory for the connections of this fracture
  fractureStencil.reserve( fractureStencil.size() + embSurfEdgeManager.size() );

  // add new connectors/connections between embedded elements to the fracture stencil
  for( localIndex ke = 0; ke <  embSurfEdgeManager.size(); ke++ )
  {
    // for now there is no generation of new elements so we add all edges.
    localIndex const numElems = 2;  // hardcoded for now but unless there is an intersection it should always be 2.
    if( edgeToEmbSurfacesMap.sizeOfSet( ke ) > 1 ) // to be a connector it need to be attached to at least 2 elements.
    {

      GEOS_ERROR_IF( numElems > maxElems, "Max stencil size exceeded by fracture-fracture connector " << ke );

      stackArray1d< localIndex, maxElems > stencilCellsRegionIndex( numElems );
      stackArray1d< localIndex, maxElems > stencilCellsSubRegionIndex( numElems );
      stackArray1d< localIndex, maxElems > stencilCellsIndex( numElems );
      stackArray1d< real64, maxElems > stencilWeights( numElems );

      stackArray1d< R1Tensor, maxElems > stencilCellCenterToEdgeCenters( numElems );
      stackArray1d< integer, maxElems > isGhostConnectors( numElems );

      // TODO get edge geometry
      real64 edgeCenter[3], edgeSegment[3];
      embSurfEdgeManager.calculateCenter( ke, X, edgeCenter );
      embSurfEdgeManager.calculateLength( ke, X, edgeSegment );
      real64 const edgeLength  = LvArray::tensorOps::l2Norm< 3 >( edgeSegment );


      // loop over all embedded surface elements attached to the connector and add them to the stencil
      for( localIndex kes = 0; kes < numElems; kes++ )
      {
        localIndex const fractureElementIndex = edgeToEmbSurfacesMap[ke][kes];

        // compute distance between cell centers
        real64 cellCenterToEdgeCenter[ 3 ];
        LvArray::tensorOps::copy< 3 >( cellCenterToEdgeCenter, edgeCenter );
        LvArray::tensorOps::subtract< 3 >( cellCenterToEdgeCenter, fractureElemCenter[fractureElementIndex] );

        // form the CellStencil entry
        stencilCellsRegionIndex[kes]    = fractureRegionIndex;
        stencilCellsSubRegionIndex[kes] = 0;  // there is only one subregion.
        stencilCellsIndex[kes]          = fractureElementIndex;

        stencilWeights[kes] = hydraulicAperture[fractureRegionIndex][0][fractureElementIndex] * edgeLength / LvArray::tensorOps::l2Norm< 3 >( cellCenterToEdgeCenter );

        LvArray::tensorOps::copy< 3 >( stencilCellCenterToEdgeCenters[kes], cellCenterToEdgeCenter );
      }

      // add/overwrite the stencil for index fci
      fractureStencil.add( numElems,
                           stencilCellsRegionIndex.data(),
                           stencilCellsSubRegionIndex.data(),
                           stencilCellsIndex.data(),
                           stencilWeights.data(),
                           connectorIndex );

      fractureStencil.add( numElems,
                           stencilCellCenterToEdgeCenters.data(),
                           connectorIndex );

      connectorIndex++;
    }
  }
}

void TwoPointFluxApproximation::addEmbeddedFracturesToStencils( MeshLevel & mesh,
                                                                string const & embeddedSurfaceRegionName ) const
{

  addFractureFractureConnectionsEDFM( mesh, embeddedSurfaceRegionName );

  addFractureMatrixConnectionsEDFM( mesh, embeddedSurfaceRegionName );

  if( m_useProjectionEmbeddedFractureMethod )
  {
    EmbeddedSurfaceToCellStencil & edfmStencil = getStencil< EmbeddedSurfaceToCellStencil >( mesh, viewKeyStruct::edfmStencilString() );
    edfmStencil.move( hostMemorySpace );

    CellElementStencilTPFA & cellStencil = getStencil< CellElementStencilTPFA >( mesh, viewKeyStruct::cellStencilString() );
    cellStencil.move( hostMemorySpace );

    ProjectionEDFMHelper pedfmHelper( mesh, cellStencil, edfmStencil, embeddedSurfaceRegionName );
    pedfmHelper.addNonNeighboringConnections();
  }
}

void TwoPointFluxApproximation::registerBoundaryStencil( Group & stencilGroup, string const & setName ) const
{
  if( !stencilGroup.hasWrapper( setName ) )
  {
    // if not there yet, let's register the set name as a wrapper
    stencilGroup.registerWrapper< BoundaryStencil >( setName ).
      setRestartFlags( RestartFlags::NO_WRITE );
  }
}

void TwoPointFluxApproximation::computeBoundaryStencil( MeshLevel & mesh,
                                                        string const & setName,
                                                        SortedArrayView< localIndex const > const & faceSet ) const
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  BoundaryStencil & stencil = getStencil< BoundaryStencil >( mesh, setName );

  arrayView2d< localIndex const > const & elemRegionList     = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList  = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList           = faceManager.elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( m_coeffName + viewKeyStruct::transMultiplierString() );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenter =
    elemManager.constructArrayViewAccessor< real64, 2 >( CellElementSubRegion::viewKeyStruct::elementCenterString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  // TODO: can we look this up better?
  string const & meshBodyName = mesh.getParent().getParent().getName();
  arrayView1d< string const > const targetRegions = m_targetRegions.at( meshBodyName );

  ArrayOfArraysView< localIndex const > const faceToNodes = faceManager.nodeList().toViewConst();

  // make a list of region indices to be included
  SortedArray< localIndex > regionFilter;
  for( string const & regionName : targetRegions )
  {
    regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
  }

  constexpr localIndex numPts = BoundaryStencil::maxNumPointsInFlux;

  stackArray1d< localIndex, numPts > stencilRegionIndices( numPts );
  stackArray1d< localIndex, numPts > stencilSubRegionIndices( numPts );
  stackArray1d< localIndex, numPts > stencilElemOrFaceIndices( numPts );
  stackArray1d< real64, numPts > stencilWeights( numPts );

  real64 const lengthTolerance = m_lengthScale * m_areaRelTol;
  real64 const areaTolerance = lengthTolerance * lengthTolerance;

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve( faceSet.size() );
  for( localIndex kf : faceSet )
  {
    real64 faceCenter[ 3 ];
    real64 faceNormal[ 3 ];
    real64 const faceArea = computationalGeometry::centroid_3DPolygon( faceToNodes[kf],
                                                                       nodePosition,
                                                                       faceCenter,
                                                                       faceNormal,
                                                                       areaTolerance );

    for( localIndex ke = 0; ke < numPts; ++ke )
    {
      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];
      localIndex const ei  = elemList[kf][ke];

      // Filter out elements not locally present
      if( er < 0 )
      {
        continue;
      }

      // Filter out elements not in target regions
      if( !regionFilter.contains( er ) )
      {
        continue;
      }

      // Filter out ghosted elements - to be handled by the owning rank
      if( elemGhostRank[er][esr][ei] >= 0 )
      {
        continue;
      }

      real64 cellToFaceVec[ 3 ];
      LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
      LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter[ er ][ esr ][ ei ] );

      if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
      }

      real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );
      real64 const faceWeight = faceArea / c2fDistance;

      stencilRegionIndices[BoundaryStencil::Order::ELEM] = er;
      stencilSubRegionIndices[BoundaryStencil::Order::ELEM] = esr;
      stencilElemOrFaceIndices[BoundaryStencil::Order::ELEM] = ei;
      stencilWeights[BoundaryStencil::Order::ELEM] = faceWeight;

      stencilRegionIndices[BoundaryStencil::Order::FACE] = -1;
      stencilSubRegionIndices[BoundaryStencil::Order::FACE] = -1;
      stencilElemOrFaceIndices[BoundaryStencil::Order::FACE] = kf;
      stencilWeights[BoundaryStencil::Order::FACE] = -faceWeight;

      stencil.add( stencilRegionIndices.size(),
                   stencilRegionIndices.data(),
                   stencilSubRegionIndices.data(),
                   stencilElemOrFaceIndices.data(),
                   stencilWeights.data(),
                   kf );

      stencil.addVectors( transMultiplier[kf], faceNormal, cellToFaceVec );
    }
  }
}

void TwoPointFluxApproximation::registerAquiferStencil( Group & stencilGroup, string const & setName ) const
{
  registerBoundaryStencil( stencilGroup, setName );
}


void TwoPointFluxApproximation::computeAquiferStencil( DomainPartition & domain, MeshLevel & mesh ) const
{
  // The computation of the aquifer stencil weights requires three passes:
  //  - In the first pass, we count the number of individual aquifers
  //  - In the second pass, we compute the sum of the face areas for each aquifer
  //  - In the third pass, we compute the area fraction for each face

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  arrayView1d< real64 const > const & faceArea              = faceManager.faceArea();
  arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList          = faceManager.elementList();

  constexpr localIndex numPts = BoundaryStencil::maxNumPointsInFlux;

  // make a list of region indices to be included
  SortedArray< localIndex > regionFilter;
  for( string const & regionName : m_targetRegions.at( mesh.getParent().getParent().getName() ) )
  {
    regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
  }
  SortedArrayView< localIndex const > const regionFilterView = regionFilter.toViewConst();

  // Step 1: count individual aquifers

  std::map< string, localIndex > aquiferNameToAquiferId;
  localIndex aquiferCounter = 0;

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
  {
    aquiferNameToAquiferId[bc.getName()] = aquiferCounter;
    aquiferCounter++;
  } );

  // Step 2: sum the face areas for each individual aquifer

  array1d< real64 > globalSumFaceAreas( aquiferNameToAquiferId.size() );
  array1d< real64 > localSumFaceAreas( aquiferNameToAquiferId.size() );
  arrayView1d< real64 > const localSumFaceAreasView = localSumFaceAreas.toView();

  fsManager.apply< FaceManager,
                   AquiferBoundaryCondition >( 0.0,
                                               domain.getMeshBody( 0 ).getBaseDiscretization(),
                                               AquiferBoundaryCondition::catalogName(),
                                               [&] ( AquiferBoundaryCondition const & bc,
                                                     string const &,
                                                     SortedArrayView< localIndex const > const & targetSet,
                                                     FaceManager const &,
                                                     string const & )
  {
    RAJA::ReduceSum< parallelHostReduce, real64 > targetSetSumFaceAreas( 0.0 );
    forAll< parallelHostPolicy >( targetSet.size(), [=] ( localIndex const i )
    {
      localIndex const iface = targetSet[i];

      bool isOwnedAquiferFaceInTarget = false;
      for( localIndex ke = 0; ke < numPts; ++ke )
      {
        // Filter out elements not locally present
        if( elemRegionList[iface][ke] < 0 )
        {
          continue;
        }

        // Filter out elements not in target regions
        if( !regionFilterView.contains( elemRegionList[iface][ke] ))
        {
          continue;
        }

        localIndex const er  = elemRegionList[iface][ke];
        localIndex const esr = elemSubRegionList[iface][ke];
        localIndex const ei  = elemList[iface][ke];

        // Filter out ghosted elements
        if( elemGhostRank[er][esr][ei] >= 0 )
        {
          continue;
        }

        isOwnedAquiferFaceInTarget = true;
      }

      if( isOwnedAquiferFaceInTarget )
      {
        targetSetSumFaceAreas += faceArea[iface];
      }
    } );
    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );
    localSumFaceAreasView[aquiferIndex] += targetSetSumFaceAreas.get();
  } );

  MpiWrapper::allReduce( localSumFaceAreas.data(),
                         globalSumFaceAreas.data(),
                         localSumFaceAreas.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOS );

  // Step 3: compute the face area fraction for each connection, and insert into boundary stencil

  fsManager.apply< FaceManager,
                   AquiferBoundaryCondition >( 0.0,
                                               domain.getMeshBody( 0 ).getBaseDiscretization(),
                                               AquiferBoundaryCondition::catalogName(),
                                               [&] ( AquiferBoundaryCondition const & bc,
                                                     string const & setName,
                                                     SortedArrayView< localIndex const > const & targetSet,
                                                     FaceManager const &,
                                                     string const & )
  {
    BoundaryStencil & stencil = getStencil< BoundaryStencil >( mesh, setName );

    stackArray1d< localIndex, numPts > stencilRegionIndices( numPts );
    stackArray1d< localIndex, numPts > stencilSubRegionIndices( numPts );
    stackArray1d< localIndex, numPts > stencilElemOrFaceIndices( numPts );
    stackArray1d< real64, numPts > stencilWeights( numPts );

    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );

    stencil.reserve( targetSet.size() );
    for( localIndex iface : targetSet )
    {

      for( localIndex ke = 0; ke < numPts; ++ke )
      {
        // Filter out elements not locally present
        if( elemRegionList[iface][ke] < 0 )
        {
          continue;
        }

        // Filter out elements not in target regions
        if( !regionFilterView.contains( elemRegionList[iface][ke] ))
        {
          continue;
        }

        localIndex const er  = elemRegionList[iface][ke];
        localIndex const esr = elemSubRegionList[iface][ke];
        localIndex const ei  = elemList[iface][ke];

        // Filter out ghosted elements
        if( elemGhostRank[er][esr][ei] >= 0 )
        {
          continue;
        }

        stencilRegionIndices[BoundaryStencil::Order::ELEM] = er;
        stencilSubRegionIndices[BoundaryStencil::Order::ELEM] = esr;
        stencilElemOrFaceIndices[BoundaryStencil::Order::ELEM] = ei;
        stencilWeights[BoundaryStencil::Order::ELEM] = faceArea[iface] / globalSumFaceAreas[aquiferIndex];

        stencilRegionIndices[BoundaryStencil::Order::FACE] = -1;
        stencilSubRegionIndices[BoundaryStencil::Order::FACE] = -1;
        stencilElemOrFaceIndices[BoundaryStencil::Order::FACE] = iface;
        stencilWeights[BoundaryStencil::Order::FACE] = -faceArea[iface] / globalSumFaceAreas[aquiferIndex]; // likely unused for aquifers

        stencil.add( stencilRegionIndices.size(),
                     stencilRegionIndices.data(),
                     stencilSubRegionIndices.data(),
                     stencilElemOrFaceIndices.data(),
                     stencilWeights.data(),
                     iface );
      }
    }
  } );
}


REGISTER_CATALOG_ENTRY( FluxApproximationBase, TwoPointFluxApproximation, string const &, Group * const )

}
