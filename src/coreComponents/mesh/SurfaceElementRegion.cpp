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
 * @file SurfaceElementRegion.cpp
 */

#include "MeshFields.hpp"
#include "EdgeManager.hpp"
#include "SurfaceElementRegion.hpp"


namespace geos
{
using namespace dataRepository;

SurfaceElementRegion::SurfaceElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent ),
  m_subRegionType( SurfaceSubRegionType::faceElement )
{
  registerWrapper( viewKeyStruct::subRegionTypeString(), &m_subRegionType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_subRegionType ).
    setDescription( "Type of surface element subregion. Valid options: {" + EnumStrings< SurfaceSubRegionType >::concat( ", " ) + "}." );

  registerWrapper( viewKeyStruct::faceBlockString(), &m_faceBlockName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( "FractureSubRegion" ).
    setDescription( "The name of the face block in the mesh, or the embedded surface." );

  registerWrapper( viewKeyStruct::defaultApertureString(), &m_defaultAperture ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The default aperture of newly formed surface elements." );
}

SurfaceElementRegion::~SurfaceElementRegion()
{}


void SurfaceElementRegion::generateMesh( Group const & faceBlocks )
{
  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );

  if( m_subRegionType == SurfaceSubRegionType::embeddedElement )
  {
    elementSubRegions.registerGroup< EmbeddedSurfaceSubRegion >( m_faceBlockName );
  }
  else if( m_subRegionType == SurfaceSubRegionType::faceElement )
  {
    FaceElementSubRegion & subRegion = elementSubRegions.registerGroup< FaceElementSubRegion >( m_faceBlockName );
    if( faceBlocks.hasGroup( m_faceBlockName ) )
    {
      FaceBlockABC const & source = faceBlocks.getGroup< FaceBlockABC >( m_faceBlockName );
      subRegion.copyFromCellBlock( source );
    }
    else
    {
      GEOS_LOG_RANK_0( "No face block \"" << m_faceBlockName << "\" was found in the mesh. Empty surface region was created." );
    }
  }
}

void SurfaceElementRegion::initializePreSubGroups()
{
  GEOS_ERROR_IF_LE_MSG( m_defaultAperture, 0.0,
                        getWrapperDataContext( viewKeyStruct::defaultApertureString() ) <<
                        ": default aperture must be larger than 0.0" );

  this->forElementSubRegions< SurfaceElementSubRegion >( [&] ( SurfaceElementSubRegion & subRegion )
  {
    subRegion.getWrapper< array1d< real64 > >( fields::elementAperture::key() ).
      setApplyDefaultValue( m_defaultAperture );
  } );
}

localIndex SurfaceElementRegion::addToFractureMesh( real64 const time_np1,
                                                    FaceManager const * const faceManager,
                                                    ArrayOfArraysView< localIndex const >  const & originalFaceToEdgeMap,
                                                    localIndex const faceIndices[2] )
{
  localIndex rval = -1;

  SortedArray< localIndex > connectedEdges;

  arrayView2d< localIndex const > const faceToElementRegion = faceManager->elementRegionList();
  arrayView2d< localIndex const > const faceToElementSubRegion = faceManager->elementSubRegionList();
  arrayView2d< localIndex const > const faceToElementIndex = faceManager->elementList();

  FaceElementSubRegion & subRegion = this->getUniqueSubRegion< FaceElementSubRegion >();
  subRegion.resize( subRegion.size() + 1 );
  rval = subRegion.size() - 1;


  arrayView1d< real64 > const ruptureTime = subRegion.getField< fields::ruptureTime >();

  arrayView2d< real64 const > const faceCenter = faceManager->faceCenter();
  arrayView2d< real64 > const elemCenter = subRegion.getElementCenter();

  arrayView1d< integer > const subRegionGhostRank = subRegion.ghostRank();

  arrayView1d< integer const > const faceGhostRank = faceManager->ghostRank();

  SurfaceElementSubRegion::NodeMapType & nodeMap = subRegion.nodeList();
  SurfaceElementSubRegion::EdgeMapType & edgeMap = subRegion.edgeList();
  FaceElementSubRegion::FaceMapType & faceMap = subRegion.faceList();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager->nodeList().toViewConst();

  localIndex const kfe = subRegion.size() - 1;
  ruptureTime( kfe ) = time_np1;

  LvArray::tensorOps::copy< 3 >( elemCenter[ kfe ], faceCenter[ faceIndices[ 0 ] ] );

//  faceMap.resize( kfe, 2 );
  faceMap[kfe][0] = faceIndices[0];
  faceMap[kfe][1] = faceIndices[1];

  globalIndex const gi = faceManager->localToGlobalMap()[faceIndices[0]];
  subRegion.localToGlobalMap()[kfe] = gi;
  subRegion.updateGlobalToLocalMap( kfe );
  subRegionGhostRank[kfe] = faceGhostRank[faceIndices[0]];

  // Add the nodes that compose the new FaceElement to the nodeList
  localIndex const numNodesInFace0 = faceToNodeMap.sizeOfArray( faceIndices[ 0 ] );
  localIndex const numNodesInFace1 = faceToNodeMap.sizeOfArray( faceIndices[ 1 ] );

  //Temporarily set the map size 8 for both quadrangle and triangle faces. TODO: need to fix for arbitrary face sizes.
  nodeMap.resizeArray( kfe, 8 );

  for( localIndex a = 0; a < numNodesInFace0; ++a )
  {
    localIndex const aa = a < 2 ? a : numNodesInFace0 - a + 1;
    localIndex const bb = aa == 0 ? aa : numNodesInFace1 - aa;

    // TODO HACK need to generalize to something other than quads
    nodeMap[ kfe ][ a ] = faceToNodeMap( faceIndices[ 0 ], aa );
    nodeMap[ kfe ][ a + numNodesInFace0 ] = faceToNodeMap( faceIndices[ 1 ], bb );
  }

  if( numNodesInFace0 == 3 )
  {
    nodeMap[ kfe ][ 6 ] = faceToNodeMap( faceIndices[ 0 ], 2 );
    nodeMap[ kfe ][ 7 ] = faceToNodeMap( faceIndices[ 1 ], 2 );
  }

  // Add the edges that compose the faceElement to the edge map. This is essentially a copy of
  // the facesToEdges entry.
  localIndex const faceIndex = faceIndices[0];
  localIndex const numEdges = originalFaceToEdgeMap.sizeOfArray( faceIndex );
  edgeMap.resizeArray( kfe, numEdges );
  for( localIndex a=0; a<numEdges; ++a )
  {
    edgeMap[kfe][a] = originalFaceToEdgeMap( faceIndex, a );
    connectedEdges.insert( originalFaceToEdgeMap( faceIndex, a ) );
  }

  // Add the cell region/subregion/index to the faceElementToCells map
  FixedToManyElementRelation & faceElementsToCells = subRegion.getToCellRelation();

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    localIndex const er = faceToElementRegion[faceIndices[ke]][ke];
    localIndex const esr = faceToElementSubRegion[faceIndices[ke]][ke];
    localIndex const ei = faceToElementIndex[faceIndices[ke]][ke];

    faceElementsToCells.m_toElementRegion[kfe][ke] = faceToElementRegion[faceIndices[ke]][ke];
    faceElementsToCells.m_toElementSubRegion[kfe][ke] = faceToElementSubRegion[faceIndices[ke]][ke];
    faceElementsToCells.m_toElementIndex[kfe][ke] = faceToElementIndex[faceIndices[ke]][ke];
    
  }

  // Fill the connectivity between FaceElement entries. This is essentially a copy of the
  // edgesToFaces map, but with differing offsets.
  for( auto const & edge : connectedEdges )
  {
    // check to see if the edgesToFractureConnectors already have an entry
    if( subRegion.m_edgesTo2dFaces.count( edge )==0 )
    {
      // if not, then fill increase the size of the fractureConnectors to face elements map and
      // fill the fractureConnectorsToEdges map with the current edge....and the inverse map too.
      subRegion.m_2dFaceTo2dElems.appendArray( 0 );
      subRegion.m_2dFaceToEdge.emplace_back( edge );
      subRegion.m_edgesTo2dFaces[edge] = subRegion.m_2dFaceToEdge.size()-1;
    }
    // now fill the fractureConnectorsToFaceElements map. This is analogous to the edge to face map
    localIndex const connectorIndex = subRegion.m_edgesTo2dFaces[edge];
    localIndex const numCells = subRegion.m_2dFaceTo2dElems.sizeOfArray( connectorIndex ) + 1;
    subRegion.m_2dFaceTo2dElems.resizeArray( connectorIndex, numCells );
    subRegion.m_2dFaceTo2dElems[connectorIndex][ numCells-1 ] = kfe;

    // And fill the list of connectors that will need stencil modifications
    subRegion.m_recalculateConnectionsFor2dFaces.insert( connectorIndex );
  }

  subRegion.calculateSingleElementGeometricQuantities( kfe, faceManager->faceArea() );

  // update the sets
  for( auto const & setIter : faceManager->sets().wrappers() )
  {
    SortedArrayView< localIndex const > const & faceSet = faceManager->sets().getReference< SortedArray< localIndex > >( setIter.first );
    SortedArray< localIndex > & faceElementSet = subRegion.sets().registerWrapper< SortedArray< localIndex > >( setIter.first ).reference();
    for( localIndex a = 0; a < faceMap.size(0); ++a )
    {
      if( faceSet.count( faceMap[a][0] ) )
      {
        faceElementSet.insert( a );
      }
    }
  }

  return rval;
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, SurfaceElementRegion, string const &, Group * const )

} /* namespace geos */
