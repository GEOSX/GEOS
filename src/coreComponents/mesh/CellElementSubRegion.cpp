/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "CellElementSubRegion.hpp"

#include "common/TypeDispatch.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/generators/CellBlockUtilities.hpp"

namespace geos
{
using namespace dataRepository;
using namespace constitutive;

CellElementSubRegion::CellElementSubRegion( string const & name, Group * const parent ):
  ElementSubRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation );
  registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation );

  registerWrapper( viewKeyStruct::constitutiveGroupingString(), &m_constitutiveGrouping ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::constitutivePointVolumeFractionString(), &m_constitutivePointVolumeFraction );

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();

  registerWrapper( viewKeyStruct::toEmbSurfString(), &m_toEmbeddedSurfaces ).setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::fracturedCellsString(), &m_fracturedCells ).setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::bubbleCellsString(), &m_bubbleCells ).setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::toFaceElementsString(), &m_toFaceElements ).setSizedFromParent( 0 );

  excludeWrappersFromPacking( { viewKeyStruct::nodeListString(),
                                viewKeyStruct::edgeListString(),
                                viewKeyStruct::faceListString(),
                                viewKeyStruct::fracturedCellsString(),
                                viewKeyStruct::toEmbSurfString() } );
}


void CellElementSubRegion::resizePerElementValues( localIndex const newNumNodesPerElement,
                                                   localIndex const newNumEdgesPerElement,
                                                   localIndex const newNumFacesPerElement )
{
  ElementSubRegionBase::resizePerElementValues( newNumNodesPerElement,
                                                newNumEdgesPerElement,
                                                newNumFacesPerElement );

  m_toNodesRelation.resize( size(), m_numNodesPerElement );
  m_toEdgesRelation.resize( size(), m_numEdgesPerElement );
  m_toFacesRelation.resize( size(), m_numFacesPerElement );
}


void CellElementSubRegion::copyFromCellBlock( CellBlockABC const & cellBlock )
{
  // Defines the (unique) element type of this cell element region,
  // and its associated number of nodes, edges, faces.
  m_elementType = cellBlock.getElementType();
  m_numNodesPerElement = cellBlock.numNodesPerElement();
  m_numEdgesPerElement = cellBlock.numEdgesPerElement();
  m_numFacesPerElement = cellBlock.numFacesPerElement();

  // We call the `resize` member function of the cell to (nodes, edges, faces) relations,
  // before calling the `CellElementSubRegion::resize` in order to keep the first dimension.
  // Be careful when refactoring.
  m_toNodesRelation.resize( this->size(), m_numNodesPerElement );
  m_toEdgesRelation.resize( this->size(), m_numEdgesPerElement );
  m_toFacesRelation.resize( this->size(), m_numFacesPerElement );
  this->resize( cellBlock.numElements() );

  this->nodeList() = cellBlock.getElemToNodes();
  this->edgeList() = cellBlock.getElemToEdges();
  this->faceList() = cellBlock.getElemToFaces();

  this->m_localToGlobalMap = cellBlock.localToGlobalMap();

  this->constructGlobalToLocalMap();
  cellBlock.forExternalProperties( [&]( WrapperBase const & wrapper )
  {
    types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;
      auto const src = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();
      this->registerWrapper( wrapper.getName(), std::make_unique< ArrayType >( &src ) );
    }, wrapper );
  } );
}

void CellElementSubRegion::addFracturedElement( localIndex const cellElemIndex,
                                                localIndex const embSurfIndex )
{
  // add the connection between the element and the embedded surface to the map
  m_toEmbeddedSurfaces.emplaceBack( cellElemIndex, embSurfIndex );
  // add the element to the fractured elements list
  m_fracturedCells.insert( cellElemIndex );
}


localIndex CellElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}


localIndex CellElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex CellElementSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > edgeLocalToGlobal = edgeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();


  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer,
                                                         nodeList().base().toViewConst(),
                                                         m_unmappedGlobalIndicesInNodelist,
                                                         packList,
                                                         localToGlobal,
                                                         nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               edgeList().base().toViewConst(),
                                               m_unmappedGlobalIndicesInEdgelist,
                                               packList,
                                               localToGlobal,
                                               edgeLocalToGlobal );


  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               faceList().base().toViewConst(),
                                               m_unmappedGlobalIndicesInFacelist,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  return packedSize;
}

localIndex CellElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const GEOS_UNUSED_PARAM( overwriteUpMaps ),
                                                   bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;
  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().base().toView(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->globalToLocalMap(),
                                     nodeList().relatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     edgeList().base(),
                                     packList,
                                     m_unmappedGlobalIndicesInEdgelist,
                                     this->globalToLocalMap(),
                                     edgeList().relatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     faceList().base(),
                                     packList,
                                     m_unmappedGlobalIndicesInFacelist,
                                     this->globalToLocalMap(),
                                     faceList().relatedObjectGlobalToLocal() );

  return unPackedSize;
}

localIndex CellElementSubRegion::packFracturedElementsSize( arrayView1d< localIndex const > const & packList,
                                                            arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  buffer_unit_type * junk = nullptr;
  return packFracturedElementsImpl< false >( junk, packList, embeddedSurfacesLocalToGlobal );
}


localIndex CellElementSubRegion::packFracturedElements( buffer_unit_type * & buffer,
                                                        arrayView1d< localIndex const > const & packList,
                                                        arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  return packFracturedElementsImpl< true >( buffer, packList, embeddedSurfacesLocalToGlobal );
}

template< bool DO_PACKING >
localIndex CellElementSubRegion::packFracturedElementsImpl( buffer_unit_type * & buffer,
                                                            arrayView1d< localIndex const > const & packList,
                                                            arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  localIndex packedSize = 0;

  // only here to use that packing function
  map< localIndex, array1d< globalIndex > > unmappedGlobalIndices;

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::toEmbSurfString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               embeddedSurfacesList().base().toViewConst(),
                                               unmappedGlobalIndices,
                                               packList,
                                               localToGlobal,
                                               embeddedSurfacesLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::fracturedCellsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_fracturedCells.toViewConst(),
                                               packList,
                                               localToGlobal );

  return packedSize;
}


localIndex CellElementSubRegion::unpackFracturedElements( buffer_unit_type const * & buffer,
                                                          localIndex_array & packList,
                                                          unordered_map< globalIndex, localIndex > const & embeddedSurfacesGlobalToLocal )
{
  localIndex unPackedSize = 0;

  string toEmbSurfString;
  unPackedSize += bufferOps::Unpack( buffer, toEmbSurfString );
  GEOS_ERROR_IF_NE( toEmbSurfString, viewKeyStruct::toEmbSurfString() );

  // only here to use that packing function
  map< localIndex, array1d< globalIndex > > unmappedGlobalIndices;

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEmbeddedSurfaces,
                                     packList,
                                     unmappedGlobalIndices,
                                     this->globalToLocalMap(),
                                     embeddedSurfacesGlobalToLocal );

  string fracturedCellsString;
  unPackedSize += bufferOps::Unpack( buffer, fracturedCellsString );
  GEOS_ERROR_IF_NE( fracturedCellsString, viewKeyStruct::fracturedCellsString() );

  SortedArray< globalIndex > junk;
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_fracturedCells,
                                     junk,
                                     this->globalToLocalMap(),
                                     false );

  return unPackedSize;
}



void CellElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( edgeList(),
                                    m_unmappedGlobalIndicesInEdgelist,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( faceList(),
                                    m_unmappedGlobalIndicesInFacelist,
                                    clearIfUnmapped );
}

localIndex CellElementSubRegion::getFaceNodes( localIndex const elementIndex,
                                               localIndex const localFaceIndex,
                                               Span< localIndex > const nodeIndices ) const
{
  return geos::getFaceNodes( m_elementType, elementIndex, localFaceIndex, m_toNodesRelation, nodeIndices );
}

void CellElementSubRegion::
  calculateElementCenterAndVolume( localIndex const k,
                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
{
  auto getElementCoordinatesaAndComputeElementCenter = [k, X, this]( auto & XLocal )
  {
    LvArray::tensorOps::fill< 3 >( m_elementCenter[k], 0 );
    for( localIndex a = 0; a < m_numNodesPerElement; ++a )
    {
      LvArray::tensorOps::copy< 3 >( XLocal[a], X[m_toNodesRelation( k, a )] );
      LvArray::tensorOps::add< 3 >( m_elementCenter[k], XLocal[a] );
    }
    LvArray::tensorOps::scale< 3 >( m_elementCenter[k], 1.0 / m_numNodesPerElement );
  };

  switch( m_elementType )
  {
    case ElementType::Hexahedron:
    {
      real64 Xlocal[8][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::hexahedronVolume( Xlocal );
      break;
    }
    case ElementType::Tetrahedron:
    {
      real64 Xlocal[4][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::tetrahedronVolume( Xlocal );
      break;
    }
    case ElementType::Wedge:
    {
      real64 Xlocal[6][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::wedgeVolume( Xlocal );
      break;
    }
    case ElementType::Pyramid:
    {
      real64 Xlocal[5][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::pyramidVolume( Xlocal );
      break;
    }
    case ElementType::Prism5:
    {
      real64 Xlocal[10][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::prismVolume< 5 >( Xlocal );
      break;
    }
    case ElementType::Prism6:
    {
      real64 Xlocal[12][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::prismVolume< 6 >( Xlocal );
      break;
    }
    case ElementType::Prism7:
    {
      real64 Xlocal[14][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::prismVolume< 7 >( Xlocal );
      break;
    }
    case ElementType::Prism8:
    {
      real64 Xlocal[16][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::prismVolume< 8 >( Xlocal );
      break;
    }
    case ElementType::Prism9:
    {
      real64 Xlocal[18][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::prismVolume< 9 >( Xlocal );
      break;
    }
    case ElementType::Prism10:
    {
      real64 Xlocal[20][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::prismVolume< 10 >( Xlocal );
      break;
    }
    case ElementType::Prism11:
    {
      real64 Xlocal[22][3];
      getElementCoordinatesaAndComputeElementCenter( Xlocal );
      m_elementVolume[k] = computationalGeometry::prismVolume< 11 >( Xlocal );
      break;
    }
    default:
    {
      GEOS_ERROR( GEOS_FMT( "Volume calculation not supported for element type {} in subregion {}",
                            m_elementType, getDataContext() ) );
    }
  }

  GEOS_ERROR_IF( m_elementVolume[k] <= 0.0,
                 GEOS_FMT( "Negative volume for element {} type {} in subregion {}",
                           k, m_elementType, getDataContext() ) );
}

void CellElementSubRegion::calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                                FaceManager const & GEOS_UNUSED_PARAM( faceManager ) )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();

  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateElementCenterAndVolume( k, X );
  } );
}

void CellElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementSubRegion, string const &, Group * const )

} /* namespace geos */
