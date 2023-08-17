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

/**
 * @file FaceElementSubRegion.cpp
 */

#include "FaceElementSubRegion.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{
using namespace dataRepository;


FaceElementSubRegion::FaceElementSubRegion( string const & name,
                                            dataRepository::Group * const parent ):
  SurfaceElementSubRegion( name, parent ),
  m_unmappedGlobalIndicesInToEdges(),
  m_unmappedGlobalIndicesInToFaces(),
  m_newFaceElements(),
  m_toFacesRelation()
{
  m_elementType = ElementType::Hexahedron;

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();

  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation ).
    setDescription( "Map to the faces attached to each FaceElement." ).
    reference().resize( 0, 2 );

  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString(), &m_edgesTo2dFaces ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of edge local indices to the fracture connector local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString(), &m_2dFaceToEdge ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices to edge local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString(), &m_2dFaceTo2dElems ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices face element local indices" ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::elem2dToCollocatedNodesString(), &m_2dElemToCollocatedNodes ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Dummy" ).
    setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::elem2dToCollocatedNodesBucketsString(), &m_2dElemToCollocatedNodesBuckets ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Dummy" ).
    setSizedFromParent( 1 );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  registerWrapper( viewKeyStruct::separationCoeffString(), &m_separationCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
    setDescription( "Scalar indicator of level of separation for a fracturing face." );
#endif

  excludeWrappersFromPacking( { viewKeyStruct::faceListString() } );

  m_2dElemToElems.resize( 0, 2 );

  m_numNodesPerElement = 8;
}

void FaceElementSubRegion::copyFromCellBlock( FaceBlockABC const & faceBlock )
{
  localIndex const num2dElements = faceBlock.num2dElements();
  resize( faceBlock.num2dElements() );

  m_toNodesRelation.base() = faceBlock.get2dElemToNodes();
  m_toEdgesRelation.base() = faceBlock.get2dElemToEdges();

  // `FaceBlockABC` is designed to be heterogeneous.
  // `FaceElementSubRegion` inherits from `ElementSubRegionBase` which is meant to be homogeneous.
  // But `FaceElementSubRegion` is sometimes used as an heterogeneous sub region,
  // which emphasizes the need of a refactoring.
  // In the meantime, we try to fill the face block into the sub region and hope for the best...
  {
    auto const hack = []( integer allSizes ) -> ElementType
    {
      GEOS_LOG_RANK( "All sizes : " << allSizes );
      switch( allSizes )
      {
        case 3:
        case 6:
          return ElementType::Wedge;
        case 4:
        case 8:
        case 0:
          return ElementType::Hexahedron;
        default:
          GEOS_ERROR( "Unsupported type of elements during the face element sub region creation." );
          return {};
      }
    };

    // Checking if all the 2d elements are homogeneous.
    // We rely on the number of nodes for each element to find out.
    std::vector< integer > sizes( num2dElements );
    for( int i = 0; i < num2dElements; ++i )
    {
      sizes[i] = m_toNodesRelation[i].size();
    }
    std::set< integer > const s( sizes.cbegin(), sizes.cend() );

    if( s.size() > 1 )
    {
      // If we have found that the input face block contains 2d elements of different types,
      // we inform the used that the situation may be at risk.
      // (We're storing the face block in a homogeneous container while it's actually heterogeneous).
      GEOS_WARNING( "Heterogeneous face element sub region found and stored as homogeneous. Use at your own risk." );
    }

    auto const it = std::max_element( s.cbegin(), s.cend() );
    integer const maxSize = *it;
    m_elementType = hack( maxSize );
    m_numNodesPerElement = maxSize;
  }

  // The `m_2dElemToElems` mappings involves element, sub regions and regions indices.
  // We store the element indices that are correct.
  // But we only have access to the cell block indices, not the sub regions indices.
  // Temporarily, and also because they share the same dimensions,
  // we store the cell block mapping at the sub region mapping location.
  // It will later be transformed into a sub regions mapping.  // Last, we fill the regions mapping with dummy -1 values that should all be replaced eventually.
  auto const elem2dToElems = faceBlock.get2dElemToElems();
  m_2dElemToElems.resize( num2dElements, 2 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( localIndex const & j: elem2dToElems.toCellIndex[i] )
    {
      m_2dElemToElems.m_toElementIndex.emplaceBack( i, j );
    }
    for( localIndex const & j: elem2dToElems.toBlockIndex[i] )
    {
      m_2dElemToElems.m_toElementSubRegion.emplaceBack( i, j );
    }
  }

  m_toFacesRelation.base() = faceBlock.get2dElemToFaces();

  m_2dFaceToEdge = faceBlock.get2dFaceToEdge();
  m_2dFaceTo2dElems = faceBlock.get2dFaceTo2dElems();

  m_localToGlobalMap = faceBlock.localToGlobalMap();
  this->constructGlobalToLocalMap();

  for( int i = 0; i < faceBlock.num2dFaces(); ++i )
  {
    m_recalculateConnectionsFor2dFaces.insert( i );
  }

  for( localIndex i = 0; i < faceBlock.num2dElements(); ++i )
  {
    m_newFaceElements.insert( i );
  }

  m_collocatedNodes = faceBlock.getCollocatedNodes();
  m_2dElemToCollocatedNodes = faceBlock.getCollocatedNodesOf2dElems();
  m_2dElemToCollocatedNodesBuckets = faceBlock.get2dElemsToCollocatedNodesBuckets();

  // TODO We still need to be able to import fields on the FaceElementSubRegion.
}

void FaceElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

void FaceElementSubRegion::calculateSingleElementGeometricQuantities( localIndex const k,
                                                                      arrayView1d< real64 const > const & faceArea )
{
  m_elementArea[k] = faceArea[ m_toFacesRelation[k][0] ];
  m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
}

void FaceElementSubRegion::calculateElementGeometricQuantities( NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                                                                FaceManager const & faceManager )
{
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();

  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateSingleElementGeometricQuantities( k, faceArea );
  } );
}



localIndex FaceElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex FaceElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}


// line 885
template< bool DO_PACKING, int USD >
localIndex Pack2( buffer_unit_type * & buffer,
                 arraySlice1d< globalIndex const, USD > const & var,
                 globalIndex const * const unmappedGlobalIndices,
                 localIndex const length )
{
  localIndex sizeOfPackedChars = bufferOps::Pack< DO_PACKING >( buffer, length );
  sizeOfPackedChars += length*sizeof(globalIndex);

  if( DO_PACKING )
  {
    globalIndex * const buffer_GI = reinterpret_cast< globalIndex * >(buffer);
    for( localIndex a=0; a<length; ++a )
    {
      if( var[a] != unmappedLocalIndexValue )
      {
//        GEOS_ERROR("HERE!");
//        buffer_GI[a] = localToGlobalMap[var[a]];
        buffer_GI[a] = var[a];
      }
      else
      {
        buffer_GI[a] = unmappedGlobalIndices[a];
      }
    }

    buffer += length * sizeof(globalIndex);
  }

  return sizeOfPackedChars;
}

// Copied from above
template< bool DO_PACKING, int USD >
localIndex Pack3( buffer_unit_type * & buffer,
                 arraySlice1d< array1d< globalIndex > const, USD > const & var,
                 globalIndex const * const unmappedGlobalIndices,
                 localIndex const length )
{
//  localIndex sizeOfPackedChars = 0;
  localIndex sizeOfPackedChars = bufferOps::Pack< DO_PACKING >( buffer, length );
  GEOS_ERROR_IF_NE_MSG(var.size(), length, "Internal error");
//  sizeOfPackedChars += length * sizeof( globalIndex );

//  if( DO_PACKING )
  {
//    globalIndex * const buffer_GI = reinterpret_cast< globalIndex * >(buffer);
    for( localIndex a = 0; a < length; ++a )
    {
      array1d< globalIndex > const & tmp = var[a];
      sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, tmp );
//      buffer += tmp.size() * sizeof( globalIndex );
//      if( var[a] != unmappedLocalIndexValue )
//      {
////        GEOS_ERROR("HERE!");
////        buffer_GI[a] = localToGlobalMap[var[a]];
//        buffer_GI[a] = var[a];
//      }
//      else
//      {
//        buffer_GI[a] = unmappedGlobalIndices[a];
//      }
    }

//    buffer += length * sizeof(globalIndex);
  }

  return sizeOfPackedChars;
}



// line 1269 in BufferOps_inline.hpp
//template< bool DO_PACKING, typename SORTED >
template< bool DO_PACKING >
localIndex
Pack2( buffer_unit_type *& buffer,
       ArrayOfArraysView< globalIndex const > const & var,
       map< localIndex, array1d< globalIndex > > const & unmappedGlobalIndices,
       arrayView1d< localIndex const > const & indices,
       arrayView1d< globalIndex const > const & localToGlobalMap )
//      arrayView1d< globalIndex const > const & localToGlobalMap,
//      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;
  array1d< globalIndex > junk;

  sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    auto iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > const & unmappedGI = iterUnmappedGI == unmappedGlobalIndices.end() ?
                                                junk :
                                                iterUnmappedGI->second;

//    sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, var[li] );
    sizeOfPackedChars += Pack2< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.data(),
                                             var.sizeOfArray( li ) );
//    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
//                                             var[li],
//                                             unmappedGI.data(),
//                                             var.sizeOfArray( li ),
//                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

// Copied from above
template< bool DO_PACKING >
localIndex
Pack3( buffer_unit_type *& buffer,
       ArrayOfArraysView< array1d< globalIndex > const > const & var,
       map< localIndex, array1d< globalIndex > > const & unmappedGlobalIndices,
       arrayView1d< localIndex const > const & indices,
       arrayView1d< globalIndex const > const & localToGlobalMap )
//      arrayView1d< globalIndex const > const & localToGlobalMap,
//      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap )
{
  localIndex sizeOfPackedChars=0;
  array1d< globalIndex > junk;

  sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, indices.size() );
  for( localIndex a = 0; a < indices.size(); ++a )
  {
    localIndex const li = indices[a];
    sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, localToGlobalMap[li] );

    auto iterUnmappedGI = unmappedGlobalIndices.find( li );

    array1d< globalIndex > const & unmappedGI = iterUnmappedGI == unmappedGlobalIndices.end() ?
                                                junk :
                                                iterUnmappedGI->second;

//    sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, var[li] );
    sizeOfPackedChars += Pack3< DO_PACKING >( buffer,
                                             var[li],
                                             unmappedGI.data(),
                                             var.sizeOfArray( li ) );
//    sizeOfPackedChars += Pack< DO_PACKING >( buffer,
//                                             var[li],
//                                             unmappedGI.data(),
//                                             var.sizeOfArray( li ),
//                                             relatedObjectLocalToGlobalMap );
  }

  return sizeOfPackedChars;
}

// copied from 977
inline
localIndex
Unpack2( buffer_unit_type const * & buffer,
        ArrayOfArrays< globalIndex > & var,
        localIndex const subArrayIndex,
        array1d< globalIndex > & unmappedGlobalIndices )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = bufferOps::Unpack( buffer, length );

  var.resizeArray( subArrayIndex, length );
  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices.setValues< serialPolicy >( unmappedLocalIndexValue );

  bool unpackedGlobalFlag = false;
  for( localIndex a = 0; a < length; ++a )
  {
    globalIndex unpackedGlobalIndex;
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, unpackedGlobalIndex );

    var( subArrayIndex, a ) = unpackedGlobalIndex;

////    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
//    auto iter = globalToLocalMap.find( unpackedGlobalIndex );
//    if( iter == globalToLocalMap.end() )
//    {
//      var( subArrayIndex, a ) = unmappedLocalIndexValue;
//      unmappedGlobalIndices[a] = unpackedGlobalIndex;
//      unpackedGlobalFlag = true;
//    }
//    else
//    {
//      var( subArrayIndex, a ) = iter->second;
//    }
  }
//  if( !unpackedGlobalFlag )
//  {
//    unmappedGlobalIndices.clear();
//  }

  return sizeOfUnpackedChars;
}

// Copied from above
inline
localIndex
Unpack3( buffer_unit_type const * & buffer,
        ArrayOfArrays< array1d< globalIndex > > & var,
        localIndex const subArrayIndex,
        array1d< globalIndex > & unmappedGlobalIndices )
{
  localIndex length;
  localIndex sizeOfUnpackedChars = bufferOps::Unpack( buffer, length );

  var.resizeArray( subArrayIndex, length );
  unmappedGlobalIndices.resize( length );
  unmappedGlobalIndices.setValues< serialPolicy >( unmappedLocalIndexValue );

  bool unpackedGlobalFlag = false;
  for( localIndex a = 0; a < length; ++a )
  {
//    globalIndex unpackedGlobalIndex;
//    sizeOfUnpackedChars += bufferOps::Unpack( buffer, unpackedGlobalIndex );
//
//    var( subArrayIndex, a ) = unpackedGlobalIndex;



//    array1d< globalIndex > tmp;
//    sizeOfUnpackedChars += bufferOps::Unpack( buffer, tmp );
//    var( subArrayIndex, a ) = tmp;
    array1d< globalIndex > & tmp = var( subArrayIndex, a );
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, tmp );

////    typename mapBase< globalIndex, localIndex, SORTED >::const_iterator
//    auto iter = globalToLocalMap.find( unpackedGlobalIndex );
//    if( iter == globalToLocalMap.end() )
//    {
//      var( subArrayIndex, a ) = unmappedLocalIndexValue;
//      unmappedGlobalIndices[a] = unpackedGlobalIndex;
//      unpackedGlobalFlag = true;
//    }
//    else
//    {
//      var( subArrayIndex, a ) = iter->second;
//    }
  }
//  if( !unpackedGlobalFlag )
//  {
//    unmappedGlobalIndices.clear();
//  }

  return sizeOfUnpackedChars;
}


// copied from BufferOps_inline 1304
inline
localIndex
Unpack2( buffer_unit_type const * & buffer,
        ArrayOfArrays< globalIndex > & var,
        array1d< localIndex > & indices,
        map< localIndex, array1d< globalIndex > > & unmappedGlobalIndices,
        unordered_map< globalIndex, localIndex > const & globalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = bufferOps::Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );
  array1d< globalIndex > unmappedIndices;

  for( localIndex a=0; a<indices.size(); ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
//      GEOS_ERROR("We should not be there");
      li = globalToLocalMap.at( gi );
    }

//    sizeOfUnpackedChars += Unpack2( buffer,
//                                    var,
//                                    li );

    unmappedIndices.resize( 0 );
    sizeOfUnpackedChars += Unpack2( buffer,
                                   var,
                                   li,
                                   unmappedIndices );

    if( unmappedIndices.size() > 0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }
  }
  return sizeOfUnpackedChars;
}

// Copied from above
inline
localIndex
Unpack3( buffer_unit_type const *& buffer,
         ArrayOfArrays< array1d< globalIndex > > & var,
         array1d< localIndex > & indices,
         map< localIndex, array1d< globalIndex > > & unmappedGlobalIndices,
         unordered_map< globalIndex, localIndex > const & globalToLocalMap )
{
  localIndex numIndicesUnpacked;
  localIndex const sizeOfIndicesPassedIn = indices.size();

  localIndex sizeOfUnpackedChars = bufferOps::Unpack( buffer, numIndicesUnpacked );

  GEOS_ERROR_IF( sizeOfIndicesPassedIn!=0 && numIndicesUnpacked!=indices.size(),
                 "number of unpacked indices("<<numIndicesUnpacked<<") does not equal size of "
                                                                    "indices passed into Unpack function("<<sizeOfIndicesPassedIn );

  indices.resize( numIndicesUnpacked );
  array1d< globalIndex > unmappedIndices;

  for( localIndex a=0; a<indices.size(); ++a )
  {
    globalIndex gi;
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, gi );

    localIndex & li = indices[a];
    if( sizeOfIndicesPassedIn > 0 )
    {
      GEOS_ERROR_IF( li!=globalToLocalMap.at( gi ),
                     "global index "<<gi<<" unpacked from buffer does not equal the lookup "
                                    <<li<<" for localIndex "<<li<<" on this rank" );
    }
    else
    {
//      GEOS_ERROR("We should not be there");
      li = globalToLocalMap.at( gi );
    }

//    sizeOfUnpackedChars += Unpack2( buffer,
//                                    var,
//                                    li );

    unmappedIndices.resize( 0 );
    sizeOfUnpackedChars += Unpack3( buffer,
                                   var,
                                   li,
                                   unmappedIndices );

    if( unmappedIndices.size() > 0 )
    {
      unmappedGlobalIndices[li] = unmappedIndices;
    }
  }
  return sizeOfUnpackedChars;
}



template< bool DO_PACKING >
localIndex FaceElementSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > const nodeLocalToGlobal = m_toNodesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const edgeLocalToGlobal = m_toEdgesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const faceLocalToGlobal = m_toFacesRelation.relatedObjectLocalToGlobal();

  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toNodesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToNodes,
                                               packList,
                                               localToGlobal,
                                               nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toEdgesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToEdges,
                                               packList,
                                               localToGlobal,
                                               edgeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toFacesRelation.toViewConst(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::surfaceElementsToCellRegionsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_2dElemToElems,
                                               packList,
                                               m_2dElemToElems.getElementRegionManager() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elem2dToCollocatedNodesString() ) );
  packedSize += Pack2< DO_PACKING >( buffer,
                                     m_2dElemToCollocatedNodes.toViewConst(),
                                     m_unmappedGlobalIndicesInToCollocatedNodes,
                                     packList,
                                     localToGlobal );

  GEOS_LOG_RANK( "before packing m_2dElemToCollocatedNodesBuckets = " << m_2dElemToCollocatedNodesBuckets );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elem2dToCollocatedNodesBucketsString() ) );
  packedSize += Pack3< DO_PACKING >( buffer,
                                     m_2dElemToCollocatedNodesBuckets.toViewConst(),
                                     m_unmappedGlobalIndicesInToCollocatedNodesBucket,
                                     packList,
                                     localToGlobal );
  GEOS_LOG_RANK( "after packing m_2dElemToCollocatedNodesBuckets" );

  GEOS_LOG_RANK("DONE pack");

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::collocatedNodesString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_collocatedNodes );
  GEOS_LOG_RANK("DONE pack2");

//  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "missingNodes" ) );
//  packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_missingNodes );

  return packedSize;
}



localIndex FaceElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const overwriteUpMaps,
                                                   bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );


  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::surfaceElementsToCellRegionsString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToElems,
                                     packList.toViewConst(),
                                     m_2dElemToElems.getElementRegionManager(),
                                     overwriteUpMaps );

//  GEOS_LOG_RANK( "m_2dElemToCollocatedNodes before = " << m_2dElemToCollocatedNodes );
  string elem2dToCollocatedNodesString;
  unPackedSize += bufferOps::Unpack( buffer, elem2dToCollocatedNodesString );
  GEOS_ERROR_IF_NE( elem2dToCollocatedNodesString, viewKeyStruct::elem2dToCollocatedNodesString() );
  unPackedSize += Unpack2( buffer,
                           m_2dElemToCollocatedNodes,
                           packList,
                           m_unmappedGlobalIndicesInToCollocatedNodes,
                           this->globalToLocalMap() );
//  GEOS_LOG_RANK( "m_2dElemToCollocatedNodes after = " << m_2dElemToCollocatedNodes );

  GEOS_LOG_RANK( "m_2dElemToCollocatedNodesBuckets before unpack = " << m_2dElemToCollocatedNodesBuckets );
  string elem2dToCollocatedNodesBucketsString;
  unPackedSize += bufferOps::Unpack( buffer, elem2dToCollocatedNodesBucketsString );
  GEOS_ERROR_IF_NE( elem2dToCollocatedNodesBucketsString, viewKeyStruct::elem2dToCollocatedNodesBucketsString() );
  unPackedSize += Unpack3( buffer,
                           m_2dElemToCollocatedNodesBuckets,
                           packList,
                           m_unmappedGlobalIndicesInToCollocatedNodesBucket,
                           this->globalToLocalMap() );
  GEOS_LOG_RANK( "m_2dElemToCollocatedNodesBuckets after unpack = " << m_2dElemToCollocatedNodesBuckets );

  GEOS_LOG_RANK("DONE unpack");

  string collocatedNodesString;
  unPackedSize += bufferOps::Unpack( buffer, collocatedNodesString );
  GEOS_ERROR_IF_NE( collocatedNodesString, viewKeyStruct::collocatedNodesString() );
  m_otherCollocatedNodes.push_back( ArrayOfArrays< globalIndex >{} );
  unPackedSize += bufferOps::Unpack( buffer, m_otherCollocatedNodes.back() );

  return unPackedSize;
}

void FaceElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}

void FaceElementSubRegion::fixSecondaryMappings( NodeManager const & nodeManager,
                                                 EdgeManager const & edgeManager,
                                                 FaceManager const & faceManager,
                                                 ElementRegionManager const & elemManager )
{
  // Here I can fix the other mappings which are not properly defined...
  localIndex const num2dElems = this->size();

  // For each collocated node, the `referenceCollocatedNodes` returns the lowest id among all the collocated nodes sharing the same position.
  // That way, it's possible to know if two nodes are collocated of each other by checking if they share the same lowest id.
  std::map< globalIndex, globalIndex > referenceCollocatedNodes;
  {
    std::set< std::set< globalIndex > > mergedCollocatedNodes;
    m_otherCollocatedNodes.push_back( m_collocatedNodes );
    for( ArrayOfArrays< globalIndex > const & dns: m_otherCollocatedNodes )
    {
      for( int i = 0; i < dns.size(); ++i )
      {
        std::set< globalIndex > tmp( dns[i].begin(), dns[i].end() );
        mergedCollocatedNodes.insert( tmp );
      }
    }

    for( std::set< globalIndex > const & collocatedNodes: mergedCollocatedNodes )
    {
      globalIndex const & ref = *std::min_element( collocatedNodes.cbegin(), collocatedNodes.cend() );
      for( globalIndex const & n: collocatedNodes )
      {
        referenceCollocatedNodes[n] = ref;
      }
    }
  }

  // The concept of "2d face" is deeply local to the `FaceElementSubRegion`.
  // Therefore, it's not exchanged during the ghosting process.
  // Furthermore, all the information is there to reconstruct the mappings related to the 2d faces.

  // To recreate the mappings, the first step is to associate a `2d face index` based on the global indices.
  // But let's be precise: some the global indices of the edges are different while the
  // Then we can resize all the mappings related to the `2d faces`.
  ArrayOfArraysView< localIndex const > const elem2dToEdges = m_toEdgesRelation.base().toViewConst();
  std::set< localIndex > edges;  // Will include twin edges.
  for( int i = 0; i < elem2dToEdges.size(); ++i )
  {
    for( localIndex const & j: elem2dToEdges[i] )
    {
      edges.insert( j );
    }
  }

  arrayView1d< globalIndex const > const nl2g = nodeManager.localToGlobalMap();

  auto const & edgeToNodes = edgeManager.nodeList();
  std::map< std::pair< globalIndex, globalIndex >, localIndex > edgesIds;
  for( localIndex const & edge: edges )
  {
    auto const nodes = edgeToNodes[edge];
    GEOS_ASSERT_EQ( nodes.size(), 2 );
    std::pair< globalIndex, globalIndex > const pg{ nl2g[nodes[0]], nl2g[nodes[1]] };
    edgesIds[pg] = edge;
  }

  // The `collocatedEdgeIds` map gathers all the collocated edges together.
  // The key of the map (`std::pair< globalIndex, globalIndex >`) represents the global indices of two nodes.
  // Those two nodes are the lowest index of collocated nodes. As such, those two nodes may not form a real edge.
  // But this trick lets us define some kind of hash that allows to compare the location of the edges:
  // edges sharing the same hash lie in the same position.
  std::map< std::pair< globalIndex, globalIndex >, std::set< localIndex > > collocatedEdgeBuckets;
  for( auto const & p: edgesIds )
  {
    std::pair< globalIndex, globalIndex > const & nodes = p.first;
    localIndex const & edge = p.second;

    auto it0 = referenceCollocatedNodes.find( nodes.first );
    globalIndex const n0 = it0 != referenceCollocatedNodes.cend() ? it0->second : nodes.first;

    auto it1 = referenceCollocatedNodes.find( nodes.second );
    globalIndex const n1 = it1 != referenceCollocatedNodes.cend() ? it1->second : nodes.second;

    std::pair< globalIndex, globalIndex > const edgeHash = std::minmax( n0, n1 );
    collocatedEdgeBuckets[edgeHash].insert( edge );
  }

  std::size_t const num2dFaces = collocatedEdgeBuckets.size();
  arrayView1d< integer const > edgeGhostRanks = edgeManager.ghostRank().toViewConst();
  m_2dFaceToEdge.clear();
  m_2dFaceToEdge.reserve( num2dFaces );
  auto cmp = [=]( int e,
                  int f ) -> bool
  {
    int const re = edgeGhostRanks[e] < 0 ? 0 : 1;
    int const rf = edgeGhostRanks[f] < 0 ? 0 : 1;
    return std::tie( re, e ) < std::tie( rf, f );
  };
  std::map< localIndex, localIndex > collocatedEdges;  // Mapping from the removed edge to its replacement.
  for( auto const & p: collocatedEdgeBuckets )
  {
    std::set< localIndex > const & input = p.second;
    localIndex const edge = *std::min_element( input.cbegin(), input.cend(), cmp );
    m_2dFaceToEdge.emplace_back( edge );
    for( localIndex const & collocatedEdge: input )  // Including `edge`
    {
      collocatedEdges[collocatedEdge] = edge;
    }
  }

  // `m_edgesTo2dFaces` is computed by the simple inversion of `m_2dFaceToEdge`
  m_edgesTo2dFaces.clear();
  for( std::size_t i = 0; i < num2dFaces; ++i )
  {
    m_edgesTo2dFaces[m_2dFaceToEdge[i]] = i;
  }

  m_newFaceElements.clear();
  m_newFaceElements.reserve( num2dElems );
  for( localIndex i = 0; i < num2dElems; ++i )
  {
    m_newFaceElements.insert( i );
  }
  m_recalculateConnectionsFor2dFaces.clear();
  m_recalculateConnectionsFor2dFaces.reserve( num2dFaces );
  for( std::size_t i = 0; i < num2dFaces; ++i )
  {
    m_recalculateConnectionsFor2dFaces.insert( i );
  }

  // Filling the 2d face to 2d elements mappings
  {
    // `tmp` contains the 2d face to 2d elements mappings as a `std` container.
    // Eventually, it's copied into an `LvArray` container.
    std::vector< std::vector< localIndex > > tmp( num2dFaces );
    for( auto i = 0; i < num2dElems; ++i )
    {
      for( auto const & e: elem2dToEdges[i] )
      {
        tmp[m_edgesTo2dFaces.at( collocatedEdges.at( e ) )].push_back( i );
      }
    }
    std::vector< localIndex > sizes;
    sizes.reserve( tmp.size() );
    for( std::vector< localIndex > const & t: tmp )
    {
      sizes.push_back( t.size() );
    }
    for( auto i = 0; i < m_2dFaceTo2dElems.size(); ++i )
    {
      m_2dFaceTo2dElems.clearArray( i );
    }
    m_2dFaceTo2dElems.resizeFromCapacities< serialPolicy >( sizes.size(), sizes.data() );

    for( std::size_t i = 0; i < tmp.size(); ++i )
    {
      for( std::size_t j = 0; j < tmp[i].size(); ++j )
      {
        m_2dFaceTo2dElems.emplaceBack( i, tmp[i][j] );
      }
    }
  }

  // When a fracture element has only one (or less!) neighbor, let's try to find the other one.
  // The `ElemPath` provides all the information of a given face: obviously its face index,
  // but also the element index that touch the face, and the node indices of the face as well.
  struct ElemPath
  {
    localIndex er;
    localIndex esr;
    localIndex ei;
    localIndex face;
    std::vector< localIndex > nodes;

    bool operator<( ElemPath const & other ) const
    {
      return std::tie( er, esr, ei, face, nodes ) < std::tie( other.er, other.esr, other.ei, other.face, other.nodes );
    }
  };

  // We are building the mapping that connect all the reference (duplicated) nodes of any face to the elements they are touching.
  // Using this nodal information will let use reconnect the the fracture 2d element to its 3d neighbor.
  std::map< std::set< globalIndex >, std::set< ElemPath > > faceRefNodesToElems;
  auto const buildFaceNodesToElems = [&]( localIndex const er,
                                          localIndex const esr,
                                          ElementRegionBase const & region,
                                          CellElementSubRegion const & subRegion )
  {
    auto const & elemToFaces = subRegion.faceList().base();
    auto const & faceToNodes = faceManager.nodeList();
    for( localIndex ei = 0; ei < elemToFaces.size( 0 ); ++ei )
    {
      for( auto const & face: elemToFaces[ei] )
      {
        // A set of the global indices of the nodes of the face is used as the "signature" of the face nodes.
        std::set< globalIndex > nodesOfFace;
        for( localIndex const & n: faceToNodes[face] )
        {
          auto const it = referenceCollocatedNodes.find( nl2g[n] );
          if( it != referenceCollocatedNodes.cend() )
          {
            nodesOfFace.insert( it->second );
          }
          else
          {
            // If all the nodes are not found, it makes no sense to keep on considering the current face as a candidate.
            // So we can exit the loop.
            break;
          }
        }
        // We still double check that all the nodes of the face were found.
        if( nodesOfFace.size() == LvArray::integerConversion< std::size_t >( faceToNodes[face].size() ) )
        {
          std::vector< localIndex > const nodes( faceToNodes[face].begin(), faceToNodes[face].end() );
          faceRefNodesToElems[nodesOfFace].insert( ElemPath{ er, esr, ei, face, nodes } );
        }
      }
    }
  };
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( buildFaceNodesToElems );

  // The `misMatches` array will contain fracture element that did not find all of its neighbors.
  // This is used to display a more precise error message.
  std::vector< localIndex > misMatches;
  // Here we loop over all the elements of the fracture.
  // When there's neighbor missing, we search for a face that would lie on the collocated nodes of the fracture element.
  for( int e2d = 0; e2d < num2dElems; ++e2d )
  {
    if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) >= 2 )  // All the neighbors are known.
    { continue; }

    std::set< globalIndex > refNodes;
    if( m_toNodesRelation[e2d].size() != 0 )
    {
      for( localIndex const & n: m_toNodesRelation[e2d] )
      {
        globalIndex const & gn = nl2g[n];
        auto const it = referenceCollocatedNodes.find( gn );
        if( it != referenceCollocatedNodes.cend() )
        {
          refNodes.insert( it->second );
        }
      }
    }
    else if( m_ghostRank[e2d] < 0 )
    {
      for( globalIndex const & gn: m_2dElemToCollocatedNodes[e2d] )
      {
        auto const it = referenceCollocatedNodes.find( gn );
        if( it != referenceCollocatedNodes.cend() )
        {
          refNodes.insert( it->second );
        }
      }
    }

    auto const match = faceRefNodesToElems.find( refNodes );
    if( match != faceRefNodesToElems.cend() )
    {
      for( ElemPath const & path: match->second )
      {
        // This `if` prevents from storing the same data twice.
        if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) == 0 || m_2dElemToElems.m_toElementIndex[e2d][0] != path.ei )
        {
          m_2dElemToElems.m_toElementRegion.emplaceBack( e2d, path.er );
          m_2dElemToElems.m_toElementSubRegion.emplaceBack( e2d, path.esr );
          m_2dElemToElems.m_toElementIndex.emplaceBack( e2d, path.ei );
          m_toFacesRelation.emplaceBack( e2d, path.face );
          for( localIndex const & n: path.nodes )
          {
            m_toNodesRelation.emplaceBack( e2d, n );
          }
        }
      }
    }
  }

  GEOS_ERROR_IF( !misMatches.empty(),
                 "Fracture " << this->getName() << " has elements {" << stringutilities::join( misMatches, ", " ) << "} without two neighbors." );

  // Checking that each face has two neighboring elements.
  // If not, we pop up an error.
  {
    std::vector< localIndex > isolatedFractureElements;
    for( int e2d = 0; e2d < num2dElems; ++e2d )
    {
      if( m_2dElemToElems.m_toElementIndex.sizeOfArray( e2d ) < 2 && m_ghostRank[e2d] < 0 )
      {
        isolatedFractureElements.push_back( e2d );
      }
    }
    GEOS_ERROR_IF( !isolatedFractureElements.empty(),
                   "Fracture " << this->getName() << " has elements {" << stringutilities::join( isolatedFractureElements, ", " ) << "} with less than two neighbors." );
  }

  // TODO clear all useless mappings!
}

void FaceElementSubRegion::inheritGhostRankFromParentFace( FaceManager const & faceManager,
                                                           std::set< localIndex > const & indices )
{
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  for( localIndex const & index: indices )
  {
    m_ghostRank[index] = faceGhostRank[m_toFacesRelation[index][0]];
  }
}

std::set< globalIndex > FaceElementSubRegion::getMissingNodes( unordered_map< globalIndex, localIndex > const & g2l ) const
{
  std::set< globalIndex > missingNodes;
  for( localIndex e2d = 0; e2d < m_2dElemToCollocatedNodesBuckets.size(); ++e2d )
  {
    for( integer ni = 0; ni < m_2dElemToCollocatedNodesBuckets[e2d].size(); ++ni )
    {
      for( globalIndex const gni: m_2dElemToCollocatedNodesBuckets( e2d, ni ) )
      {
        auto const it = g2l.find( gni );
        if( it == g2l.cend() )
        {
          missingNodes.insert( gni );
        }
      }
    }
  }

//  for( int i = 0; i < m_collocatedNodes.size(); ++i )
//  {
//    for( globalIndex const & n:  m_collocatedNodes[i] )
//    {
//      auto const it = g2l.find( n );
//      if( it == g2l.cend() )
//      {
//        missingNodes.insert( n );
//      }
//    }
//  }

  return missingNodes;
}

} /* namespace geos */
