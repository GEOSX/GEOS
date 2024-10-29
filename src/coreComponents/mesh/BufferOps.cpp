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


#include "BufferOps.hpp"
#include "dataRepository/BufferOps.hpp"
#include "ToElementRelation.hpp"
#include "codingUtilities/Utilities.hpp"
#include "ElementRegionManager.hpp"

namespace geos
{
namespace bufferOps
{


template< bool DO_PACKING >
localIndex Pack( buffer_unit_type * & buffer,
                 OrderedVariableToManyElementRelation const & var,
                 arrayView1d< localIndex const > const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, packList.size() );
  for( localIndex a=0; a<packList.size(); ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, var.m_toElementRegion.sizeOfArray( index ) );
    for( localIndex b=0; b<var.m_toElementRegion.sizeOfArray( index ); ++b )
    {
      localIndex elemRegionIndex    = var.m_toElementRegion[index][b];
      localIndex elemSubRegionIndex = var.m_toElementSubRegion[index][b];
      localIndex elemIndex          = var.m_toElementIndex[index][b];

      sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, elemRegionIndex );
      sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, elemSubRegionIndex );

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 && elemIndex!=-1 )
      {
        ElementRegionBase const & elemRegion = elementRegionManager->getRegion( elemRegionIndex );
        ElementSubRegionBase const & elemSubRegion = elemRegion.getSubRegion( elemSubRegionIndex );
        sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, elemSubRegion.localToGlobalMap()[elemIndex] );
      }
      else
      {
        sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, globalIndex( elemIndex ) );
      }
    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack< true >( buffer_unit_type * &,
                                  OrderedVariableToManyElementRelation const &,
                                  arrayView1d< localIndex const > const &,
                                  ElementRegionManager const * const );
template localIndex Pack< false >( buffer_unit_type * &,
                                   OrderedVariableToManyElementRelation const &,
                                   arrayView1d< localIndex const > const &,
                                   ElementRegionManager const * const );


localIndex Unpack( buffer_unit_type const * & buffer,
                   OrderedVariableToManyElementRelation & var,
                   arrayView1d< localIndex const > const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( numIndicesUnpacked != packList.size(), "" );

  using ElementID = std::array< localIndex, 3 >;

  // Allocate some memory to store map entries that don't fit in existing capacity
  array1d< localIndex > indicesToReplace;
  ArrayOfArrays< ElementID > valuesToReplace;
  indicesToReplace.reserve( numIndicesUnpacked );
  valuesToReplace.reserve( numIndicesUnpacked );
  valuesToReplace.reserveValues( numIndicesUnpacked * 12 ); // guesstimate

  std::vector< ElementID > values;
  for( localIndex a=0; a<packList.size(); ++a )
  {
    values.clear();

    localIndex const index = packList[a];
    sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );

    if( !clearFlag )
    {
      for( localIndex b=0; b<var.m_toElementRegion.sizeOfArray( index ); ++b )
      {
        values.push_back( { var.m_toElementRegion[index][b],
                            var.m_toElementSubRegion[index][b],
                            var.m_toElementIndex[index][b] } );
      }
    }

    for( localIndex b=0; b<numIndicesUnpacked; ++b )
    {

      localIndex elemRegionIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, elemRegionIndex );

      localIndex elemSubRegionIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, elemSubRegionIndex );

      globalIndex globalElementIndex;
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 )
      {
        ElementRegionBase const & elemRegion = elementRegionManager->getRegion( elemRegionIndex );
        ElementSubRegionBase const & elemSubRegion = elemRegion.getSubRegion( elemSubRegionIndex );

        localIndex localElementIndex = softMapLookup( elemSubRegion.globalToLocalMap(),
                                                      globalElementIndex,
                                                      localIndex( -1 ) );
        if( localElementIndex!=-1 )
        {
          values.push_back( { elemRegionIndex,
                              elemSubRegionIndex,
                              localElementIndex } );
        }
      }
      else
      {
        values.push_back( { -1, -1, -1 } );
      }
    }

    localIndex const numUniqueValues = LvArray::sortedArrayManipulation::makeSortedUnique( values.begin(), values.end() );

    if( numUniqueValues <= var.m_toElementIndex.capacityOfArray( index ) )
    {
      var.m_toElementRegion.resizeArray( index, numUniqueValues );
      var.m_toElementSubRegion.resizeArray( index, numUniqueValues );
      var.m_toElementIndex.resizeArray( index, numUniqueValues );

      for( localIndex b = 0; b < numUniqueValues; ++b )
      {
        var.m_toElementRegion( index, b ) = values[b][0];
        var.m_toElementSubRegion( index, b ) = values[b][1];
        var.m_toElementIndex( index, b ) = values[b][2];
      }
    }
    else
    {
      localIndex const k = indicesToReplace.size();
      indicesToReplace.emplace_back( index );
      valuesToReplace.appendArray( numUniqueValues );
      std::copy( values.begin(), values.end(), valuesToReplace[k].begin() );
    }
  }

  // If there were element lists that didn't fit in the map, rebuild the whole thing
  if( !indicesToReplace.empty() )
  {
    ArrayOfArraysView< localIndex const > const toRegion = var.m_toElementRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const toSubRegion = var.m_toElementSubRegion.toViewConst();
    ArrayOfArraysView< localIndex const > const toIndex = var.m_toElementIndex.toViewConst();
    localIndex const numEntries = toRegion.size();

    // Copy old capacities (no way to direct copy, must kernel launch)
    array1d< localIndex > newCapacities( numEntries );
    forAll< parallelHostPolicy >( numEntries, [toRegion, newCapacities = newCapacities.toView()]( localIndex const k )
    {
      newCapacities[k] = toRegion.capacityOfArray( k );
    } );

    // Replace with new capacities where needed
    for( localIndex i = 0; i < indicesToReplace.size(); ++i )
    {
      newCapacities[indicesToReplace[i]] = valuesToReplace.sizeOfArray( i );
    }

    // Scan to generate new offsets
    array1d< localIndex > newOffsets( newCapacities.size() + 1 );
    newOffsets[ 0 ] = 0;
    RAJA::inclusive_scan< parallelHostPolicy >( RAJA::make_span( newCapacities.data(), numEntries ),
                                                RAJA::make_span( newOffsets.data() + 1, numEntries ) );

    // Allocate new maps
    ArrayOfArrays< localIndex > toRegionNew, toSubRegionNew, toIndexNew;
    toRegionNew.resizeFromOffsets( numEntries, newOffsets.data() );
    toSubRegionNew.resizeFromOffsets( numEntries, newOffsets.data() );
    toIndexNew.resizeFromOffsets( numEntries, newOffsets.data() );

    // Fill new maps with old values
    forAll< parallelHostPolicy >( numEntries, [toRegion, toSubRegion, toIndex,
                                               toRegionNew = toRegionNew.toView(),
                                               toSubRegionNew = toSubRegionNew.toView(),
                                               toIndexNew = toIndexNew.toView() ]( localIndex const k )
    {
      toRegionNew.appendToArray( k, toRegion[k].begin(), toRegion[k].end() );
      toSubRegionNew.appendToArray( k, toSubRegion[k].begin(), toSubRegion[k].end() );
      toIndexNew.appendToArray( k, toIndex[k].begin(), toIndex[k].end() );
    } );

    // Replace with new values
    for( localIndex i = 0; i < indicesToReplace.size(); ++i )
    {
      localIndex const k = indicesToReplace[i];
      localIndex const numValues = valuesToReplace.sizeOfArray( i );
      toRegionNew.resizeArray( k, numValues );
      toSubRegionNew.resizeArray( k, numValues );
      toIndexNew.resizeArray( k, numValues );
      for( localIndex j = 0; j < numValues; ++j )
      {
        toRegionNew[k][j] = valuesToReplace[i][j][0];
        toSubRegionNew[k][j] = valuesToReplace[i][j][1];
        toIndexNew[k][j] = valuesToReplace[i][j][2];
      }
    }

    // Move into place
    var.m_toElementRegion = std::move( toRegionNew );
    var.m_toElementSubRegion = std::move( toSubRegionNew );
    var.m_toElementIndex = std::move( toIndexNew );
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex Pack( buffer_unit_type * & buffer,
                 FixedToManyElementRelation const & var,
                 arrayView1d< localIndex const > const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, packList.size() );
  for( localIndex a=0; a<packList.size(); ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, var.m_toElementRegion.size( 1 ) );

    for( localIndex b=0; b<var.m_toElementRegion.size( 1 ); ++b )
    {
      localIndex elemRegionIndex    = var.m_toElementRegion[index][b];
      localIndex elemSubRegionIndex = var.m_toElementSubRegion[index][b];
      localIndex elemIndex          = var.m_toElementIndex[index][b];

      sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, elemRegionIndex );
      sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, elemSubRegionIndex );

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 && elemIndex!=-1 )
      {
        ElementRegionBase const & elemRegion = elementRegionManager->getRegion( elemRegionIndex );
        ElementSubRegionBase const & elemSubRegion = elemRegion.getSubRegion( elemSubRegionIndex );
        sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, elemSubRegion.localToGlobalMap()[elemIndex] );
      }
      else
      {
        sizeOfPackedChars += bufferOps::Pack< DO_PACKING >( buffer, globalIndex( elemIndex ) );
      }
    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack< true >( buffer_unit_type * &,
                                  FixedToManyElementRelation const &,
                                  arrayView1d< localIndex const > const &,
                                  ElementRegionManager const * const );
template localIndex Pack< false >( buffer_unit_type * &,
                                   FixedToManyElementRelation const &,
                                   arrayView1d< localIndex const > const &,
                                   ElementRegionManager const * const );


localIndex Unpack( buffer_unit_type const * & buffer,
                   FixedToManyElementRelation & var,
                   arrayView1d< localIndex const > const & packList,
                   ElementRegionManager const * const elementRegionManager,
                   bool const clearFlag )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += bufferOps::Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( numIndicesUnpacked != packList.size(), "" );

  for( localIndex a=0; a<packList.size(); ++a )
  {
    localIndex index = packList[a];
    localIndex numSubIndicesUnpacked;

    sizeOfUnpackedChars += bufferOps::Unpack( buffer, numSubIndicesUnpacked );
    GEOS_ERROR_IF( numSubIndicesUnpacked != var.m_toElementRegion.size( 1 ), "" );

    for( localIndex b=0; b<numSubIndicesUnpacked; ++b )
    {
      localIndex recvElemRegionIndex;
      localIndex recvElemSubRegionIndex;
      globalIndex globalElementIndex;

      sizeOfUnpackedChars += bufferOps::Unpack( buffer, recvElemRegionIndex );
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, recvElemSubRegionIndex );
      sizeOfUnpackedChars += bufferOps::Unpack( buffer, globalElementIndex );

      if( recvElemRegionIndex!=-1 && recvElemSubRegionIndex!=-1 && globalElementIndex!=-1 )
      {
        ElementRegionBase const & elemRegion = elementRegionManager->getRegion( recvElemRegionIndex );

        ElementSubRegionBase const & elemSubRegion = elemRegion.getSubRegion( recvElemSubRegionIndex );

        localIndex const recvElemIndex = softMapLookup( elemSubRegion.globalToLocalMap(),
                                                        globalElementIndex,
                                                        localIndex( -1 ) );

        for( localIndex c=0; c<var.m_toElementRegion.size( 1 ); ++c )
        {
          localIndex & elemRegionIndex = var.m_toElementRegion[index][c];
          localIndex & elemSubRegionIndex = var.m_toElementSubRegion[index][c];
          localIndex & elemIndex = var.m_toElementIndex[index][c];
          if( ( elemRegionIndex==recvElemRegionIndex &&
                elemSubRegionIndex==recvElemSubRegionIndex &&
                elemIndex==recvElemIndex ) )
          {
            break;
          }
          else if( ( elemRegionIndex==-1 || elemSubRegionIndex==-1 || elemIndex==-1 ) )
          {
            elemRegionIndex = recvElemRegionIndex;
            elemSubRegionIndex = recvElemSubRegionIndex;
            elemIndex = recvElemIndex;
            break;
          }
          else
          {
            //TODO need a better criteria and an error check here
          }
        }
      }
      else if( clearFlag )
      {
        var.m_toElementRegion[index][b] = -1;
        var.m_toElementSubRegion[index][b] = -1;
        var.m_toElementIndex[index][b] = -1;
      }
    }
  }

  return sizeOfUnpackedChars;
}

}
}
