/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


#include "BufferOps.hpp"
#include "dataRepository/Packing.hpp"
#include "ToElementRelation.hpp"
#include "ElementRegionManager.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{
namespace bufferOps
{


template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 OrderedVariableToManyElementRelation const & var,
                 arrayView1d<localIndex> const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, packList.size() );
  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, var.m_toElementRegion[index].size() );
    for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
    {
      localIndex elemRegionIndex             = var.m_toElementRegion[index][b];
      ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

      localIndex elemSubRegionIndex                  = var.m_toElementSubRegion[index][b];
      CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

      localIndex elemIndex = var.m_toElementIndex[index][b];

      sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, elemRegionIndex );
      sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, elemSubRegionIndex );
      sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, elemSubRegion->m_localToGlobalMap[elemIndex] );

    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack<true>( char*&,
                                OrderedVariableToManyElementRelation const &,
                                arrayView1d<localIndex> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( char*&,
                                 OrderedVariableToManyElementRelation const &,
                                 arrayView1d<localIndex> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( char const * & buffer,
                   OrderedVariableToManyElementRelation & var,
                   arrayView1d<localIndex> const & packList,
                   ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Packing::Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( numIndicesUnpacked != packList.size(), "");

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfUnpackedChars += Packing::Unpack( buffer, numIndicesUnpacked );
    var.m_toElementRegion[index].resize( numIndicesUnpacked );
    var.m_toElementSubRegion[index].resize( numIndicesUnpacked );
    var.m_toElementIndex[index].resize( numIndicesUnpacked );
//    GEOS_ERROR_IF( numIndicesUnpacked != var.m_toElementRegion[index].size(), "")

    for( localIndex b=0 ; b<var.m_toElementRegion[index].size() ; ++b )
    {
      localIndex elemRegionIndex             = var.m_toElementRegion[index][b];
      ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

      localIndex elemSubRegionIndex                  = var.m_toElementSubRegion[index][b];
      CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);


      sizeOfUnpackedChars += Packing::Unpack( buffer, var.m_toElementRegion[index][b] );
      sizeOfUnpackedChars += Packing::Unpack( buffer, var.m_toElementSubRegion[index][b] );

      globalIndex globalElementIndex;
      sizeOfUnpackedChars += Packing::Unpack( buffer, globalElementIndex );

      var.m_toElementIndex[index][b] = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                      globalElementIndex,
                                                      localIndex(-1) );
//      var.m_toElementIndex[index][b] = elemSubRegion->m_globalToLocalMap.at(globalElementIndex);

    }
  ;}

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING >
localIndex Pack( char*& buffer,
                 FixedToManyElementRelation const & var,
                 arrayView1d<localIndex> const & packList,
                 ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, packList.size() );
  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, var.m_toElementRegion.size(1) );
    for( localIndex b=0 ; b<var.m_toElementRegion.size(1) ; ++b )
    {
      localIndex elemRegionIndex             = var.m_toElementRegion[index][b];

      if( elemRegionIndex == -1 )
      {
        sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, localIndex(-1) );
        sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, localIndex(-1) );
        sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, localIndex(-1) );
      }
      else
      {
        ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

        localIndex elemSubRegionIndex                  = var.m_toElementSubRegion[index][b];
        CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

        localIndex elemIndex = var.m_toElementIndex[index][b];

        sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, elemRegionIndex );
        sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, elemSubRegionIndex );
        sizeOfPackedChars += Packing::Pack<DO_PACKING>( buffer, elemSubRegion->m_localToGlobalMap[elemIndex] );
      }
    }
  }

  return sizeOfPackedChars;
}
template localIndex Pack<true>( char*&,
                                FixedToManyElementRelation const &,
                                arrayView1d<localIndex> const &,
                                ElementRegionManager const * const );
template localIndex Pack<false>( char*&,
                                 FixedToManyElementRelation const &,
                                 arrayView1d<localIndex> const &,
                                 ElementRegionManager const * const );


localIndex Unpack( char const * & buffer,
                   FixedToManyElementRelation & var,
                   arrayView1d<localIndex> const & packList,
                   ElementRegionManager const * const elementRegionManager )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex numIndicesUnpacked;
  sizeOfUnpackedChars += Packing::Unpack( buffer, numIndicesUnpacked );
  GEOS_ERROR_IF( numIndicesUnpacked != packList.size(), "");

  for( localIndex a=0 ; a<packList.size() ; ++a )
  {
    localIndex index = packList[a];
    sizeOfUnpackedChars += Packing::Unpack( buffer, numIndicesUnpacked );
    GEOS_ERROR_IF( numIndicesUnpacked != var.m_toElementRegion.size(1), "");

    for( localIndex b=0 ; b<var.m_toElementRegion.size(1) ; ++b )
    {
      localIndex & elemRegionIndex = var.m_toElementRegion[index][b];
      sizeOfUnpackedChars += Packing::Unpack( buffer, elemRegionIndex );

      if( elemRegionIndex==-1 )
      {
        sizeOfUnpackedChars += Packing::Unpack( buffer, var.m_toElementSubRegion[index][b] );
        sizeOfUnpackedChars += Packing::Unpack( buffer, var.m_toElementIndex[index][b] );
      }
      else
      {
        ElementRegion const * const elemRegion = elementRegionManager->GetRegion(elemRegionIndex);

        localIndex & elemSubRegionIndex = var.m_toElementSubRegion[index][b];
        sizeOfUnpackedChars += Packing::Unpack( buffer, elemSubRegionIndex );

        CellBlockSubRegion const * const elemSubRegion = elemRegion->GetSubRegion(elemSubRegionIndex);

        globalIndex globalElementIndex;
        sizeOfUnpackedChars += Packing::Unpack( buffer, globalElementIndex );
        var.m_toElementIndex[index][b] = softMapLookup( elemSubRegion->m_globalToLocalMap,
                                                        globalElementIndex,
                                                        localIndex(-1) );
      }
    }
  }

  return sizeOfUnpackedChars;
}

}
}
