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
 * @file EmbeddedSurfaceToCellStencil.cpp
 */

#include "EmbeddedSurfaceToCellStencil.hpp"

namespace geos
{

void EmbeddedSurfaceToCellStencil::move( LvArray::MemorySpace const space )
{
  StencilBase::move( space );
}

void EmbeddedSurfaceToCellStencil::add( localIndex const numPts,
                                        localIndex const * const elementRegionIndices,
                                        localIndex const * const elementSubRegionIndices,
                                        localIndex const * const elementIndices,
                                        real64 const * const weights,
                                        localIndex const connectorIndex )
{
  GEOS_ERROR_IF_NE_MSG( numPts, 2, "Number of cells in TPFA stencil should be 2" );

  localIndex const oldSize = m_elementRegionIndices.size( 0 );
  localIndex const newSize = oldSize + 1;
  m_elementRegionIndices.resize( newSize, numPts );
  m_elementSubRegionIndices.resize( newSize, numPts );
  m_elementIndices.resize( newSize, numPts );
  m_weights.resize( newSize, numPts );

  for( localIndex a=0; a<numPts; ++a )
  {
    m_elementRegionIndices( oldSize, a ) = elementRegionIndices[a];
    m_elementSubRegionIndices( oldSize, a ) = elementSubRegionIndices[a];
    m_elementIndices( oldSize, a ) = elementIndices[a];
    m_weights( oldSize, a ) = weights[a];
  }
  m_connectorIndices[connectorIndex] = oldSize;
}

EmbeddedSurfaceToCellStencil::KernelWrapper
EmbeddedSurfaceToCellStencil::createKernelWrapper() const
{
  return { m_elementRegionIndices,
           m_elementSubRegionIndices,
           m_elementIndices,
           m_weights };
}

EmbeddedSurfaceToCellStencilWrapper::
  EmbeddedSurfaceToCellStencilWrapper( IndexContainerType const & elementRegionIndices,
                                       IndexContainerType const & elementSubRegionIndices,
                                       IndexContainerType const & elementIndices,
                                       WeightContainerType const & weights )
  : StencilWrapperBase( elementRegionIndices,
                        elementSubRegionIndices,
                        elementIndices,
                        weights )
{}

} /* namespace geos */
