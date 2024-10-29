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
 * @file SurfaceElementStencil.cpp
 */

#include "SurfaceElementStencil.hpp"

namespace geos
{

void SurfaceElementStencil::move( LvArray::MemorySpace const space )
{
  StencilBase::move( space );
  m_cellCenterToEdgeCenters.move( space, true );
}

void SurfaceElementStencil::add( localIndex const numPts,
                                 localIndex const * const elementRegionIndices,
                                 localIndex const * const elementSubRegionIndices,
                                 localIndex const * const elementIndices,
                                 real64 const * const weights,
                                 localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= maxStencilSize, "Maximum stencil size exceeded" );

  auto const iter = m_connectorIndices.find( connectorIndex );
  if( iter == m_connectorIndices.end() )
  {
    m_elementRegionIndices.appendArray( elementRegionIndices, elementRegionIndices + numPts );
    m_elementSubRegionIndices.appendArray( elementSubRegionIndices, elementSubRegionIndices + numPts );
    m_elementIndices.appendArray( elementIndices, elementIndices + numPts );
    m_weights.appendArray( weights, weights + numPts );

    m_connectorIndices[connectorIndex] = m_weights.size() - 1;
  }
  else
  {
    localIndex const stencilIndex = iter->second;
    m_elementRegionIndices.clearArray( stencilIndex );
    m_elementSubRegionIndices.clearArray( stencilIndex );
    m_elementIndices.clearArray( stencilIndex );
    m_weights.clearArray( stencilIndex );

    m_elementRegionIndices.appendToArray( stencilIndex, elementRegionIndices, elementRegionIndices + numPts );
    m_elementSubRegionIndices.appendToArray( stencilIndex, elementSubRegionIndices, elementSubRegionIndices + numPts );
    m_elementIndices.appendToArray( stencilIndex, elementIndices, elementIndices + numPts );
    m_weights.appendToArray( stencilIndex, weights, weights + numPts );
  }
}

void SurfaceElementStencil::add( localIndex const numPts,
                                 R1Tensor const * const cellCenterToEdgeCenter,
                                 localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= maxStencilSize, "Maximum stencil size exceeded" );

  auto const iter = m_connectorIndices.find( connectorIndex );
  if( iter == m_connectorIndices.end() )
  {
    GEOS_ERROR( "Wrong connectorIndex" );
  }
  else
  {
    localIndex const stencilIndex = iter->second;
    if( stencilIndex < m_cellCenterToEdgeCenters.size())
    {
      m_cellCenterToEdgeCenters.clearArray( stencilIndex );
      m_cellCenterToEdgeCenters.appendToArray( stencilIndex, cellCenterToEdgeCenter, cellCenterToEdgeCenter + numPts );
    }
    else
    {
      m_cellCenterToEdgeCenters.appendArray( cellCenterToEdgeCenter, cellCenterToEdgeCenter + numPts );
    }
  }
}

SurfaceElementStencil::KernelWrapper
SurfaceElementStencil::createKernelWrapper() const
{
  return { m_elementRegionIndices,
           m_elementSubRegionIndices,
           m_elementIndices,
           m_weights,
           m_cellCenterToEdgeCenters,
           m_meanPermCoefficient };
}

SurfaceElementStencilWrapper::
  SurfaceElementStencilWrapper( IndexContainerType const & elementRegionIndices,
                                IndexContainerType const & elementSubRegionIndices,
                                IndexContainerType const & elementIndices,
                                WeightContainerType const & weights,
                                ArrayOfArrays< R1Tensor > const & cellCenterToEdgeCenters,
                                real64 const meanPermCoefficient )

  : StencilWrapperBase( elementRegionIndices,
                        elementSubRegionIndices,
                        elementIndices,
                        weights ),
  m_cellCenterToEdgeCenters( cellCenterToEdgeCenters.toView() ),
  m_meanPermCoefficient( meanPermCoefficient )
{}

} /* namespace geos */
