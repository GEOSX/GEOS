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
 * @file BoundaryStencil.cpp
 */

#include "BoundaryStencil.hpp"

namespace geos
{

BoundaryStencil::BoundaryStencil()
  : StencilBase()
{
  m_faceNormal.resize( 0, 3 );
  m_cellToFaceVec.resize( 0, 3 );
}

void BoundaryStencil::add( localIndex const numPts,
                           localIndex const * const elementRegionIndices,
                           localIndex const * const elementSubRegionIndices,
                           localIndex const * const elementIndices,
                           real64 const * const weights,
                           localIndex const connectorIndex )
{
  GEOS_ERROR_IF_NE_MSG( numPts, 2, "Number of points in boundary stencil should be 2" );

  localIndex const oldSize = m_elementRegionIndices.size( 0 );
  localIndex const newSize = oldSize + 1;
  m_elementRegionIndices.resize( newSize, numPts );
  m_elementSubRegionIndices.resize( newSize, numPts );
  m_elementIndices.resize( newSize, numPts );
  m_weights.resize( newSize, numPts );

  for( localIndex a = 0; a < numPts; ++a )
  {
    m_elementRegionIndices( oldSize, a ) = elementRegionIndices[a];
    m_elementSubRegionIndices( oldSize, a ) = elementSubRegionIndices[a];
    m_elementIndices( oldSize, a ) = elementIndices[a];
    m_weights( oldSize, a ) = weights[a];
  }
  m_connectorIndices[connectorIndex] = oldSize;
}

void BoundaryStencil::addVectors( real64 const & transMultiplier,
                                  real64 const (&faceNormal)[3],
                                  real64 const (&cellToFaceVec)[3] )
{
  localIndex const size = m_faceNormal.size( 0 );
  m_faceNormal.resize( size + 1 );
  m_cellToFaceVec.resize( size + 1 );
  m_weightMultiplier.emplace_back( transMultiplier );
  LvArray::tensorOps::copy< 3 >( m_faceNormal[size], faceNormal );
  LvArray::tensorOps::copy< 3 >( m_cellToFaceVec[size], cellToFaceVec );
}

BoundaryStencil::KernelWrapper BoundaryStencil::createKernelWrapper() const
{
  return { m_elementRegionIndices,
           m_elementSubRegionIndices,
           m_elementIndices,
           m_weights,
           m_faceNormal,
           m_cellToFaceVec,
           m_weightMultiplier };
}

BoundaryStencilWrapper::
  BoundaryStencilWrapper( IndexContainerType const & elementRegionIndices,
                          IndexContainerType const & elementSubRegionIndices,
                          IndexContainerType const & elementIndices,
                          WeightContainerType const & weights,
                          arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & faceNormal,
                          arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & cellToFaceVec,
                          arrayView1d< real64 > const & weightMultiplier )
  : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ),
  m_faceNormal( faceNormal ),
  m_cellToFaceVec( cellToFaceVec ),
  m_weightMultiplier( weightMultiplier )
{}

} // namespace geos
