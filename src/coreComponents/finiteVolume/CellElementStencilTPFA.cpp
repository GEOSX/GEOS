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
 * @file CellElementStencilTPFA.cpp
 */


#include "CellElementStencilTPFA.hpp"

namespace geos
{

CellElementStencilTPFA::CellElementStencilTPFA()
  : StencilBase()
{
  m_faceNormal.resize( 0, 3 );
  m_cellToFaceVec.resize( 0, 2, 3 );
}

void CellElementStencilTPFA::reserve( localIndex const size )
{
  StencilBase::reserve( size );

  m_faceNormal.reserve( 3 * size );
  m_cellToFaceVec.reserve( 6 * size );
  m_transMultiplier.reserve( size );
  m_geometricStabilizationCoef.reserve( size );
}

void CellElementStencilTPFA::add( localIndex const numPts,
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

void CellElementStencilTPFA::addVectors( real64 const & transMultiplier,
                                         real64 const & geometricStabilizationCoef,
                                         real64 const (&faceNormal)[3],
                                         real64 const (&cellToFaceVec)[2][3] )
{
  localIndex const oldSize = m_faceNormal.size( 0 );
  localIndex const newSize = oldSize + 1;
  m_faceNormal.resize( newSize );
  m_cellToFaceVec.resize( newSize );
  m_transMultiplier.resize( newSize );
  m_geometricStabilizationCoef.resize( newSize );

  m_transMultiplier[oldSize] = transMultiplier;
  m_geometricStabilizationCoef[oldSize] = geometricStabilizationCoef;

  LvArray::tensorOps::copy< 3 >( m_faceNormal[oldSize], faceNormal );
  for( localIndex a=0; a<2; a++ )
  {
    LvArray::tensorOps::copy< 3 >( m_cellToFaceVec[oldSize][a], cellToFaceVec[a] );
  }
}

CellElementStencilTPFA::KernelWrapper
CellElementStencilTPFA::createKernelWrapper() const
{
  return { m_elementRegionIndices,
           m_elementSubRegionIndices,
           m_elementIndices,
           m_weights,
           m_faceNormal,
           m_cellToFaceVec,
           m_transMultiplier,
           m_geometricStabilizationCoef };
}

CellElementStencilTPFAWrapper::
  CellElementStencilTPFAWrapper( IndexContainerType const & elementRegionIndices,
                                 IndexContainerType const & elementSubRegionIndices,
                                 IndexContainerType const & elementIndices,
                                 WeightContainerType const & weights,
                                 arrayView2d< real64 > const & faceNormal,
                                 arrayView3d< real64 > const & cellToFaceVec,
                                 arrayView1d< real64 > const & transMultiplier,
                                 arrayView1d< real64 > const & geometricStabilizationCoef )
  : StencilWrapperBase( elementRegionIndices,
                        elementSubRegionIndices,
                        elementIndices,
                        weights ),
  m_faceNormal( faceNormal ),
  m_cellToFaceVec( cellToFaceVec ),
  m_transMultiplier( transMultiplier ),
  m_geometricStabilizationCoef( geometricStabilizationCoef )
{}

void
StencilUtils::computeVelocity( const geos::CellElementStencilTPFAWrapper & stencil,
                               localIndex iconn, localIndex ip,
                               const real64 (&phaseFlux),
                               arraySlice1d< real64 const > const (&cellCartDim)[2],
                               localIndex const (&ghostRank)[2],
                               ElementRegionManager::ElementView< arrayView3d< real64 > > const & phaseVelocity )
{

  real64 surface[2];

  for( localIndex i = 0; i < 2; i++ )
  {
    localIndex const er = stencil.m_elementRegionIndices[iconn][i];
    localIndex const esr = stencil.m_elementSubRegionIndices[iconn][i];
    localIndex const ei = stencil.m_elementIndices[iconn][i];

    if( ghostRank[i] < 0 )
    {
      //halfWeight is d(Cell,Face)/area, we want area
      real64 const halfWeight = stencil.m_weights[iconn][i];
      real64 c2fVec[3];
      LvArray::tensorOps::copy< 3 >( c2fVec, stencil.m_cellToFaceVec[iconn][i] );
      real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( c2fVec );

      surface[i] = c2fDistance * halfWeight;

      real64 faceNormal[3], invDist[3], velocityNorm[3], phaseVel[3];
      LvArray::tensorOps::copy< 3 >( faceNormal, stencil.m_faceNormal[iconn] );

      LvArray::tensorOps::scale< 3 >( faceNormal, phaseFlux / surface[i] );
      //change sign
      if( LvArray::tensorOps::AiBi< 3 >( stencil.m_cellToFaceVec[iconn][i], faceNormal ) < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
      }

      LvArray::tensorOps::hadamardProduct< 3 >( velocityNorm, faceNormal, stencil.m_cellToFaceVec[iconn][i] );
      for( int dir = 0; dir < 3; ++dir )
      {
        invDist[dir] = (LvArray::math::abs( cellCartDim[i][dir] ) > LvArray::NumericLimits< real64 >::epsilon) ?
                       1. / cellCartDim[i][dir] : LvArray::NumericLimits< real64 >::epsilon;
      }
      LvArray::tensorOps::hadamardProduct< 3 >( phaseVel, velocityNorm, invDist );
      LvArray::tensorOps::add< 3 >( phaseVelocity[er][esr][ei][ip], phaseVel );
    }
  }
}

void
StencilUtils::initVelocity( const CellElementStencilTPFAWrapper & stencil, localIndex iconn,
                            ElementRegionManager::ElementView< arrayView3d< real64 > > const & phaseVelocity )
{
  for( localIndex i = 0; i < 2; i++ )
  {
    localIndex const er = stencil.m_elementRegionIndices[iconn][i];
    localIndex const esr = stencil.m_elementSubRegionIndices[iconn][i];

    phaseVelocity[er][esr].zero();

  }
}


} /* namespace geos */
