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
 * @file singlePhaseProppantFluxKernels.cpp
 */

#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"
#include "SinglePhaseProppantFluxKernels.hpp"

namespace geos
{

namespace singlePhaseProppantFluxKernels
{

using namespace fluxKernelsHelper;


void FaceElementFluxKernel::
  launch( SurfaceElementStencilWrapper const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPerm_dDispJump,
          ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
          R1Tensor const & gravityVector,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  constexpr localIndex maxNumFluxElems = SurfaceElementStencilWrapper::maxNumPointsInFlux;
  constexpr localIndex maxStencilSize  = SurfaceElementStencilWrapper::maxStencilSize;
  constexpr localIndex maxNumConnections  = SurfaceElementStencilWrapper::maxNumConnections;

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );
    localIndex const numDofs = stencilSize; // pressures

    // For now, we have to filter out connections for which numElems == 1 in this function and not early on in
    // TwoPointFluxApproximation.cpp.
    // The reason for keeping the connections numElems == 1 is that the ProppantTransport solver needs these connections to produce correct
    // results.
    if( numFluxElems > 1 )
    {

      // working arrays
      stackArray1d< globalIndex, maxNumFluxElems > dofColIndices( numDofs );
      stackArray1d< localIndex, maxNumFluxElems > localColIndices( numFluxElems );

      stackArray1d< real64, maxNumFluxElems > localFlux( numFluxElems );
      stackArray2d< real64, maxNumFluxElems * maxStencilSize > localFluxJacobian( numFluxElems, numDofs );

      // need to store this for later use in determining the dFlux_dU terms when using better permeabilty approximations.
      stackArray2d< real64, maxNumFluxElems * maxStencilSize > dFlux_dAper( numFluxElems, stencilSize );

      // compute transmissibility
      real64 transmissibility[maxNumConnections][2]{};
      real64 dTrans_dPres[maxNumConnections][2]{};
      real64 dTrans_dDispJump[maxNumConnections][2][3]{};

      GEOS_UNUSED_VAR( dPerm_dPres, dPerm_dDispJump );
      stencilWrapper.computeWeights( iconn,
                                     permeability,
                                     permeabilityMultiplier,
                                     gravityVector,
                                     transmissibility );

      compute( stencilSize,
               seri[iconn],
               sesri[iconn],
               sei[iconn],
               transmissibility,
               dTrans_dPres,
               dTrans_dDispJump,
               pres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               dt,
               localFlux,
               localFluxJacobian,
               dFlux_dAper );

      // extract DOF numbers
      for( localIndex i = 0; i < numDofs; ++i )
      {
        dofColIndices[i] = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localColIndices[i] = sei( iconn, i );
      }

      for( localIndex i = 0; i < numFluxElems; ++i )
      {
        if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
        {
          globalIndex const globalRow = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );
          GEOS_ASSERT_GT( localMatrix.numRows(), localRow );

          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i].dataIfContiguous(),
                                                                            stencilSize );

        }
      }
    }
  } );

}


template< localIndex MAX_NUM_CONNECTIONS >
GEOS_HOST_DEVICE
void
FaceElementFluxKernel::compute( localIndex const numFluxElems,
                                arraySlice1d< localIndex const > const & seri,
                                arraySlice1d< localIndex const > const & sesri,
                                arraySlice1d< localIndex const > const & sei,
                                real64 const (&transmissibility)[MAX_NUM_CONNECTIONS][2],
                                real64 const (&dTrans_dPres)[MAX_NUM_CONNECTIONS][2],
                                real64 const (&dTrans_dDispJump)[MAX_NUM_CONNECTIONS][2][3],
                                ElementViewConst< arrayView1d< real64 const > > const & pres,
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const > > const & dens,
                                ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                ElementViewConst< arrayView1d< real64 const > > const & mob,
                                ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                real64 const dt,
                                arraySlice1d< real64 > const & flux,
                                arraySlice2d< real64 > const & fluxJacobian,
                                arraySlice2d< real64 > const & dFlux_dAperture )
{

  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numFluxElems; ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numFluxElems; ++k[1] )
    {
      real64 fluxVal = 0.0;
      real64 dFlux_dTrans = 0.0;
      real64 alpha = 0.0;
      real64 mobility = 0.0;
      real64 potGrad = 0.0;
      real64 trans[2] = {transmissibility[connectionIndex][0], transmissibility[connectionIndex][1]};
      real64 dTrans[2] = { dTrans_dPres[connectionIndex][0], dTrans_dPres[connectionIndex][1] };
      real64 dFlux_dP[2] = {0.0, 0.0};
      localIndex const regionIndex[2]    = {seri[k[0]], seri[k[1]]};
      localIndex const subRegionIndex[2] = {sesri[k[0]], sesri[k[1]]};
      localIndex const elementIndex[2]   = {sei[k[0]], sei[k[1]]};

      computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                              trans,
                              dTrans,
                              pres,
                              gravCoef,
                              dens,
                              dDens_dPres,
                              mob,
                              dMob_dPres,
                              alpha,
                              mobility,
                              potGrad,
                              fluxVal,
                              dFlux_dP,
                              dFlux_dTrans );

      // populate local flux vector and derivatives
      flux[k[0]] += dt * fluxVal;
      flux[k[1]] -= dt * fluxVal;

      real64 dFlux_dAper[2] = {0.0, 0.0};
      dFlux_dAper[0] = dt * dFlux_dTrans * dTrans_dDispJump[connectionIndex][0][0];
      dFlux_dAper[1] = -dt * dFlux_dTrans * dTrans_dDispJump[connectionIndex][1][0];

      fluxJacobian[k[0]][k[0]] += dFlux_dP[0] * dt;
      fluxJacobian[k[0]][k[1]] += dFlux_dP[1] * dt;
      fluxJacobian[k[1]][k[0]] -= dFlux_dP[0] * dt;
      fluxJacobian[k[1]][k[1]] -= dFlux_dP[1] * dt;

      dFlux_dAperture[k[0]][k[0]] += dFlux_dAper[0];
      dFlux_dAperture[k[0]][k[1]] += dFlux_dAper[1];
      dFlux_dAperture[k[1]][k[0]] -= dFlux_dAper[0];
      dFlux_dAperture[k[1]][k[1]] -= dFlux_dAper[1];

      connectionIndex++;
    }
  }
}


}// namespace singlePhaseProppantFluxKernels

} // namespace geos
