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
 * @file ProppantTransportKernels.cpp
 */

#include "ProppantTransportKernels.hpp"

#include "constitutive/fluid/singlefluid/ParticleFluidBase.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

namespace proppantTransportKernels
{

GEOS_HOST_DEVICE
inline
void
AccumulationKernel::
  compute( localIndex const numComps,
           real64 const proppantConc_n,
           real64 const proppantConcNew,
           arraySlice1d< real64 const > const & componentDens_n,
           arraySlice1d< real64 const > const & componentDensNew,
           arraySlice1d< real64 const > const & GEOS_UNUSED_PARAM( dCompDens_dPres ),
           arraySlice2d< real64 const > const & dCompDens_dCompConc,
           real64 const volume,
           real64 const packPoreVolume,
           real64 const proppantLiftVolume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian )
{

  // proppant mass conservation
  localAccum[0] = (proppantConcNew - proppantConc_n) * volume - proppantLiftVolume;

  for( localIndex c1 = 0; c1 < numComps; ++c1 )
  {
    for( localIndex c2 = 0; c2 < numComps; ++c2 )
    {
      localAccumJacobian[c1][c2] = 0.0;
    }
  }

  localAccumJacobian[0][0] = volume;

  // component mass conservation
  for( localIndex c1 = 0; c1 < numComps; ++c1 )
  {

    localAccum[c1+1] = ( componentDensNew[c1] * (1.0 - proppantConcNew) - componentDens_n[c1] * (1.0 - proppantConc_n) ) * volume
                       + (componentDensNew[c1] - componentDens_n[c1]) * packPoreVolume;

    for( localIndex c2 = 0; c2 < numComps; ++c2 )
    {
      localAccumJacobian[c1+1][c2+1] = dCompDens_dCompConc[c1][c2] * ( 1.0 - proppantConcNew ) * volume
                                       + dCompDens_dCompConc[c1][c2] * packPoreVolume;
    }

    localAccumJacobian[c1+1][0] = -componentDensNew[c1] * volume;
  }
}

void
AccumulationKernel::
  launch( localIndex const size,
          localIndex const numComps,
          localIndex const nDofs,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & proppantConc_n,
          arrayView1d< real64 const > const & proppantConc,
          arrayView2d< real64 const > const & componentDens_n,
          arrayView3d< real64 const > const & componentDens,
          arrayView3d< real64 const > const & dCompDens_dPres,
          arrayView4d< real64 const > const & dCompDens_dCompConc,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & proppantPackVolFrac,
          arrayView1d< real64 const > const & proppantLiftFlux,
          real64 const dt,
          real64 const maxProppantConcentration,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    if( elemGhostRank[ei] < 0 )
    {
      localIndex constexpr MAX_NC = constitutive::ParticleFluidBase::MAX_NUM_COMPONENTS;
      stackArray1d< globalIndex, MAX_NC > localAccumDOF( nDofs );
      stackArray1d< real64, MAX_NC > localAccum( nDofs );
      stackArray2d< real64, MAX_NC * MAX_NC > localAccumJacobian( nDofs, nDofs );

      real64 effectiveVolume = volume[ei];
      real64 packPoreVolume = 0.0;

      if( proppantPackVolFrac[ei] < 1.0 )
      {
        effectiveVolume = volume[ei] * ( 1.0 - proppantPackVolFrac[ei] );
        packPoreVolume = volume[ei] * proppantPackVolFrac[ei] * ( 1.0 - maxProppantConcentration );
      }

      real64 const proppantLiftVolume = proppantLiftFlux[ei] * dt;

      compute( numComps,
               proppantConc_n[ei],
               proppantConc[ei],
               componentDens_n[ei],
               componentDens[ei][0],
               dCompDens_dPres[ei][0],
               dCompDens_dCompConc[ei][0],
               effectiveVolume,
               packPoreVolume,
               proppantLiftVolume,
               localAccum,
               localAccumJacobian );

      globalIndex const elemDOF = dofNumber[ei];

      for( localIndex idof = 0; idof < nDofs; ++idof )
      {
        localAccumDOF[idof] = elemDOF + idof;
      }

      // add contribution to global residual and dRdP
      localIndex const localRow = dofNumber[ei] - rankOffset;
      for( localIndex idof = 0; idof < nDofs; ++idof )
      {
        localRhs[localRow + idof] += localAccum[idof];
        localMatrix.addToRow< serialAtomic >( localRow + idof,
                                              localAccumDOF.data(),
                                              localAccumJacobian[idof].dataIfContiguous(),
                                              nDofs );
      }
    }
  } );
}


template< localIndex MAX_NUM_FLUX_ELEMS >
GEOS_HOST_DEVICE
void
FluxKernel::
  computeJunction( localIndex const numElems,
                   localIndex const numDofPerCell,
                   arraySlice1d< localIndex const > const & stencilElementIndices,
                   arrayView1d< real64 const > const & pres,
                   arrayView1d< real64 const > const & proppantConc,
                   arrayView3d< real64 const > const & componentDens,
                   arrayView3d< real64 const > const & dComponentDens_dPres,
                   arrayView4d< real64 const > const & dComponentDens_dComponentConc,
                   arrayView1d< real64 const > const & gravDepth,
                   arrayView2d< real64 const > const & dens,
                   arrayView2d< real64 const > const & dDens_dPres,
                   arrayView2d< real64 const > const & dDens_dProppantConc,
                   arrayView3d< real64 const > const & dDens_dComponentConc,
                   arrayView2d< real64 const > const & visc,
                   arrayView2d< real64 const > const & dVisc_dPres,
                   arrayView2d< real64 const > const & dVisc_dProppantConc,
                   arrayView3d< real64 const > const & dVisc_dComponentConc,
                   arrayView2d< real64 const > const & fluidDensity,
                   arrayView2d< real64 const > const & GEOS_UNUSED_PARAM( dFluidDens_dPres ),
                   arrayView3d< real64 const > const & GEOS_UNUSED_PARAM( dFluidDens_dComponentConc ),
                   arrayView1d< real64 const > const & settlingFactor,
                   arrayView1d< real64 const > const & GEOS_UNUSED_PARAM( dSettlingFactor_dPres ),
                   arrayView1d< real64 const > const & dSettlingFactor_dProppantConc,
                   arrayView2d< real64 const > const & GEOS_UNUSED_PARAM( dSettlingFactor_dComponentConc ),
                   arrayView1d< real64 const > const & collisionFactor,
                   arrayView1d< real64 const > const & dCollisionFactor_dProppantConc,
                   arrayView1d< integer const > const & isProppantMobile,
                   real64 const (&transmissibility)[MAX_NUM_FLUX_ELEMS],
                   real64 const (&apertureWeight)[MAX_NUM_FLUX_ELEMS],
                   real64 const (&geometricWeight)[MAX_NUM_FLUX_ELEMS],
                   real64 const dt,
                   arraySlice1d< real64 > const & localFlux,
                   arraySlice2d< real64 > const & localFluxJacobian )
{

  // We assume numElems == stencilSize;

  constexpr real64 TINY = 1e-10;

  constexpr localIndex MAX_NUM_COMPONENTS = constitutive::ParticleFluidBase::MAX_NUM_COMPONENTS;

  localIndex const numComps = numDofPerCell - 1;

  // mixture density and fluid density in each face
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > mixDens( numElems );
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > fluidDens( numElems );

  // related to slip velocity calculation
  real64 edgeDensity = 0.0;
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dEdgeDens_dP( numElems );
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dEdgeDens_dProppantC( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dEdgeDens_dComponentC( numElems, numComps );

  real64 edgeViscosity = 0.0;
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dEdgeVisc_dP( numElems );
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dEdgeVisc_dProppantC( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dEdgeVisc_dComponentC( numElems, numComps );

  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dPe_dP( numElems );
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dPe_dProppantC( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dPe_dComponentC( numElems, numComps );

  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > edgeToFaceFlux( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS > dEdgeToFaceFlux_dP( numElems, numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS > dEdgeToFaceFlux_dProppantC( numElems, numElems );
  stackArray3d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dEdgeToFaceFlux_dComponentC( numElems, numElems, numComps );

  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > edgeToFaceProppantFlux( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS > dEdgeToFaceProppantFlux_dP( numElems, numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS > dEdgeToFaceProppantFlux_dProppantC( numElems, numElems );
  stackArray3d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dEdgeToFaceProppantFlux_dComponentC( numElems, numElems, numComps );

  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > edgeToFaceFluidFlux( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS > dEdgeToFaceFluidFlux_dP( numElems, numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS > dEdgeToFaceFluidFlux_dProppantC( numElems, numElems );
  stackArray3d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dEdgeToFaceFluidFlux_dComponentC( numElems, numElems, numComps );

  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > proppantC( numElems );

  localIndex numberOfMobileProppantElems = 0;

  //get averaged edgeDensity and edgeViscosity
  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    edgeDensity += geometricWeight[i] * dens[ei][0];
    dEdgeDens_dP[i] = geometricWeight[i] * dDens_dPres[ei][0];
    dEdgeDens_dProppantC[i] = geometricWeight[i] * dDens_dProppantConc[ei][0];

    edgeViscosity += geometricWeight[i] * visc[ei][0];
    dEdgeVisc_dP[i] = geometricWeight[i] * dVisc_dPres[ei][0];
    dEdgeVisc_dProppantC[i] = geometricWeight[i] * dVisc_dProppantConc[ei][0];

    proppantC[i] = proppantConc[ei];

    mixDens[i] = dens[ei][0];
    fluidDens[i] = fluidDensity[ei][0];

    if( isProppantMobile[ei] == 1 )
    {
      numberOfMobileProppantElems++;
    }

    for( localIndex c = 0; c < numComps; ++c )
    {
      dEdgeDens_dComponentC[i][c] = geometricWeight[i] * dDens_dComponentConc[ei][0][c];
      dEdgeVisc_dComponentC[i][c] = geometricWeight[i] * dVisc_dComponentConc[ei][0][c];
    }
  }

  real64 const proppantFluxCoef = ( numberOfMobileProppantElems > 1 ) ? 1.0 : 0.0;

  real64 transmissibilitySum = 0.0;
  real64 Pe = 0.0;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    Pe += transmissibility[i] * (pres[ei] - gravTerm);
    transmissibilitySum += transmissibility[i];
    dPe_dP[i] += transmissibility[i];

    for( localIndex j = 0; j < numElems; ++j )
    {
      dPe_dP[j] += -transmissibility[i] * gravD * dEdgeDens_dP[j];
      dPe_dProppantC[j] += -transmissibility[i] * gravD * dEdgeDens_dProppantC[j];

      for( localIndex c = 0; c < numComps; ++c )
      {
        dPe_dComponentC[j][c] += -transmissibility[i] * gravD * dEdgeDens_dComponentC[j][c];
      }
    }
  }

  for( localIndex i = 0; i < numElems; ++i )
  {
    dPe_dP[i] /= transmissibilitySum;
    dPe_dProppantC[i] /= transmissibilitySum;

    for( localIndex c = 0; c < numComps; ++c )
    {
      dPe_dComponentC[i][c] /= transmissibilitySum;
    }
  }

  Pe /= transmissibilitySum;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    real64 const fluxTerm = Pe - (pres[ei] - gravTerm);

    edgeToFaceFlux[i] = transmissibility[i] * fluxTerm / edgeViscosity;
    dEdgeToFaceFlux_dP[i][i] += -transmissibility[i] / edgeViscosity;

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFlux_dP[i][j] += -transmissibility[i] * fluxTerm * dEdgeVisc_dP[j] / (edgeViscosity * edgeViscosity)
                                  + transmissibility[i] * (dPe_dP[j] + dEdgeDens_dP[j] * gravD) / edgeViscosity;

      dEdgeToFaceFlux_dProppantC[i][j] += -transmissibility[i] * fluxTerm * dEdgeVisc_dProppantC[j] / (edgeViscosity * edgeViscosity)
                                          + transmissibility[i] * (dPe_dProppantC[j] + dEdgeDens_dProppantC[j] * gravD) / edgeViscosity;

      for( localIndex c = 0; c < numComps; ++c )
      {
        dEdgeToFaceFlux_dComponentC[i][j][c] += -transmissibility[i] * fluxTerm * dEdgeVisc_dComponentC[j][c] / (edgeViscosity * edgeViscosity)
                                                + transmissibility[i] * (dPe_dComponentC[j][c] + dEdgeDens_dComponentC[j][c] * gravD) / edgeViscosity;
      }
    }

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceProppantFlux_dProppantC[i][j] = 0.0;
      for( localIndex c = 0; c < numComps; ++c )
      {
        dEdgeToFaceProppantFlux_dComponentC[i][j][c] = 0.0;
      }
    }

    if( LvArray::math::abs( apertureWeight[i] ) > TINY )
    {
      // vertical
      edgeToFaceProppantFlux[i] = (1.0 - proppantC[i]) * settlingFactor[ei] * apertureWeight[i] * fluidDens[i] / mixDens[i];

      dEdgeToFaceProppantFlux_dProppantC[i][i] = (-settlingFactor[ei] + (1 - proppantC[i]) * dSettlingFactor_dProppantConc[ei]) *
                                                 apertureWeight[i] * fluidDens[i] / mixDens[i];

      edgeToFaceProppantFlux[i] += proppantFluxCoef * edgeToFaceFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dEdgeToFaceProppantFlux_dProppantC[i][j] += proppantFluxCoef * dEdgeToFaceFlux_dProppantC[i][j];

        for( localIndex c = 0; c < numComps; ++c )
        {
          dEdgeToFaceProppantFlux_dComponentC[i][j][c] += proppantFluxCoef * dEdgeToFaceFlux_dComponentC[i][j][c];
        }
      }
    }
    else
    {
      // horizontal
      edgeToFaceProppantFlux[i] = (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFactor[ei]) * proppantFluxCoef * edgeToFaceFlux[i];

      dEdgeToFaceProppantFlux_dProppantC[i][i] = -fluidDens[i] / mixDens[i] * (collisionFactor[ei] - (1.0 - proppantC[i]) * dCollisionFactor_dProppantConc[ei]) *
                                                 proppantFluxCoef * edgeToFaceFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dEdgeToFaceProppantFlux_dProppantC[i][j] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFactor[ei]) *
                                                    proppantFluxCoef * dEdgeToFaceFlux_dProppantC[i][j];

        for( localIndex c = 0; c < numComps; ++c )
        {
          dEdgeToFaceProppantFlux_dComponentC[i][j][c] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFactor[ei]) *
                                                          proppantFluxCoef * dEdgeToFaceFlux_dComponentC[i][j][c];
        }
      }
    }

    // fluid flux for component transport
    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFluidFlux_dP[i][j] = 0.0;
      dEdgeToFaceFluidFlux_dProppantC[i][j] = 0.0;

      for( localIndex c = 0; c < numComps; ++c )
      {
        dEdgeToFaceFluidFlux_dComponentC[i][j][c] = 0.0;
      }
    }

    // note that all the fluid properties are from previous time step
    real64 const fluidFluxCoef = ( isProppantMobile[ei] == 0 || numElems == 1 ) ? 0.0 : 1.0;

    edgeToFaceFluidFlux[i] = mixDens[i] / fluidDens[i] * edgeToFaceFlux[i] - fluidFluxCoef * (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) /
                             fluidDens[i] * edgeToFaceProppantFlux[i];

    dEdgeToFaceFluidFlux_dProppantC[i][i] = -fluidFluxCoef * edgeToFaceProppantFlux[i];

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFluidFlux_dProppantC[i][j] += mixDens[i] / fluidDens[i] * dEdgeToFaceFlux_dProppantC[i][j] - fluidFluxCoef *
                                               (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) / fluidDens[i] * dEdgeToFaceProppantFlux_dProppantC[i][j];
      for( localIndex c = 0; c < numComps; ++c )
      {
        dEdgeToFaceFluidFlux_dComponentC[i][j][c] += mixDens[i] / fluidDens[i] * dEdgeToFaceFlux_dComponentC[i][j][c] - fluidFluxCoef *
                                                     (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) / fluidDens[i] *
                                                     dEdgeToFaceProppantFlux_dComponentC[i][j][c];
      }
    }

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFluidFlux_dProppantC[i][j] = 0.0;
      for( localIndex c = 0; c < numComps; ++c )
      {
        dEdgeToFaceFluidFlux_dComponentC[i][j][c] = 0.0;
      }
    }
  }

  // get proppantCe

  real64 proppantCe = 0.0;
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dProppantCe_dP( numElems );
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dProppantCe_dProppantC( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dProppantCe_dComponentC( numElems, numComps );

  real64 downStreamFlux = 0.0;
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dDownStreamFlux_dP( numElems );
  stackArray1d< real64, MAX_NUM_FLUX_ELEMS > dDownStreamFlux_dProppantC( numElems );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dDownStreamFlux_dComponentC( numElems, numComps );

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei = stencilElementIndices[i];

    if( isProppantMobile[ei] == 0 )
    {
      continue;
    }

    if( edgeToFaceProppantFlux[i] >= 0.0 )
    {
      // downstream
      downStreamFlux += edgeToFaceProppantFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dDownStreamFlux_dP[j] += dEdgeToFaceProppantFlux_dP[i][j];
        dDownStreamFlux_dProppantC[j] += dEdgeToFaceProppantFlux_dProppantC[i][j];

        for( localIndex c = 0; c < numComps; ++c )
        {
          dDownStreamFlux_dComponentC[j][c] += dEdgeToFaceProppantFlux_dComponentC[i][j][c];
        }
      }
    }
    else
    {
      // upstream
      proppantCe += -edgeToFaceProppantFlux[i] * proppantC[i];
      dProppantCe_dProppantC[i] += -edgeToFaceProppantFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dProppantCe_dP[j] += -dEdgeToFaceProppantFlux_dP[i][j] * proppantC[i];
        dProppantCe_dProppantC[j] += -dEdgeToFaceProppantFlux_dProppantC[i][j] * proppantC[i];

        for( localIndex c = 0; c < numComps; ++c )
        {
          dProppantCe_dComponentC[j][c] += -dEdgeToFaceProppantFlux_dComponentC[i][j][c] * proppantC[i];
        }
      }
    }
  }

  if( downStreamFlux > 0.0 )
  {
    for( localIndex i = 0; i < numElems; ++i )
    {
      localIndex const ei = stencilElementIndices[i];

      if( isProppantMobile[ei] == 0 )
      {
        continue;
      }

      dProppantCe_dP[i] = dProppantCe_dP[i] / downStreamFlux - proppantCe * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
      dProppantCe_dProppantC[i] = dProppantCe_dProppantC[i] / downStreamFlux - proppantCe * dDownStreamFlux_dProppantC[i] /
                                  (downStreamFlux * downStreamFlux);;

      for( localIndex c = 0; c < numComps; ++c )
      {
        dProppantCe_dComponentC[i][c] = dProppantCe_dComponentC[i][c] / downStreamFlux - proppantCe * dDownStreamFlux_dComponentC[i][c] /
                                        (downStreamFlux * downStreamFlux);
      }
    }

    proppantCe = proppantCe / downStreamFlux;
  }
  else
  {
    proppantCe = 0.0;
    for( localIndex i = 0; i < numElems; ++i )
    {
      localIndex const ei = stencilElementIndices[i];

      if( isProppantMobile[ei] == 0 )
      {
        continue;
      }

      dProppantCe_dP[i] = 0.0;
      dProppantCe_dProppantC[i] =  geometricWeight[i];
      proppantCe += proppantC[i] * geometricWeight[i];
    }
  }

  // get componentCe

  stackArray1d< real64, MAX_NUM_COMPONENTS > componentCe( numComps );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dComponentCe_dP( numElems, numComps );
  stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS > dComponentCe_dProppantC( numElems, numComps );
  stackArray3d< real64, MAX_NUM_FLUX_ELEMS * MAX_NUM_COMPONENTS * MAX_NUM_COMPONENTS > dComponentCe_dComponentC( numElems, numComps, numComps );

  downStreamFlux = 0.0;
  for( localIndex i = 0; i < numElems; ++i )
  {
    dDownStreamFlux_dP[i] = 0.0;
    dDownStreamFlux_dProppantC[i] = 0.0;
    for( localIndex c = 0; c < numComps; ++c )
    {
      dDownStreamFlux_dComponentC[i][c] = 0.0;
    }
  }

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei = stencilElementIndices[i];

    if( edgeToFaceFluidFlux[i] >= 0.0 )
    {
      // downstream
      downStreamFlux += edgeToFaceFluidFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dDownStreamFlux_dP[j] += dEdgeToFaceFluidFlux_dP[i][j];
        dDownStreamFlux_dProppantC[j] += dEdgeToFaceFluidFlux_dProppantC[i][j];

        for( localIndex c = 0; c < numComps; ++c )
        {
          dDownStreamFlux_dComponentC[j][c] += dEdgeToFaceFluidFlux_dComponentC[i][j][c];
        }
      }
    }
    else
    {
      // upstream
      for( localIndex c1 = 0; c1 < numComps; ++c1 )
      {
        componentCe[c1] += -edgeToFaceFluidFlux[i] * componentDens[ei][0][c1];
        dComponentCe_dP[i][c1] += -edgeToFaceFluidFlux[i] * dComponentDens_dPres[ei][0][c1];

        for( localIndex c2 = 0; c2 < numComps; ++c2 )
        {
          dComponentCe_dComponentC[i][c1][c2] += -edgeToFaceFluidFlux[i] * dComponentDens_dComponentConc[ei][0][c1][c2];
        }

        for( localIndex j = 0; j < numElems; ++j )
        {
          dComponentCe_dP[j][c1] += -dEdgeToFaceFluidFlux_dP[i][j] * componentDens[ei][0][c1];
          dComponentCe_dProppantC[j][c1] += -dEdgeToFaceFluidFlux_dProppantC[i][j] * componentDens[ei][0][c1];

          for( localIndex c2 = 0; c2 < numComps; ++c2 )
          {
            dComponentCe_dComponentC[j][c1][c2] += -dEdgeToFaceFluidFlux_dComponentC[i][j][c2] * componentDens[ei][0][c1];
          }
        }
      }
    }
  }

  if( downStreamFlux > 0.0 )
  {
    for( localIndex c1 = 0; c1 < numComps; ++c1 )
    {
      for( localIndex i = 0; i < numElems; ++i )
      {
        dComponentCe_dP[i][c1] = dComponentCe_dP[i][c1] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
        dComponentCe_dProppantC[i][c1] =
          dComponentCe_dProppantC[i][c1] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dProppantC[i] / (downStreamFlux * downStreamFlux);

        for( localIndex c2 = 0; c2 < numComps; ++c2 )
        {
          dComponentCe_dComponentC[i][c1][c2] =
            dComponentCe_dComponentC[i][c1][c2] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dComponentC[i][c2] / (downStreamFlux * downStreamFlux);

        }
      }

      componentCe[c1] = componentCe[c1] / downStreamFlux;
    }
  }
  else
  {
    for( localIndex c = 0; c < numComps; ++c )
    {
      componentCe[c] = 0.0;

      for( localIndex i = 0; i < numElems; ++i )
      {
        localIndex const ei = stencilElementIndices[i];

        componentCe[c] += componentDens[ei][0][c] * geometricWeight[i];
        dComponentCe_dP[i][c] = dComponentDens_dPres[ei][0][c] * geometricWeight[i];

        for( localIndex c2 = 0; c2 < numComps; ++c2 )
        {
          dComponentCe_dComponentC[i][c][c2] = dComponentDens_dComponentConc[ei][0][c][c2] * geometricWeight[i];
        }
      }
    }
  }

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei = stencilElementIndices[i];

    localIndex idx1 = i * numDofPerCell; // proppant

    if( isProppantMobile[ei] == 1 && !(numElems == 1 && apertureWeight[i] > TINY) )
    {
      if( edgeToFaceProppantFlux[i] >= 0.0 )
      {
        localFlux[idx1] = -proppantCe * edgeToFaceProppantFlux[i] * dt;
      }
      else
      {
        localFlux[idx1] = -proppantC[i] * edgeToFaceProppantFlux[i] * dt;
      }

      for( localIndex j = 0; j < numElems; ++j )
      {
        localIndex idx2 = j * numDofPerCell;

        if( edgeToFaceProppantFlux[i] >= 0.0 )
        {
          localFluxJacobian[idx1][idx2] =
            -(dProppantCe_dProppantC[j] * edgeToFaceProppantFlux[i] + proppantCe * dEdgeToFaceProppantFlux_dProppantC[i][j]) * dt;
        }
        else
        {
          localFluxJacobian[idx1][idx2] = -proppantC[i] * dEdgeToFaceProppantFlux_dProppantC[i][j] * dt;

          if( i == j )
          {
            localFluxJacobian[idx1][idx2] += -edgeToFaceProppantFlux[i] * dt;
          }
        }
      }
    }
    else
    {
      localFlux[idx1] = 0.0;

      for( localIndex j = 0; j < numElems; ++j )
      {
        for( localIndex c = 0; c < numDofPerCell; ++c )
        {
          localIndex idx2 = j * numDofPerCell + c;
          localFluxJacobian[idx1][idx2] = 0.0;
        }
      }
    }

    // component

    if( numElems > 1 )
    {
      for( localIndex c1 = 0; c1 < numComps; ++c1 )
      {
        idx1 = i * numDofPerCell + 1 + c1;

        if( edgeToFaceFluidFlux[i] >= 0.0 )
        {
          localFlux[idx1] = -componentCe[c1] * edgeToFaceFluidFlux[i] * dt;
        }
        else
        {
          localFlux[idx1] = -componentDens[ei][0][c1] * edgeToFaceFluidFlux[i] * dt;
        }

        for( localIndex j = 0; j < numElems; ++j )
        {
          localIndex idx2 = j * numDofPerCell;

          if( edgeToFaceFluidFlux[i] >= 0.0 )
          {
            for( localIndex c2 = 0; c2 < numComps; ++c2 )
            {
              localFluxJacobian[idx1][idx2+1+c2] = -( dComponentCe_dComponentC[j][c1][c2] * edgeToFaceFluidFlux[i] ) * dt;
            }
          }
          else
          {
            if( i == j )
            {
              for( localIndex c2 = 0; c2 < numComps; ++c2 )
              {
                localFluxJacobian[idx1][idx2 + 1 + c2] += -dComponentDens_dComponentConc[ei][0][c1][c2] * edgeToFaceFluidFlux[i] * dt;
              }
            }
          }
        }
      }
    }
  }
}


void FluxKernel::
  launch( SurfaceElementStencilWrapper const & stencilWrapper,
          localIndex const numDofPerCell,
          real64 const dt,
          globalIndex const rankOffset,
          integer const updateProppantPacking,
          R1Tensor const & unitGravityVector,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & proppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & componentDens,
          ElementViewConst< arrayView3d< real64 const > > const & dComponentDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dComponentDens_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dDens_dComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & visc,
          ElementViewConst< arrayView2d< real64 const > > const & dVisc_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & dVisc_dProppantConc,
          ElementViewConst< arrayView3d< real64 const > > const & dVisc_dComponentConc,
          ElementViewConst< arrayView2d< real64 const > > const & fluidDensity,
          ElementViewConst< arrayView2d< real64 const > > const & dFluidDens_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dFluidDens_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & settlingFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dProppantConc,
          ElementViewConst< arrayView2d< real64 const > > const & dSettlingFactor_dComponentConc,
          ElementViewConst< arrayView1d< real64 const > > const & collisionFactor,
          ElementViewConst< arrayView1d< real64 const > > const & dCollisionFactor_dProppantConc,
          ElementViewConst< arrayView1d< integer const > > const & isProppantMobile,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
          ElementViewConst< arrayView1d< real64 const > > const & aperture,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  constexpr localIndex maxNumFluxElems = SurfaceElementStencilWrapper::maxNumPointsInFlux;
  constexpr localIndex maxStencilSize = SurfaceElementStencilWrapper::maxStencilSize;

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  constexpr localIndex DOF1 = maxNumFluxElems * constitutive::ParticleFluidBase::MAX_NUM_COMPONENTS;
  constexpr localIndex DOF2 = maxStencilSize * constitutive::ParticleFluidBase::MAX_NUM_COMPONENTS;

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );

    if( ( numFluxElems > 1 || updateProppantPacking != 0 ) )
    {
      localIndex const stencilSize  = numFluxElems;
      localIndex const nDofs = numFluxElems * numDofPerCell;

      // working arrays
      stackArray1d< globalIndex, DOF2 > dofColIndices( nDofs );

      stackArray1d< real64, DOF1 > localFlux( nDofs );
      stackArray2d< real64, DOF1 * DOF2 > localFluxJacobian( nDofs, nDofs );

      localIndex const er = seri[iconn][0];
      localIndex const esr = sesri[iconn][0];

      // compute transmissibility
      real64 transmissibility[maxNumFluxElems]{};
      real64 apertureWeight[maxNumFluxElems]{};
      real64 geometricWeight[maxNumFluxElems]{};
      stencilWrapper.computeWeights( iconn,
                                     permeability,
                                     permeabilityMultiplier,
                                     aperture,
                                     unitGravityVector,
                                     transmissibility,
                                     apertureWeight,
                                     geometricWeight );

      computeJunction( numFluxElems,
                       numDofPerCell,
                       sei[iconn],
                       pres[er][esr],
                       proppantConc[er][esr],
                       componentDens[er][esr],
                       dComponentDens_dPres[er][esr],
                       dComponentDens_dComponentConc[er][esr],
                       gravDepth[er][esr],
                       dens[er][esr],
                       dDens_dPres[er][esr],
                       dDens_dProppantConc[er][esr],
                       dDens_dComponentConc[er][esr],
                       visc[er][esr],
                       dVisc_dPres[er][esr],
                       dVisc_dProppantConc[er][esr],
                       dVisc_dComponentConc[er][esr],
                       fluidDensity[er][esr],
                       dFluidDens_dPres[er][esr],
                       dFluidDens_dComponentConc[er][esr],
                       settlingFactor[er][esr],
                       dSettlingFactor_dPres[er][esr],
                       dSettlingFactor_dProppantConc[er][esr],
                       dSettlingFactor_dComponentConc[er][esr],
                       collisionFactor[er][esr],
                       dCollisionFactor_dProppantConc[er][esr],
                       isProppantMobile[er][esr],
                       transmissibility,
                       apertureWeight,
                       geometricWeight,
                       dt,
                       localFlux,
                       localFluxJacobian );

      for( localIndex i = 0; i < stencilSize; ++i )
      {
        for( localIndex j = 0; j < numDofPerCell; ++j )
        {
          dofColIndices[i * numDofPerCell + j] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + j;
        }
      }

      for( localIndex i = 0; i < numFluxElems; ++i )
      {
        if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
        {
          globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );
          GEOS_ASSERT_GE( localMatrix.numRows(), localRow + numDofPerCell );

          for( localIndex idof = 0; idof < numDofPerCell; ++idof )
          {
            RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + idof], localFlux[i * numDofPerCell + idof] );
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + idof,
                                                                              dofColIndices.data(),
                                                                              localFluxJacobian[i * numDofPerCell + idof].dataIfContiguous(),
                                                                              stencilSize * numDofPerCell );
          }
        }
      }
    }
  } );
}


template< localIndex MAX_NUM_FLUX_ELEMS >
GEOS_HOST_DEVICE
void
FluxKernel::
  computeCellBasedFlux( localIndex const numElems,
                        arraySlice1d< localIndex const > const & stencilElementIndices,
                        arrayView1d< real64 const > const & pres,
                        arrayView1d< real64 const > const & gravDepth,
                        arrayView2d< real64 const > const & dens,
                        arrayView2d< real64 const > const & visc,
                        arraySlice1d< R1Tensor const > const & cellCenterToEdgeCenters,
                        real64 const (&transmissibility)[MAX_NUM_FLUX_ELEMS],
                        real64 const (&apertureWeight)[MAX_NUM_FLUX_ELEMS],
                        real64 const (&geometricWeight)[MAX_NUM_FLUX_ELEMS],
                        arrayView2d< real64 > const & cellBasedFlux )
{
  GEOS_UNUSED_VAR( apertureWeight );

  // get averaged edgeDensity and edgeViscosity
  real64 edgeDensity = 0.0;
  real64 edgeViscosity = 0.0;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];
    edgeDensity += geometricWeight[i] * dens[ei][0];
    edgeViscosity += geometricWeight[i] * visc[ei][0];
  }

  real64 transmissibilitySum = 0.0;
  real64 Pe = 0.0;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    Pe += transmissibility[i] * (pres[ei] - gravTerm);
    transmissibilitySum += transmissibility[i];
  }

  Pe /= transmissibilitySum;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    real64 const fluxTerm = Pe - (pres[ei] - gravTerm);
    real64 const edgeToFaceFlux = transmissibility[i] * fluxTerm / edgeViscosity;

    for( localIndex k = 0; k < 3; ++k )
    {
      RAJA::atomicAdd( parallelDeviceAtomic{}, &cellBasedFlux[ei][k], -edgeToFaceFlux * cellCenterToEdgeCenters[i][k] );
    }
  }
}

void FluxKernel::
  launchCellBasedFluxCalculation( SurfaceElementStencilWrapper const & stencilWrapper,
                                  R1Tensor const & unitGravityVector,
                                  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dens,
                                  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & visc,
                                  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
                                  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                  FluxKernel::ElementView< arrayView2d< real64 > > const & cellBasedFlux )
{

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  constexpr localIndex maxNumFluxElems = SurfaceElementStencilWrapper::maxNumPointsInFlux;

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencilWrapper.getCellCenterToEdgeCenters();

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    // compute transmissibility
    real64 transmissibility[maxNumFluxElems]{};
    real64 apertureWeight[maxNumFluxElems]{};
    real64 geometricWeight[maxNumFluxElems]{};
    stencilWrapper.computeWeights( iconn,
                                   permeability,
                                   permeabilityMultiplier,
                                   aperture,
                                   unitGravityVector,
                                   transmissibility,
                                   apertureWeight,
                                   geometricWeight );

    computeCellBasedFlux( numFluxElems,
                          sei[iconn],
                          pres[er][esr],
                          gravDepth[er][esr],
                          dens[er][esr],
                          visc[er][esr],
                          cellCenterToEdgeCenters[iconn],
                          transmissibility,
                          apertureWeight,
                          geometricWeight,
                          cellBasedFlux[er][esr] );
  } );
}


GEOS_HOST_DEVICE
inline
void
ProppantPackVolumeKernel::
  computeProppantPackVolume( localIndex const numElems,
                             real64 const dt,
                             real64 const proppantDensity,
                             real64 const proppantDiameter,
                             real64 const maxProppantConcentration,
                             R1Tensor const & unitGravityVector,
                             real64 const criticalShieldsNumber,
                             real64 const frictionCoefficient,
                             arraySlice1d< localIndex const > const & stencilElementIndices,
                             arraySlice1d< real64 const > const & stencilWeights,
                             arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                             arrayView1d< real64 const > const & settlingFactor,
                             arrayView2d< real64 const > const & density,
                             arrayView2d< real64 const > const & fluidDensity,
                             arrayView2d< real64 const > const &,
                             arrayView1d< real64 const > const & volume,
                             arrayView1d< real64 const > const & aperture,
                             arrayView1d< integer const > const & elemGhostRank,
                             arrayView1d< integer const > const & isProppantBoundaryElement,
                             arrayView1d< integer const > const & isProppantMobile,
                             arrayView2d< real64 const > const & cellBasedFlux,
                             arrayView1d< real64 > const & conc,
                             arrayView1d< real64 > const & proppantPackVolFrac,
                             arrayView1d< real64 > const & proppantExcessPackVolume,
                             arrayView1d< real64 > const & proppantLiftFlux )
{

  integer faceIndex = -1;

  real64 constexpr TINY = 1e-10;

  real64 const stencilCellToEdgeDistance0 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[0] );
  real64 const edgeLength = stencilWeights[0] * stencilCellToEdgeDistance0;

  if( numElems == 1 )
  {
    localIndex const ei  = stencilElementIndices[0];

    real64 const stencilEdgeToFaceDownDistance =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance0;

    if( stencilEdgeToFaceDownDistance < -TINY && isProppantMobile[ei] == 1 )
    {
      faceIndex = 0;
    }
  }
  else if( numElems == 2 )
  {
    localIndex const ei0  = stencilElementIndices[0];
    localIndex const ei1  = stencilElementIndices[1];

    real64 const stencilCellToEdgeDistance1 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[1] );

    real64 const stencilEdgeToFaceDownDistance0 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance0;

    real64 const stencilEdgeToFaceDownDistance1 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[1], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance1;

    //0: top  1: bottom

    if( stencilEdgeToFaceDownDistance0 < -TINY && stencilEdgeToFaceDownDistance1 > TINY && isProppantMobile[ei0] == 1 && isProppantMobile[ei1] == 0 )
    {
      faceIndex = 0;
    }

    //0: bottom  1: top

    if( stencilEdgeToFaceDownDistance0 > TINY && stencilEdgeToFaceDownDistance1 < -TINY && isProppantMobile[ei1] == 1 && isProppantMobile[ei0] == 0 )
    {
      faceIndex = 1;
    }
  }

  if( faceIndex >= 0 )
  {
    localIndex const ei = stencilElementIndices[faceIndex];

    if( elemGhostRank[ei] < 0 && isProppantBoundaryElement[ei] == 0 )
    {
      real64 const L = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[faceIndex] ) * 2.0;

      real64 velocity[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( cellBasedFlux[ei] );

      real64 const downVelocity = LvArray::tensorOps::AiBi< 3 >( cellBasedFlux[ei], unitGravityVector );

      LvArray::tensorOps::scaledAdd< 3 >( velocity, unitGravityVector, -downVelocity );

      real64 const velocityMag = LvArray::tensorOps::l2Norm< 3 >( velocity ) / volume[ei];

      real64 dH = fluidDensity[ei][0] / density[ei][0] * (1.0 - conc[ei]) * settlingFactor[ei] * conc[ei] / maxProppantConcentration * dt;

      real64 const tau = 1.0/8.0 * frictionCoefficient * fluidDensity[ei][0] * velocityMag * velocityMag;

      real64 const ShieldsNumber = tau / (proppantDensity - fluidDensity[ei][0]) / 9.81 / proppantDiameter;

      proppantLiftFlux[ei] = 0.0;

      if( ShieldsNumber > criticalShieldsNumber )
      {
        proppantLiftFlux[ei] = aperture[ei] *
                               (proppantDiameter * sqrt( 9.81 * proppantDiameter * (proppantDensity - fluidDensity[ei][0]) / fluidDensity[ei][0] )) *
                               (9.64 * pow( ShieldsNumber, 0.166 )) * pow( ShieldsNumber - criticalShieldsNumber, 1.5 );
      }

      if( proppantLiftFlux[ei] < 0.0 )
      {
        proppantLiftFlux[ei] = 0.0;
      }

      real64 liftH =  proppantLiftFlux[ei] / edgeLength / aperture[ei] / maxProppantConcentration * dt;

      real64 volFrac = proppantPackVolFrac[ei] + (dH - liftH) / L;

      if( volFrac < 0.0 )
      {
        volFrac = 0.0;
        liftH = dH + proppantPackVolFrac[ei] * L;
        proppantLiftFlux[ei] = liftH * edgeLength * aperture[ei] * maxProppantConcentration / dt;
      }

      if( volFrac >= 1.0 )
      {
        proppantExcessPackVolume[ei] = (volFrac - 1.0) * L;
        proppantPackVolFrac[ei] = 1.0;
        conc[ei] = maxProppantConcentration;
      }
      else
      {
        proppantPackVolFrac[ei] = volFrac;
      }
    }
  }
}

void ProppantPackVolumeKernel::
  launchProppantPackVolumeCalculation( SurfaceElementStencil const & stencil,
                                       real64 const dt,
                                       real64 const proppantDensity,
                                       real64 const proppantDiameter,
                                       real64 const maxProppantConcentration,
                                       R1Tensor const & unitGravityVector,
                                       real64 const criticalShieldsNumber,
                                       real64 const frictionCoefficient,
                                       ElementView< arrayView1d< real64 const > > const & settlingFactor,
                                       ElementView< arrayView2d< real64 const > > const & density,
                                       ElementView< arrayView2d< real64 const > > const & fluidDensity,
                                       ElementView< arrayView2d< real64 const > > const & fluidViscosity,
                                       ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                       ElementView< arrayView1d< integer const > > const & isProppantBoundaryElement,
                                       ElementView< arrayView1d< real64 const > > const & aperture,
                                       ElementView< arrayView1d< real64 const > > const & volume,
                                       ElementView< arrayView1d< integer const > > const & elemGhostRank,
                                       ElementView< arrayView2d< real64 const > > const & cellBasedFlux,
                                       ElementView< arrayView1d< real64 > > const & conc,
                                       ElementView< arrayView1d< real64 > > const & proppantPackVolFrac,
                                       ElementView< arrayView1d< real64 > > const & proppantExcessPackVolume,
                                       ElementView< arrayView1d< real64 > > const & proppantLiftFlux )
{

  typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename SurfaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = seri.sizeOfArray( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    computeProppantPackVolume( numFluxElems,
                               dt,
                               proppantDensity,
                               proppantDiameter,
                               maxProppantConcentration,
                               unitGravityVector,
                               criticalShieldsNumber,
                               frictionCoefficient,
                               sei[iconn],
                               weights[iconn],
                               cellCenterToEdgeCenters[iconn],
                               settlingFactor[er][esr],
                               density[er][esr],
                               fluidDensity[er][esr],
                               fluidViscosity[er][esr],
                               volume[er][esr],
                               aperture[er][esr],
                               elemGhostRank[er][esr],
                               isProppantBoundaryElement[er][esr],
                               isProppantMobile[er][esr],
                               cellBasedFlux[er][esr],
                               conc[er][esr],
                               proppantPackVolFrac[er][esr],
                               proppantExcessPackVolume[er][esr],
                               proppantLiftFlux[er][esr] );
  } );
}

GEOS_HOST_DEVICE
inline
void
ProppantPackVolumeKernel::
  updateProppantPackVolume( localIndex const numElems,
                            arraySlice1d< localIndex const > const & stencilElementIndices,
                            arraySlice1d< real64 const > const & stencilWeights,
                            arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                            R1Tensor const & unitGravityVector,
                            real64 const maxProppantConcentration,
                            arrayView1d< integer const > const & isProppantMobile,
                            arrayView1d< real64 const > const & proppantExcessPackVolume,
                            arrayView1d< real64 > const & conc,
                            arrayView1d< real64 > const & proppantPackVolFrac )
{
  integer faceIndex = -1;
  real64 excessVolume = 0.0;

  real64 constexpr TINY = 1e-10;

  real64 const stencilCellToEdgeDistance0 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[0] );
  real64 const edgeLength = stencilWeights[0] * stencilCellToEdgeDistance0;

  if( numElems == 2 )
  {
    localIndex const ei0  = stencilElementIndices[0];
    localIndex const ei1  = stencilElementIndices[1];

    real64 const stencilCellToEdgeDistance1 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[1] );

    real64 const stencilEdgeToFaceDownDistance0 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance0;

    real64 const stencilEdgeToFaceDownDistance1 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[1], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance1;

    //0: top  1: bottom

    if( stencilEdgeToFaceDownDistance0 < -TINY && stencilEdgeToFaceDownDistance1 > TINY && isProppantMobile[ei0] == 1 && isProppantMobile[ei1] == 0 &&
        proppantExcessPackVolume[ei1] > 0.0 )
    {
      faceIndex = 0;
      excessVolume = proppantExcessPackVolume[ei1];
    }

    //0: bottom  1: top

    if( stencilEdgeToFaceDownDistance0 > TINY && stencilEdgeToFaceDownDistance1 < -TINY && isProppantMobile[ei1] == 1 && isProppantMobile[ei0] == 0 &&
        proppantExcessPackVolume[ei0] > 0.0 )
    {
      faceIndex = 1;
      excessVolume = proppantExcessPackVolume[ei0];
    }
  }

  if( faceIndex >= 0 )
  {
    localIndex const ei = stencilElementIndices[faceIndex];
    real64 const L = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[faceIndex] ) * 2.0;
    real64 const volFrac = proppantPackVolFrac[ei] + excessVolume / L;

    if( volFrac >= 1.0 )
    {
      proppantPackVolFrac[ei] = 1.0;
      conc[ei] = maxProppantConcentration;
    }
    else
    {
      proppantPackVolFrac[ei] = volFrac;
    }
  }
}


void ProppantPackVolumeKernel::
  launchProppantPackVolumeUpdate( SurfaceElementStencil const & stencil,
                                  R1Tensor const & unitGravityVector,
                                  real64 const maxProppantConcentration,
                                  ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                  ElementView< arrayView1d< real64 const > > const & proppantExcessPackVolume,
                                  ElementView< arrayView1d< real64 > > const & conc,
                                  ElementView< arrayView1d< real64 > > const & proppantPackVolFrac )
{

  typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename SurfaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = seri.sizeOfArray( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    updateProppantPackVolume( numFluxElems,
                              sei[iconn],
                              weights[iconn],
                              cellCenterToEdgeCenters[iconn],
                              unitGravityVector,
                              maxProppantConcentration,
                              isProppantMobile[er][esr],
                              proppantExcessPackVolume[er][esr],
                              conc[er][esr],
                              proppantPackVolFrac[er][esr] );
  } );
}

} // namespace proppantTransportKernels

} // namespace geos
