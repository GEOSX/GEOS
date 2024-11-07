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
 * @file IsothermalCompositionalMultiphaseFVMKernels.cpp
 */

#include "IsothermalCompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"

#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/SurfaceElementStencil.hpp"
#include "finiteVolume/EmbeddedSurfaceToCellStencil.hpp"
#include "finiteVolume/FaceElementToCellStencil.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"

namespace geos
{
using namespace constitutive;

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** FaceBasedAssemblyKernel ********************************/

FaceBasedAssemblyKernelBase::FaceBasedAssemblyKernelBase( integer const numPhases,
                                                          globalIndex const rankOffset,
                                                          DofNumberAccessor const & dofNumberAccessor,
                                                          CompFlowAccessors const & compFlowAccessors,
                                                          MultiFluidAccessors const & multiFluidAccessors,
                                                          real64 const dt,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs,
                                                          BitFlags< FaceBasedAssemblyKernelFlags > kernelFlags )
  : m_numPhases( numPhases ),
  m_rankOffset( rankOffset ),
  m_dt( dt ),
  m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
  m_ghostRank( compFlowAccessors.get( fields::ghostRank {} ) ),
  m_gravCoef( compFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
  m_pres( compFlowAccessors.get( fields::flow::pressure {} ) ),
  m_dCompFrac_dCompDens( compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity {} ) ),
  m_dPhaseVolFrac( compFlowAccessors.get( fields::flow::dPhaseVolumeFraction {} ) ),
  m_phaseCompFrac( multiFluidAccessors.get( fields::multifluid::phaseCompFraction {} ) ),
  m_dPhaseCompFrac( multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction {} ) ),
  m_localMatrix( localMatrix ),
  m_localRhs( localRhs ),
  m_kernelFlags( kernelFlags )
{}

/******************************** CFLFluxKernel ********************************/

template< integer NC, localIndex NUM_ELEMS, localIndex maxStencilSize >
GEOS_HOST_DEVICE
inline
void
CFLFluxKernel::
  compute( integer const numPhases,
           localIndex const stencilSize,
           real64 const dt,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           real64 const (&transmissibility)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & phaseRelPerm,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseVisc,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
           ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux )
{
  // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
  for( integer ip = 0; ip < numPhases; ++ip )
  {

    // clear working arrays
    real64 densMean{};

    // create local work arrays
    real64 presGrad{};
    real64 gravHead{};

    // calculate quantities on primary connected cells
    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      // average density across the face
      densMean += 0.5 * phaseMassDens[er][esr][ei][0][ip];
    }

    //***** calculation of phase volumetric flux *****

    // compute potential difference MPFA-style
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      presGrad += transmissibility[i] * pres[er][esr][ei];
      gravHead += transmissibility[i] * densMean * gravCoef[er][esr][ei];
    }

    // *** upwinding ***

    // compute phase potential gradient
    real64 const potGrad = presGrad - gravHead;

    // choose upstream cell
    localIndex const k_up = (potGrad >= 0) ? 0 : 1;

    localIndex const er_up  = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up  = sei[k_up];

    // compute the phase flux only if the phase is present
    bool const phaseExists = (phaseVolFrac[er_up][esr_up][ei_up][ip] > 0);
    if( !phaseExists )
    {
      continue;
    }

    real64 const mobility = phaseRelPerm[er_up][esr_up][ei_up][0][ip] / phaseVisc[er_up][esr_up][ei_up][0][ip];

    // increment the phase (volumetric) outflux of the upstream cell
    real64 const absPhaseFlux = LvArray::math::abs( dt * mobility * potGrad );
    RAJA::atomicAdd( parallelDeviceAtomic{}, &phaseOutflux[er_up][esr_up][ei_up][ip], absPhaseFlux );

    // increment the component (mass/molar) outflux of the upstream cell
    for( integer ic = 0; ic < NC; ++ic )
    {
      real64 const absCompFlux = phaseCompFrac[er_up][esr_up][ei_up][0][ip][ic]
                                 * phaseDens[er_up][esr_up][ei_up][0][ip]
                                 * absPhaseFlux;
      RAJA::atomicAdd( parallelDeviceAtomic{}, &compOutflux[er_up][esr_up][ei_up][ic], absCompFlux );
    }
  }
}

template< integer NC, typename STENCILWRAPPER_TYPE >
void
CFLFluxKernel::
  launch( integer const numPhases,
          real64 const dt,
          STENCILWRAPPER_TYPE const & stencilWrapper,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & phaseRelPerm,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseVisc,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux,
          ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux )
{
  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  localIndex constexpr numElems = STENCILWRAPPER_TYPE::maxNumPointsInFlux;
  localIndex constexpr maxStencilSize = STENCILWRAPPER_TYPE::maxStencilSize;

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    // compute transmissibility
    real64 transmissibility[STENCILWRAPPER_TYPE::maxNumConnections][2];
    real64 dTrans_dPres[STENCILWRAPPER_TYPE::maxNumConnections][2];

    stencilWrapper.computeWeights( iconn,
                                   permeability,
                                   dPerm_dPres,
                                   transmissibility,
                                   dTrans_dPres );

    CFLFluxKernel::compute< NC, numElems, maxStencilSize >( numPhases,
                                                            sei[iconn].size(),
                                                            dt,
                                                            seri[iconn],
                                                            sesri[iconn],
                                                            sei[iconn],
                                                            transmissibility[0],
                                                            pres,
                                                            gravCoef,
                                                            phaseVolFrac,
                                                            phaseRelPerm,
                                                            phaseVisc,
                                                            phaseDens,
                                                            phaseMassDens,
                                                            phaseCompFrac,
                                                            phaseOutflux,
                                                            compOutflux );
  } );
}

#define INST_CFLFluxKernel( NC, STENCILWRAPPER_TYPE ) \
  template \
  void CFLFluxKernel:: \
    launch< NC, STENCILWRAPPER_TYPE >( integer const numPhases, \
                                       real64 const dt, \
                                       STENCILWRAPPER_TYPE const & stencil, \
                                       ElementViewConst< arrayView1d< real64 const > > const & pres, \
                                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef, \
                                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac, \
                                       ElementViewConst< arrayView3d< real64 const > > const & permeability, \
                                       ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres, \
                                       ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & phaseRelPerm, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseVisc, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                       ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                       ElementView< arrayView2d< real64, compflow::USD_PHASE > > const & phaseOutflux, \
                                       ElementView< arrayView2d< real64, compflow::USD_COMP > > const & compOutflux )

INST_CFLFluxKernel( 1, CellElementStencilTPFAWrapper );
INST_CFLFluxKernel( 2, CellElementStencilTPFAWrapper );
INST_CFLFluxKernel( 3, CellElementStencilTPFAWrapper );
INST_CFLFluxKernel( 4, CellElementStencilTPFAWrapper );
INST_CFLFluxKernel( 5, CellElementStencilTPFAWrapper );

INST_CFLFluxKernel( 1, SurfaceElementStencilWrapper );
INST_CFLFluxKernel( 2, SurfaceElementStencilWrapper );
INST_CFLFluxKernel( 3, SurfaceElementStencilWrapper );
INST_CFLFluxKernel( 4, SurfaceElementStencilWrapper );
INST_CFLFluxKernel( 5, SurfaceElementStencilWrapper );

INST_CFLFluxKernel( 1, EmbeddedSurfaceToCellStencilWrapper );
INST_CFLFluxKernel( 2, EmbeddedSurfaceToCellStencilWrapper );
INST_CFLFluxKernel( 3, EmbeddedSurfaceToCellStencilWrapper );
INST_CFLFluxKernel( 4, EmbeddedSurfaceToCellStencilWrapper );
INST_CFLFluxKernel( 5, EmbeddedSurfaceToCellStencilWrapper );

INST_CFLFluxKernel( 1, FaceElementToCellStencilWrapper );
INST_CFLFluxKernel( 2, FaceElementToCellStencilWrapper );
INST_CFLFluxKernel( 3, FaceElementToCellStencilWrapper );
INST_CFLFluxKernel( 4, FaceElementToCellStencilWrapper );
INST_CFLFluxKernel( 5, FaceElementToCellStencilWrapper );

#undef INST_CFLFluxKernel

/******************************** CFLKernel ********************************/

template< integer NP >
GEOS_HOST_DEVICE
void
CFLKernel::
  computePhaseCFL( real64 const poreVol,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
                   arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > phaseRelPerm,
                   arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > dPhaseRelPerm_dPhaseVolFrac,
                   arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseVisc,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseOutflux,
                   real64 & phaseCFLNumber )
{
  // first, check which phases are mobile in the cell
  real64 mob[NP]{};
  localIndex mobilePhases[NP]{};
  localIndex numMobilePhases = 0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    if( phaseVolFrac[ip] > 0 )
    {
      mob[ip] = phaseRelPerm[ip] / phaseVisc[ip];
      if( mob[ip] > minPhaseMobility )
      {
        mobilePhases[numMobilePhases] = ip;
        numMobilePhases++;
      }
    }
  }

  // then, depending on the regime, apply the appropriate CFL formula
  phaseCFLNumber = 0;

  // single-phase flow regime
  if( numMobilePhases == 1 )
  {
    phaseCFLNumber = phaseOutflux[mobilePhases[0]] / poreVol;
  }
  // two-phase flow regime
  else if( numMobilePhases == 2 )
  {
    // from Hui Cao's PhD thesis
    localIndex const ip0 = mobilePhases[0];
    localIndex const ip1 = mobilePhases[1];
    real64 const dMob_dVolFrac[2] = { dPhaseRelPerm_dPhaseVolFrac[ip0][ip0] / phaseVisc[ip0],
                                      -dPhaseRelPerm_dPhaseVolFrac[ip1][ip1] / phaseVisc[ip1] }; // using S0 = 1 - S1
    real64 const denom = 1. / ( poreVol * ( mob[ip0] + mob[ip1] ) );
    real64 const coef0 = denom * mob[ip1] / mob[ip0] * dMob_dVolFrac[ip0];
    real64 const coef1 = -denom * mob[ip0] / mob[ip1] * dMob_dVolFrac[ip1];

    phaseCFLNumber = LvArray::math::abs( coef0*phaseOutflux[ip0] + coef1*phaseOutflux[ip1] );
  }
  // three-phase flow regime
  else if( numMobilePhases == 3 )
  {
    // from Keith Coats, IMPES stability: Selection of stable timesteps (2003)
    real64 totalMob = 0.0;
    for( integer ip = 0; ip < numMobilePhases; ++ip )
    {
      totalMob += mob[ip];
    }

    real64 f[2][2]{};
    for( integer i = 0; i < 2; ++i )
    {
      for( integer j = 0; j < 2; ++j )
      {
        f[i][j]  = ( i == j )*totalMob - mob[i];
        f[i][j] /= (totalMob * mob[j]);
        real64 sum = 0;
        for( integer k = 0; k < 3; ++k )
        {
          sum += dPhaseRelPerm_dPhaseVolFrac[k][j] / phaseVisc[k]
                 * phaseOutflux[j];
        }
        f[i][j] *= sum;
      }
    }
    phaseCFLNumber = f[0][0] + f[1][1];
    phaseCFLNumber += sqrt( phaseCFLNumber*phaseCFLNumber - 4 * ( f[0][0]*f[1][1] - f[1][0]*f[0][1] ) );
    phaseCFLNumber = 0.5 * LvArray::math::abs( phaseCFLNumber ) / poreVol;
  }
}


template< integer NC >
GEOS_HOST_DEVICE
void
CFLKernel::
  computeCompCFL( real64 const poreVol,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compFrac,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compOutflux,
                  real64 & compCFLNumber )
{


  compCFLNumber = 0.0;
  for( integer ic = 0; ic < NC; ++ic )
  {
    if( compFrac[ic] > minComponentFraction )
    {
      real64 const compMoles = compDens[ic] * poreVol;
      real64 const CFL = compOutflux[ic] / compMoles;
      if( CFL > compCFLNumber )
      {
        compCFLNumber = CFL;
      }
    }
  }
}

template< integer NC, integer NP >
void
CFLKernel::
  launch( localIndex const size,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux,
          arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux,
          arrayView1d< real64 > const & phaseCFLNumber,
          arrayView1d< real64 > const & compCFLNumber,
          real64 & maxPhaseCFLNumber,
          real64 & maxCompCFLNumber )
{
  RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionPhaseCFLNumber( 0.0 );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionCompCFLNumber( 0.0 );

  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const poreVol = volume[ei] * porosity[ei][0];

    // phase CFL number
    real64 cellPhaseCFLNumber = 0.0;
    computePhaseCFL< NP >( poreVol,
                           phaseVolFrac[ei],
                           phaseRelPerm[ei][0],
                           dPhaseRelPerm_dPhaseVolFrac[ei][0],
                           phaseVisc[ei][0],
                           phaseOutflux[ei],
                           cellPhaseCFLNumber );
    subRegionPhaseCFLNumber.max( cellPhaseCFLNumber );
    phaseCFLNumber[ei] = cellPhaseCFLNumber;

    // component CFL number
    real64 cellCompCFLNumber = 0.0;
    computeCompCFL< NC >( poreVol,
                          compDens[ei],
                          compFrac[ei],
                          compOutflux[ei],
                          cellCompCFLNumber );
    subRegionCompCFLNumber.max( cellCompCFLNumber );
    compCFLNumber[ei] = cellCompCFLNumber;
  } );

  maxPhaseCFLNumber = subRegionPhaseCFLNumber.get();
  maxCompCFLNumber = subRegionCompCFLNumber.get();
}

#define INST_CFLKernel( NC, NP ) \
  template \
  void CFLKernel:: \
    launch< NC, NP >( localIndex const size, \
                      arrayView1d< real64 const > const & volume, \
                      arrayView2d< real64 const > const & porosity, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & compDens, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & compFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac, \
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm, \
                      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux, \
                      arrayView1d< real64 > const & phaseCFLNumber, \
                      arrayView1d< real64 > const & compCFLNumber, \
                      real64 & maxPhaseCFLNumber, \
                      real64 & maxCompCFLNumber )
INST_CFLKernel( 1, 2 );
INST_CFLKernel( 2, 2 );
INST_CFLKernel( 3, 2 );
INST_CFLKernel( 4, 2 );
INST_CFLKernel( 5, 2 );

INST_CFLKernel( 1, 3 );
INST_CFLKernel( 2, 3 );
INST_CFLKernel( 3, 3 );
INST_CFLKernel( 4, 3 );
INST_CFLKernel( 5, 3 );

#undef INST_CFLKernel

/******************************** AquiferBCKernel ********************************/

template< integer NC >
GEOS_HOST_DEVICE
void
AquiferBCKernel::
  compute( integer const numPhases,
           integer const ipWater,
           bool const allowAllPhasesIntoAquifer,
           real64 const aquiferVolFlux,
           real64 const dAquiferVolFlux_dPres,
           real64 const aquiferWaterPhaseDens,
           arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac,
           arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens,
           real64 const dt,
           real64 (& localFlux)[NC],
           real64 (& localFluxJacobian)[NC][NC+1] )
{
  using Deriv = multifluid::DerivativeOffset;

  real64 dProp_dC[NC]{};
  real64 dPhaseFlux_dCompDens[NC]{};

  if( aquiferVolFlux > 0 ) // aquifer is upstream
  {
    // in this case, we assume that:
    //    - only the water phase is present in the aquifer
    //    - the aquifer water phase composition is constant

    for( integer ic = 0; ic < NC; ++ic )
    {
      real64 const phaseFlux = aquiferVolFlux * aquiferWaterPhaseDens;
      localFlux[ic] -= dt * phaseFlux * aquiferWaterPhaseCompFrac[ic];
      localFluxJacobian[ic][0] -= dt * dAquiferVolFlux_dPres * aquiferWaterPhaseDens * aquiferWaterPhaseCompFrac[ic];
    }
  }
  else // reservoir is upstream
  {
    for( integer ip = 0; ip < numPhases; ++ip )
    {

      // Why two options below:
      //   - The aquifer model assumes single-phase water flow, so ideally, we should only allow water phase flow from the reservoir to the
      // aquifer
      //   - But, if/when the CO2 plume reaches the reservoir cell connected to the aquifer and saturates it, the aquifer flux becomes zero
      //     if we don't let some CO2 go into the aquifer

      if( ip == ipWater || allowAllPhasesIntoAquifer )
      {
        real64 const phaseDensVolFrac = phaseDens[ip] * phaseVolFrac[ip];
        real64 const phaseFlux = aquiferVolFlux * phaseDensVolFrac;
        real64 const dPhaseFlux_dPres = dAquiferVolFlux_dPres * phaseDensVolFrac
                                        + aquiferVolFlux * ( dPhaseDens[ip][Deriv::dP] * phaseVolFrac[ip] + phaseDens[ip] * dPhaseVolFrac[ip][Deriv::dP] );

        applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens[ip], dProp_dC, Deriv::dC );
        for( integer ic = 0; ic < NC; ++ic )
        {
          dPhaseFlux_dCompDens[ic] = aquiferVolFlux * ( dProp_dC[ic] * phaseVolFrac[ip] + phaseDens[ip] * dPhaseVolFrac[ip][Deriv::dC+ic] );
        }

        for( integer ic = 0; ic < NC; ++ic )
        {
          localFlux[ic] -= dt * phaseFlux * phaseCompFrac[ip][ic];
          localFluxJacobian[ic][0] -= dt * ( dPhaseFlux_dPres * phaseCompFrac[ip][ic] + phaseFlux * dPhaseCompFrac[ip][ic][Deriv::dP] );

          applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac[ip][ic], dProp_dC, Deriv::dC );
          for( integer jc = 0; jc < NC; ++jc )
          {
            localFluxJacobian[ic][jc+1] -= dt * ( dPhaseFlux_dCompDens[jc] * phaseCompFrac[ip][ic] + phaseFlux * dProp_dC[jc] );
          }
        }
      }
    }
  }
}

template< integer NC >
void
AquiferBCKernel::
  launch( integer const numPhases,
          integer const ipWater,
          bool const allowAllPhasesIntoAquifer,
          integer const useTotalMassEquation,
          BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const aquiferWaterPhaseDens,
          arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & presOld,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
          real64 const timeAtBeginningOfStep,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{

  using namespace compositionalMultiphaseUtilities;
  using Order = BoundaryStencil::Order;

  BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
  BoundaryStencil::WeightContainerViewConstType const & weight = stencil.getWeights();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
  {
    constexpr integer NDOF = NC + 1;

    // working arrays
    globalIndex dofColIndices[NDOF]{};
    real64 localFlux[NC]{};
    real64 localFluxJacobian[NC][NDOF]{};

    localIndex const er  = seri( iconn, Order::ELEM );
    localIndex const esr = sesri( iconn, Order::ELEM );
    localIndex const ei  = sefi( iconn, Order::ELEM );
    real64 const areaFraction = weight( iconn, Order::ELEM );

    // compute the aquifer influx rate using the pressure influence function and the aquifer props
    real64 dAquiferVolFlux_dPres = 0.0;
    real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                            dt,
                                                            pres[er][esr][ei],
                                                            presOld[er][esr][ei],
                                                            gravCoef[er][esr][ei],
                                                            areaFraction,
                                                            dAquiferVolFlux_dPres );

    // compute the phase/component aquifer flux
    AquiferBCKernel::compute< NC >( numPhases,
                                    ipWater,
                                    allowAllPhasesIntoAquifer,
                                    aquiferVolFlux,
                                    dAquiferVolFlux_dPres,
                                    aquiferWaterPhaseDens,
                                    aquiferWaterPhaseCompFrac,
                                    phaseDens[er][esr][ei][0],
                                    dPhaseDens[er][esr][ei][0],
                                    phaseVolFrac[er][esr][ei],
                                    dPhaseVolFrac[er][esr][ei],
                                    phaseCompFrac[er][esr][ei][0],
                                    dPhaseCompFrac[er][esr][ei][0],
                                    dCompFrac_dCompDens[er][esr][ei],
                                    dt,
                                    localFlux,
                                    localFluxJacobian );

    // populate dof indices
    globalIndex const offset = dofNumber[er][esr][ei];
    for( integer jdof = 0; jdof < NDOF; ++jdof )
    {
      dofColIndices[jdof] = offset + jdof;
    }

    if( useTotalMassEquation > 0 )
    {
      // Apply equation/variable change transformation(s)
      real64 work[NDOF];
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NDOF, localFluxJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, localFlux );
    }

    // Add to residual/jacobian
    if( ghostRank[er][esr][ei] < 0 )
    {
      globalIndex const globalRow = dofNumber[er][esr][ei];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
      GEOS_ASSERT_GE( localRow, 0 );
      GEOS_ASSERT_GT( localMatrix.numRows(), localRow + NC );

      for( integer ic = 0; ic < NC; ++ic )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + ic], localFlux[ic] );
        localMatrix.addToRow< parallelDeviceAtomic >( localRow + ic,
                                                      dofColIndices,
                                                      localFluxJacobian[ic],
                                                      NDOF );
      }
    }
  } );
}

#define INST_AquiferBCKernel( NC ) \
  template \
  void AquiferBCKernel:: \
    launch< NC >( integer const numPhases, \
                  integer const ipWater, \
                  bool const allowAllPhasesIntoAquifer, \
                  integer const useTotalMassEquation, \
                  BoundaryStencil const & stencil, \
                  globalIndex const rankOffset, \
                  ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber, \
                  AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper, \
                  real64 const aquiferWaterPhaseDens, \
                  arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac, \
                  ElementViewConst< arrayView1d< integer const > > const & ghostRank, \
                  ElementViewConst< arrayView1d< real64 const > > const & pres, \
                  ElementViewConst< arrayView1d< real64 const > > const & dPres, \
                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef, \
                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac, \
                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac, \
                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                  ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                  real64 const timeAtBeginningOfStep, \
                  real64 const dt, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_AquiferBCKernel( 1 );
INST_AquiferBCKernel( 2 );
INST_AquiferBCKernel( 3 );
INST_AquiferBCKernel( 4 );
INST_AquiferBCKernel( 5 );

#undef INST_AquiferBCKernel

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos
