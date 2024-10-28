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
 * @file IHU2PhaseFlux.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHU2PHASEFLUX_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHU2PHASEFLUX_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernelUtilities
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

/*** HU 2 phase simplified version ***/

struct HU2PhaseFlux
{

  using UPWIND_SCHEME = HybridUpwind;

  /**
   * @brief Form the Implicit Hybrid Upwind from pressure gradient and gravitational head
   * @tparam numComp number of components
   * @tparam numFluxSupportPoints number of flux support points
   * @param numPhase number of phases
   * @param ip phase index
   * @param hasCapPressure flag indicating if there is capillary pressure
   * @param seri arraySlice of the stencil-implied element region index
   * @param sesri arraySlice of the stencil-implied element subregion index
   * @param sei arraySlice of the stencil-implied element index
   * @param trans transmissibility at the connection
   * @param dTrans_dPres derivative of transmissibility wrt pressure
   * @param pres pressure
   * @param gravCoef gravitational coefficient
   * @param phaseMob phase mobility
   * @param dPhaseMob derivative of phase mobility wrt pressure, temperature, comp density
   * @param dPhaseVolFrac derivative of phase volume fraction wrt pressure, temperature, comp density
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseMassDens phase mass density
   * @param dPhaseMassDens derivative of phase mass density wrt pressure, temperature, comp fraction
   * @param phaseCapPressure phase capillary pressure
   * @param dPhaseCapPressure_dPhaseVolFrac derivative of phase capillary pressure wrt phase volume fraction
   * @param k_up uptream index for this phase
   * @param potGrad potential gradient for this phase
   * @param phaseFlux phase flux
   * @param dPhaseFlux_dP derivative of phase flux wrt pressure
   * @param dPhaseFlux_dC derivative of phase flux wrt comp density
   */
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( integer const numPhase,
           integer const ip,
           integer const hasCapPressure,
           localIndex const ( &seri )[numFluxSupportPoints],
           localIndex const ( &sesri )[numFluxSupportPoints],
           localIndex const ( &sei )[numFluxSupportPoints],
           real64 const ( &trans )[2],
           real64 const ( &dTrans_dPres )[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
           real64 & GEOS_UNUSED_PARAM( potGrad ),
           real64 ( &phaseFlux ),
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp] )
  {
    // viscous part
    computeViscousFlux< numComp, numFluxSupportPoints >( ip, numPhase, hasCapPressure,
                                                         seri, sesri, sei,
                                                         trans, dTrans_dPres,
                                                         pres, gravCoef,
                                                         dCompFrac_dCompDens,
                                                         phaseMassDens, dPhaseMassDens,
                                                         phaseMob, dPhaseMob,
                                                         dPhaseVolFrac,
                                                         phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                                                         phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );
    // gravity part
    computeGravityFlux< numComp, numFluxSupportPoints >( ip, numPhase,
                                                         seri, sesri, sei,
                                                         trans, dTrans_dPres, gravCoef,
                                                         phaseMob, dPhaseMob,
                                                         dCompFrac_dCompDens,
                                                         phaseMassDens, dPhaseMassDens,
                                                         phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    // capillary part
    if( hasCapPressure )
    {
      computeCapillaryFlux< numComp, numFluxSupportPoints >( ip, numPhase, hasCapPressure,
                                                             seri, sesri, sei,
                                                             trans, dTrans_dPres,
                                                             pres, gravCoef,
                                                             phaseMob, dPhaseMob,
                                                             dPhaseVolFrac,
                                                             dCompFrac_dCompDens,
                                                             phaseMassDens, dPhaseMassDens,
                                                             phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                                                             phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );
    }
  }

protected:

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  computeViscousFlux( const integer & ip, const integer & numPhase, const integer & hasCapPressure,
                      const localIndex (& seri)[numFluxSupportPoints],
                      const localIndex (& sesri)[numFluxSupportPoints],
                      const localIndex (& sei)[numFluxSupportPoints],
                      const real64 (& trans)[2], const real64 (& dTrans_dPres)[2],
                      const ElementViewConst< arrayView1d< const real64 > > & pres,
                      const ElementViewConst< arrayView1d< const real64 > > & gravCoef,
                      ElementViewConst< arrayView3d< const real64 > > const & dCompFrac_dCompDens,
                      const ElementViewConst< arrayView3d< const real64 > > & phaseMassDens,
                      const ElementViewConst< arrayView4d< const real64 > > & dPhaseMassDens,
                      const ElementViewConst< arrayView2d< const real64, 1 > > & phaseMob,
                      const ElementViewConst< arrayView3d< const real64 > > & dPhaseMob,
                      const ElementViewConst< arrayView3d< const real64 > > & dPhaseVolFrac,
                      const ElementViewConst< arrayView3d< const real64 > > & phaseCapPressure,
                      const ElementViewConst< arrayView4d< const real64 > > & dPhaseCapPressure_dPhaseVolFrac,
                      real64 ( &phaseFlux ),
                      real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
                      real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp] )
  {
    // form total velocity and derivatives (TODO move it OUT!)
    real64 totFlux{};
    real64 dTotFlux_dP[numFluxSupportPoints]{};
    real64 dTotFlux_dC[numFluxSupportPoints][numComp]{};

    computeTotalFlux( numPhase, hasCapPressure,
                      seri, sesri, sei,
                      trans, dTrans_dPres,
                      pres, gravCoef,
                      phaseMob, dPhaseMob,
                      dPhaseVolFrac,
                      dCompFrac_dCompDens,
                      phaseMassDens, dPhaseMassDens,
                      phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                      totFlux, dTotFlux_dP, dTotFlux_dC );

    localIndex k_up = -1;
    real64 mob{};
    real64 dMob_dP{};
    real64 dMob_dC[numComp]{};

    real64 totMob{};
    real64 dTotMob_dP{};
    real64 dTotMob_dC[numComp]{};

    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( jp == ip ) // ip will be computed later
        continue;

      // upwind based on totFlux sign
      upwindMobility< numComp, numFluxSupportPoints >( jp,
                                                       seri,
                                                       sesri,
                                                       sei,
                                                       totFlux,
                                                       phaseMob,
                                                       dPhaseMob,
                                                       k_up,
                                                       mob,
                                                       dMob_dP,
                                                       dMob_dC );
      totMob += mob;
      dTotMob_dP += dMob_dP;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dTotMob_dC[jc] += dMob_dC[jc];
      }
    }

    // upwind based on totFlux sign
    upwindMobility< numComp, numFluxSupportPoints >( ip,
                                                     seri,
                                                     sesri,
                                                     sei,
                                                     totFlux,
                                                     phaseMob,
                                                     dPhaseMob,
                                                     k_up,
                                                     mob,
                                                     dMob_dP,
                                                     dMob_dC );
    totMob += mob;
    dTotMob_dP += dMob_dP;
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dTotMob_dC[jc] += dMob_dC[jc];
    }

    totMob = LvArray::math::max( totMob, minTotMob );
    real64 const invTotMob = 1 / totMob;

    // fractional flow for viscous part as \lambda_i^{up}/\sum_{NP}(\lambda_j^{up})

    real64 const fractionalFlow = mob * invTotMob;
    real64 const dFractionalFlow_dP = (dMob_dP - fractionalFlow * dTotMob_dP) * invTotMob;
    real64 dFractionalFlow_dC[numComp]{};
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dFractionalFlow_dC[jc] = (dMob_dC[jc] - fractionalFlow * dTotMob_dC[jc]) * invTotMob;
    }

    /// Assembling the viscous flux (and derivatives) from fractional flow and total velocity as \phi_{\mu} = f_i^{up,\mu} uT

    real64 const viscousPhaseFlux = fractionalFlow * totFlux;
    real64 dViscousPhaseFlux_dP[numFluxSupportPoints]{};
    real64 dViscousPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    // fractionalFlow derivatives
    dViscousPhaseFlux_dP[k_up] += dFractionalFlow_dP * totFlux;
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dViscousPhaseFlux_dC[k_up][jc] += dFractionalFlow_dC[jc] * totFlux;
    }

    // Ut derivatives
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dViscousPhaseFlux_dP[ke] += fractionalFlow * dTotFlux_dP[ke];
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dViscousPhaseFlux_dC[ke][jc] += fractionalFlow * dTotFlux_dC[ke][jc];
      }
    }

    // accumulate in the flux and its derivatives (need to be very careful doing that)
    phaseFlux += viscousPhaseFlux;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += dViscousPhaseFlux_dP[ke];
      for( localIndex ic = 0; ic < numComp; ++ic )
        dPhaseFlux_dC[ke][ic] += dViscousPhaseFlux_dC[ke][ic];
    }
  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  computeGravityFlux( const integer & ip, const integer & numPhase,
                      const localIndex (& seri)[numFluxSupportPoints],
                      const localIndex (& sesri)[numFluxSupportPoints],
                      const localIndex (& sei)[numFluxSupportPoints],
                      const real64 (& trans)[2], const real64 (& dTrans_dPres)[2],
                      const ElementViewConst< arrayView1d< const real64 > > & gravCoef,
                      const ElementViewConst< arrayView2d< const real64, 1 > > & phaseMob,
                      const ElementViewConst< arrayView3d< const real64 > > & dPhaseMob,
                      ElementViewConst< arrayView3d< const real64 > > const & dCompFrac_dCompDens,
                      const ElementViewConst< arrayView3d< const real64 > > & phaseMassDens,
                      const ElementViewConst< arrayView4d< const real64 > > & dPhaseMassDens,
                      real64 & phaseFlux,
                      real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                      real64 (& dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
  {
    /// Assembling the gravitational flux (and derivatives)
    real64 gravPhaseFlux{};
    real64 dGravPhaseFlux_dP[numFluxSupportPoints]{};
    real64 dGravPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    real64 pot_i{};
    real64 dPot_i_dP[numFluxSupportPoints]{};
    real64 dPot_i_dC[numFluxSupportPoints][numComp]{};
    computeGravityPotential< numComp, numFluxSupportPoints >( ip,
                                                              seri,
                                                              sesri,
                                                              sei,
                                                              trans,
                                                              dTrans_dPres,
                                                              gravCoef,
                                                              dCompFrac_dCompDens,
                                                              phaseMassDens,
                                                              dPhaseMassDens,
                                                              pot_i,
                                                              dPot_i_dP,
                                                              dPot_i_dC );

    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( ip != jp )
      {
        real64 pot_j{};
        real64 dPot_j_dP[numFluxSupportPoints]{};
        real64 dPot_j_dC[numFluxSupportPoints][numComp]{};
        computeGravityPotential< numComp, numFluxSupportPoints >( jp,
                                                                  seri,
                                                                  sesri,
                                                                  sei,
                                                                  trans,
                                                                  dTrans_dPres,
                                                                  gravCoef,
                                                                  dCompFrac_dCompDens,
                                                                  phaseMassDens,
                                                                  dPhaseMassDens,
                                                                  pot_j,
                                                                  dPot_j_dP,
                                                                  dPot_j_dC );

        // upwind based on pot diff sign
        real64 const potDiff = pot_j - pot_i;
        real64 dPotDiff_dP[numFluxSupportPoints]{};
        real64 dPotDiff_dC[numFluxSupportPoints][numComp]{};
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dPotDiff_dP[ke] += dPot_j_dP[ke] - dPot_i_dP[ke];
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPotDiff_dC[ke][jc] += dPot_j_dC[ke][jc] - dPot_i_dC[ke][jc];
          }
        }

        localIndex k_up_i = -1;
        real64 mob_i{};
        real64 dMob_i_dP{};
        real64 dMob_i_dC[numComp]{};
        upwindMobility< numComp, numFluxSupportPoints >( ip,
                                                         seri,
                                                         sesri,
                                                         sei,
                                                         potDiff,
                                                         phaseMob,
                                                         dPhaseMob,
                                                         k_up_i,
                                                         mob_i,
                                                         dMob_i_dP,
                                                         dMob_i_dC );
        localIndex k_up_j = -1;
        real64 mob_j{};
        real64 dMob_j_dP{};
        real64 dMob_j_dC[numComp]{};
        upwindMobility< numComp, numFluxSupportPoints >( jp,
                                                         seri,
                                                         sesri,
                                                         sei,
                                                         -potDiff,
                                                         phaseMob,
                                                         dPhaseMob,
                                                         k_up_j,
                                                         mob_j,
                                                         dMob_j_dP,
                                                         dMob_j_dC );

        real64 const mobTot = LvArray::math::max( mob_i + mob_j, minTotMob );
        real64 const mobTotInv = 1 / mobTot;
        real64 dMobTot_dP[numFluxSupportPoints]{};
        real64 dMobTot_dC[numFluxSupportPoints][numComp]{};
        dMobTot_dP[k_up_i] += dMob_i_dP;
        dMobTot_dP[k_up_j] += dMob_j_dP;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dMobTot_dC[k_up_i][jc] += dMob_i_dC[jc];
          dMobTot_dC[k_up_j][jc] += dMob_j_dC[jc];
        }

        // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
        gravPhaseFlux += mob_i * mob_j * mobTotInv * potDiff;

        // mob_i derivatives
        dGravPhaseFlux_dP[k_up_i] += dMob_i_dP * mob_j * mobTotInv * potDiff;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dGravPhaseFlux_dC[k_up_i][jc] += dMob_i_dC[jc] * mob_j * mobTotInv * potDiff;
        }

        // mob_j derivatives
        dGravPhaseFlux_dP[k_up_j] += mob_i * dMob_j_dP * mobTotInv * potDiff;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dGravPhaseFlux_dC[k_up_j][jc] += mob_i * dMob_j_dC[jc] * mobTotInv * potDiff;
        }

        // mobTot derivatives
        real64 const mobTotInv2 = mobTotInv * mobTotInv;
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dGravPhaseFlux_dP[ke] -= mob_i * mob_j * dMobTot_dP[ke] * mobTotInv2 * potDiff;
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dGravPhaseFlux_dC[ke][jc] -= mob_i * mob_j * dMobTot_dC[ke][jc] * mobTotInv2 * potDiff;
          }
        }

        // potDiff derivatives
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dGravPhaseFlux_dP[ke] += mob_i * mob_j * mobTotInv * dPotDiff_dP[ke];
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dGravPhaseFlux_dC[ke][jc] += mob_i * mob_j * mobTotInv * dPotDiff_dC[ke][jc];
          }
        }
      }
    }

    // update phaseFlux from gravitational
    phaseFlux += gravPhaseFlux;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += dGravPhaseFlux_dP[ke];
      for( localIndex ic = 0; ic < numComp; ++ic )
        dPhaseFlux_dC[ke][ic] += dGravPhaseFlux_dC[ke][ic];
    }
  }

   template< localIndex numComp, localIndex numFluxSupportPoints >
   GEOS_HOST_DEVICE
   static void
   computeCapillaryFlux( const integer & ip, const integer & numPhase,
                        const localIndex (& seri)[numFluxSupportPoints],
                        const localIndex (& sesri)[numFluxSupportPoints],
                        const localIndex (& sei)[numFluxSupportPoints],
                        const real64 (& trans)[2], const real64 (& dTrans_dPres)[2],
                        const ElementViewConst< arrayView1d< const real64 > > & pres,
                        const ElementViewConst< arrayView1d< const real64 > > & gravCoef,
                        const ElementViewConst< arrayView2d< const real64, 1 > > & phaseMob,
                        const ElementViewConst< arrayView3d< const real64 > > & dPhaseMob,
                        const ElementViewConst< arrayView3d< const real64 > > & dPhaseVolFrac,
                        ElementViewConst< arrayView4d< const real64 > > const & phaseCompFrac,
                        ElementViewConst< arrayView5d< const real64, 4 > > const & dPhaseCompFrac,
                        ElementViewConst< arrayView3d< const real64 > > const & dCompFrac_dCompDens,
                        const ElementViewConst< arrayView3d< const real64 > > & phaseMassDens,
                        const ElementViewConst< arrayView4d< const real64 > > & dPhaseMassDens,
                        const ElementViewConst< arrayView3d< const real64 > > & phaseCapPressure,
                        const ElementViewConst< arrayView4d< const real64 > > & dPhaseCapPressure_dPhaseVolFrac,
                        real64 & phaseFlux,
                        real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                        real64 (& dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
   {
    /// Assembling the capillary flux (and derivatives)
    real64 capPhaseFlux{};
    real64 dCapPhaseFlux_dP[numFluxSupportPoints]{};
    real64 dCapPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    real64 pot_i{};
    real64 dPot_i_dP[numFluxSupportPoints]{};
    real64 dPot_i_dC[numFluxSupportPoints][numComp]{};

  static void computeCapillaryPotential( localIndex const ip,
                       localIndex const numPhase,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                       real64 & capPot,
                       real64 ( & dCapPot_dPres )[numFluxSupportPoints],
                       real64 ( & dCapPot_dComp )[numFluxSupportPoints][numComp] )

    computeCapillaryPotential< numComp, numFluxSupportPoints >( ip,
                                                              numPhase,
                                                              seri,
                                                              sesri,
                                                              sei,
                                                              trans,
                                                              dTrans_dPres,
                                                              gravCoef,
                                                              dCompFrac_dCompDens,
                                                              phaseMassDens,
                                                              dPhaseMassDens,
                                                              pot_i,
                                                              dPot_i_dP,
                                                              dPot_i_dC );

    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( ip != jp )
      {
        real64 pot_j{};
        real64 dPot_j_dP[numFluxSupportPoints]{};
        real64 dPot_j_dC[numFluxSupportPoints][numComp]{};
        computeGravityPotential< numComp, numFluxSupportPoints >( jp,
                                                                  seri,
                                                                  sesri,
                                                                  sei,
                                                                  trans,
                                                                  dTrans_dPres,
                                                                  gravCoef,
                                                                  dCompFrac_dCompDens,
                                                                  phaseMassDens,
                                                                  dPhaseMassDens,
                                                                  pot_j,
                                                                  dPot_j_dP,
                                                                  dPot_j_dC );

        // upwind based on pot diff sign
        real64 const potDiff = pot_j - pot_i;
        real64 dPotDiff_dP[numFluxSupportPoints]{};
        real64 dPotDiff_dC[numFluxSupportPoints][numComp]{};
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dPotDiff_dP[ke] += dPot_j_dP[ke] - dPot_i_dP[ke];
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dPotDiff_dC[ke][jc] += dPot_j_dC[ke][jc] - dPot_i_dC[ke][jc];
          }
        }

        localIndex k_up_i = -1;
        real64 mob_i{};
        real64 dMob_i_dP{};
        real64 dMob_i_dC[numComp]{};
        upwindMobility< numComp, numFluxSupportPoints >( ip,
                                                         seri,
                                                         sesri,
                                                         sei,
                                                         potDiff,
                                                         phaseMob,
                                                         dPhaseMob,
                                                         k_up_i,
                                                         mob_i,
                                                         dMob_i_dP,
                                                         dMob_i_dC );
        localIndex k_up_j = -1;
        real64 mob_j{};
        real64 dMob_j_dP{};
        real64 dMob_j_dC[numComp]{};
        upwindMobility< numComp, numFluxSupportPoints >( jp,
                                                         seri,
                                                         sesri,
                                                         sei,
                                                         -potDiff,
                                                         phaseMob,
                                                         dPhaseMob,
                                                         k_up_j,
                                                         mob_j,
                                                         dMob_j_dP,
                                                         dMob_j_dC );

        real64 const mobTot = LvArray::math::max( mob_i + mob_j, minTotMob );
        real64 const mobTotInv = 1 / mobTot;
        real64 dMobTot_dP[numFluxSupportPoints]{};
        real64 dMobTot_dC[numFluxSupportPoints][numComp]{};
        dMobTot_dP[k_up_i] += dMob_i_dP;
        dMobTot_dP[k_up_j] += dMob_j_dP;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dMobTot_dC[k_up_i][jc] += dMob_i_dC[jc];
          dMobTot_dC[k_up_j][jc] += dMob_j_dC[jc];
        }

        // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
        gravPhaseFlux += mob_i * mob_j * mobTotInv * potDiff;

        // mob_i derivatives
        dGravPhaseFlux_dP[k_up_i] += dMob_i_dP * mob_j * mobTotInv * potDiff;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dGravPhaseFlux_dC[k_up_i][jc] += dMob_i_dC[jc] * mob_j * mobTotInv * potDiff;
        }

        // mob_j derivatives
        dGravPhaseFlux_dP[k_up_j] += mob_i * dMob_j_dP * mobTotInv * potDiff;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dGravPhaseFlux_dC[k_up_j][jc] += mob_i * dMob_j_dC[jc] * mobTotInv * potDiff;
        }

        // mobTot derivatives
        real64 const mobTotInv2 = mobTotInv * mobTotInv;
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dGravPhaseFlux_dP[ke] -= mob_i * mob_j * dMobTot_dP[ke] * mobTotInv2 * potDiff;
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dGravPhaseFlux_dC[ke][jc] -= mob_i * mob_j * dMobTot_dC[ke][jc] * mobTotInv2 * potDiff;
          }
        }

        // potDiff derivatives
        for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dGravPhaseFlux_dP[ke] += mob_i * mob_j * mobTotInv * dPotDiff_dP[ke];
          for( localIndex jc = 0; jc < numComp; ++jc )
          {
            dGravPhaseFlux_dC[ke][jc] += mob_i * mob_j * mobTotInv * dPotDiff_dC[ke][jc];
          }
        }
      }
    }

    // update phaseFlux from gravitational
    phaseFlux += gravPhaseFlux;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += dGravPhaseFlux_dP[ke];
      for( localIndex ic = 0; ic < numComp; ++ic )
        dPhaseFlux_dC[ke][ic] += dGravPhaseFlux_dC[ke][ic];
    }
  }
 
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  computeTotalFlux( integer const & numPhase, const integer & hasCapPressure,
                    localIndex const (&seri)[numFluxSupportPoints],
                    localIndex const (&sesri)[numFluxSupportPoints],
                    localIndex const (&sei)[numFluxSupportPoints],
                    real64 const (&trans)[2], real64 const (&dTrans_dPres)[2],
                    ElementViewConst< arrayView1d< const real64 > > const & pres,
                    ElementViewConst< arrayView1d< const real64 > > const & gravCoef,
                    ElementViewConst< arrayView2d< const real64, 1 > > const & phaseMob,
                    ElementViewConst< arrayView3d< const real64 > > const & dPhaseMob,
                    ElementViewConst< arrayView3d< const real64 > > const & dPhaseVolFrac,
                    ElementViewConst< arrayView3d< const real64 > > const & dCompFrac_dCompDens,
                    ElementViewConst< arrayView3d< const real64 > > const & phaseMassDens,
                    ElementViewConst< arrayView4d< const real64 > > const & dPhaseMassDens,
                    ElementViewConst< arrayView3d< const real64 > > const & phaseCapPressure,
                    ElementViewConst< arrayView4d< const real64 > > const & dPhaseCapPressure_dPhaseVolFrac,
                    real64 & totFlux, real64 (& dTotFlux_dP)[numFluxSupportPoints], real64 (& dTotFlux_dC)[numFluxSupportPoints][numComp] )
  {
    for( integer jp = 0; jp < numPhase; ++jp )
    {
      // working arrays for phase flux
      real64 potGrad{};
      real64 phaseFlux{};
      real64 dPhaseFlux_dP[numFluxSupportPoints]{};
      real64 dPhaseFlux_dC[numFluxSupportPoints][numComp]{};
      PPUPhaseFlux::compute( numPhase, jp, hasCapPressure,
                             seri, sesri, sei,
                             trans, dTrans_dPres,
                             pres, gravCoef,
                             phaseMob, dPhaseMob,
                             dPhaseVolFrac,
                             dCompFrac_dCompDens,
                             phaseMassDens, dPhaseMassDens,
                             phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                             potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

      totFlux += phaseFlux;
      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dTotFlux_dP[ke] += dPhaseFlux_dP[ke];
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dTotFlux_dC[ke][jc] += dPhaseFlux_dC[ke][jc];
        }
      }
    }
  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  upwindMobility( localIndex const ip,
                  localIndex const (&seri)[numFluxSupportPoints],
                  localIndex const (&sesri)[numFluxSupportPoints],
                  localIndex const (&sei)[numFluxSupportPoints],
                  real64 const pot,
                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                  localIndex & upwindDir,
                  real64 & mobility,
                  real64 & dMobility_dP,
                  real64 ( & dMobility_dC)[numComp] )
  {
    upwindDir = (pot > 0) ? 0 : 1;

    localIndex const er_up = seri[upwindDir];
    localIndex const esr_up = sesri[upwindDir];
    localIndex const ei_up = sei[upwindDir];

    mobility = phaseMob[er_up][esr_up][ei_up][ip];
    dMobility_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
    }
  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void computeGravityPotential( localIndex const ip,
                                       localIndex const (&seri)[numFluxSupportPoints],
                                       localIndex const (&sesri)[numFluxSupportPoints],
                                       localIndex const (&sei)[numFluxSupportPoints],
                                       real64 const (&trans)[2],
                                       real64 const (&dTrans_dPres)[2],
                                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                       real64 & gravPot,
                                       real64 ( & dGravPot_dP )[numFluxSupportPoints],
                                       real64 ( & dGravPot_dC )[numFluxSupportPoints][numComp] )
  {
    // init
    gravPot = 0.0;
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      dGravPot_dP[i] = 0.0;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dGravPot_dC[i][jc] = 0.0;
      }
    }

    // get average density

    real64 densMean{};
    real64 dDensMean_dP[numFluxSupportPoints]{};
    real64 dDensMean_dC[numFluxSupportPoints][numComp]{};

    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      // density
      real64 const density = phaseMassDens[er][esr][ei][0][ip];
      real64 const dDens_dPres = dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];
      real64 dDens_dC[numComp]{};
      applyChainRule( numComp,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens[er][esr][ei][0][ip],
                      dDens_dC,
                      Deriv::dC );

      // average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dP[i] = 0.5 * dDens_dPres;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dC[i][jc] = 0.5 * dDens_dC[jc];
      }
    }

    // compute potential difference MPFA-style
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      real64 const gravD = trans[i] * gravCoef[er][esr][ei];
      real64 const dGravD_dP = dTrans_dPres[i] * gravCoef[er][esr][ei];
      gravPot += densMean * gravD;

      // need to add contributions from both cells the mean density depends on
      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dGravPot_dP[ke] += dDensMean_dP[ke] * gravD + densMean * dGravD_dP;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dGravPot_dC[ke][jc] += dDensMean_dC[ke][jc] * gravD;
        }
      }
    }

  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void computeCapillaryPotential( localIndex const ip,
                       localIndex const numPhase,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                       real64 & capPot,
                       real64 ( & dCapPot_dPres )[numFluxSupportPoints],
                       real64 ( & dCapPot_dComp )[numFluxSupportPoints][numComp] )
  {
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      capPot += transmissibility[i] * phaseCapPressure[er][esr][ei][0][ip];
      // need to add contributions from both cells
      for( localIndex jp = 0; jp < numPhase; ++jp )
      {
        real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
        dPot_dPres[i] +=
          transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP]
          + dTrans_dPres[i] * phaseCapPressure[er][esr][ei][0][jp];

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPot_dComp[i][jc] += transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dC + jc];
        }
      }
    }
  }

  static constexpr double minTotMob = 1e-12;

};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHU2PHASEFLUX_HPP
