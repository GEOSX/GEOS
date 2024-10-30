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

  static constexpr double minTotMob = 1e-12;

  /**
   * @brief Simplified 2-phase version of hybrid upwinding
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
      computeCapillaryFlux< numComp, numFluxSupportPoints >( ip, numPhase,
                                                             seri, sesri, sei,
                                                             trans, dTrans_dPres,
                                                             phaseMob, dPhaseMob,
                                                             dPhaseVolFrac,
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
      // accumulate total mobility
      UpwindHelpers::addValueAndDerivatives( totMob, dTotMob_dP, dTotMob_dC, mob, dMob_dP, dMob_dC );
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
    // accumulate total mobility
    UpwindHelpers::addValueAndDerivatives( totMob, dTotMob_dP, dTotMob_dC, mob, dMob_dP, dMob_dC );

    // safeguard
    totMob = LvArray::math::max( totMob, minTotMob );
    real64 const invTotMob = 1 / totMob;

    // fractional flow for viscous part as \lambda_i^{up}/\sum_{NP}(\lambda_j^{up})

    real64 const fractionalFlow = mob * invTotMob;
    real64 dFractionalFlow_dP{};
    real64 dFractionalFlow_dC[numComp]{};
    UpwindHelpers::addDerivativesScaled( dFractionalFlow_dP, dFractionalFlow_dC, dMob_dP, dMob_dC, invTotMob );
    UpwindHelpers::addDerivativesScaled( dFractionalFlow_dP, dFractionalFlow_dC, dTotMob_dP, dTotMob_dC, -fractionalFlow * invTotMob );

    /// Assembling the viscous flux (and derivatives) from fractional flow and total velocity as \phi_{\mu} = f_i^{up,\mu} uT

    real64 const viscousPhaseFlux = fractionalFlow * totFlux;
    real64 dViscousPhaseFlux_dP[numFluxSupportPoints]{};
    real64 dViscousPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    // fractionalFlow derivatives
    UpwindHelpers::addDerivativesScaled( dViscousPhaseFlux_dP[k_up], dViscousPhaseFlux_dC[k_up], dFractionalFlow_dP, dFractionalFlow_dC, totFlux );

    // Ut derivatives
    UpwindHelpers::addDerivativesScaled( dViscousPhaseFlux_dP, dViscousPhaseFlux_dC, dTotFlux_dP, dTotFlux_dC, fractionalFlow );

    // accumulate in the flux and its derivatives (need to be very careful doing that)
    UpwindHelpers::addValueAndDerivatives( phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC,
                                           viscousPhaseFlux, dViscousPhaseFlux_dP, dViscousPhaseFlux_dC );
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

        computePotDiffFlux( ip, jp,
                            seri, sesri, sei,
                            pot_i, dPot_i_dP, dPot_i_dC,
                            pot_j, dPot_j_dP, dPot_j_dC,
                            phaseMob, dPhaseMob,
                            gravPhaseFlux, dGravPhaseFlux_dP, dGravPhaseFlux_dC );
      }
    }

    // update phaseFlux from gravitational
    UpwindHelpers::addValueAndDerivatives( phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC,
                                           gravPhaseFlux, dGravPhaseFlux_dP, dGravPhaseFlux_dC );
  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  computeCapillaryFlux( const integer & ip, const integer & numPhase,
                        const localIndex (& seri)[numFluxSupportPoints],
                        const localIndex (& sesri)[numFluxSupportPoints],
                        const localIndex (& sei)[numFluxSupportPoints],
                        const real64 (& trans)[2], const real64 (& dTrans_dPres)[2],
                        const ElementViewConst< arrayView2d< const real64, 1 > > & phaseMob,
                        const ElementViewConst< arrayView3d< const real64 > > & dPhaseMob,
                        const ElementViewConst< arrayView3d< const real64 > > & dPhaseVolFrac,
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

    computeCapillaryPotential< numComp, numFluxSupportPoints >( ip,
                                                                numPhase,
                                                                seri,
                                                                sesri,
                                                                sei,
                                                                trans,
                                                                dTrans_dPres,
                                                                dPhaseVolFrac,
                                                                phaseCapPressure,
                                                                dPhaseCapPressure_dPhaseVolFrac,
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
        computeCapillaryPotential< numComp, numFluxSupportPoints >( jp,
                                                                    numPhase,
                                                                    seri,
                                                                    sesri,
                                                                    sei,
                                                                    trans,
                                                                    dTrans_dPres,
                                                                    dPhaseVolFrac,
                                                                    phaseCapPressure,
                                                                    dPhaseCapPressure_dPhaseVolFrac,
                                                                    pot_j,
                                                                    dPot_j_dP,
                                                                    dPot_j_dC );

        computePotDiffFlux( ip, jp,
                            seri, sesri, sei,
                            pot_i, dPot_i_dP, dPot_i_dC,
                            pot_j, dPot_j_dP, dPot_j_dC,
                            phaseMob, dPhaseMob,
                            capPhaseFlux, dCapPhaseFlux_dP, dCapPhaseFlux_dC );
      }
    }

    // update phaseFlux from capillary flux
    UpwindHelpers::addValueAndDerivatives( phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC,
                                           capPhaseFlux, dCapPhaseFlux_dP, dCapPhaseFlux_dC );
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

      UpwindHelpers::addValueAndDerivatives( totFlux, dTotFlux_dP, dTotFlux_dC,
                                             phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );
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
    UpwindHelpers::assignToZero( mobility, dMobility_dP, dMobility_dC );
    UpwindHelpers::assignMobilityAndDerivatives( ip, upwindDir, seri, sesri, sei, phaseMob, dPhaseMob, mobility, dMobility_dP, dMobility_dC );
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
    UpwindHelpers::assignToZero( gravPot, dGravPot_dP, dGravPot_dC );

    // get average density TODO change after #3337 is merged

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
      dGravPot_dP[i] += densMean * dGravD_dP;

      // need to add contributions from both cells the mean density depends on
      UpwindHelpers::addDerivativesScaled( dGravPot_dP, dGravPot_dC, dDensMean_dP, dDensMean_dC, gravD );
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
                                         real64 ( & dCapPot_dP )[numFluxSupportPoints],
                                         real64 ( & dCapPot_dC )[numFluxSupportPoints][numComp] )
  {
    // init
    UpwindHelpers::assignToZero( capPot, dCapPot_dP, dCapPot_dC );

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
        dCapPot_dP[i] +=
          transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP]
          + dTrans_dPres[i] * phaseCapPressure[er][esr][ei][0][jp];

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dCapPot_dC[i][jc] += transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dC + jc];
        }
      }
    }
  }

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  computePotDiffFlux( integer const & ip, integer const & jp,
                      localIndex const (&seri)[numFluxSupportPoints],
                      localIndex const (&sesri)[numFluxSupportPoints],
                      localIndex const (&sei)[numFluxSupportPoints],
                      real64 const & pot_i, real64 const ( & dPot_i_dP )[numFluxSupportPoints], real64 const (&dPot_i_dC )[numFluxSupportPoints][numComp],
                      real64 const & pot_j, real64 const ( & dPot_j_dP )[numFluxSupportPoints], real64 const ( &dPot_j_dC)[numFluxSupportPoints][numComp],
                      ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                      ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                      real64 & phaseFlux, real64 ( & dPhaseFlux_dP )[numFluxSupportPoints], real64 ( & dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
  {
    // upwind based on pot diff sign
    real64 const potDiff = pot_j - pot_i;
    real64 dPotDiff_dP[numFluxSupportPoints]{};
    real64 dPotDiff_dC[numFluxSupportPoints][numComp]{};
    UpwindHelpers::addDerivativesScaled( dPotDiff_dP, dPotDiff_dC, dPot_j_dP, dPot_j_dC, 1.0 );
    UpwindHelpers::addDerivativesScaled( dPotDiff_dP, dPotDiff_dC, dPot_i_dP, dPot_i_dC, -1.0 );

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

    // safeguard
    real64 const mobTot = LvArray::math::max( mob_i + mob_j, minTotMob );
    real64 const mobTotInv = 1 / mobTot;
    real64 dMobTot_dP[numFluxSupportPoints]{};
    real64 dMobTot_dC[numFluxSupportPoints][numComp]{};
    UpwindHelpers::addDerivatives( dMobTot_dP[k_up_i], dMobTot_dC[k_up_i], dMob_i_dP, dMob_i_dC );
    UpwindHelpers::addDerivatives( dMobTot_dP[k_up_j], dMobTot_dC[k_up_j], dMob_j_dP, dMob_j_dC );

    // Assembling flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (Pot_i - Pot_k)
    phaseFlux += mob_i * mob_j * mobTotInv * potDiff;

    // mob_i derivatives
    UpwindHelpers::addDerivativesScaled( dPhaseFlux_dP[k_up_i], dPhaseFlux_dC[k_up_i], dMob_i_dP, dMob_i_dC, mob_j * mobTotInv * potDiff );

    // mob_j derivatives
    UpwindHelpers::addDerivativesScaled( dPhaseFlux_dP[k_up_j], dPhaseFlux_dC[k_up_j], dMob_j_dP, dMob_j_dC, mob_i * mobTotInv * potDiff );

    // mobTot derivatives
    real64 const mobTotInv2 = mobTotInv * mobTotInv;
    UpwindHelpers::addDerivativesScaled( dPhaseFlux_dP, dPhaseFlux_dC, dMobTot_dP, dMobTot_dC, -mob_i * mob_j * mobTotInv2 * potDiff );

    // potDiff derivatives
    UpwindHelpers::addDerivativesScaled( dPhaseFlux_dP, dPhaseFlux_dC, dPotDiff_dP, dPotDiff_dC, mob_i * mob_j * mobTotInv );
  }

};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHU2PHASEFLUX_HPP
