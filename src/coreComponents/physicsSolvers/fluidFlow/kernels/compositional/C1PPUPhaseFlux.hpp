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
 * @file C1PPUPhaseFlux.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_C1PPUPHASEFLUX_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_C1PPUPHASEFLUX_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/PotGrad.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/PhaseComponentFlux.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernelUtilities
{

// TODO make input parameter
static constexpr real64 epsC1PPU = 5000;

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

struct C1PPUPhaseFlux
{
  /**
   * @brief Form the PhasePotentialUpwind from pressure gradient and gravitational head
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
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
           ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
           localIndex & k_up,
           real64 & potGrad,
           real64 ( &phaseFlux ),
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp],
           real64 ( & compFlux )[numComp],
           real64 ( & dCompFlux_dP )[numFluxSupportPoints][numComp],
           real64 ( & dCompFlux_dC )[numFluxSupportPoints][numComp][numComp] )
  {
    real64 dPotGrad_dTrans = 0.0;
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};
    PotGrad::compute< numComp, numFluxSupportPoints >( numPhase, ip, hasCapPressure, seri, sesri, sei, trans, dTrans_dPres, pres,
                                                       gravCoef, phaseVolFrac, dPhaseVolFrac, dCompFrac_dCompDens, phaseMassDens, dPhaseMassDens,
                                                       phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac, potGrad, dPotGrad_dTrans,
                                                       dPresGrad_dP, dPresGrad_dC, dGravHead_dP, dGravHead_dC );

    // gravity head
    real64 gravHead = 0.0;
    for( integer i = 0; i < numFluxSupportPoints; i++ )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      real64 const gravD = trans[i] * gravCoef[er][esr][ei];

      gravHead += gravD;
    }

    // *** upwinding ***

    // phase flux and derivatives

    // assuming TPFA in the code below

    real64 Ttrans = fabs( trans[0] );
    potGrad = potGrad / Ttrans;

    real64 const mobility_i = phaseMob[seri[0]][sesri[0]][sei[0]][ip];
    real64 const mobility_j = phaseMob[seri[1]][sesri[1]][sei[1]][ip];

    // compute phase flux, see Eqs. (66) and (69) from the reference above
    real64 smoEps = epsC1PPU;
    if( fabs( gravHead ) <= 1e-20 )
      smoEps = 1000;
    real64 const tmpSqrt = sqrt( potGrad * potGrad + smoEps * smoEps );
    real64 const smoMax = 0.5 * (-potGrad + tmpSqrt);

    phaseFlux = Ttrans * ( potGrad * mobility_i - smoMax * (mobility_j - mobility_i) );

    // derivativess

    // first part, mobility derivative

    // dP
    {
      real64 const dMob_dP = dPhaseMob[seri[0]][sesri[0]][sei[0]][ip][Deriv::dP];
      dPhaseFlux_dP[0] += Ttrans * potGrad * dMob_dP;
    }

    // dC
    {
      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 >
      dPhaseMobSub = dPhaseMob[seri[0]][sesri[0]][sei[0]][ip];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[0][jc] += Ttrans * potGrad * dPhaseMobSub[Deriv::dC + jc];
      }
    }

    real64 const tmpInv = 1.0 / tmpSqrt;
    real64 const dSmoMax_x = 0.5 * (1.0 - potGrad * tmpInv);

    // pressure gradient and mobility difference depend on all points in the stencil
    real64 const dMobDiff_sign[numFluxSupportPoints] = {-1.0, 1.0};
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      // dP

      real64 const dPotGrad_dP = dPresGrad_dP[ke] - dGravHead_dP[ke];

      // first part
      dPhaseFlux_dP[ke] += dPotGrad_dP * mobility_i;

      // second part
      real64 const dSmoMax_dP = -dPotGrad_dP * dSmoMax_x;
      dPhaseFlux_dP[ke] += -dSmoMax_dP * (mobility_j - mobility_i);

      real64 const dMob_dP = dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][ip][Deriv::dP];
      dPhaseFlux_dP[ke] += -Ttrans * smoMax * dMobDiff_sign[ke] * dMob_dP;

      // dC

      arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 >
      dPhaseMobSub = dPhaseMob[seri[ke]][sesri[ke]][sei[ke]][ip];

      for( integer jc = 0; jc < numComp; ++jc )
      {
        real64 const dPotGrad_dC = dPresGrad_dC[ke][jc] - dGravHead_dC[ke][jc];

        // first part
        dPhaseFlux_dC[ke][jc] += dPotGrad_dC * mobility_i;

        // second part
        real64 const dSmoMax_dC = -dPotGrad_dC * dSmoMax_x;
        dPhaseFlux_dC[ke][jc] += -dSmoMax_dC * (mobility_j - mobility_i);
        dPhaseFlux_dC[ke][jc] += -Ttrans * smoMax * dMobDiff_sign[ke] * dPhaseMobSub[Deriv::dC + jc];
      }
    }

    potGrad = potGrad * Ttrans;

    // choose upstream cell for composition upwinding
    k_up = (phaseFlux >= 0) ? 0 : 1;

    //distribute on phaseComponentFlux here
    PhaseComponentFlux::compute( ip, k_up, seri, sesri, sei, phaseCompFrac, dPhaseCompFrac, dCompFrac_dCompDens, phaseFlux
                                 , dPhaseFlux_dP, dPhaseFlux_dC, compFlux, dCompFlux_dP, dCompFlux_dC );
  }
};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_C1PPUPHASEFLUX_HPP
