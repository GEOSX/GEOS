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
 * @file PPUPhaseFluxZFormulation.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUXZFORMULATION_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUXZFORMULATION_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/PotGradZFormulation.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/PhaseComponentFluxZFormulation.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernelUtilities
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

struct PPUPhaseFluxZFormulation
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
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
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
    real64 dPresGrad_dP[numFluxSupportPoints]{};
    real64 dPresGrad_dC[numFluxSupportPoints][numComp]{};
    real64 dGravHead_dP[numFluxSupportPoints]{};
    real64 dGravHead_dC[numFluxSupportPoints][numComp]{};
    PotGradZFormulation::compute< numComp, numFluxSupportPoints >( numPhase, ip, hasCapPressure, seri, sesri, sei, trans, dTrans_dPres, pres,
                                                                   gravCoef, dPhaseVolFrac, phaseMassDens, dPhaseMassDens,
                                                                   phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac, potGrad, dPresGrad_dP,
                                                                   dPresGrad_dC, dGravHead_dP, dGravHead_dC );

    // *** upwinding ***

    // choose upstream cell
    k_up = (potGrad >= 0) ? 0 : 1;

    localIndex const er_up  = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up  = sei[k_up];

    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

    // pressure gradient depends on all points in the stencil
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += dPresGrad_dP[ke] - dGravHead_dP[ke];
      dPhaseFlux_dP[ke] *= mobility;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc] - dGravHead_dC[ke][jc];
        dPhaseFlux_dC[ke][jc] *= mobility;
      }
    }
    // compute phase flux using upwind mobility.
    phaseFlux = mobility * potGrad;

    real64 const dMob_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dPhaseMobSub =
      dPhaseMob[er_up][esr_up][ei_up][ip];

    // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
    for( integer jc = 0; jc < numComp; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dPhaseMobSub[Deriv::dC+jc] * potGrad;
    }

    //distribute on phaseComponentFlux here
    PhaseComponentFluxZFormulation::compute( ip, k_up, seri, sesri, sei, phaseCompFrac, dPhaseCompFrac, phaseFlux
                                             , dPhaseFlux_dP, dPhaseFlux_dC, compFlux, dCompFlux_dP, dCompFlux_dC );

  }
};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PPUPHASEFLUXZFORMULATION_HPP
