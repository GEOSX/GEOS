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
 * @file PhaseComponentFlux.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASECOMPONENTFLUX_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASECOMPONENTFLUX_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernelUtilities
{

template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

using Deriv = constitutive::multifluid::DerivativeOffset;

struct PhaseComponentFlux
{
  /**
   * @brief Compute the component flux for a given phase
   * @tparam numComp number of components
   * @tparam numFluxSupportPoints number of flux support points
   * @param ip phase index
   * @param k_up uptream index for this phase
   * @param seri arraySlice of the stencil-implied element region index
   * @param sesri arraySlice of the stencil-implied element subregion index
   * @param sei arraySlice of the stencil-implied element index
   * @param phaseCompFrac phase component fraction
   * @param dPhaseCompFrac derivative of phase component fraction wrt pressure, temperature, component fraction
   * @param dCompFrac_dCompDens derivative of component fraction wrt component density
   * @param phaseFlux phase flux
   * @param dPhaseFlux_dP derivative of phase flux wrt pressure
   * @param dPhaseFlux_dC derivative of phase flux wrt comp density
   * @param compFlux component flux
   * @param dCompFlux_dP derivative of phase flux wrt pressure
   * @param dCompFlux_dC derivative of phase flux wrt comp density
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute( localIndex const ip,
           localIndex const k_up,
           localIndex const ( &seri )[numFluxSupportPoints],
           localIndex const ( &sesri )[numFluxSupportPoints],
           localIndex const ( &sei )[numFluxSupportPoints],
           ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           real64 const & phaseFlux,
           real64 const ( &dPhaseFlux_dP )[numFluxSupportPoints],
           real64 const ( &dPhaseFlux_dC )[numFluxSupportPoints][numComp],
           real64 ( & compFlux )[numComp],
           real64 ( & dCompFlux_dP )[numFluxSupportPoints][numComp],
           real64 ( & dCompFlux_dC )[numFluxSupportPoints][numComp][numComp] )
  {
    localIndex const er_up = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up = sei[k_up];

    real64 dProp_dC[numComp]{};

    // slice some constitutive arrays to avoid too much indexing in component loop
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP-3 > phaseCompFracSub =
      phaseCompFrac[er_up][esr_up][ei_up][0][ip];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFracSub =
      dPhaseCompFrac[er_up][esr_up][ei_up][0][ip];

    // compute component fluxes and derivatives using upstream cell composition
    for( integer ic = 0; ic < numComp; ++ic )
    {
      real64 const ycp = phaseCompFracSub[ic];
      compFlux[ic] += phaseFlux * ycp;

      // derivatives stemming from phase flux
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
        }
      }

      // additional derivatives stemming from upstream cell phase composition
      dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFracSub[ic][Deriv::dP];

      // convert derivatives of comp fraction w.r.t. comp fractions to derivatives w.r.t. comp densities
      applyChainRule( numComp,
                      dCompFrac_dCompDens[er_up][esr_up][ei_up],
                      dPhaseCompFracSub[ic],
                      dProp_dC,
                      Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
      }
    }
  }
};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASECOMPONENTFLUX_HPP
