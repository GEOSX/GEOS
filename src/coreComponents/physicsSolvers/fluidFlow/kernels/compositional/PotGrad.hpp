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
 * @file PotGrad.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_POTGRAD_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_POTGRAD_HPP

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

struct PotGrad
{
  template< integer numComp, integer numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void
  compute ( integer const numPhase,
            integer const ip,
            integer const hasCapPressure,
            localIndex const ( &seri )[numFluxSupportPoints],
            localIndex const ( &sesri )[numFluxSupportPoints],
            localIndex const ( &sei )[numFluxSupportPoints],
            real64 const ( &trans )[numFluxSupportPoints],
            real64 const ( &dTrans_dPres )[numFluxSupportPoints],
            ElementViewConst< arrayView1d< real64 const > > const & pres,
            ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
            ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
            ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
            ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
            ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
            ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
            real64 & potGrad,
            real64 & dPotGrad_dTrans,
            real64 ( & dPresGrad_dP )[numFluxSupportPoints],
            real64 ( & dPresGrad_dC )[numFluxSupportPoints][numComp],
            real64 ( & dGravHead_dP )[numFluxSupportPoints],
            real64 ( & dGravHead_dC )[numFluxSupportPoints][numComp] )
  {
    // assign derivatives arrays to zero
    for( integer i = 0; i < numFluxSupportPoints; ++i )
    {
      dPresGrad_dP[i] = 0.0;
      dGravHead_dP[i] = 0.0;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPresGrad_dC[i][jc] = 0.0;
        dGravHead_dC[i][jc] = 0.0;
      }
    }

    // create local work arrays
    real64 densMean = 0.0;
    real64 dDensMean_dP[numFluxSupportPoints]{};
    real64 dDensMean_dC[numFluxSupportPoints][numComp]{};

    real64 presGrad = 0.0;
    real64 dPresGrad_dTrans = 0.0;
    real64 gravHead = 0.0;
    real64 dGravHead_dTrans = 0.0;
    real64 dCapPressure_dC[numComp]{};

    real64 dProp_dC[numComp]{};

    // calculate quantities on primary connected cells
    integer denom = 0;
    for( integer i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      bool const phaseExists = (phaseVolFrac[er][esr][ei][ip] > 0);
      if( !phaseExists )
      {
        continue;
      }

      // density
      real64 const density  = phaseMassDens[er][esr][ei][0][ip];
      real64 const dDens_dP = dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];

      applyChainRule( numComp,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens[er][esr][ei][0][ip],
                      dProp_dC,
                      Deriv::dC );

      // average density and derivatives
      densMean += density;
      dDensMean_dP[i] = dDens_dP;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dC[i][jc] = dProp_dC[jc];
      }
      denom++;
    }
    if( denom > 1 )
    {
      densMean /= denom;
      for( integer i = 0; i < numFluxSupportPoints; ++i )
      {
        dDensMean_dP[i] /= denom;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dDensMean_dC[i][jc] /= denom;
        }
      }
    }

    /// compute the TPFA potential difference
    for( integer i = 0; i < numFluxSupportPoints; i++ )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      // capillary pressure
      real64 capPressure     = 0.0;
      real64 dCapPressure_dP = 0.0;

      for( integer ic = 0; ic < numComp; ++ic )
      {
        dCapPressure_dC[ic] = 0.0;
      }

      if( hasCapPressure )
      {
        capPressure = phaseCapPressure[er][esr][ei][0][ip];

        for( integer jp = 0; jp < numPhase; ++jp )
        {
          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP];

          for( integer jc = 0; jc < numComp; ++jc )
          {
            dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dC+jc];
          }
        }
      }

      real64 const dP = pres[er][esr][ei] - capPressure;
      presGrad += trans[i] * dP;
      dPresGrad_dTrans += dP;
      dPresGrad_dP[i] += trans[i] * (1 - dCapPressure_dP) + dTrans_dPres[i] * dP;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPresGrad_dC[i][jc] += -trans[i] * dCapPressure_dC[jc];
      }

      real64 const gC = gravCoef[er][esr][ei];
      real64 const gravD     = trans[i] * gC;
      real64 const dGravD_dTrans = gC;
      real64 const dGravD_dP = dTrans_dPres[i] * gC;

      // the density used in the potential difference is always a mass density
      // unlike the density used in the phase mobility, which is a mass density
      // if useMass == 1 and a molar density otherwise
      gravHead += densMean * gravD;
      dGravHead_dTrans += densMean * dGravD_dTrans;
      // need to add contributions from both cells the mean density depends on
      for( integer j = 0; j < numFluxSupportPoints; ++j )
      {
        dGravHead_dP[j] += dDensMean_dP[j] * gravD + dGravD_dP * densMean;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
        }
      }
    }

    // compute phase potential gradient
    potGrad = presGrad - gravHead;
    dPotGrad_dTrans = dPresGrad_dTrans - dGravHead_dTrans;
  }

};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_POTGRAD_HPP
