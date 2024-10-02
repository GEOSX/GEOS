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
 * @file CompositionalMultiphaseHybridFVMKernels.cpp
 */

#include "CompositionalMultiphaseHybridFVMKernels.hpp"
#include "CompositionalMultiphaseUtilities.hpp"

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"

#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/HybridFVMHelperKernels.hpp"

namespace geos
{
using namespace constitutive;

namespace compositionalMultiphaseHybridFVMKernels
{

/******************************** UpwindingHelper ********************************/

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  upwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
                            localIndex const (&neighborIds)[ 3 ],
                            ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                            ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                            ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                            ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                            ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                            ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                            real64 const & oneSidedVolFlux,
                            real64 ( & upwPhaseViscCoef )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dPres )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ],
                            globalIndex & upwViscDofNumber )
{
  using Deriv = multifluid::DerivativeOffset;

  real64 dUpwMobRatio_dCompDens[ NC ]{};
  real64 dUpwDensMobRatio_dCompDens[ NC ]{};
  real64 dPhaseDens_dC[ NC ]{};
  real64 dPhaseCompFrac_dC[ NC ]{};

  // 1) Upwind
  localIndex const er  = ( oneSidedVolFlux > 0 ) ? localIds[0] : neighborIds[0];
  localIndex const esr = ( oneSidedVolFlux > 0 ) ? localIds[1] : neighborIds[1];
  localIndex const ei  = ( oneSidedVolFlux > 0 ) ? localIds[2] : neighborIds[2];

  // 2) Compute total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0;
  real64 dTotalMob_dPres = 0;
  real64 dTotalMob_dCompDens[ NC ]{};
  for( integer ip = 0; ip < NP; ++ip )
  {
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres = dTotalMob_dPres + dPhaseMob[er][esr][ei][ip][Deriv::dP];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[ic] = dTotalMob_dCompDens[ic] + dPhaseMob[er][esr][ei][ip][Deriv::dC+ic];
    }
  }
  real64 const totalMobInv = 1.0 / totalMob;

  for( integer ip = 0; ip < NP; ++ip )
  {
    // 3) Compute viscous mobility ratio: \frac{\lambda_{\ell}}{\lambda_T}
    real64 const upwMobRatio = phaseMob[er][esr][ei][ip] * totalMobInv;
    real64 const dUpwMobRatio_dPres = ( dPhaseMob[er][esr][ei][ip][Deriv::dP] - upwMobRatio * dTotalMob_dPres )
                                      * totalMobInv;
    for( integer ic = 0; ic < NC; ++ic )
    {
      dUpwMobRatio_dCompDens[ic] = ( dPhaseMob[er][esr][ei][ip][Deriv::dC+ic] - upwMobRatio * dTotalMob_dCompDens[ic] )
                                   * totalMobInv;
    }

    // 4) Multiply mobility ratio by phase density: \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseDens[er][esr][ei][0][ip],
                    dPhaseDens_dC,
                    Deriv::dC );
    real64 const upwDensMobRatio = phaseDens[er][esr][ei][0][ip] * upwMobRatio;
    real64 const dUpwDensMobRatio_dPres = dPhaseDens[er][esr][ei][0][ip][Deriv::dP] * upwMobRatio
                                          + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dPres;
    for( integer ic = 0; ic < NC; ++ic )
    {
      dUpwDensMobRatio_dCompDens[ic] = dPhaseDens_dC[ic] * upwMobRatio
                                       + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dCompDens[ic];
    }

    // 5) Multiply density mobility ratio by phase comp fraction: x_{c,\ell} \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    for( integer ic = 0; ic < NC; ++ic )
    {
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseCompFrac[er][esr][ei][0][ip][ic],
                      dPhaseCompFrac_dC,
                      Deriv::dC );
      upwPhaseViscCoef[ip][ic] = phaseCompFrac[er][esr][ei][0][ip][ic] * upwDensMobRatio;
      dUpwPhaseViscCoef_dPres[ip][ic] = dPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dP] * upwDensMobRatio
                                        + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dPres;
      for( integer jc = 0; jc < NC; ++jc )
      {
        dUpwPhaseViscCoef_dCompDens[ip][ic][jc] = dPhaseCompFrac_dC[jc] * upwDensMobRatio
                                                  + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dCompDens[jc];
      }
    }
  }
  // 6) Save the dof number of the upwind cell
  upwViscDofNumber = elemDofNumber[er][esr][ei];

}

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  upwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                             localIndex const (&neighborIds)[ 3 ],
                             real64 const & transGravCoef,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                             ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                             ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                             ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                             ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                             real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                             real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                             real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ],
                             real64 ( & upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ],
                             real64 ( & dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ],
                             real64 ( & dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  // 1) Compute the driving force: T ( \rho^{avg}_{\ell} - \rho^{avg}_m ) g \Delta z
  computePhaseGravTerm( localIds,
                        neighborIds,
                        transGravCoef,
                        phaseMassDens,
                        dPhaseMassDens,
                        dCompFrac_dCompDens,
                        phaseGravTerm,
                        dPhaseGravTerm_dPres,
                        dPhaseGravTerm_dCompDens );

  // 2) Compute the total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0.0;
  real64 dTotalMob_dPres[ 2 ]{};
  real64 dTotalMob_dCompDens[ 2 ][ NC ]{};

  computeUpwindedTotalMobility( localIds,
                                neighborIds,
                                phaseMob,
                                dPhaseMob,
                                phaseGravTerm,
                                totalMob,
                                dTotalMob_dPres,
                                dTotalMob_dCompDens );
  real64 const totalMobInv = 1.0 / totalMob;

  // 3) Compute the quantities \x_{up}_{c,p} \rho_p \frac{\lambda_p \lambda_m}{\lambda_T}
  real64 dMobRatio_dPres[ 2 ]{};
  real64 dMobRatio_dCompDens[ 2 ][ NC ]{};
  real64 dDensMobRatio_dPres[ 2 ]{};
  real64 dDensMobRatio_dCompDens[ 2 ][ NC ]{};

  real64 dPhaseDens_dC[ NC ]{};
  real64 dPhaseCompFrac_dC[ NC ]{};

  for( integer ip = 0; ip < NP; ++ip )
  {
    localIndex k = 0;
    for( integer jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      // 3.a) Upwinding using the gravity term
      localIndex eru, esru, eiu, posu; // upwind
      localIndex erd, esrd, eid, posd; // downwind
      setIndicesForMobilityRatioUpwinding( localIds, neighborIds,
                                           phaseGravTerm[ip][k],
                                           eru, esru, eiu, posu,
                                           erd, esrd, eid, posd );

      // 3.b) Compute mobility ratio \frac{\lambda_l \lambda_m}{\lambda_T}
      real64 const mobRatio = phaseMob[eru][esru][eiu][ip] * phaseMob[erd][esrd][eid][jp]
                              * totalMobInv;
      dMobRatio_dPres[posu] = ( dPhaseMob[eru][esru][eiu][ip][Deriv::dP] * phaseMob[erd][esrd][eid][jp]
                                - mobRatio * dTotalMob_dPres[posu] ) * totalMobInv;
      dMobRatio_dPres[posd] = ( dPhaseMob[erd][esrd][eid][jp][Deriv::dP] * phaseMob[eru][esru][eiu][ip]
                                - mobRatio * dTotalMob_dPres[posd] ) * totalMobInv;

      for( integer ic = 0; ic < NC; ++ic )
      {
        dMobRatio_dCompDens[posu][ic] = ( dPhaseMob[eru][esru][eiu][ip][Deriv::dC+ic] * phaseMob[erd][esrd][eid][jp]
                                          - mobRatio * dTotalMob_dCompDens[posu][ic] ) * totalMobInv;
        dMobRatio_dCompDens[posd][ic] = ( dPhaseMob[erd][esrd][eid][jp][Deriv::dC+ic] * phaseMob[eru][esru][eiu][ip]
                                          - mobRatio * dTotalMob_dCompDens[posd][ic] ) * totalMobInv;
      }

      // 3.c) Compute mobility ratio multiplied by upwinded phase density \rho_l \frac{\lambda_l \lambda_m}{\lambda_T}
      applyChainRule( NC,
                      dCompFrac_dCompDens[eru][esru][eiu],
                      dPhaseDens[eru][esru][eiu][0][ip],
                      dPhaseDens_dC,
                      Deriv::dC );
      real64 const densMobRatio = phaseDens[eru][esru][eiu][0][ip] * mobRatio;
      dDensMobRatio_dPres[posu] = dPhaseDens[eru][esru][eiu][0][ip][Deriv::dP] * mobRatio
                                  + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posu];
      dDensMobRatio_dPres[posd] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posd];
      for( integer ic = 0; ic < NC; ++ic )
      {
        dDensMobRatio_dCompDens[posu][ic] = dPhaseDens_dC[ic] * mobRatio
                                            + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posu][ic];
        dDensMobRatio_dCompDens[posd][ic] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posd][ic];
      }

      // 3.d) Compute the final gravity coefficient \x_{up}_{c,p} \rho_p \frac{\lambda_l \lambda_m}{\lambda_T}
      for( integer ic = 0; ic < NC; ++ic )
      {
        applyChainRule( NC,
                        dCompFrac_dCompDens[eru][esru][eiu],
                        dPhaseCompFrac[eru][esru][eiu][0][ip][ic],
                        dPhaseCompFrac_dC,
                        Deriv::dC );
        upwPhaseGravCoef[ip][k][ic] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * densMobRatio;
        dUpwPhaseGravCoef_dPres[ip][k][ic][posu] = dPhaseCompFrac[eru][esru][eiu][0][ip][ic][Deriv::dP] * densMobRatio
                                                   + phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posu];
        dUpwPhaseGravCoef_dPres[ip][k][ic][posd] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posd];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dUpwPhaseGravCoef_dCompDens[ip][k][ic][posu][jc] = dPhaseCompFrac_dC[jc] * densMobRatio
                                                             + phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posu][jc];
          dUpwPhaseGravCoef_dCompDens[ip][k][ic][posd][jc] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posd][jc];
        }
      }
      ++k;
    }
  }

}

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  computePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                        localIndex const (&neighborIds)[ 3 ],
                        real64 const & transGravCoef,
                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                        ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                        real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                        real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                        real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  localIndex const er   = localIds[0];
  localIndex const esr  = localIds[1];
  localIndex const ei   = localIds[2];
  localIndex const ern  = neighborIds[0];
  localIndex const esrn = neighborIds[1];
  localIndex const ein  = neighborIds[2];

  real64 dPhaseMassDens_dCLoc[ NC ]{};
  real64 dPhaseMassDens_dCNeighbor[ NC ]{};
  real64 dPhaseMassDens_dC[ NC ]{};

  for( integer ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseMassDens[er][esr][ei][0][ip],
                    dPhaseMassDens_dCLoc,
                    Deriv::dC );
    applyChainRule( NC,
                    dCompFrac_dCompDens[ern][esrn][ein],
                    dPhaseMassDens[ern][esrn][ein][0][ip],
                    dPhaseMassDens_dCNeighbor,
                    Deriv::dC );

    localIndex k = 0;
    for( integer jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      phaseGravTerm[ip][k] = -( phaseMassDens[er][esr][ei][0][ip] + phaseMassDens[ern][esrn][ein][0][ip] );
      phaseGravTerm[ip][k] += ( phaseMassDens[er][esr][ei][0][jp] + phaseMassDens[ern][esrn][ein][0][jp] );
      phaseGravTerm[ip][k] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] = ( -dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP] + dPhaseMassDens[er][esr][ei][0][jp][Deriv::dP] );
      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] = ( -dPhaseMassDens[ern][esrn][ein][0][ip][Deriv::dP] + dPhaseMassDens[ern][esrn][ein][0][jp][Deriv::dP] );
      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] *= 0.5 * transGravCoef;

      for( integer ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCLoc[ic];
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCNeighbor[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens[er][esr][ei][0][jp],
                      dPhaseMassDens_dC,
                      Deriv::dC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[ern][esrn][ein],
                      dPhaseMassDens[ern][esrn][ein][0][jp],
                      dPhaseMassDens_dC,
                      Deriv::dC );
      for( integer ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
      }
      ++k;
    }
  }
}

template< integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  computeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                localIndex const (&neighborIds)[ 3 ],
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                real64 & totalMob,
                                real64 ( & dTotalMob_dPres )[ 2 ],
                                real64 ( & dTotalMob_dCompDens )[ 2 ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  localIndex totalMobIds[ NP ][ 3 ]{};
  localIndex totalMobPos[ NP ]{};
  setIndicesForTotalMobilityUpwinding< NP >( localIds,
                                             neighborIds,
                                             phaseGravTerm,
                                             totalMobIds,
                                             totalMobPos );
  for( integer ip = 0; ip < NP; ++ip )
  {
    localIndex const er  = totalMobIds[ip][0];
    localIndex const esr = totalMobIds[ip][1];
    localIndex const ei  = totalMobIds[ip][2];
    localIndex const pos = totalMobPos[ip];
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres[pos] = dTotalMob_dPres[pos] + dPhaseMob[er][esr][ei][pos][Deriv::dP];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[pos][ic] = dTotalMob_dCompDens[pos][ic] + dPhaseMob[er][esr][ei][ip][Deriv::dC+ic];
    }
  }
  if( totalMob < 1e-12 )
  {
    totalMob = 1e-12;
  }
}

GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  setIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const & gravTerm,
                                       localIndex & eru, localIndex & esru, localIndex & eiu, localIndex & posu,
                                       localIndex & erd, localIndex & esrd, localIndex & eid, localIndex & posd )
{
  if( gravTerm > 0 )
  {
    eru  = localIds[0];
    esru = localIds[1];
    eiu  = localIds[2];
    posu = Pos::LOCAL;
    erd  = neighborIds[0];
    esrd = neighborIds[1];
    eid  = neighborIds[2];
    posd = Pos::NEIGHBOR;
  }
  else
  {
    eru  = neighborIds[0];
    esru = neighborIds[1];
    eiu  = neighborIds[2];
    posu = Pos::NEIGHBOR;
    erd  = localIds[0];
    esrd = localIds[1];
    eid  = localIds[2];
    posd = Pos::LOCAL;
  }
}

template< integer NP >
GEOS_HOST_DEVICE
inline
void
UpwindingHelper::
  setIndicesForTotalMobilityUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const (&gravTerm)[ NP ][ NP-1 ],
                                       localIndex ( & totalMobIds )[ NP ][ 3 ],
                                       localIndex ( & totalMobPos )[ NP ] )
{
  if( NP == 2 )
  {
    if( gravTerm[0][0] > 0 )
    {
      totalMobIds[0][0] = localIds[0];
      totalMobIds[0][1] = localIds[1];
      totalMobIds[0][2] = localIds[2];
      totalMobPos[0] = Pos::LOCAL;
      totalMobIds[1][0] = neighborIds[0];
      totalMobIds[1][1] = neighborIds[1];
      totalMobIds[1][2] = neighborIds[2];
      totalMobPos[1] = Pos::NEIGHBOR;
    }
    else
    {
      totalMobIds[0][0] = neighborIds[0];
      totalMobIds[0][1] = neighborIds[1];
      totalMobIds[0][2] = neighborIds[2];
      totalMobPos[0] = Pos::NEIGHBOR;
      totalMobIds[1][0] = localIds[0];
      totalMobIds[1][1] = localIds[1];
      totalMobIds[1][2] = localIds[2];
      totalMobPos[1] = Pos::LOCAL;
    }
  }
  else if( NP == 3 )
  {
    // TODO Francois: this should be improved
    // currently this implements the algorithm proposed by SH Lee
    for( integer ip = 0; ip < NP; ++ip )
    {
      if( ( gravTerm[ip][0] >= 0 && gravTerm[ip][1] >= 0 ) || // includes the no-buoyancy case
          ( ( LvArray::math::abs( gravTerm[ip][0] ) >= LvArray::math::abs( gravTerm[ip][1] ) ) && gravTerm[ip][1] >= 0 ) ||
          ( ( LvArray::math::abs( gravTerm[ip][1] ) >= LvArray::math::abs( gravTerm[ip][0] ) ) && gravTerm[ip][0] >= 0 ) )
      {
        totalMobIds[ip][0] = localIds[0];
        totalMobIds[ip][1] = localIds[1];
        totalMobIds[ip][2] = localIds[2];
        totalMobPos[ip] = Pos::LOCAL;
      }
      else
      {
        totalMobIds[ip][0] = neighborIds[0];
        totalMobIds[ip][1] = neighborIds[1];
        totalMobIds[ip][2] = neighborIds[2];
        totalMobPos[ip] = Pos::NEIGHBOR;
      }
    }
  }
}


#define INST_UpwindingHelperNCNP( NC, NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    upwindViscousCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                        localIndex const (&neighborIds)[ 3 ], \
                                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                        ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                        ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                        ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                        ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                        ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                        real64 const & oneSidedVolFlux, \
                                        real64 ( &upwPhaseViscCoef )[ NP ][ NC ], \
                                        real64 ( &dUpwPhaseViscCoef_dPres )[ NP ][ NC ], \
                                        real64 ( &dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ], \
                                        globalIndex & upwViscDofNumber ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    upwindBuoyancyCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                         localIndex const (&neighborIds)[ 3 ], \
                                         real64 const & transGravCoef,  \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                         ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                         real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                         real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                         real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ], \
                                         real64 ( &upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ], \
                                         real64 ( &dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ], \
                                         real64 ( &dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    computePhaseGravTerm< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                    localIndex const (&neighborIds)[ 3 ], \
                                    real64 const & transGravCoef, \
                                    ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                    ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                    real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                    real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                    real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    computeUpwindedTotalMobility< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                            localIndex const (&neighborIds)[ 3 ], \
                                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                            real64 const (&phaseGravTerm)[ NP ][ NP-1 ], \
                                            real64 & totalMob,    \
                                            real64 ( &dTotalMob_dPres )[ 2 ], \
                                            real64 ( &dTotalMob_dCompDens )[ 2 ][ NC ] )

INST_UpwindingHelperNCNP( 1, 2 );
INST_UpwindingHelperNCNP( 2, 2 );
INST_UpwindingHelperNCNP( 3, 2 );
INST_UpwindingHelperNCNP( 4, 2 );
INST_UpwindingHelperNCNP( 5, 2 );

INST_UpwindingHelperNCNP( 1, 3 );
INST_UpwindingHelperNCNP( 2, 3 );
INST_UpwindingHelperNCNP( 3, 3 );
INST_UpwindingHelperNCNP( 4, 3 );
INST_UpwindingHelperNCNP( 5, 3 );

#undef INST_UpwindingHelperNCNP

#define INST_UpwindingHelperNP( NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  UpwindingHelper:: \
    setIndicesForTotalMobilityUpwinding< NP >( localIndex const (&localIds)[ 3 ], \
                                               localIndex const (&neighborIds)[ 3 ], \
                                               real64 const (&gravTerm)[ NP ][ NP-1 ], \
                                               localIndex ( &totalMobIds )[ NP ][ 3 ], \
                                               localIndex ( &totalMobPos )[ NP ] )

INST_UpwindingHelperNP( 2 );
INST_UpwindingHelperNP( 3 );

#undef INST_UpwindingHelperNP

/******************************** AssemblerKernelHelper ********************************/

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  applyGradient( arrayView1d< real64 const > const & facePres,
                 arrayView1d< real64 const > const & faceGravCoef,
                 arraySlice1d< localIndex const > const & elemToFaces,
                 real64 const & elemPres,
                 real64 const & elemGravCoef,
                 arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & elemPhaseMassDens,
                 arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dElemPhaseMassDens,
                 arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & elemPhaseMob,
                 arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dElemPhaseMob,
                 arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dElemCompFrac_dCompDens,
                 arraySlice2d< real64 const > const & transMatrix,
                 real64 ( & oneSidedVolFlux )[ NF ],
                 real64 ( & dOneSidedVolFlux_dPres )[ NF ],
                 real64 ( & dOneSidedVolFlux_dFacePres )[ NF ][ NF ],
                 real64 ( & dOneSidedVolFlux_dCompDens )[ NF ][ NC ] )
{
  using Deriv = multifluid::DerivativeOffset;

  real64 dPhaseMassDens_dC[ NP ][ NC ]{};
  real64 dPresDif_dCompDens[ NC ]{};
  real64 dPhaseGravDif_dCompDens[ NC ]{};
  real64 dPhaseMobPotDif_dCompDens[ NC ]{};

  // 0) precompute dPhaseDens_dC since it is always computed at the element center
  for( integer ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dElemCompFrac_dCompDens,
                    dElemPhaseMassDens[ip],
                    dPhaseMassDens_dC[ip],
                    Deriv::dC );
  }

  for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {

      // depth difference between element center and face center
      real64 const ccGravCoef = elemGravCoef;
      real64 const fGravCoef = faceGravCoef[elemToFaces[jfaceLoc]];
      real64 const gravCoefDif = ccGravCoef - fGravCoef;

      for( integer ip = 0; ip < NP; ++ip )
      {

        // 1) compute the potential diff between the cell center and the face center
        real64 const ccPres = elemPres;
        real64 const fPres  = facePres[elemToFaces[jfaceLoc]];

        // pressure difference
        real64 const presDif = ccPres - fPres;
        real64 const dPresDif_dPres = 1;
        real64 const dPresDif_dFacePres = -1;
        for( integer ic = 0; ic < NC; ++ic )
        {
          dPresDif_dCompDens[ic] = 0.0; // no capillary pressure
        }

        // gravity term
        real64 const phaseGravDif = elemPhaseMassDens[ip] * gravCoefDif;
        real64 const dPhaseGravDif_dPres = dElemPhaseMassDens[ip][Deriv::dP] * gravCoefDif;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseGravDif_dCompDens[ic] = dPhaseMassDens_dC[ip][ic] * gravCoefDif;
        }
        // no density evaluated at the face center

        // potential difference
        real64 const phasePotDif = presDif - phaseGravDif;
        real64 const phaseMobPotDif = elemPhaseMob[ip] * phasePotDif;
        real64 const dPhaseMobPotDif_dPres = dElemPhaseMob[ip][Deriv::dP] * phasePotDif
                                             + elemPhaseMob[ip] * (dPresDif_dPres - dPhaseGravDif_dPres);
        real64 const dPhaseMobPotDif_dFacePres = elemPhaseMob[ip] * dPresDif_dFacePres;
        for( integer ic = 0; ic < NC; ++ic )
        {
          dPhaseMobPotDif_dCompDens[ic] = dElemPhaseMob[ip][Deriv::dC+ic] * phasePotDif
                                          + elemPhaseMob[ip] * (dPresDif_dCompDens[ic] - dPhaseGravDif_dCompDens[ic]);
        }

        // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
        oneSidedVolFlux[ifaceLoc] = oneSidedVolFlux[ifaceLoc]
                                    + transMatrix[ifaceLoc][jfaceLoc] * phaseMobPotDif;
        dOneSidedVolFlux_dPres[ifaceLoc] = dOneSidedVolFlux_dPres[ifaceLoc]
                                           + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dPres;
        dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc] = dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc]
                                                         + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dFacePres;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dOneSidedVolFlux_dCompDens[ifaceLoc][ic] = dOneSidedVolFlux_dCompDens[ifaceLoc][ic]
                                                     + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dCompDens[ic];
        }
      }
    }
  }
}

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleFluxDivergence( localIndex const (&localIds)[ 3 ],
                          globalIndex const rankOffset,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< real64 const > const & mimFaceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemGravCoef,
                          integer const useTotalMassEquation,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
{
  using namespace compositionalMultiphaseUtilities;
  integer constexpr NDOF = NC+1;

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF*(NF+1) ]{};
  globalIndex dofColIndicesFaceVars[ NF ]{};
  for( integer idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber[localIds[0]][localIds[1]][localIds[2]] + idof;
  }

  // divergence of fluxes
  real64 divMassFluxes[ NC ]{};
  real64 dDivMassFluxes_dElemVars[ NC ][ NDOF*(NF+1) ]{};
  real64 dDivMassFluxes_dFaceVars[ NC ][ NF ]{};

  // auxiliary variables for upwinding

  // upwinding phase buoyancy transport coefficients
  real64 upwPhaseViscCoef[ NP ][ NC ]{};
  real64 dUpwPhaseViscCoef_dPres[ NP ][ NC ]{};
  real64 dUpwPhaseViscCoef_dCompDens[ NP ][ NC ][ NC ]{};
  globalIndex upwViscDofNumber = 0;

  // gravity term: ( \rho_l - \rho_m ) g \Delta z
  real64 phaseGravTerm[ NP ][ NP-1 ]{};
  real64 dPhaseGravTerm_dPres[ NP ][ NP-1 ][ 2 ]{};
  real64 dPhaseGravTerm_dCompDens[ NP ][ NP-1 ][ 2 ][ NC ]{};

  // upwinding phase buoyancy transport coefficients
  real64 upwPhaseGravCoef[ NP ][ NP-1 ][ NC ]{};
  real64 dUpwPhaseGravCoef_dPres[ NP ][ NP-1 ][ NC ][ 2 ]{};
  real64 dUpwPhaseGravCoef_dCompDens[ NP ][ NP-1 ][ NC ][ 2 ][ NC ]{};

  // for each element, loop over the one-sided faces
  for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {

    // 1) Find if there is a neighbor, and if there is, grab the indices of the neighbor element

    localIndex neighborIds[ 3 ] = { localIds[0], localIds[1], localIds[2] };
    hybridFVMKernels::CellConnectivity::isNeighborFound( localIds,
                                                         ifaceLoc,
                                                         elemRegionList,
                                                         elemSubRegionList,
                                                         elemList,
                                                         regionFilter,
                                                         elemToFaces,
                                                         neighborIds );
    localIndex const neighborDofNumber = elemDofNumber[neighborIds[0]][neighborIds[1]][neighborIds[2]];

    // 2) *************** Assemble viscous terms ******************

    // 2.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} for each phase at this face
    UpwindingHelper::upwindViscousCoefficient< NC, NP >( localIds,
                                                         neighborIds,
                                                         phaseDens,
                                                         dPhaseDens,
                                                         phaseMob,
                                                         dPhaseMob,
                                                         dCompFrac_dCompDens,
                                                         phaseCompFrac,
                                                         dPhaseCompFrac,
                                                         elemDofNumber,
                                                         oneSidedVolFlux[ifaceLoc],
                                                         upwPhaseViscCoef,
                                                         dUpwPhaseViscCoef_dPres,
                                                         dUpwPhaseViscCoef_dCompDens,
                                                         upwViscDofNumber );

    // 2.b) Add the \x_{c,\ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} q_T of this face to the divergence of the flux in this cell
    assembleViscousFlux< NF, NC, NP >( ifaceLoc,
                                       oneSidedVolFlux,
                                       dOneSidedVolFlux_dPres,
                                       dOneSidedVolFlux_dFacePres,
                                       dOneSidedVolFlux_dCompDens,
                                       upwPhaseViscCoef,
                                       dUpwPhaseViscCoef_dPres,
                                       dUpwPhaseViscCoef_dCompDens,
                                       elemDofNumber[localIds[0]][localIds[1]][localIds[2]],
                                       neighborDofNumber,
                                       upwViscDofNumber,
                                       faceDofNumber[elemToFaces[ifaceLoc]],
                                       dt,
                                       divMassFluxes,
                                       dDivMassFluxes_dElemVars,
                                       dDivMassFluxes_dFaceVars,
                                       dofColIndicesElemVars,
                                       dofColIndicesFaceVars );

    // 3) *************** Assemble buoyancy terms ******************

    real64 const transGravCoef = (localIds[0] != neighborIds[0] || localIds[1] != neighborIds[1] || localIds[2] != neighborIds[2])
                                 * transMatrixGrav[ifaceLoc][ifaceLoc] * (elemGravCoef - mimFaceGravCoef[elemToFaces[ifaceLoc]]);

    // 3.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}\lambda_m}{\lambda_T}
    //      and (\rho_{\ell} - \rho_m) g \Delta z for each phase at this face
    UpwindingHelper::upwindBuoyancyCoefficient< NC, NP >( localIds,
                                                          neighborIds,
                                                          transGravCoef,
                                                          phaseDens,
                                                          dPhaseDens,
                                                          phaseMassDens,
                                                          dPhaseMassDens,
                                                          phaseMob,
                                                          dPhaseMob,
                                                          dCompFrac_dCompDens,
                                                          phaseCompFrac,
                                                          dPhaseCompFrac,
                                                          phaseGravTerm,
                                                          dPhaseGravTerm_dPres,
                                                          dPhaseGravTerm_dCompDens,
                                                          upwPhaseGravCoef,
                                                          dUpwPhaseGravCoef_dPres,
                                                          dUpwPhaseGravCoef_dCompDens );

    // 3.b) Add the buoyancy term of this face to the divergence of the flux in this cell
    assembleBuoyancyFlux< NF, NC, NP >( ifaceLoc,
                                        phaseGravTerm,
                                        dPhaseGravTerm_dPres,
                                        dPhaseGravTerm_dCompDens,
                                        upwPhaseGravCoef,
                                        dUpwPhaseGravCoef_dPres,
                                        dUpwPhaseGravCoef_dCompDens,
                                        dt,
                                        divMassFluxes,
                                        dDivMassFluxes_dElemVars );

  }

  if( useTotalMassEquation > 0 )
  {
    // Apply equation/variable change transformation(s)
    real64 work[NDOF * ( NF + 1 )];
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NDOF * ( NF + 1 ), dDivMassFluxes_dElemVars, work );
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NF, dDivMassFluxes_dFaceVars, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, divMassFluxes );
  }

  // we are ready to assemble the local flux and its derivatives
  // no need for atomic adds - each row is assembled by a single thread

  for( integer ic = 0; ic < NC; ++ic )
  {
    localIndex const eqnRowLocalIndex =
      LvArray::integerConversion< localIndex >( elemDofNumber[localIds[0]][localIds[1]][localIds[2]] + ic - rankOffset );

    GEOS_ASSERT_GE( eqnRowLocalIndex, 0 );
    GEOS_ASSERT_GT( localMatrix.numRows(), eqnRowLocalIndex );

    // residual
    localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + divMassFluxes[ic];

    // jacobian -- derivative wrt elem centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &dofColIndicesElemVars[0],
                                                              &dDivMassFluxes_dElemVars[0][0] + ic * NDOF * (NF+1),
                                                              NDOF * (NF+1) );

    // jacobian -- derivatives wrt face centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &dofColIndicesFaceVars[0],
                                                              &dDivMassFluxes_dFaceVars[0][0] + ic * NF,
                                                              NF );
  }
}

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleViscousFlux( localIndex const ifaceLoc,
                       real64 const (&oneSidedVolFlux)[ NF ],
                       real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                       real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                       real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                       real64 const (&upwPhaseViscCoef)[ NP ][ NC ],
                       real64 const (&dUpwPhaseViscCoef_dPres)[ NP ][ NC ],
                       real64 const (&dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ],
                       globalIndex const elemDofNumber,
                       globalIndex const neighborDofNumber,
                       globalIndex const upwViscDofNumber,
                       globalIndex const faceDofNumber,
                       real64 const & dt,
                       real64 ( & divMassFluxes )[ NC ],
                       real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ],
                       real64 ( & dDivMassFluxes_dFaceVars )[ NC ][ NF ],
                       globalIndex ( & dofColIndicesElemVars )[ (NC+1)*(NF+1) ],
                       globalIndex ( & dofColIndicesFaceVars )[ NF ] )
{
  integer constexpr NDOF = NC+1;
  localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

  for( integer ip = 0; ip < NP; ++ip )
  {
    for( integer ic = 0; ic < NC; ++ic )
    {
      // compute the mass flux at the one-sided face plus its derivatives
      // add the newly computed flux to the sum

      real64 const dt_upwPhaseViscCoef = dt * upwPhaseViscCoef[ip][ic];

      // residual
      divMassFluxes[ic] = divMassFluxes[ic] + dt_upwPhaseViscCoef * oneSidedVolFlux[ifaceLoc];

      // local derivatives
      dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                        + dt_upwPhaseViscCoef * dOneSidedVolFlux_dPres[ifaceLoc];
      dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                        + ( elemDofNumber == upwViscDofNumber )
                                        * dt * dUpwPhaseViscCoef_dPres[ip][ic] * oneSidedVolFlux[ifaceLoc];
      for( integer jc = 0; jc < NC; ++jc )
      {
        dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                             + dt_upwPhaseViscCoef * dOneSidedVolFlux_dCompDens[ifaceLoc][jc];
        dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                             + ( elemDofNumber == upwViscDofNumber )
                                             * dt * dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
      }

      // neighbor derivatives
      dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                     + ( elemDofNumber != upwViscDofNumber )
                                                     * dt * dUpwPhaseViscCoef_dPres[ip][ic] * oneSidedVolFlux[ifaceLoc];

      for( integer jc = 0; jc < NC; ++jc )
      {
        dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                            + ( elemDofNumber != upwViscDofNumber )
                                                            * dt * dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
      }

      for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        dDivMassFluxes_dFaceVars[ic][jfaceLoc] = dDivMassFluxes_dFaceVars[ic][jfaceLoc]
                                                 + dt_upwPhaseViscCoef * dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
      }
    }
  }

  // collect the relevant dof numbers
  for( integer idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[elemVarsOffset+idof] = neighborDofNumber + idof;
  }
  dofColIndicesFaceVars[ifaceLoc] = faceDofNumber;
}

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleBuoyancyFlux( localIndex const ifaceLoc,
                        real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                        real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ],
                        real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ],
                        real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ],
                        real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ],
                        real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ],
                        real64 const & dt,
                        real64 ( & divMassFluxes )[ NC ],
                        real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] )
{
  integer constexpr NDOF = NC+1;
  localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

  for( integer ip = 0; ip < NP; ++ip )
  {
    for( integer jp = 0; jp < NP - 1; ++jp )
    {
      for( integer ic = 0; ic < NC; ++ic )
      {
        real64 const dt_upwPhaseGravCoef = dt * upwPhaseGravCoef[ip][jp][ic];

        // residual
        divMassFluxes[ic] = divMassFluxes[ic] + dt_upwPhaseGravCoef * phaseGravTerm[ip][jp];

        // local derivatives
        dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                          + dt * dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::LOCAL] * phaseGravTerm[ip][jp];
        dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                          + dt_upwPhaseGravCoef * dPhaseGravTerm_dPres[ip][jp][Pos::LOCAL];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                               + dt * dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::LOCAL][jc] * phaseGravTerm[ip][jp];
          dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                               + dt_upwPhaseGravCoef * dPhaseGravTerm_dCompDens[ip][jp][Pos::LOCAL][jc];
        }

        // neighbor derivatives
        dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                       + dt * dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::NEIGHBOR] * phaseGravTerm[ip][jp];
        dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                       + dt_upwPhaseGravCoef * dPhaseGravTerm_dPres[ip][jp][Pos::NEIGHBOR];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                              + dt * dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::NEIGHBOR][jc] * phaseGravTerm[ip][jp];
          dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                              + dt_upwPhaseGravCoef * dPhaseGravTerm_dCompDens[ip][jp][Pos::NEIGHBOR][jc];
        }
      }
    }
  }
}

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernelHelper::
  assembleFaceConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                           arrayView1d< integer const > const & faceGhostRank,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           globalIndex const elemDofNumber,
                           globalIndex const rankOffset,
                           real64 const (&oneSidedVolFlux)[ NF ],
                           real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                           real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                           real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
{
  integer constexpr NDOF = NC+1;

  // fluxes
  real64 dFlux_dElemVars[ NDOF ]{};
  real64 dFlux_dFaceVars[ NF ]{};

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF ]{};
  globalIndex dofColIndicesFaceVars[ NF ]{};
  for( integer idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber + idof;
  }

  // for each element, loop over the local (one-sided) faces
  for( integer ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    if( faceGhostRank[elemToFaces[ifaceLoc]] >= 0 )
    {
      continue;
    }

    // flux at this face
    real64 const flux = oneSidedVolFlux[ifaceLoc];
    dFlux_dElemVars[0] = dOneSidedVolFlux_dPres[ifaceLoc];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dFlux_dElemVars[ic+1] = dOneSidedVolFlux_dCompDens[ifaceLoc][ic];
    }

    for( integer jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dFlux_dFaceVars[jfaceLoc] = dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
      dofColIndicesFaceVars[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
    }

    // dof number of this face constraint
    localIndex const eqnLocalRowIndex = LvArray::integerConversion< localIndex >( faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset );

    GEOS_ASSERT_GE( eqnLocalRowIndex, 0 );
    GEOS_ASSERT_GT( localMatrix.numRows(), eqnLocalRowIndex );

    // residual
    RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnLocalRowIndex], flux );

    // jacobian -- derivatives wrt elem-centered terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesElemVars[0],
                                                                      &dFlux_dElemVars[0],
                                                                      NDOF );

    // jacobian -- derivatives wrt face pressure terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesFaceVars[0],
                                                                      &dFlux_dFaceVars[0],
                                                                      NF );
  }
}


#define INST_AssemblerKernelHelper( NF, NC, NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    applyGradient< NF, NC, NP >( arrayView1d< real64 const > const & facePres, \
                                 arrayView1d< real64 const > const & faceGravCoef, \
                                 arraySlice1d< localIndex const > const & elemToFaces, \
                                 real64 const & elemPres, \
                                 real64 const & elemGravCoef, \
                                 arraySlice1d< real64 const, multifluid::USD_PHASE-2 > const & elemPhaseMassDens, \
                                 arraySlice2d< real64 const, multifluid::USD_PHASE_DC-2 > const & dElemPhaseMassDens_dCompFrac, \
                                 arraySlice1d< real64 const, compflow::USD_PHASE-1 > const & elemPhaseMob, \
                                 arraySlice2d< real64 const, compflow::USD_PHASE_DC-1 > const & dElemPhaseMob, \
                                 arraySlice2d< real64 const, compflow::USD_COMP_DC-1 > const & dElemCompFrac_dCompDens, \
                                 arraySlice2d< real64 const > const & transMatrix, \
                                 real64 ( &oneSidedVolFlux )[ NF ], \
                                 real64 ( &dOneSidedVolFlux_dPres )[ NF ], \
                                 real64 ( &dOneSidedVolFlux_dFacePres )[ NF ][ NF ], \
                                 real64 ( &dOneSidedVolFlux_dCompDens )[ NF ][ NC ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleFluxDivergence< NF, NC, NP >( localIndex const (&localIds)[ 3 ], \
                                          globalIndex const rankOffset, \
                                          arrayView2d< localIndex const > const & elemRegionList, \
                                          arrayView2d< localIndex const > const & elemSubRegionList, \
                                          arrayView2d< localIndex const > const & elemList, \
                                          SortedArrayView< localIndex const > const & regionFilter, \
                                          arrayView1d< globalIndex const > const & faceDofNumber, \
                                          arrayView1d< real64 const > const & mimFaceGravCoef, \
                                          arraySlice1d< localIndex const > const & elemToFaces, \
                                          real64 const & elemGravCoef, \
                                          integer const useTotalMassEquation, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                          arraySlice2d< real64 const > const & transMatrixGrav, \
                                          real64 const (&oneSidedVolFlux)[ NF ], \
                                          real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                          real64 const & dt, \
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                          arrayView1d< real64 > const & localRhs ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleViscousFlux< NF, NC, NP >( localIndex const ifaceLoc, \
                                       real64 const (&oneSidedVolFlux)[ NF ], \
                                       real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                       real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                       real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                       real64 const (&upwPhaseViscCoef)[ NP ][ NC ], \
                                       real64 const (&dUpwPhaseViscCoef_dPres)[ NP ][ NC ], \
                                       real64 const (&dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ], \
                                       globalIndex const elemDofNumber, \
                                       globalIndex const neighborDofNumber, \
                                       globalIndex const upwViscDofNumber, \
                                       globalIndex const faceDofNumber, \
                                       real64 const & dt, \
                                       real64 ( &divMassFluxes )[ NC ], \
                                       real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ], \
                                       real64 ( &dDivMassFluxes_dFaceVars )[ NC ][ NF ], \
                                       globalIndex ( &dofColIndicesElemVars )[ (NC+1)*(NF+1) ], \
                                       globalIndex ( &dofColIndicesFaceVars )[ NF ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleBuoyancyFlux< NF, NC, NP >( localIndex const ifaceLoc, \
                                        real64 const (&phaseGravTerm)[ NP ][ NP-1 ], \
                                        real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ], \
                                        real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ], \
                                        real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ], \
                                        real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ], \
                                        real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ], \
                                        real64 const & dt, \
                                        real64 ( &divMassFluxes )[ NC ], \
                                        real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] ); \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernelHelper:: \
    assembleFaceConstraints< NF, NC, NP >( arrayView1d< globalIndex const > const & faceDofNumber, \
                                           arrayView1d< integer const > const & faceGhostRank, \
                                           arraySlice1d< localIndex const > const & elemToFaces, \
                                           globalIndex const elemDofNumber, \
                                           globalIndex const rankOffset,  \
                                           real64 const (&oneSidedVolFlux)[ NF ], \
                                           real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                           real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                           real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                           arrayView1d< real64 > const & localRhs )

INST_AssemblerKernelHelper( 4, 1, 2 );
INST_AssemblerKernelHelper( 4, 2, 2 );
INST_AssemblerKernelHelper( 4, 3, 2 );
INST_AssemblerKernelHelper( 4, 4, 2 );
INST_AssemblerKernelHelper( 4, 5, 2 );

INST_AssemblerKernelHelper( 4, 1, 3 );
INST_AssemblerKernelHelper( 4, 2, 3 );
INST_AssemblerKernelHelper( 4, 3, 3 );
INST_AssemblerKernelHelper( 4, 4, 3 );
INST_AssemblerKernelHelper( 4, 5, 3 );

INST_AssemblerKernelHelper( 5, 1, 2 );
INST_AssemblerKernelHelper( 5, 2, 2 );
INST_AssemblerKernelHelper( 5, 3, 2 );
INST_AssemblerKernelHelper( 5, 4, 2 );
INST_AssemblerKernelHelper( 5, 5, 2 );

INST_AssemblerKernelHelper( 5, 1, 3 );
INST_AssemblerKernelHelper( 5, 2, 3 );
INST_AssemblerKernelHelper( 5, 3, 3 );
INST_AssemblerKernelHelper( 5, 4, 3 );
INST_AssemblerKernelHelper( 5, 5, 3 );

INST_AssemblerKernelHelper( 6, 1, 2 );
INST_AssemblerKernelHelper( 6, 2, 2 );
INST_AssemblerKernelHelper( 6, 3, 2 );
INST_AssemblerKernelHelper( 6, 4, 2 );
INST_AssemblerKernelHelper( 6, 5, 2 );

INST_AssemblerKernelHelper( 6, 1, 3 );
INST_AssemblerKernelHelper( 6, 2, 3 );
INST_AssemblerKernelHelper( 6, 3, 3 );
INST_AssemblerKernelHelper( 6, 4, 3 );
INST_AssemblerKernelHelper( 6, 5, 3 );

#undef INST_AssemblerKernelHelper

/******************************** AssemblerKernel ********************************/

template< integer NF, integer NC, integer NP >
GEOS_HOST_DEVICE
inline
void
AssemblerKernel::
  compute( localIndex const er, localIndex const esr, localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arrayView1d< real64 const > const & mimFaceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
           real64 const & elemGravCoef,
           integer const useTotalMassEquation,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           arraySlice2d< real64 const > const & transMatrixGrav,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  // one sided flux
  real64 oneSidedVolFlux[ NF ]{};
  real64 dOneSidedVolFlux_dPres[ NF ]{};
  real64 dOneSidedVolFlux_dFacePres[ NF ][ NF ]{};
  real64 dOneSidedVolFlux_dCompDens[ NF ][ NC ]{};

  localIndex const localIds[3] = { er, esr, ei };

  /*
   * compute auxiliary quantities at the one sided faces of this element:
   * 1) One-sided volumetric fluxes
   * 2) Upwinded mobilities
   */

  // for each one-sided face of the elem,
  // compute the volumetric flux using transMatrix
  AssemblerKernelHelper::applyGradient< NF, NC, NP >( facePres,
                                                      faceGravCoef,
                                                      elemToFaces,
                                                      elemPres,
                                                      elemGravCoef,
                                                      phaseMassDens[er][esr][ei][0],
                                                      dPhaseMassDens[er][esr][ei][0],
                                                      phaseMob[er][esr][ei],
                                                      dPhaseMob[er][esr][ei],
                                                      dCompFrac_dCompDens[er][esr][ei],
                                                      transMatrix,
                                                      oneSidedVolFlux,
                                                      dOneSidedVolFlux_dPres,
                                                      dOneSidedVolFlux_dFacePres,
                                                      dOneSidedVolFlux_dCompDens );

  // at this point, we know the local flow direction in the element
  // so we can upwind the transport coefficients (mobilities) at the one sided faces
  // ** this function needs non-local information **
  if( elemGhostRank < 0 )
  {
    /*
     * perform assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // use the computed one sided vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
    AssemblerKernelHelper::assembleFluxDivergence< NF, NC, NP >( localIds,
                                                                 rankOffset,
                                                                 elemRegionList,
                                                                 elemSubRegionList,
                                                                 elemList,
                                                                 regionFilter,
                                                                 faceDofNumber,
                                                                 mimFaceGravCoef,
                                                                 elemToFaces,
                                                                 elemGravCoef,
                                                                 useTotalMassEquation,
                                                                 phaseDens,
                                                                 dPhaseDens,
                                                                 phaseMassDens,
                                                                 dPhaseMassDens,
                                                                 phaseMob,
                                                                 dPhaseMob,
                                                                 dCompFrac_dCompDens,
                                                                 phaseCompFrac,
                                                                 dPhaseCompFrac,
                                                                 elemDofNumber,
                                                                 transMatrixGrav,
                                                                 oneSidedVolFlux,
                                                                 dOneSidedVolFlux_dPres,
                                                                 dOneSidedVolFlux_dFacePres,
                                                                 dOneSidedVolFlux_dCompDens,
                                                                 dt,
                                                                 localMatrix,
                                                                 localRhs );
  }

  // use the computed one sided vol fluxes to assemble the constraints
  // enforcing flux continuity at this element's faces
  AssemblerKernelHelper::assembleFaceConstraints< NF, NC, NP >( faceDofNumber,
                                                                faceGhostRank,
                                                                elemToFaces,
                                                                elemDofNumber[er][esr][ei],
                                                                rankOffset,
                                                                oneSidedVolFlux,
                                                                dOneSidedVolFlux_dPres,
                                                                dOneSidedVolFlux_dFacePres,
                                                                dOneSidedVolFlux_dCompDens,
                                                                localMatrix,
                                                                localRhs );
}

#define INST_AssemblerKernel( NF, NC, NP ) \
  template \
  GEOS_HOST_DEVICE \
  void \
  AssemblerKernel:: \
    compute< NF, NC, NP >( localIndex const er, localIndex const esr, localIndex const ei, \
                           SortedArrayView< localIndex const > const & regionFilter, \
                           arrayView2d< localIndex const > const & elemRegionList, \
                           arrayView2d< localIndex const > const & elemSubRegionList, \
                           arrayView2d< localIndex const > const & elemList, \
                           arrayView1d< globalIndex const > const & faceDofNumber, \
                           arrayView1d< integer const > const & faceGhostRank, \
                           arrayView1d< real64 const > const & facePres, \
                           arrayView1d< real64 const > const & faceGravCoef, \
                           arrayView1d< real64 const > const & mimFaceGravCoef, \
                           arraySlice1d< localIndex const > const & elemToFaces, \
                           real64 const & elemPres, \
                           real64 const & elemGravCoef, \
                           integer const useTotalMassEquation, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                           integer const elemGhostRank, \
                           globalIndex const rankOffset, \
                           real64 const & dt, \
                           arraySlice2d< real64 const > const & transMatrix, \
                           arraySlice2d< real64 const > const & transMatrixGrav, \
                           CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                           arrayView1d< real64 > const & localRhs )

INST_AssemblerKernel( 4, 1, 2 );
INST_AssemblerKernel( 4, 2, 2 );
INST_AssemblerKernel( 4, 3, 2 );
INST_AssemblerKernel( 4, 4, 2 );
INST_AssemblerKernel( 4, 5, 2 );

INST_AssemblerKernel( 4, 1, 3 );
INST_AssemblerKernel( 4, 2, 3 );
INST_AssemblerKernel( 4, 3, 3 );
INST_AssemblerKernel( 4, 4, 3 );
INST_AssemblerKernel( 4, 5, 3 );

INST_AssemblerKernel( 5, 1, 2 );
INST_AssemblerKernel( 5, 2, 2 );
INST_AssemblerKernel( 5, 3, 2 );
INST_AssemblerKernel( 5, 4, 2 );
INST_AssemblerKernel( 5, 5, 2 );

INST_AssemblerKernel( 5, 1, 3 );
INST_AssemblerKernel( 5, 2, 3 );
INST_AssemblerKernel( 5, 3, 3 );
INST_AssemblerKernel( 5, 4, 3 );
INST_AssemblerKernel( 5, 5, 3 );

INST_AssemblerKernel( 6, 1, 2 );
INST_AssemblerKernel( 6, 2, 2 );
INST_AssemblerKernel( 6, 3, 2 );
INST_AssemblerKernel( 6, 4, 2 );
INST_AssemblerKernel( 6, 5, 2 );

INST_AssemblerKernel( 6, 1, 3 );
INST_AssemblerKernel( 6, 2, 3 );
INST_AssemblerKernel( 6, 3, 3 );
INST_AssemblerKernel( 6, 4, 3 );
INST_AssemblerKernel( 6, 5, 3 );

#undef INST_AssemblerKernel



/******************************** FluxKernel ********************************/

template< integer NF, integer NC, integer NP, typename IP_TYPE >
void
FluxKernel::
  launch( localIndex er, localIndex esr,
          CellElementSubRegion const & subRegion,
          constitutive::PermeabilityBase const & permeabilityModel,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView1d< real64 const > const & mimFaceGravCoef,
          arrayView1d< real64 const > const & transMultiplier,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          globalIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          integer const useTotalMassEquation,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  // get the cell-centered DOF numbers and ghost rank for the assembly
  arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

  // get the map from elem to faces
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

  // get the cell-centered pressures
  arrayView1d< real64 const > const & elemPres  =
    subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );

  // get the element data needed for transmissibility computation
  arrayView2d< real64 const > const & elemCenter =
    subRegion.getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );

  // TODO add this dependency to the compute function
  //arrayView3d< real64 const > const elemdPermdPres = permeabilityModel.dPerm_dPressure();

  arrayView3d< real64 const > const & elemPerm = permeabilityModel.permeability();

  // get the cell-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.getReference< array1d< real64 > >( fields::flow::gravityCoefficient::key() );

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const ei )
  {
    // transmissibility matrix
    stackArray2d< real64, NF *NF > transMatrix( NF, NF );
    stackArray2d< real64, NF *NF > transMatrixGrav( NF, NF );

    real64 const perm[ 3 ] = { elemPerm[ei][0][0], elemPerm[ei][0][1], elemPerm[ei][0][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    IP_TYPE::template compute< NF >( nodePosition,
                                     transMultiplier,
                                     faceToNodes,
                                     elemToFaces[ei],
                                     elemCenter[ei],
                                     elemVolume[ei],
                                     perm,
                                     lengthTolerance,
                                     transMatrix );

    // currently the gravity term in the transport scheme is treated as in MRST, that is, always with TPFA
    // this is why below we have to recompute the TPFA transmissibility in addition to the transmissibility matrix above
    // TODO: treat the gravity term with a consistent inner product
    mimeticInnerProduct::TPFAInnerProduct::compute< NF >( nodePosition,
                                                          transMultiplier,
                                                          faceToNodes,
                                                          elemToFaces[ei],
                                                          elemCenter[ei],
                                                          elemVolume[ei],
                                                          perm,
                                                          lengthTolerance,
                                                          transMatrixGrav );

    // perform flux assembly in this element
    compositionalMultiphaseHybridFVMKernels::AssemblerKernel::compute< NF, NC, NP >( er, esr, ei,
                                                                                     regionFilter,
                                                                                     elemRegionList,
                                                                                     elemSubRegionList,
                                                                                     elemList,
                                                                                     faceDofNumber,
                                                                                     faceGhostRank,
                                                                                     facePres,
                                                                                     faceGravCoef,
                                                                                     mimFaceGravCoef,
                                                                                     elemToFaces[ei],
                                                                                     elemPres[ei],
                                                                                     elemGravCoef[ei],
                                                                                     useTotalMassEquation,
                                                                                     phaseDens,
                                                                                     dPhaseDens,
                                                                                     phaseMassDens,
                                                                                     dPhaseMassDens,
                                                                                     phaseMob,
                                                                                     dPhaseMob,
                                                                                     dCompFrac_dCompDens,
                                                                                     phaseCompFrac,
                                                                                     dPhaseCompFrac,
                                                                                     elemDofNumber,
                                                                                     elemGhostRank[ei],
                                                                                     rankOffset,
                                                                                     dt,
                                                                                     transMatrix,
                                                                                     transMatrixGrav,
                                                                                     localMatrix,
                                                                                     localRhs );
  } );
}

#define INST_FluxKernel( NF, NC, NP, IP_TYPE ) \
  template \
  void \
  FluxKernel:: \
    launch< NF, NC, NP, IP_TYPE >( localIndex er, localIndex esr, \
                                   CellElementSubRegion const & subRegion, \
                                   constitutive::PermeabilityBase const & permeabilityModel, \
                                   SortedArrayView< localIndex const > const & regionFilter, \
                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                   arrayView2d< localIndex const > const & elemRegionList, \
                                   arrayView2d< localIndex const > const & elemSubRegionList, \
                                   arrayView2d< localIndex const > const & elemList, \
                                   ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                   arrayView1d< globalIndex const > const & faceDofNumber, \
                                   arrayView1d< integer const > const & faceGhostRank, \
                                   arrayView1d< real64 const > const & facePres, \
                                   arrayView1d< real64 const > const & faceGravCoef, \
                                   arrayView1d< real64 const > const & mimFaceGravCoef, \
                                   arrayView1d< real64 const > const & transMultiplier, \
                                   ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                   ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
                                   ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                   globalIndex const rankOffset, \
                                   real64 const lengthTolerance, \
                                   real64 const dt, \
                                   integer const useTotalMassEquation, \
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                   arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 4, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 4, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 5, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 5, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 6, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 6, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );


INST_FluxKernel( 4, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 4, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 5, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 5, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 6, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 6, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );


#undef INST_FluxKernel

} // namespace compositionalMultiphaseHybridFVMKernels

} // namespace geos
