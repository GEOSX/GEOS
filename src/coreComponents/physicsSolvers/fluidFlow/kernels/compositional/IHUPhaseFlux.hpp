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
 * @file IHUPhaseFlux.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHUPHASEFLUX_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHUPHASEFLUX_HPP

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

/************************* HELPERS ******************/
namespace UpwindHelpers
{

static constexpr double minTotMob = 1e-12;

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
upwindMobilityViscous( localIndex const ip,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const totFlux,
                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                       localIndex & upwindDir,
                       real64 & mobility,
                       real64 ( & dMobility_dP ),
                       real64 ( & dMobility_dC )[numComp] )
{
  // reinit
  mobility = 0.0;
  dMobility_dP = 0.0;
  for( localIndex ic = 0; ic < numComp; ++ic )
  {
    dMobility_dC[ic] = 0.0;
  }

  upwindDir = (totFlux > 0) ? 0 : 1;

  localIndex const er_up = seri[upwindDir];
  localIndex const esr_up = sesri[upwindDir];
  localIndex const ei_up = sei[upwindDir];

  if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
  {
    mobility = phaseMob[er_up][esr_up][ei_up][ip];
    dMobility_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
    }
  }
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
upwindMobilityGravity( localIndex const numPhase,
                       localIndex const ip,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       real64 const totFlux,            //in fine should be a ElemnetViewConst once seq form are in place
                       ElementViewConst< arrayView1d< real64 const > > const & pres,
                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                       integer const capPressureFlag,
                       localIndex & upwindDir,
                       real64 & mobility,
                       real64( &dMobility_dP),
                       real64 ( & dMobility_dC)[numComp]
                       )
{

  //reinit
  mobility = 0.0;
  dMobility_dP = 0.0;
  for( localIndex ic = 0; ic < numComp; ++ic )
  {
    dMobility_dC[ic] = 0.0;
  }

  UPWIND scheme;
  scheme.template getUpwindDirectionGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                      ip,
                                                                                      seri,
                                                                                      sesri,
                                                                                      sei,
                                                                                      transmissibility,
                                                                                      dTrans_dPres,
                                                                                      totFlux,
                                                                                      pres,
                                                                                      gravCoef,
                                                                                      phaseMob,
                                                                                      dCompFrac_dCompDens,
                                                                                      phaseMassDens,
                                                                                      dPhaseMassDens,
                                                                                      dPhaseVolFrac,
                                                                                      phaseCapPressure,
                                                                                      dPhaseCapPressure_dPhaseVolFrac,
                                                                                      capPressureFlag,
                                                                                      upwindDir );

  localIndex const er_up = seri[upwindDir];
  localIndex const esr_up = sesri[upwindDir];
  localIndex const ei_up = sei[upwindDir];

  if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
  {
    mobility = phaseMob[er_up][esr_up][ei_up][ip];
    dMobility_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
    }
  }
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
upwindMobilityCapillary( localIndex const numPhase,
                         localIndex const ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         real64 const (&transmissibility)[2],
                         real64 const (&dTrans_dPres)[2],
                         real64 const totFlux,          //in fine should be a ElemnetViewConst once seq form are in place
                         ElementViewConst< arrayView1d< real64 const > > const & pres,
                         ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                         ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                         ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                         ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                         ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                         integer const capPressureFlag,
                         localIndex & upwindDir,
                         real64 & mobility,
                         real64( &dMobility_dP),
                         real64 ( & dMobility_dC)[numComp]
                         )
{

  //reinit
  mobility = 0.0;
  dMobility_dP = 0.0;
  for( localIndex ic = 0; ic < numComp; ++ic )
  {
    dMobility_dC[ic] = 0.0;
  }

  UPWIND scheme;
  scheme.template getUpwindDirectionCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                        ip,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        pres,
                                                                                        gravCoef,
                                                                                        phaseMob,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        capPressureFlag,
                                                                                        upwindDir );

  localIndex const er_up = seri[upwindDir];
  localIndex const esr_up = sesri[upwindDir];
  localIndex const ei_up = sei[upwindDir];

  if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
  {
    mobility = phaseMob[er_up][esr_up][ei_up][ip];
    dMobility_dP = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dP];
    for( localIndex ic = 0; ic < numComp; ++ic )
    {
      dMobility_dC[ic] = dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dC + ic];
    }
  }
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
computeFractionalFlowViscous( localIndex const ip,
                              localIndex const (&seri)[numFluxSupportPoints],
                              localIndex const (&sesri)[numFluxSupportPoints],
                              localIndex const (&sei)[numFluxSupportPoints],
                              real64 const totFlux,
                              real64 const totMob,
                              real64 const (&dTotMob_dP)[numFluxSupportPoints],
                              real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                              ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                              real64 & fractionalFlow,
                              real64 ( & dFractionalFlow_dP )[numFluxSupportPoints],
                              real64 ( & dFractionalFlow_dC )[numFluxSupportPoints][numComp] )
{
  // reinit
  fractionalFlow = 0;
  for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
  {
    dFractionalFlow_dP[ke] = 0;
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dFractionalFlow_dC[ke][jc] = 0;
    }
  }

  localIndex k_up;
  real64 mob{};
  real64 dMob_dP{};
  real64 dMob_dC[numComp]{};

  upwindMobilityViscous< numComp, numFluxSupportPoints, UPWIND >( ip,
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

  // fractional flow too low to let the upstream phase flow
  // guard against no flow region
  if( std::fabs( mob ) > 1e-20 )
  {
    fractionalFlow = mob / LvArray::math::max( totMob, minTotMob );
    dFractionalFlow_dP[k_up] = dMob_dP / LvArray::math::max( totMob, minTotMob );
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dFractionalFlow_dC[k_up][jc] = dMob_dC[jc] / totMob;
    }

    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dFractionalFlow_dP[ke] -= fractionalFlow * dTotMob_dP[ke] / LvArray::math::max( totMob, minTotMob );

      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[ke][jc] -= fractionalFlow * dTotMob_dC[ke][jc] / LvArray::math::max( totMob, minTotMob );
      }
    }
  }
}


template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
computeFractionalFlowGravity( localIndex const numPhase,
                              localIndex const ip,
                              localIndex const (&seri)[numFluxSupportPoints],
                              localIndex const (&sesri)[numFluxSupportPoints],
                              localIndex const (&sei)[numFluxSupportPoints],
                              real64 const (&transmissibility)[2],
                              real64 const (&dTrans_dPres)[2],
                              real64 const totFlux,
                              real64 const totMob,
                              real64 const (&dTotMob_dP)[numFluxSupportPoints],
                              real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                              ElementViewConst< arrayView1d< real64 const > > const & pres,
                              ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                              ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                              ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                              ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                              ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                              ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                              ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                              integer const capPressureFlag,
                              real64 & fractionalFlow,
                              real64 ( & dFractionalFlow_dP)[numFluxSupportPoints],
                              real64 ( & dFractionalFlow_dC)[numFluxSupportPoints][numComp]
                              )
{
  // get var to memorized the numerator mobility properly upwinded
  real64 mainMob{};
  real64 dMMob_dP{};
  real64 dMMob_dC[numComp]{};

  //reinit
  //fractional flow too low to let the upstream phase flow
  localIndex k_up_main = -1;                           //to throw error if unmodified
  fractionalFlow = 0;
  for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
  {
    dFractionalFlow_dP[ke] = 0;
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dFractionalFlow_dC[ke][jc] = 0;
    }
  }


  localIndex k_up;
  real64 mob{};
  real64 dMob_dP{};
  real64 dMob_dC[numComp]{};

  upwindMobilityGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                  ip,
                                                                  seri,
                                                                  sesri,
                                                                  sei,
                                                                  transmissibility,
                                                                  dTrans_dPres,
                                                                  totFlux,
                                                                  pres,
                                                                  gravCoef,
                                                                  dCompFrac_dCompDens,
                                                                  phaseMassDens,
                                                                  dPhaseMassDens,
                                                                  phaseMob,
                                                                  dPhaseMob,
                                                                  dPhaseVolFrac,
                                                                  phaseCapPressure,
                                                                  dPhaseCapPressure_dPhaseVolFrac,
                                                                  capPressureFlag,
                                                                  k_up,
                                                                  mob,
                                                                  dMob_dP,
                                                                  dMob_dC );

  k_up_main = k_up;
  mainMob = mob;
  dMMob_dP = dMob_dP;
  for( localIndex ic = 0; ic < numComp; ++ic )
  {
    dMMob_dC[ic] = dMob_dC[ic];
  }

  //guard against no flow region
  if( std::fabs( mainMob ) > 1e-20 )
  {
    fractionalFlow = mainMob / LvArray::math::max( totMob, minTotMob );
    dFractionalFlow_dP[k_up_main] = dMMob_dP / LvArray::math::max( totMob, minTotMob );
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dFractionalFlow_dC[k_up_main][jc] = dMMob_dC[jc] / totMob;

    }

    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dFractionalFlow_dP[ke] -= fractionalFlow * dTotMob_dP[ke] / LvArray::math::max( totMob, minTotMob );

      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[ke][jc] -= fractionalFlow * dTotMob_dC[ke][jc] / LvArray::math::max( totMob, minTotMob );
      }
    }
  }
}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void
computeFractionalFlowCapillary( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totFlux,
                                real64 const totMob,
                                real64 const (&dTotMob_dP)[numFluxSupportPoints],
                                real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
                                ElementViewConst< arrayView1d< real64 const > > const & pres,
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const capPressureFlag,
                                localIndex & k_up_main,
                                real64 & fractionalFlow,
                                real64 ( & dFractionalFlow_dP)[numFluxSupportPoints],
                                real64 ( & dFractionalFlow_dC)[numFluxSupportPoints][numComp]
                                )
{
  // get var to memorized the numerator mobility properly upwinded
  real64 mainMob{};
  real64 dMMob_dP{};
  real64 dMMob_dC[numComp]{};

  //reinit
  //fractional flow too low to let the upstream phase flow
  k_up_main = -1;                           //to throw error if unmodified
  fractionalFlow = 0;
  for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
  {
    dFractionalFlow_dP[ke] = 0;
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dFractionalFlow_dC[ke][jc] = 0;
    }
  }

  localIndex k_up;
  real64 mob{};
  real64 dMob_dP{};
  real64 dMob_dC[numComp]{};

  upwindMobilityCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                    ip,
                                                                    seri,
                                                                    sesri,
                                                                    sei,
                                                                    transmissibility,
                                                                    dTrans_dPres,
                                                                    totFlux,
                                                                    pres,
                                                                    gravCoef,
                                                                    dCompFrac_dCompDens,
                                                                    phaseMassDens,
                                                                    dPhaseMassDens,
                                                                    phaseMob,
                                                                    dPhaseMob,
                                                                    dPhaseVolFrac,
                                                                    phaseCapPressure,
                                                                    dPhaseCapPressure_dPhaseVolFrac,
                                                                    capPressureFlag,
                                                                    k_up,
                                                                    mob,
                                                                    dMob_dP,
                                                                    dMob_dC );

  k_up_main = k_up;
  mainMob = mob;
  dMMob_dP = dMob_dP;
  for( localIndex ic = 0; ic < numComp; ++ic )
  {
    dMMob_dC[ic] = dMob_dC[ic];
  }

  //guard against no flow region
  if( std::fabs( mainMob ) > 1e-20 )
  {
    fractionalFlow = mainMob / LvArray::math::max( totMob, minTotMob );
    dFractionalFlow_dP[k_up_main] = dMMob_dP / LvArray::math::max( totMob, minTotMob );
    for( localIndex jc = 0; jc < numComp; ++jc )
    {
      dFractionalFlow_dC[k_up_main][jc] = dMMob_dC[jc] / LvArray::math::max( totMob, minTotMob );

    }

    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dFractionalFlow_dP[ke] -= fractionalFlow * dTotMob_dP[ke] / LvArray::math::max( totMob, minTotMob );

      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dFractionalFlow_dC[ke][jc] -= fractionalFlow * dTotMob_dC[ke][jc] / LvArray::math::max( totMob, minTotMob );
      }
    }
  }
}

/**
 * @brief  Struct defining formation of potential from different Physics (flagged by enum type T) to be used
 *            in Upwind discretization schemes
 * @tparam numComp number of component
 * @tparam T the concerned physics (Viscou,Gravity or Capillary)
 * @tparam numFluxSupportPoints number of point in the stencil
 */

/*! @copydoc computePotential
 */
struct computePotentialGravity
{
  /*! @copydoc computePotential::compute
   *
   * @brief specialization for gravitational driving forces which only relies on total flux
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void compute( localIndex const GEOS_UNUSED_PARAM( numPhase ),
                       localIndex const ip,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       real64 const GEOS_UNUSED_PARAM( totFlux ),
                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const &
                       GEOS_UNUSED_PARAM( dPhaseVolFrac ),
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const &
                       GEOS_UNUSED_PARAM( phaseCapPressure ),
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const &
                       GEOS_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ),
                       real64 & pot,
                       real64 ( & dPot_dPres )[numFluxSupportPoints],
                       real64 (& dPot_dComp )[numFluxSupportPoints][numComp],
                       real64 ( & dProp_dComp )[numComp] )
  {
    //working arrays
    real64 densMean{};
    real64 dDensMean_dPres[numFluxSupportPoints]{};
    real64 dDensMean_dComp[numFluxSupportPoints][numComp]{};

    //init
    pot = 0.0;
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      dPot_dPres[i] = 0.0;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dPot_dComp[i][jc] = 0.0;
        dProp_dComp[jc] = 0.0;
      }
    }

    //inner loop to get average density
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      // density
      real64 const density = phaseMassDens[er][esr][ei][0][ip];
      real64 const dDens_dPres = dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];

      applyChainRule( numComp,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens[er][esr][ei][0][ip],
                      dProp_dComp,
                      Deriv::dC );

      // average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dPres[i] = 0.5 * dDens_dPres;
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dComp[i][jc] = 0.5 * dProp_dComp[jc];
      }
    }

    // compute potential difference MPFA-style
    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      real64 const gravD = transmissibility[i] * gravCoef[er][esr][ei];
      real64 const dGravD_dP = dTrans_dPres[i] * gravCoef[er][esr][ei];
      pot += densMean * gravD;

      // need to add contributions from both cells the mean density depends on
      for( localIndex j = 0; j < numFluxSupportPoints; ++j )
      {
        dPot_dPres[j] += dDensMean_dPres[j] * gravD + densMean * dGravD_dP;
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPot_dComp[j][jc] += dDensMean_dComp[j][jc] * gravD;
        }
      }
    }

  }
};

/*! @copydoc computePotential
 */
struct computePotentialCapillary
{
  /*! @copydoc computePotential::compute
   *
   * @brief specialization for capillary driving forces which only relies on total flux
   */
  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static void compute( localIndex const numPhase,
                       localIndex const ip,
                       localIndex const (&seri)[numFluxSupportPoints],
                       localIndex const (&sesri)[numFluxSupportPoints],
                       localIndex const (&sei)[numFluxSupportPoints],
                       real64 const (&transmissibility)[2],
                       real64 const (&dTrans_dPres)[2],
                       real64 const GEOS_UNUSED_PARAM( totFlux ),
                       ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( gravCoef ),
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const &
                       GEOS_UNUSED_PARAM( dCompFrac_dCompDens ),
                       ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const &
                       GEOS_UNUSED_PARAM( phaseMassDens ),
                       ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const &
                       GEOS_UNUSED_PARAM( dPhaseMassDens ),
                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                       ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                       ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                       real64 & pot,
                       real64 ( & dPot_dPres)[numFluxSupportPoints],
                       real64 (& dPot_dComp)[numFluxSupportPoints][numComp],
                       real64( &GEOS_UNUSED_PARAM( dProp_dComp ))[numComp] )
  {


    for( localIndex i = 0; i < numFluxSupportPoints; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      pot += transmissibility[i] * phaseCapPressure[er][esr][ei][0][ip];
      // need to add contributions from both cells
      for( localIndex jp = 0; jp < numPhase; ++jp )
      {

        real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
        dPot_dPres[i] +=
          transmissibility[i] * dCapPressure_dS * dPhaseVolFrac[er][esr][ei][jp][Deriv::dP]
          + dTrans_dPres[i] * phaseCapPressure[er][esr][ei][0][jp];

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPot_dComp[i][jc] += transmissibility[i] * dCapPressure_dS *
                               dPhaseVolFrac[er][esr][ei][jp][Deriv::dC + jc];
        }
      }
    }
  }
};

/// Form potential-related parts of fluxes

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void computePotentialFluxesGravity( localIndex const numPhase,
                                           localIndex const ip,
                                           localIndex const (&seri)[numFluxSupportPoints],
                                           localIndex const (&sesri)[numFluxSupportPoints],
                                           localIndex const (&sei)[numFluxSupportPoints],
                                           real64 const (&transmissibility)[2],
                                           real64 const (&dTrans_dPres)[2],
                                           real64 const totFlux,
                                           real64 const totMob,
                                           real64 const (&dTotMob_dP)[numFluxSupportPoints],
                                           real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
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
                                           localIndex const capPressureFlag,
                                           localIndex (&k_up_o),
                                           real64 & phaseFlux,
                                           real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                                           real64 ( & dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
{

  real64 fflow{};
  real64 dFflow_dP[numFluxSupportPoints]{};
  real64 dFflow_dC[numFluxSupportPoints][numComp]{};

  real64 pot{};
  real64 dPot_dP[numFluxSupportPoints]{};
  real64 dPot_dC[numFluxSupportPoints][numComp]{};
  real64 dProp_dC[numComp]{};

  //
  UpwindHelpers::computePotentialGravity::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                                    ip,
                                                                                    seri,
                                                                                    sesri,
                                                                                    sei,
                                                                                    transmissibility,
                                                                                    dTrans_dPres,
                                                                                    totFlux,
                                                                                    gravCoef,
                                                                                    dCompFrac_dCompDens,
                                                                                    phaseMassDens,
                                                                                    dPhaseMassDens,
                                                                                    dPhaseVolFrac,
                                                                                    phaseCapPressure,
                                                                                    dPhaseCapPressure_dPhaseVolFrac,
                                                                                    pot,
                                                                                    dPot_dP,
                                                                                    dPot_dC,
                                                                                    dProp_dC );

  // and the fractional flow for gravitational part as \lambda_i^{up}/\sum_{numPhase}(\lambda_k^{up}) with up decided upon
  // the Upwind strategy
  UpwindHelpers::computeFractionalFlowGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                        ip,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        totMob,
                                                                                        dTotMob_dP,
                                                                                        dTotMob_dC,
                                                                                        pres,
                                                                                        gravCoef,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        phaseMob,
                                                                                        dPhaseMob,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        capPressureFlag,
                                                                                        fflow,
                                                                                        dFflow_dP,
                                                                                        dFflow_dC );


  for( localIndex jp = 0; jp < numPhase; ++jp )
  {
    if( ip != jp )
    {

      real64 potOther{};
      real64 dPotOther_dP[numFluxSupportPoints]{};
      real64 dPotOther_dC[numFluxSupportPoints][numComp]{};
      real64 dPropOther_dC[numComp]{};

      //Fetch pot for phase j!=i defined as \rho_j g dz/dx
      UpwindHelpers::computePotentialGravity::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                                        jp,
                                                                                        seri,
                                                                                        sesri,
                                                                                        sei,
                                                                                        transmissibility,
                                                                                        dTrans_dPres,
                                                                                        totFlux,
                                                                                        gravCoef,
                                                                                        dCompFrac_dCompDens,
                                                                                        phaseMassDens,
                                                                                        dPhaseMassDens,
                                                                                        dPhaseVolFrac,
                                                                                        phaseCapPressure,
                                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                                        potOther,
                                                                                        dPotOther_dP,
                                                                                        dPotOther_dC,
                                                                                        dPropOther_dC );

      //Eventually get the mobility of the second phase
      real64 mobOther{};
      real64 dMobOther_dP{};
      real64 dMobOther_dC[numComp]{};

      // and the other mobility for gravitational part as \lambda_j^{up} with up decided upon
      // the Upwind strategy - Note that it should be the same as the gravitational fractional flow

      UpwindHelpers::upwindMobilityGravity< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                     jp,
                                                                                     seri,
                                                                                     sesri,
                                                                                     sei,
                                                                                     transmissibility,
                                                                                     dTrans_dPres,
                                                                                     totFlux,
                                                                                     pres,
                                                                                     gravCoef,
                                                                                     dCompFrac_dCompDens,
                                                                                     phaseMassDens,
                                                                                     dPhaseMassDens,
                                                                                     phaseMob,
                                                                                     dPhaseMob,
                                                                                     dPhaseVolFrac,
                                                                                     phaseCapPressure,
                                                                                     dPhaseCapPressure_dPhaseVolFrac,
                                                                                     capPressureFlag,
                                                                                     k_up_o,
                                                                                     mobOther,
                                                                                     dMobOther_dP,
                                                                                     dMobOther_dC );


      // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
      phaseFlux -= fflow * mobOther * (pot - potOther);
      dPhaseFlux_dP[k_up_o] -= fflow * dMobOther_dP * (pot - potOther);
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[k_up_o][jc] -= fflow * dMobOther_dC[jc] * (pot - potOther);
      }

      //mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dP[ke] -= dFflow_dP[ke] * mobOther * (pot - potOther);

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] -= dFflow_dC[ke][jc] * mobOther * (pot - potOther);
        }
      }

      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dP[ke] -= fflow * mobOther * (dPot_dP[ke] - dPotOther_dP[ke]);
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] -= fflow * mobOther * (dPot_dC[ke][jc] - dPotOther_dC[ke][jc]);
        }
      }
    }
  }

}

template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
GEOS_HOST_DEVICE
static void computePotentialFluxesCapillary( localIndex const numPhase,
                                             localIndex const ip,
                                             localIndex const (&seri)[numFluxSupportPoints],
                                             localIndex const (&sesri)[numFluxSupportPoints],
                                             localIndex const (&sei)[numFluxSupportPoints],
                                             real64 const (&transmissibility)[2],
                                             real64 const (&dTrans_dPres)[2],
                                             real64 const totFlux,
                                             real64 const totMob,
                                             real64 const (&dTotMob_dP)[numFluxSupportPoints],
                                             real64 const (&dTotMob_dC)[numFluxSupportPoints][numComp],
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
                                             localIndex const capPressureFlag,
                                             localIndex( &k_up),
                                             localIndex (&k_up_o),
                                             real64 & phaseFlux,
                                             real64 (& dPhaseFlux_dP)[numFluxSupportPoints],
                                             real64 ( & dPhaseFlux_dC)[numFluxSupportPoints][numComp] )
{

  real64 fflow{};
  real64 dFflow_dP[numFluxSupportPoints]{};
  real64 dFflow_dC[numFluxSupportPoints][numComp]{};

  real64 pot{};
  real64 dPot_dP[numFluxSupportPoints]{};
  real64 dPot_dC[numFluxSupportPoints][numComp]{};
  real64 dProp_dC[numComp]{};

  UpwindHelpers::computePotentialCapillary::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                                      ip,
                                                                                      seri,
                                                                                      sesri,
                                                                                      sei,
                                                                                      transmissibility,
                                                                                      dTrans_dPres,
                                                                                      totFlux,
                                                                                      gravCoef,
                                                                                      dCompFrac_dCompDens,
                                                                                      phaseMassDens,
                                                                                      dPhaseMassDens,
                                                                                      dPhaseVolFrac,
                                                                                      phaseCapPressure,
                                                                                      dPhaseCapPressure_dPhaseVolFrac,
                                                                                      pot,
                                                                                      dPot_dP,
                                                                                      dPot_dC,
                                                                                      dProp_dC );

  // and the fractional flow for gravitational part as \lambda_i^{up}/\sum_{numPhase}(\lambda_k^{up}) with up decided upon
  // the Upwind strategy
  UpwindHelpers::computeFractionalFlowCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                          ip,
                                                                                          seri,
                                                                                          sesri,
                                                                                          sei,
                                                                                          transmissibility,
                                                                                          dTrans_dPres,
                                                                                          totFlux,
                                                                                          totMob,
                                                                                          dTotMob_dP,
                                                                                          dTotMob_dC,
                                                                                          pres,
                                                                                          gravCoef,
                                                                                          dCompFrac_dCompDens,
                                                                                          phaseMassDens,
                                                                                          dPhaseMassDens,
                                                                                          phaseMob,
                                                                                          dPhaseMob,
                                                                                          dPhaseVolFrac,
                                                                                          phaseCapPressure,
                                                                                          dPhaseCapPressure_dPhaseVolFrac,
                                                                                          capPressureFlag,
                                                                                          k_up,
                                                                                          fflow,
                                                                                          dFflow_dP,
                                                                                          dFflow_dC );


  for( localIndex jp = 0; jp < numPhase; ++jp )
  {
    if( ip != jp )
    {

      real64 potOther{};
      real64 dPotOther_dP[numFluxSupportPoints]{};
      real64 dPotOther_dC[numFluxSupportPoints][numComp]{};
      real64 dPropOther_dC[numComp]{};

      //Fetch pot for phase j!=i defined as \rho_j g dz/dx
      UpwindHelpers::computePotentialCapillary::compute< numComp, numFluxSupportPoints >( numPhase,
                                                                                          jp,
                                                                                          seri,
                                                                                          sesri,
                                                                                          sei,
                                                                                          transmissibility,
                                                                                          dTrans_dPres,
                                                                                          totFlux,
                                                                                          gravCoef,
                                                                                          dCompFrac_dCompDens,
                                                                                          phaseMassDens,
                                                                                          dPhaseMassDens,
                                                                                          dPhaseVolFrac,
                                                                                          phaseCapPressure,
                                                                                          dPhaseCapPressure_dPhaseVolFrac,
                                                                                          potOther,
                                                                                          dPotOther_dP,
                                                                                          dPotOther_dC,
                                                                                          dPropOther_dC );

      //Eventually get the mobility of the second phase
      real64 mobOther{};
      real64 dMobOther_dP{};
      real64 dMobOther_dC[numComp]{};

      // and the other mobility for gravitational part as \lambda_j^{up} with up decided upon
      // the Upwind strategy - Note that it should be the same as the gravitational fractional flow

      UpwindHelpers::upwindMobilityCapillary< numComp, numFluxSupportPoints, UPWIND >( numPhase,
                                                                                       jp,
                                                                                       seri,
                                                                                       sesri,
                                                                                       sei,
                                                                                       transmissibility,
                                                                                       dTrans_dPres,
                                                                                       totFlux,
                                                                                       pres,
                                                                                       gravCoef,
                                                                                       dCompFrac_dCompDens,
                                                                                       phaseMassDens,
                                                                                       dPhaseMassDens,
                                                                                       phaseMob,
                                                                                       dPhaseMob,
                                                                                       dPhaseVolFrac,
                                                                                       phaseCapPressure,
                                                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                                                       capPressureFlag,
                                                                                       k_up_o,
                                                                                       mobOther,
                                                                                       dMobOther_dP,
                                                                                       dMobOther_dC );


      // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
      phaseFlux -= fflow * mobOther * (pot - potOther);
      dPhaseFlux_dP[k_up_o] -= fflow * dMobOther_dP * (pot - potOther);
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[k_up_o][jc] -= fflow * dMobOther_dC[jc] * (pot - potOther);
      }

      //mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dP[ke] -= dFflow_dP[ke] * mobOther * (pot - potOther);

        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] -= dFflow_dC[ke][jc] * mobOther * (pot - potOther);
        }
      }

      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dP[ke] -= fflow * mobOther * (dPot_dP[ke] - dPotOther_dP[ke]);
        for( localIndex jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[ke][jc] -= fflow * mobOther * (dPot_dC[ke][jc] - dPotOther_dC[ke][jc]);
        }
      }
    }
  }
}

}//end of struct UpwindHelpers

/************************* UPWIND ******************/

/**
 * @brief Template base class for different upwind Scheme
 * @tparam T physics concerned by the scheme if specialized
 */
class UpwindScheme
{

public:

  //default ctor
  UpwindScheme() = default;

  //usual copy ctor
  UpwindScheme( UpwindScheme const & scheme ) = default;

  //default move ctor
  UpwindScheme( UpwindScheme && ) = default;

  //deleted copy and move assignement
  UpwindScheme & operator=( UpwindScheme const & ) = delete;

  UpwindScheme & operator=( UpwindScheme && ) = delete;

  virtual ~UpwindScheme() = default;

  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOS_HOST_DEVICE
  void getUpwindDirectionGravity( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totFlux,
                                  ElementViewConst< arrayView1d< real64 const > > const & pres,
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const capPressureFlag,
                                  localIndex & upwindDir
                                  )
  {
    real64 pot{};

    /// each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    /// lamda defining how to get these potentials
    UPWIND::template computePotentialGravity< numComp, numFluxSupportPoints >( numPhase,
                                                                               ip,
                                                                               seri,
                                                                               sesri,
                                                                               sei,
                                                                               transmissibility,
                                                                               dTrans_dPres,
                                                                               totFlux,
                                                                               pres,
                                                                               gravCoef,
                                                                               phaseMob,
                                                                               dCompFrac_dCompDens,
                                                                               phaseMassDens,
                                                                               dPhaseMassDens,
                                                                               dPhaseVolFrac,
                                                                               phaseCapPressure,
                                                                               dPhaseCapPressure_dPhaseVolFrac,
                                                                               capPressureFlag,
                                                                               pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot >= 0) ? 0 : 1;
  }


  template< localIndex numComp, localIndex numFluxSupportPoints, class UPWIND >
  GEOS_HOST_DEVICE
  void getUpwindDirectionCapillary( localIndex const numPhase,
                                    localIndex const ip,
                                    localIndex const (&seri)[numFluxSupportPoints],
                                    localIndex const (&sesri)[numFluxSupportPoints],
                                    localIndex const (&sei)[numFluxSupportPoints],
                                    real64 const (&transmissibility)[2],
                                    real64 const (&dTrans_dPres)[2],
                                    real64 const totFlux,                   //in fine should be a ElemnetViewConst once seq form are in
                                                                            // place
                                    ElementViewConst< arrayView1d< real64 const > > const & pres,
                                    ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                    ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                    ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                    ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                    ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                    ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                    integer const capPressureFlag,
                                    localIndex & upwindDir
                                    )
  {
    real64 pot{};

    // each derived concrete class has to define a computePotential method that is calling UpwindScheme::potential method with a specific
    // lamda defining how to get these potentials
    UPWIND::template computePotentialCapillary< numComp, numFluxSupportPoints >( numPhase,
                                                                                 ip,
                                                                                 seri,
                                                                                 sesri,
                                                                                 sei,
                                                                                 transmissibility,
                                                                                 dTrans_dPres,
                                                                                 totFlux,
                                                                                 pres,
                                                                                 gravCoef,
                                                                                 phaseMob,
                                                                                 dCompFrac_dCompDens,
                                                                                 phaseMassDens,
                                                                                 dPhaseMassDens,
                                                                                 dPhaseVolFrac,
                                                                                 phaseCapPressure,
                                                                                 dPhaseCapPressure_dPhaseVolFrac,
                                                                                 capPressureFlag,
                                                                                 pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = (pot >= 0) ? 0 : 1;
  }



  // templated way of evaluating the potential (to the exception of viscous one) which relies on
  // up-or-downwinded mobility terms pre-multiplying potential differences
  template< localIndex numComp, localIndex numFluxSupportPoints, typename LAMBDA >
  GEOS_HOST_DEVICE
  static void potential( localIndex numPhase,
                         localIndex ip,
                         localIndex const (&seri)[numFluxSupportPoints],
                         localIndex const (&sesri)[numFluxSupportPoints],
                         localIndex const (&sei)[numFluxSupportPoints],
                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                         real64 & weightedPotential,
                         LAMBDA && fn )
  {
    //getPhase Pot
    real64 pot{};
    real64 pot_dP[numFluxSupportPoints]{};
    real64 pot_dC[numFluxSupportPoints][numComp]{};
    real64 dProp_dC[numComp]{};

    fn( ip, pot, pot_dP, pot_dC, dProp_dC );

    localIndex const k_up = 0;
    localIndex const k_dw = 1;

    //loop other other phases to form
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( jp != ip )
      {
        localIndex const er_up = seri[k_up];
        localIndex const esr_up = sesri[k_up];
        localIndex const ei_up = sei[k_up];

        localIndex const er_dw = seri[k_dw];
        localIndex const esr_dw = sesri[k_dw];
        localIndex const ei_dw = sei[k_dw];

        real64 potOther{};
        real64 potOther_dP[numFluxSupportPoints]{};
        real64 potOther_dC[numFluxSupportPoints][numComp]{};
        real64 dPropOther_dC[numComp]{};

        fn( jp, potOther, potOther_dP, potOther_dC, dPropOther_dC );

        real64 const mob_up = phaseMob[er_up][esr_up][ei_up][jp];
        real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][jp];

        weightedPotential += (pot - potOther >= 0) ? mob_dw * (potOther - pot) : mob_up * (potOther - pot);

      }
    }
  }

};

/**
 * @brief  Class describing the Hybrid Upwind scheme as defined in "Consistent upwinding for sequential fully implicit
 *         multiscale compositional simulation" (Moncorge,2020)
 */
class HybridUpwind : public UpwindScheme
{

public:

  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static
  void computePotentialGravity( localIndex const numPhase,
                                localIndex const ip,
                                localIndex const (&seri)[numFluxSupportPoints],
                                localIndex const (&sesri)[numFluxSupportPoints],
                                localIndex const (&sei)[numFluxSupportPoints],
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const totalFlux,
                                ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( pres ),
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                integer const GEOS_UNUSED_PARAM( capPressureFlag ),
                                real64 & potential
                                )
  {

    //Form total velocity
    potential = 0;

    //the arg lambda allows us to access some genericity
    UpwindScheme::template potential< numComp, numFluxSupportPoints >( numPhase, ip, seri, sesri, sei,
                                                                       phaseMob, potential,
                                                                       [&]( localIndex ipp,
                                                                            real64 & potential_,
                                                                            real64 (& dPotential_dP_)[numFluxSupportPoints],
                                                                            real64 (& dPotential_dC_)[numFluxSupportPoints][numComp],
                                                                            real64 (& dProp_dC)[numComp] ) {

      UpwindHelpers::computePotentialGravity::compute< numComp, numFluxSupportPoints >(
        numPhase,
        ipp,
        seri,
        sesri,
        sei,
        transmissibility,
        dTrans_dPres,
        totalFlux,
        gravCoef,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        dPhaseVolFrac,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        potential_,
        dPotential_dP_,
        dPotential_dC_,
        dProp_dC );

    } );
  }


  template< localIndex numComp, localIndex numFluxSupportPoints >
  GEOS_HOST_DEVICE
  static
  void computePotentialCapillary( localIndex const numPhase,
                                  localIndex const ip,
                                  localIndex const (&seri)[numFluxSupportPoints],
                                  localIndex const (&sesri)[numFluxSupportPoints],
                                  localIndex const (&sei)[numFluxSupportPoints],
                                  real64 const (&transmissibility)[2],
                                  real64 const (&dTrans_dPres)[2],
                                  real64 const totalFlux,
                                  ElementViewConst< arrayView1d< real64 const > > const & GEOS_UNUSED_PARAM( pres ),
                                  ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const &
                                  phaseMob,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                  ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                  ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                  ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                  integer const GEOS_UNUSED_PARAM( capPressureFlag ),
                                  real64 & potential )

  {

    //Form total velocity
    potential = 0;

    //the arg lambda allows us to access some genericity
    UpwindScheme::template potential< numComp, numFluxSupportPoints >( numPhase, ip, seri, sesri, sei,
                                                                       phaseMob, potential,
                                                                       [&]( localIndex ipp,
                                                                            real64 & potential_,
                                                                            real64 (& dPotential_dP_)[numFluxSupportPoints],
                                                                            real64 (& dPotential_dC_)[numFluxSupportPoints][numComp],
                                                                            real64 (& dProp_dC)[numComp] ) {

      UpwindHelpers::computePotentialCapillary::compute< numComp, numFluxSupportPoints >(
        numPhase,
        ipp,
        seri,
        sesri,
        sei,
        transmissibility,
        dTrans_dPres,
        totalFlux,
        gravCoef,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        dPhaseVolFrac,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        potential_,
        dPotential_dP_,
        dPotential_dC_,
        dProp_dC );
    } );
  }

};

/*** IHU ***/

struct IHUPhaseFlux
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
           real64 & potGrad,
           real64 ( &phaseFlux ),
           real64 ( & dPhaseFlux_dP )[numFluxSupportPoints],
           real64 ( & dPhaseFlux_dC )[numFluxSupportPoints][numComp] )
  {

    //loop over all phases to form total velocity
    real64 totFlux{};
    real64 dTotFlux_dP[numFluxSupportPoints]{};
    real64 dTotFlux_dC[numFluxSupportPoints][numComp]{};

    //store totMob upwinded by PPU for later schemes
    real64 totMob{};
    real64 dTotMob_dP[numFluxSupportPoints]{};
    real64 dTotMob_dC[numFluxSupportPoints][numComp]{};

    for( integer jp = 0; jp < numPhase; ++jp )
    {
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

      // choose upstream cell
      localIndex const k_up = (phaseFlux >= 0) ? 0 : 1;
      localIndex const er_up = seri[k_up];
      localIndex const esr_up = sesri[k_up];
      localIndex const ei_up = sei[k_up];
      totMob += phaseMob[er_up][esr_up][ei_up][jp];
      dTotMob_dP[k_up] += dPhaseMob[er_up][esr_up][ei_up][jp][Deriv::dP];
      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dTotMob_dC[k_up][jc] += dPhaseMob[er_up][esr_up][ei_up][jp][Deriv::dC + jc];
      }
    }

    // assign to zero
    phaseFlux = 0;
    for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] = 0;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseFlux_dC[ke][jc] = 0;
      }
    }

    // fractional flow loop with IHU
    // maybe needed to have density out for upwinding

    // choose upstream cell
    // create local work arrays
    real64 viscousPhaseFlux{};
    real64 dViscousPhaseFlux_dP[numFluxSupportPoints]{};
    real64 dViscousPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    real64 fractionalFlow{};
    real64 dFractionalFlow_dP[numFluxSupportPoints]{};
    real64 dFractionalFlow_dC[numFluxSupportPoints][numComp]{};

    // and the fractional flow for viscous part as \lambda_i^{up}/\sum_{NP}(\lambda_j^{up}) with up decided upon
    // the Upwind strategy
    UpwindHelpers::computeFractionalFlowViscous< numComp, numFluxSupportPoints,
                                                 UPWIND_SCHEME >( ip,
                                                                  seri,
                                                                  sesri,
                                                                  sei,
                                                                  totFlux,
                                                                  totMob,
                                                                  dTotMob_dP,
                                                                  dTotMob_dC,
                                                                  phaseMob,
                                                                  dPhaseMob,
                                                                  fractionalFlow,
                                                                  dFractionalFlow_dP,
                                                                  dFractionalFlow_dC );

    /// Assembling the viscous flux (and derivatives) from fractional flow and total velocity as \phi_{\mu} = f_i^{up,\mu} uT
    viscousPhaseFlux = fractionalFlow * totFlux;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dViscousPhaseFlux_dP[ke] += dFractionalFlow_dP[ke] * totFlux;

      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dViscousPhaseFlux_dC[ke][jc] += dFractionalFlow_dC[ke][jc] * totFlux;
      }
    }

    // NON-FIXED UT -- to be canceled out if considered fixed
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dViscousPhaseFlux_dP[ke] += fractionalFlow * dTotFlux_dP[ke];

      for( localIndex jc = 0; jc < numComp; ++jc )
      {
        dViscousPhaseFlux_dC[ke][jc] += fractionalFlow * dTotFlux_dC[ke][jc];
      }
    }

    // accumulate in the flux and its derivatives
    phaseFlux += viscousPhaseFlux;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += dViscousPhaseFlux_dP[ke];

      for( localIndex ic = 0; ic < numComp; ++ic )
        dPhaseFlux_dC[ke][ic] += dViscousPhaseFlux_dC[ke][ic];
    }

    /// Assembling the gravitational flux (and derivatives) from fractional flow and total velocity as \phi_{g} = f_i^{up,g} uT
    localIndex k_up_og = -1;

    real64 gravitationalPhaseFlux{};
    real64 gravitationalPhaseFlux_dP[numFluxSupportPoints]{};
    real64 gravitationalPhaseFlux_dC[numFluxSupportPoints][numComp]{};

    UpwindHelpers::computePotentialFluxesGravity< numComp,
                                                  numFluxSupportPoints, UPWIND_SCHEME >(
      numPhase,
      ip,
      seri,
      sesri,
      sei,
      trans,
      dTrans_dPres,
      totFlux,
      totMob,
      dTotMob_dP,
      dTotMob_dC,
      pres,
      gravCoef,
      phaseMob,
      dPhaseMob,
      dPhaseVolFrac,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      hasCapPressure,
      k_up_og,
      gravitationalPhaseFlux,
      gravitationalPhaseFlux_dP,
      gravitationalPhaseFlux_dC );

    //update phaseFlux from gravitational
    phaseFlux += gravitationalPhaseFlux;
    for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
    {
      dPhaseFlux_dP[ke] += gravitationalPhaseFlux_dP[ke];
      for( localIndex ic = 0; ic < numComp; ++ic )
        dPhaseFlux_dC[ke][ic] += gravitationalPhaseFlux_dC[ke][ic];
    }

    if( hasCapPressure )
    {
      /// Assembling the capillary flux (and derivatives) from fractional flow and total velocity as \phi_{g} = f_i^{up,g} uT
      localIndex k_up_pc = -1;
      localIndex k_up_opc = -1;

      real64 capillaryPhaseFlux{};
      real64 capillaryPhaseFlux_dP[numFluxSupportPoints]{};
      real64 capillaryPhaseFlux_dC[numFluxSupportPoints][numComp]{};

      UpwindHelpers::computePotentialFluxesCapillary< numComp,
                                                      numFluxSupportPoints, UPWIND_SCHEME >(
        numPhase,
        ip,
        seri,
        sesri,
        sei,
        trans,
        dTrans_dPres,
        totFlux,
        totMob,
        dTotMob_dP,
        dTotMob_dC,
        pres,
        gravCoef,
        phaseMob,
        dPhaseMob,
        dPhaseVolFrac,
        dCompFrac_dCompDens,
        phaseMassDens,
        dPhaseMassDens,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        hasCapPressure,
        k_up_pc,
        k_up_opc,
        capillaryPhaseFlux,
        capillaryPhaseFlux_dP,
        capillaryPhaseFlux_dC );

      //update phaseFlux from capillary
      phaseFlux += capillaryPhaseFlux;
      for( localIndex ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dP[ke] += capillaryPhaseFlux_dP[ke];
        for( localIndex ic = 0; ic < numComp; ++ic )
          dPhaseFlux_dC[ke][ic] += capillaryPhaseFlux_dC[ke][ic];

      }

    }//end if cappres

  }

};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_IHUPHASEFLUX_HPP
