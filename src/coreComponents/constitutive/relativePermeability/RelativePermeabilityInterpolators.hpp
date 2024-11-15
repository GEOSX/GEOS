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
 * @file RelativePermeabilityInterpolators.hpp
 */

#ifndef GEOS_CONSTITUTIVE_RELATIVEPERMEABILITYINTERPOLATORS_HPP
#define GEOS_CONSTITUTIVE_RELATIVEPERMEABILITYINTERPOLATORS_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

namespace geos
{

namespace constitutive
{

namespace relpermInterpolators
{

struct Baker
{

  /**
   * @brief Interpolate the two-phase relperms to compute the three-phase relperm
   * @param[in] shiftedWaterVolFrac
   * @param[in] gasVolFrac
   * @param[out] threePhaseRelPerm
   * @param[out] dThreePhaseRelPerm_dVolFrac
   * @param[in] woRelPerm
   * @param[in] dWoRelPerm_dOilVolFrac
   * @param[in] goRelPerm
   * @param[in] dGoRelPerm_dOilVolFrac
   *
   * This function interpolates the two-phase relperms to compute the three-phase relperm
   * The interpolation is based on the modified Baker method, also used as default in Eclipse
   * Reference: Eclipse technical description and PetroWiki
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  compute( real64 const & shiftedWaterVolFrac,
           real64 const & gasVolFrac,
           arraySlice1d< integer const > const & phaseOrder,
           real64 const & woRelPerm,
           real64 const & dWoRelPerm_dOilVolFrac,
           real64 const & goRelPerm,
           real64 const & dGoRelPerm_dOilVolFrac,
           real64 & threePhaseRelPerm,
           arraySlice1d< real64, constitutive::relperm::USD_RELPERM_DS - 3 > const & dThreePhaseRelPerm_dVolFrac )
  {
    using PT = RelativePermeabilityBase::PhaseType;
    integer const ipWater = phaseOrder[PT::WATER];
    integer const ipOil   = phaseOrder[PT::OIL];
    integer const ipGas   = phaseOrder[PT::GAS];

    // if water phase is immobile, then use the two-phase gas-oil data only
    if( shiftedWaterVolFrac <= 0.0 )
    {
      threePhaseRelPerm = goRelPerm;
      dThreePhaseRelPerm_dVolFrac[ipOil] = dGoRelPerm_dOilVolFrac;
    }
    // if gas phase is immobile, then use the two-phase water-oil data only
    else if( gasVolFrac <= 0.0 )
    {
      threePhaseRelPerm = woRelPerm;
      dThreePhaseRelPerm_dVolFrac[ipOil] = dWoRelPerm_dOilVolFrac;
    }
    // if both the water phase and the gas phase are mobile,
    // then use a saturation-weighted interpolation of the two-phase oil rel perms
    else
    {
      real64 const sumRelPerm = (shiftedWaterVolFrac * woRelPerm
                                 + gasVolFrac   * goRelPerm);
      real64 const dSumRelPerm_dWaterVolFrac = woRelPerm;
      real64 const dSumRelPerm_dOilVolFrac   = shiftedWaterVolFrac * dWoRelPerm_dOilVolFrac
                                               + gasVolFrac   * dGoRelPerm_dOilVolFrac;
      real64 const dSumRelPerm_dGasVolFrac   = goRelPerm;


      real64 const sumVolFrac    = shiftedWaterVolFrac + gasVolFrac;
      real64 const sumVolFracInv = 1 / sumVolFrac; // div by 0 handled by the if statement above
      real64 const dSumVolFracInv_dWaterVolFrac = -sumVolFracInv * sumVolFracInv;
      real64 const dSumVolFracInv_dGasVolFrac   = dSumVolFracInv_dWaterVolFrac;

      // three-phase oil rel perm
      threePhaseRelPerm = sumRelPerm * sumVolFracInv;
      // derivative w.r.t. Sw
      dThreePhaseRelPerm_dVolFrac[ipWater] = dSumRelPerm_dWaterVolFrac * sumVolFracInv
                                             + sumRelPerm                * dSumVolFracInv_dWaterVolFrac;
      // derivative w.r.t. So
      dThreePhaseRelPerm_dVolFrac[ipOil]   = dSumRelPerm_dOilVolFrac   * sumVolFracInv;
      // derivative w.r.t. Sg
      dThreePhaseRelPerm_dVolFrac[ipGas]   = dSumRelPerm_dGasVolFrac   * sumVolFracInv
                                             + sumRelPerm                * dSumVolFracInv_dGasVolFrac;
    }
  }

};

struct Stone2
{
  /**
   * @brief Interpolate the two-phase relperms to compute the three-phase relperm
   * @param[in] shiftedWaterVolFrac
   * @param[in] gasVolFrac
   * @param[in] woRelPerm
   * @param[in] dWoRelPerm_dOilVolFrac
   * @param[in] connatewoRelPerm
   * @param[in] goRelPerm
   * @param[in] dGoRelPerm_dOilVolFrac
   * @param[in] wRelPerm
   * @param[in] dWRelPerm_dWaterVolFrac
   * @param[in] gRelPerm
   * @param[in] dGRelPerm_dGasVolFrac
   * @param[out] threePhaseRelPerm
   * @param[out] dThreePhaseRelPerm_dVolFrac
   *
   * This function interpolates the two-phase relperms to compute the three-phase relperm
   * The interpolation is based on the modified Stone 2 method
   * Reference: Eclipse technical description
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( real64 const & shiftedWaterVolFrac,
                       real64 const & gasVolFrac,
                       arraySlice1d< integer const > const & phaseOrder,
                       real64 const & connatewoRelPerm,
                       real64 const & woRelPerm,
                       real64 const & dWoRelPerm_dOilVolFrac,
                       real64 const & goRelPerm,
                       real64 const & dGoRelPerm_dOilVolFrac,
                       real64 const & wRelPerm,
                       real64 const & dWRelPerm_dWaterVolFrac,
                       real64 const & gRelPerm,
                       real64 const & dGRelPerm_dGasVolFrac,
                       real64 & threePhaseRelPerm,
                       arraySlice1d< real64, constitutive::relperm::USD_RELPERM_DS - 3 > const & dThreePhaseRelPerm_dVolFrac )
  {

    using PT = RelativePermeabilityBase::PhaseType;
    integer const ipWater = phaseOrder[PT::WATER];
    integer const ipOil   = phaseOrder[PT::OIL];
    integer const ipGas   = phaseOrder[PT::GAS];

    if( gasVolFrac <= 0 )
    {
      threePhaseRelPerm = woRelPerm;
      dThreePhaseRelPerm_dVolFrac[ipOil] = dWoRelPerm_dOilVolFrac;
    }
    else if( shiftedWaterVolFrac <= 0 )
    {
      threePhaseRelPerm = goRelPerm;
      dThreePhaseRelPerm_dVolFrac[ipOil] = dGoRelPerm_dOilVolFrac;
    }
    else
    {
      threePhaseRelPerm = connatewoRelPerm * ((woRelPerm/connatewoRelPerm + wRelPerm)*(goRelPerm/connatewoRelPerm + gRelPerm) - (wRelPerm + gRelPerm));
      // derivative w.r.t. Sw
      dThreePhaseRelPerm_dVolFrac[ipWater] = dWRelPerm_dWaterVolFrac * (goRelPerm + gRelPerm*connatewoRelPerm - connatewoRelPerm);
      // derivative w.r.t. So
      dThreePhaseRelPerm_dVolFrac[ipOil]   = dWoRelPerm_dOilVolFrac/connatewoRelPerm*goRelPerm + woRelPerm*dGoRelPerm_dOilVolFrac/connatewoRelPerm + wRelPerm*dGoRelPerm_dOilVolFrac + gRelPerm*
                                             dWoRelPerm_dOilVolFrac;
      // derivative w.r.t. Sg
      dThreePhaseRelPerm_dVolFrac[ipGas] = dGRelPerm_dGasVolFrac * (woRelPerm + wRelPerm*connatewoRelPerm - connatewoRelPerm);
      if( threePhaseRelPerm < 0 )
      {
        threePhaseRelPerm = 0;
        dThreePhaseRelPerm_dVolFrac[ipWater] = 0;
        // derivative w.r.t. So
        dThreePhaseRelPerm_dVolFrac[ipOil]   = 0;
        // derivative w.r.t. Sg
        dThreePhaseRelPerm_dVolFrac[ipGas] = 0;
      }
    }


  }
};


} // namespace relpermInterpolators

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_RELATIVEPERMEABILITYINTERPOLATORS_HPP
