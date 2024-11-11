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
 * @file ThermalPhaseVolumeFractionKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALPHASEVOLUMEFRACTIONKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALPHASEVOLUMEFRACTIONKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/PhaseVolumeFractionKernel.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @class PhaseVolumeFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseVolumeFractionKernel : public isothermalCompositionalMultiphaseBaseKernels::PhaseVolumeFractionKernel< NUM_COMP, NUM_PHASE >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PhaseVolumeFractionKernel< NUM_COMP, NUM_PHASE >;
  using Base::m_dPhaseDens;
  using Base::m_dPhaseFrac;
  using Base::m_dPhaseVolFrac;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  PhaseVolumeFractionKernel( ObjectManagerBase & subRegion,
                             constitutive::MultiFluidBase const & fluid )
    : Base( subRegion, fluid )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @param[in] ei the element index
   */
  GEOS_HOST_DEVICE
  real64 compute( localIndex const ei ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseFrac = m_dPhaseFrac[ei][0];

    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];

    // Call the base compute the compute the phase volume fraction
    return Base::compute( ei, [&] ( localIndex const ip,
                                    real64 const & phaseVolFrac,
                                    real64 const & phaseDensInv,
                                    real64 const & totalDensity )
    {
      // when this lambda is called, we are in the phase loop
      // for each phase ip, compute the derivative of phase volume fraction wrt temperature
      dPhaseVolFrac[ip][Deriv::dT] = (dPhaseFrac[ip][Deriv::dT] - phaseVolFrac * dPhaseDens[ip][Deriv::dT]) * phaseDensInv;
      dPhaseVolFrac[ip][Deriv::dT] *= totalDensity;
    } );
  }

};

/**
 * @class PhaseVolumeFractionKernelFactory
 */
class PhaseVolumeFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static real64
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   constitutive::MultiFluidBase const & fluid )
  {
    real64 maxDeltaPhaseVolFrac = 0.0;
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    return maxDeltaPhaseVolFrac;
  }
};

} // namespace thermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALPHASEVOLUMEFRACTIONKERNEL_HPP
