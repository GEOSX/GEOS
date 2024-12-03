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
 * @file ThermalPhaseMobilityKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALPHASEMOBILITYKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALPHASEMOBILITYKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/PhaseMobilityKernel.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseFVMKernels
{

/******************************** PhaseMobilityKernel ********************************/

/**
 * @class PhaseMobilityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseMobilityKernel : public isothermalCompositionalMultiphaseFVMKernels::PhaseMobilityKernel< NUM_COMP, NUM_PHASE >
{
public:

  using Base = isothermalCompositionalMultiphaseFVMKernels::PhaseMobilityKernel< NUM_COMP, NUM_PHASE >;
  using Base::numPhase;
  using Base::m_dPhaseVolFrac;
  using Base::m_dPhaseMob;
  using Base::m_phaseDens;
  using Base::m_dPhaseDens;
  using Base::m_phaseVisc;
  using Base::m_dPhaseVisc;
  using Base::m_dPhaseRelPerm_dPhaseVolFrac;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityKernel( ObjectManagerBase & subRegion,
                       constitutive::MultiFluidBase const & fluid,
                       constitutive::RelativePermeabilityBase const & relperm )
    : Base( subRegion, fluid, relperm )
  {}

  /**
   * @brief Compute the phase mobilities in an element
   * @param[in] ei the element index
   */
  GEOS_HOST_DEVICE
  inline
  void compute( localIndex const ei ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseVisc = m_phaseVisc[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseVisc = m_dPhaseVisc[ei][0];
    arraySlice2d< real64 const, constitutive::relperm::USD_RELPERM_DS - 2 > const dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[ei][0];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];

    Base::compute( ei, [&] ( localIndex const ip,
                             real64 const & phaseMob,
                             arraySlice1d< real64, compflow::USD_PHASE_DC - 2 > const & dPhaseMob )
    {
      // Step 1: compute the derivative of relPerm[ip] wrt temperature
      real64 dRelPerm_dT = 0.0;
      for( integer jp = 0; jp < numPhase; ++jp )
      {
        dRelPerm_dT += dPhaseRelPerm_dPhaseVolFrac[ip][jp] * dPhaseVolFrac[jp][Deriv::dT];
      }

      // Step 2: compute the derivative of phaseMob[ip] wrt temperature
      dPhaseMob[Deriv::dT] = dRelPerm_dT * phaseDens[ip] / phaseVisc[ip]
                             + phaseMob * (dPhaseDens[ip][Deriv::dT] / phaseDens[ip] - dPhaseVisc[ip][Deriv::dT] / phaseVisc[ip] );
    } );
  }

};

/**
 * @class PhaseMobilityKernelFactory
 */
class PhaseMobilityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::RelativePermeabilityBase const & relperm )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 2 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityKernel< NUM_COMP, 3 > kernel( subRegion, fluid, relperm );
        PhaseMobilityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};

} // namespace thermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALPHASEMOBILITYKERNEL_HPP
