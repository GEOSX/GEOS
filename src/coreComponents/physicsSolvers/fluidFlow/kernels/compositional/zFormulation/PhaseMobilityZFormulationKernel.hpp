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
 * @file PhaseMobilityZFormulationKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASEMOBILITYZFORMULATIONKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASEMOBILITYZFORMULATIONKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/PropertyKernelBase.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** PhaseMobilityZFormulationKernel ********************************/

/**
 * @class PhaseMobilityZFormulationKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Defines the interface for the property kernel in charge of computing the phase mobilities
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseMobilityZFormulationKernel : public isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] relperm the relperm model
   */
  PhaseMobilityZFormulationKernel( ObjectManagerBase & subRegion,
                                   constitutive::MultiFluidBase const & fluid,
                                   constitutive::RelativePermeabilityBase const & relperm )
    : Base(),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseVisc( fluid.phaseViscosity() ),
    m_dPhaseVisc( fluid.dPhaseViscosity() ),
    m_phaseRelPerm( relperm.phaseRelPerm() ),
    m_dPhaseRelPerm_dPhaseVolFrac( relperm.dPhaseRelPerm_dPhaseVolFraction() ),
    m_phaseMob( subRegion.getField< fields::flow::phaseMobility >() ),
    m_dPhaseMob( subRegion.getField< fields::flow::dPhaseMobility >() )
  {}

  /**
   * @brief Compute the phase mobilities in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] phaseMobilityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                FUNC && phaseMobilityKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseVisc = m_phaseVisc[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseVisc = m_dPhaseVisc[ei][0];
    arraySlice1d< real64 const, constitutive::relperm::USD_RELPERM - 2 > const phaseRelPerm = m_phaseRelPerm[ei][0];
    arraySlice2d< real64 const, constitutive::relperm::USD_RELPERM_DS - 2 > const dPhaseRelPerm_dPhaseVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[ei][0];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice1d< real64, compflow::USD_PHASE - 1 > const phaseMob = m_phaseMob[ei];
    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseMob = m_dPhaseMob[ei];

    real64 dRelPerm_dC[numComp]{};
    real64 dDens_dC[numComp]{};
    real64 dVisc_dC[numComp]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {

      // compute the phase mobility only if the phase is present
      bool const phaseExists = (phaseVolFrac[ip] > 0);
      if( !phaseExists )
      {
        phaseMob[ip] = 0.0;
        for( integer jc = 0; jc < numComp + 2; ++jc )
        {
          dPhaseMob[ip][jc] = 0.0;
        }
        continue;
      }

      real64 const density = phaseDens[ip];
      real64 const dDens_dP = dPhaseDens[ip][Deriv::dP];
      for( integer jc = 0; jc < numComp; ++jc )
        dDens_dC[jc] = dPhaseDens[ip][Deriv::dC+jc];

      real64 const viscosity = phaseVisc[ip];
      real64 const dVisc_dP = dPhaseVisc[ip][Deriv::dP];
      for( integer jc = 0; jc < numComp; ++jc )
        dDens_dC[jc] = dPhaseVisc[ip][Deriv::dC+jc];

      real64 const relPerm = phaseRelPerm[ip];
      real64 dRelPerm_dP = 0.0;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dRelPerm_dC[ic] = 0.0;
      }

      for( integer jp = 0; jp < numPhase; ++jp )
      {
        real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dP];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac[jp][Deriv::dC+jc];
        }
      }

      real64 const mobility = relPerm * density / viscosity;

      phaseMob[ip] = mobility;
      dPhaseMob[ip][Deriv::dP] = dRelPerm_dP * density / viscosity
                                 + mobility * (dDens_dP / density - dVisc_dP / viscosity);

      // compositional derivatives
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseMob[ip][Deriv::dC+jc] = dRelPerm_dC[jc] * density / viscosity
                                      + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
      }

      // call the lambda in the phase loop to allow the reuse of the relperm, density, viscosity, and mobility
      // possible use: assemble the derivatives wrt temperature
      phaseMobilityKernelOp( ip, phaseMob[ip], dPhaseMob[ip] );
    }
  }

protected:

  // inputs

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac;

  /// Views on the phase densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseDens;

  /// Views on the phase viscosities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseVisc;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseVisc;

  /// Views on the phase relative permeabilities
  arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > m_phaseRelPerm;
  arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > m_dPhaseRelPerm_dPhaseVolFrac;

  // outputs

  /// Views on the phase mobilities
  arrayView2d< real64, compflow::USD_PHASE > m_phaseMob;
  arrayView3d< real64, compflow::USD_PHASE_DC > m_dPhaseMob;

};

/**
 * @class PhaseMobilityZFormulationKernelFactory
 */
class PhaseMobilityZFormulationKernelFactory
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
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityZFormulationKernel< NUM_COMP, 2 > kernel( subRegion, fluid, relperm );
        PhaseMobilityZFormulationKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseMobilityZFormulationKernel< NUM_COMP, 3 > kernel( subRegion, fluid, relperm );
        PhaseMobilityZFormulationKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_PHASEMOBILITYZFORMULATIONKERNEL_HPP
