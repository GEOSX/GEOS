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
 * @file ThermalDiffusionDispersionFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIFFUSIONDISPERSIONFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIFFUSIONDISPERSIONFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/DiffusionDispersionFluxComputeKernel.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseFVMKernels
{

/******************************** DiffusionDispersionFluxComputeKernel ********************************/

/**
 * @class DiffusionDispersionFluxComputeKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of diffusion/dispersion flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class DiffusionDispersionFluxComputeKernel :
  public isothermalCompositionalMultiphaseFVMKernels::DiffusionDispersionFluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using AbstractBase::m_dt;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dPhaseVolFrac;

  using Base = typename isothermalCompositionalMultiphaseFVMKernels::DiffusionDispersionFluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
  using DiffusionAccessors = typename Base::DiffusionAccessors;
  using DispersionAccessors = typename Base::DispersionAccessors;
  using PorosityAccessors = typename Base::PorosityAccessors;
  using Base::numFluxSupportPoints;
  using Base::numEqn;
  using Base::numComp;
  using Base::numDof;
  using Base::m_referencePorosity;
  using Base::m_phaseVolFrac;
  using Base::m_phaseDens;
  using Base::m_dPhaseDens;
  using Base::m_phaseDiffusivityMultiplier;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] diffusionAccessors
   * @param[in] dispersionAccessors
   * @param[in] porosityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  DiffusionDispersionFluxComputeKernel( integer const numPhases,
                                        globalIndex const rankOffset,
                                        STENCILWRAPPER const & stencilWrapper,
                                        DofNumberAccessor const & dofNumberAccessor,
                                        CompFlowAccessors const & compFlowAccessors,
                                        MultiFluidAccessors const & multiFluidAccessors,
                                        DiffusionAccessors const & diffusionAccessors,
                                        DispersionAccessors const & dispersionAccessors,
                                        PorosityAccessors const & porosityAccessors,
                                        real64 const dt,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs,
                                        BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags )
    : Base( numPhases,
            rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            diffusionAccessors,
            dispersionAccessors,
            porosityAccessors,
            dt,
            localMatrix,
            localRhs,
            kernelFlags )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    using Base::StackVariables::transmissibility;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

  };

  /**
   * @brief Compute the local diffusion flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void computeDiffusionFlux( localIndex const iconn,
                             StackVariables & stack ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    // ***********************************************
    // First, we call the base computeFlux to compute the diffusionFlux and its derivatives (including derivatives wrt temperature),
    //
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeDiffusionFlux( iconn, stack, [&] ( integer const ip,
                                                    integer const ic,
                                                    localIndex const (&k)[2],
                                                    localIndex const (&seri)[2],
                                                    localIndex const (&sesri)[2],
                                                    localIndex const (&sei)[2],
                                                    localIndex const connectionIndex,
                                                    localIndex const k_up,
                                                    localIndex const er_up,
                                                    localIndex const esr_up,
                                                    localIndex const ei_up,
                                                    real64 const compFracGrad,
                                                    real64 const upwindCoefficient )
    {
      // We are in the loop over phases and components, ip provides the current phase index.

      real64 dCompFracGrad_dT[numFluxSupportPoints]{};
      real64 dDiffusionFlux_dT[numFluxSupportPoints]{};

      /// compute the TPFA component difference
      for( integer i = 0; i < numFluxSupportPoints; i++ )
      {
        localIndex const er  = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei  = sei[i];

        dCompFracGrad_dT[i] += stack.transmissibility[connectionIndex][i] * m_dPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dT];
      }

      // add contributions of the derivatives of component fractions wrt pressure/component fractions
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dDiffusionFlux_dT[ke] += upwindCoefficient * dCompFracGrad_dT[ke];
      }

      // add contributions of the derivatives of upwind coefficient wrt temperature
      real64 const dUpwindCoefficient_dT =
        m_referencePorosity[er_up][esr_up][ei_up] *
        m_phaseDiffusivityMultiplier[er_up][esr_up][ei_up][0][ip] *
        ( m_dPhaseDens[er_up][esr_up][ei_up][0][ip][Deriv::dT] * m_phaseVolFrac[er_up][esr_up][ei_up][ip]
          + m_phaseDens[er_up][esr_up][ei_up][0][ip] * m_dPhaseVolFrac[er_up][esr_up][ei_up][ip][Deriv::dT] );
      dDiffusionFlux_dT[k_up] += dUpwindCoefficient_dT * compFracGrad;

      // finally, increment local flux and local Jacobian
      integer const eqIndex0 = k[0] * numEqn + ic;
      integer const eqIndex1 = k[1] * numEqn + ic;

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numComp + 1;
        stack.localFluxJacobian[eqIndex0][localDofIndexTemp] += m_dt * dDiffusionFlux_dT[ke];
        stack.localFluxJacobian[eqIndex1][localDofIndexTemp] -= m_dt * dDiffusionFlux_dT[ke];
      }
    } );
  }

  /**
   * @brief Compute the local dispersion flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void computeDispersionFlux( localIndex const iconn,
                              StackVariables & stack ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    // ***********************************************
    // First, we call the base computeFlux to compute the dispersionFlux and its derivatives (including derivatives wrt temperature),
    //
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeDispersionFlux( iconn, stack, [&] ( integer const ip,
                                                     integer const ic,
                                                     localIndex const (&k)[2],
                                                     localIndex const (&seri)[2],
                                                     localIndex const (&sesri)[2],
                                                     localIndex const (&sei)[2],
                                                     localIndex const connectionIndex,
                                                     localIndex const k_up,
                                                     localIndex const er_up,
                                                     localIndex const esr_up,
                                                     localIndex const ei_up,
                                                     real64 const compFracGrad )
    {
      // We are in the loop over phases and components, ip provides the current phase index.

      real64 dCompFracGrad_dT[numFluxSupportPoints]{};
      real64 dDispersionFlux_dT[numFluxSupportPoints]{};

      /// compute the TPFA component difference
      for( integer i = 0; i < numFluxSupportPoints; i++ )
      {
        localIndex const er  = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei  = sei[i];

        dCompFracGrad_dT[i] += stack.transmissibility[connectionIndex][i] * m_dPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dT];
      }

      // add contributions of the derivatives of component fractions wrt pressure/component fractions
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dDispersionFlux_dT[ke] += m_phaseDens[er_up][esr_up][ei_up][0][ip] * dCompFracGrad_dT[ke];
      }

      // add contributions of the derivatives of upwind coefficient wrt temperature
      dDispersionFlux_dT[k_up] += m_dPhaseDens[er_up][esr_up][ei_up][0][ip][Deriv::dT] * compFracGrad;

      // finally, increment local flux and local Jacobian
      integer const eqIndex0 = k[0] * numEqn + ic;
      integer const eqIndex1 = k[1] * numEqn + ic;

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numComp + 1;
        stack.localFluxJacobian[eqIndex0][localDofIndexTemp] += m_dt * dDispersionFlux_dT[ke];
        stack.localFluxJacobian[eqIndex1][localDofIndexTemp] -= m_dt * dDispersionFlux_dT[ke];
      }
    } );
  }
};

/**
 * @class DiffusionDispersionFluxComputeKernelFactory
 */
class DiffusionDispersionFluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] hasDiffusion flag specifying whether diffusion is used or not
   * @param[in] hasDispersion flag specifying whether dispersion is used or not
   * @param[in] solverName the name of the solver
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasDiffusion,
                   integer const hasDispersion,
                   integer const useTotalMassEquation,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 2;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::KernelFlags::TotalMassEquation );

      using kernelType = DiffusionDispersionFluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename kernelType::DiffusionAccessors diffusionAccessors( elemManager, solverName );
      typename kernelType::DispersionAccessors dispersionAccessors( elemManager, solverName );
      typename kernelType::PorosityAccessors porosityAccessors( elemManager, solverName );

      kernelType kernel( numPhases, rankOffset, stencilWrapper,
                         dofNumberAccessor, compFlowAccessors, multiFluidAccessors,
                         diffusionAccessors, dispersionAccessors, porosityAccessors,
                         dt, localMatrix, localRhs, kernelFlags );
      kernelType::template launch< POLICY >( stencilWrapper.size(),
                                             hasDiffusion, hasDispersion,
                                             kernel );
    } );
  }
};

} // namespace thermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIFFUSIONDISPERSIONFLUXCOMPUTEKERNEL_HPP
