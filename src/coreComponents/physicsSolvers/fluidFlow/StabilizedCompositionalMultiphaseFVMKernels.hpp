/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StabilizedCompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_STABILIZEDCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_STABILIZEDCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

namespace geos
{

namespace stabilizedCompositionalMultiphaseFVMKernels
{

using namespace constitutive;


/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using StabCompFlowAccessors =
    StencilAccessors< fields::flow::macroElementIndex,
                      fields::flow::elementStabConstant,
                      fields::flow::pressure_n >;

  using StabMultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseDensity_n,
                              fields::multifluid::phaseCompFraction_n >;

  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase, fields::relperm::phaseRelPerm_n >;

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_kernelFlags;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;
  using AbstractBase::m_pres;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::numFluxSupportPoints;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_phaseMob;
  using Base::m_dPhaseMassDens;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] compFlowAccessor accessor for wrappers registered by the solver
   * @param[in] stabCompFlowAccessor accessor for wrappers registered by the solver needed for stabilization
   * @param[in] multiFluidAccessor accessor for wrappers registered by the multifluid model
   * @param[in] stabMultiFluidAccessor accessor for wrappers registered by the multifluid model needed for stabilization
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] relPermAccessors accessor for wrappers registered by the relative permeability model needed for stabilization
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  FaceBasedAssemblyKernel( integer const numPhases,
                           globalIndex const rankOffset,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           StabCompFlowAccessors const & stabCompFlowAccessors,
                           MultiFluidAccessors const & multiFluidAccessors,
                           StabMultiFluidAccessors const & stabMultiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           RelPermAccessors const & relPermAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           BitFlags< isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags > kernelFlags )
    : Base( numPhases,
            rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs,
            kernelFlags ),
    m_pres_n( stabCompFlowAccessors.get( fields::flow::pressure_n {} ) ),
    m_phaseDens_n( stabMultiFluidAccessors.get( fields::multifluid::phaseDensity_n {} ) ),
    m_phaseCompFrac_n( stabMultiFluidAccessors.get( fields::multifluid::phaseCompFraction_n {} ) ),
    m_phaseRelPerm_n( relPermAccessors.get( fields::relperm::phaseRelPerm_n {} ) ),
    m_macroElementIndex( stabCompFlowAccessors.get( fields::flow::macroElementIndex {} ) ),
    m_elementStabConstant( stabCompFlowAccessors.get( fields::flow::elementStabConstant {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numConnectedElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    /// Stabilization transmissibility
    real64 stabTransmissibility[maxNumConns][numFluxSupportPoints]{};

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian, including stabilization
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {

    m_stencilWrapper.computeStabilizationWeights( iconn,
                                                  stack.stabTransmissibility );

    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives,
    //
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to compute stabilization terms
    Base::computeFlux( iconn, stack, [&] ( integer const ip,
                                           localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const connectionIndex,
                                           localIndex const k_up,
                                           localIndex const er_up,
                                           localIndex const esr_up,
                                           localIndex const ei_up,
                                           real64 const potGrad,
                                           real64 const phaseFlux,
                                           real64 const (&dPhaseFlux_dP)[2],
                                           real64 const (&dPhaseFlux_dC)[2][numComp] )
    {
      GEOS_UNUSED_VAR( k_up, potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, er_up, esr_up, ei_up );

      /// stabilization flux and derivatives
      real64 stabFlux[numComp]{};
      real64 dStabFlux_dP[numFluxSupportPoints][numComp]{};

      real64 const stabTrans[numFluxSupportPoints] = { stack.stabTransmissibility[connectionIndex][0],
                                                       stack.stabTransmissibility[connectionIndex][1] };

      // we are in the loop over phases, ip provides the current phase index.

      real64 dPresGradStab = 0.0;
      integer stencilMacroElements[numFluxSupportPoints]{};

      // Step 1: compute the pressure jump at the interface
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const er  = seri[ke];
        localIndex const esr = sesri[ke];
        localIndex const ei  = sei[ke];

        stencilMacroElements[ke] = m_macroElementIndex[er][esr][ei];

        // use the jump in *delta* pressure in the stabilization term
        dPresGradStab +=
          m_elementStabConstant[er][esr][ei] * stabTrans[ke] * (m_pres[er][esr][ei] - m_pres_n[er][esr][ei]);
      }

      // Step 2: compute the stabilization flux
      integer const k_up_stab = (dPresGradStab >= 0) ? 0 : 1;

      localIndex const er_up_stab  = seri[k_up_stab];
      localIndex const esr_up_stab = sesri[k_up_stab];
      localIndex const ei_up_stab  = sei[k_up_stab];

      bool const areInSameMacroElement = stencilMacroElements[0] == stencilMacroElements[1];
      bool const isStabilizationActive = stencilMacroElements[0] >= 0 && stencilMacroElements[1] >= 0;
      if( isStabilizationActive && areInSameMacroElement )
      {

        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const laggedUpwindCoef = m_phaseDens_n[er_up_stab][esr_up_stab][ei_up_stab][0][ip]
                                          * m_phaseCompFrac_n[er_up_stab][esr_up_stab][ei_up_stab][0][ip][ic]
                                          * m_phaseRelPerm_n[er_up_stab][esr_up_stab][ei_up_stab][0][ip];
          stabFlux[ic] += dPresGradStab * laggedUpwindCoef;

          for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
          {
            real64 const tauStab = m_elementStabConstant[seri[ke]][sesri[ke]][sei[ke]];
            dStabFlux_dP[ke][ic] += tauStab * stabTrans[ke] * laggedUpwindCoef;
          }
        }
      }

      // Step 3: add the stabilization flux and its derivatives to the residual and Jacobian
      for( integer ic = 0; ic < numComp; ++ic )
      {
        integer const eqIndex0 = k[0] * numEqn + ic;
        integer const eqIndex1 = k[1] * numEqn + ic;

        stack.localFlux[eqIndex0] +=  stabFlux[ic];
        stack.localFlux[eqIndex1] += -stabFlux[ic];

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          localIndex const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[eqIndex0][localDofIndexPres] +=  dStabFlux_dP[ke][ic];
          stack.localFluxJacobian[eqIndex1][localDofIndexPres] += -dStabFlux_dP[ke][ic];
        }
      }

    } ); // end call to Base::computeFlux

  }

protected:

  /// Views on flow properties at the previous converged time step
  ElementViewConst< arrayView1d< real64 const > > const m_pres_n;
  ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const m_phaseDens_n;
  ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const m_phaseCompFrac_n;
  ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const m_phaseRelPerm_n;

  /// Views on the macroelement indices and stab constant
  ElementViewConst< arrayView1d< integer const > > const m_macroElementIndex;
  ElementViewConst< arrayView1d< real64 const > > const m_elementStabConstant;

};

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class FaceBasedAssemblyKernelFactory
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
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] solverName name of the solver (to name accessors)
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
                   integer const hasCapPressure,
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
      integer constexpr NUM_DOF = NC() + 1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      BitFlags< isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags > kernelFlags;
      if( hasCapPressure )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags::CapPressure );
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags::TotalMassEquation );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::StabCompFlowAccessors stabCompFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::StabMultiFluidAccessors stabMultiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::RelPermAccessors relPermAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor,
                          compFlowAccessors, stabCompFlowAccessors, multiFluidAccessors, stabMultiFluidAccessors,
                          capPressureAccessors, permeabilityAccessors, relPermAccessors,
                          dt, localMatrix, localRhs, kernelFlags );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace stabilizedCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_STABILIZEDCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
