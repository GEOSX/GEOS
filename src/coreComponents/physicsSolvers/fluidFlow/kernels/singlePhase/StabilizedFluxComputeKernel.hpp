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
 * @file StabilizedFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_STABILIZEDFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_STABILIZEDFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernel.hpp"

namespace geos
{

namespace stabilizedSinglePhaseFVMKernels
{

/******************************** FluxComputeKernel ********************************/

/**
 * @class FluxComputeKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeKernel : public singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >
{
public:

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using AbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using StabSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::macroElementIndex,
                      fields::flow::elementStabConstant,
                      fields::flow::pressure_n >;

  using StabSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density_n >;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;

  using Base = singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] singlePhaseFlowAccessor accessor for wrappers registered by the solver
   * @param[in] stabSinglePhaseFlowAccessor accessor for wrappers registered by the solver needed for stabilization
   * @param[in] singlePhaseFluidAccessor accessor for wrappers registered by the single fluid model
   * @param[in] stabSinglePhaseFluidAccessor accessor for wrappers registered by the single fluid model needed for stabilization
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FluxComputeKernel( globalIndex const rankOffset,
                     STENCILWRAPPER const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor,
                     SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                     StabSinglePhaseFlowAccessors const & stabSinglePhaseFlowAccessors,
                     SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                     StabSinglePhaseFluidAccessors const & stabSinglePhaseFluidAccessors,
                     PermeabilityAccessors const & permeabilityAccessors,
                     real64 const & dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_pres_n( stabSinglePhaseFlowAccessors.get( fields::flow::pressure_n {} ) ),
    m_dens_n( stabSinglePhaseFluidAccessors.get( fields::singlefluid::density_n {} ) ),
    m_macroElementIndex( stabSinglePhaseFlowAccessors.get( fields::flow::macroElementIndex {} ) ),
    m_elementStabConstant( stabSinglePhaseFlowAccessors.get( fields::flow::elementStabConstant {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    /// Stabilization transmissibility
    real64 stabTransmissibility[maxNumConns][2]{};

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
    // We call the base computeFlux,
    //
    // We use the lambda below to compute stabilization terms
    Base::computeFlux( iconn, stack, [&] ( localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const connectionIndex,
                                           real64 const alpha,
                                           real64 const mobility,
                                           real64 const potGrad,
                                           real64 const fluxVal,
                                           real64 const (&dFlux_dP)[2] )
    {

      GEOS_UNUSED_VAR( alpha, mobility, potGrad, fluxVal, dFlux_dP );

      /// stabilization flux and derivatives
      real64 stabFlux{};
      real64 dStabFlux_dP[2]{};

      real64 const stabTrans[2] = { stack.stabTransmissibility[connectionIndex][0],
                                    stack.stabTransmissibility[connectionIndex][1] };


      real64 dPresGradStab = 0.0;
      integer stencilMacroElements[2]{};

      // Step 1: compute the pressure jump at the interface
      for( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        localIndex const er  = seri[ke];
        localIndex const esr = sesri[ke];
        localIndex const ei  = sei[ke];

        stencilMacroElements[ke] = m_macroElementIndex[er][esr][ei];

        // use the jump in *delta* pressure in the stabilization term
        dPresGradStab += m_elementStabConstant[er][esr][ei] * stabTrans[ke] * (m_pres[er][esr][ei] - m_pres_n[er][esr][ei]);
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

        real64 const laggedUpwindCoef = m_dens_n[er_up_stab][esr_up_stab][ei_up_stab][0];
        stabFlux += dPresGradStab * laggedUpwindCoef;

        for( integer ke = 0; ke < stack.numFluxElems; ++ke )
        {
          real64 const tauStab = m_elementStabConstant[seri[ke]][sesri[ke]][sei[ke]];
          dStabFlux_dP[ke] += tauStab * stabTrans[ke] * laggedUpwindCoef;
        }
      }

      // Step 3: add the stabilization flux and its derivatives to the residual and Jacobian
      integer const eqIndex0 = k[0] * numEqn;
      integer const eqIndex1 = k[1] * numEqn;

      stack.localFlux[eqIndex0] +=  stabFlux;
      stack.localFlux[eqIndex1] += -stabFlux;

      for( integer ke = 0; ke < stack.numFluxElems; ++ke )
      {
        localIndex const localDofIndexPres = k[ke] * numDof;
        stack.localFluxJacobian[eqIndex0][localDofIndexPres] +=  dStabFlux_dP[ke];
        stack.localFluxJacobian[eqIndex1][localDofIndexPres] += -dStabFlux_dP[ke];
      }

    } ); // end call to Base::computeFlux

  }

protected:

  /// Views on flow properties at the previous converged time step
  ElementViewConst< arrayView1d< real64 const > > const m_pres_n;
  ElementViewConst< arrayView2d< real64 const > > const m_dens_n;

  /// Views on the macroelement indices and stab constant
  ElementViewConst< arrayView1d< integer const > > const m_macroElementIndex;
  ElementViewConst< arrayView1d< real64 const > > const m_elementStabConstant;

};

/**
 * @class FluxComputeKernelFactory
 */
class FluxComputeKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {

    integer constexpr NUM_EQN = 1;
    integer constexpr NUM_DOF = 1;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using KERNEL_TYPE = FluxComputeKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename KERNEL_TYPE::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, solverName );
    typename KERNEL_TYPE::SinglePhaseFluidAccessors singlePhaseFluidAccessors( elemManager, solverName );
    typename KERNEL_TYPE::StabSinglePhaseFlowAccessors stabSinglePhaseFlowAccessors( elemManager, solverName );
    typename KERNEL_TYPE::StabSinglePhaseFluidAccessors stabSinglePhaseFluidAccessors( elemManager, solverName );
    typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

    KERNEL_TYPE kernel( rankOffset, stencilWrapper, dofNumberAccessor,
                        singlePhaseFlowAccessors, stabSinglePhaseFlowAccessors, singlePhaseFluidAccessors, stabSinglePhaseFluidAccessors,
                        permeabilityAccessors,
                        dt, localMatrix, localRhs );
    KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace stabilizedSinglePhaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_STABILIZEDFLUXCOMPUTEKERNEL_HPP
