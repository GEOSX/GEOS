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
 * @file ThermalDirichletFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/DirichletFluxComputeKernel.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/thermalConductivity/MultiPhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseFVMKernels
{

/******************************** DirichletFluxComputeKernel ********************************/

/**
 * @class DirichletFluxComputeKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFluxComputeKernel : public isothermalCompositionalMultiphaseFVMKernels::DirichletFluxComputeKernel< NUM_COMP,
                                                                                                                   NUM_DOF,
                                                                                                                   FLUIDWRAPPER >
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
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;

  using Base = isothermalCompositionalMultiphaseFVMKernels::DirichletFluxComputeKernel< NUM_COMP, NUM_DOF, FLUIDWRAPPER >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_phaseMob;
  using Base::m_dPhaseMob;
  using Base::m_dPhaseMassDens;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::m_faceTemp;
  using Base::m_faceGravCoef;


  using ThermalCompFlowAccessors =
    StencilAccessors< fields::flow::temperature >;

  using ThermalMultiFluidAccessors =
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseEnthalpy,
                              fields::multifluid::dPhaseEnthalpy >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< constitutive::MultiPhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity >;
  // for now, we treat thermal conductivity explicitly

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] compFlowAccessor accessor for wrappers registered by the solver
   * @param[in] thermalCompFlowAccessors accessor for *thermal* wrappers registered by the solver
   * @param[in] multiFluidAccessor accessor for wrappers registered by the multifluid model
   * @param[in] thermalMultiFluidAccessors accessor for *thermal* wrappers registered by the multifluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] thermalConductivityAccessors accessor for wrappers registered by the thermal conductivity model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  DirichletFluxComputeKernel( integer const numPhases,
                              globalIndex const rankOffset,
                              FaceManager const & faceManager,
                              BoundaryStencilWrapper const & stencilWrapper,
                              FLUIDWRAPPER const & fluidWrapper,
                              DofNumberAccessor const & dofNumberAccessor,
                              CompFlowAccessors const & compFlowAccessors,
                              ThermalCompFlowAccessors const & thermalCompFlowAccessors,
                              MultiFluidAccessors const & multiFluidAccessors,
                              ThermalMultiFluidAccessors const & thermalMultiFluidAccessors,
                              CapPressureAccessors const & capPressureAccessors,
                              PermeabilityAccessors const & permeabilityAccessors,
                              ThermalConductivityAccessors const & thermalConductivityAccessors,
                              real64 const dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags )
    : Base( numPhases,
            rankOffset,
            faceManager,
            stencilWrapper,
            fluidWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs,
            kernelFlags ),
    m_temp( thermalCompFlowAccessors.get( fields::flow::temperature {} ) ),
    m_phaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::phaseEnthalpy {} ) ),
    m_dPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::dPhaseEnthalpy {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    // Component fluxes and derivatives

    /// Derivatives of component fluxes wrt temperature
    real64 dCompFlux_dT[numComp]{};


    // Energy fluxes and derivatives

    /// Energy fluxes
    real64 energyFlux = 0.0;
    /// Derivative of energy fluxes wrt pressure
    real64 dEnergyFlux_dP = 0.0;
    /// Derivative of energy fluxes wrt temperature
    real64 dEnergyFlux_dT = 0.0;
    /// Derivatives of energy fluxes wrt component densities
    real64 dEnergyFlux_dC[numComp]{};

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    using Order = BoundaryStencil::Order;
    using Deriv = constitutive::multifluid::DerivativeOffset;

    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives (including derivatives wrt temperature),
    //  2) enthalpy part of energyFlux  and its derivatives (including derivatives wrt temperature)
    //
    // Computing dCompFlux_dT and the enthalpy flux requires quantities already computed in the base computeFlux,
    // such as potGrad, phaseFlux, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeFlux( iconn, stack, [&] ( integer const ip,
                                           localIndex const er,
                                           localIndex const esr,
                                           localIndex const ei,
                                           localIndex const kf,
                                           real64 const f, // potGrad times trans
                                           real64 const facePhaseMob,
                                           arraySlice1d< const real64, constitutive::multifluid::USD_PHASE - 2 > const & facePhaseEnthalpy,
                                           arraySlice2d< const real64, constitutive::multifluid::USD_PHASE_COMP-2 > const & facePhaseCompFrac,
                                           real64 const phaseFlux,
                                           real64 const dPhaseFlux_dP,
                                           real64 const (&dPhaseFlux_dC)[numComp] )
    {
      // We are in the loop over phases, ip provides the current phase index.

      // Step 1: compute the derivatives of the mean density at the interface wrt temperature

      real64 const dDensMean_dT = 0.5 * m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dT];

      // Step 2: compute the derivatives of the phase potential difference wrt temperature
      //***** calculation of flux *****

      real64 const dF_dT = -stack.transmissibility * dDensMean_dT * ( m_gravCoef[er][esr][ei] - m_faceGravCoef[kf] );

      // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
      // *** upwinding ***

      // note: the upwinding is done in the base class, which is in charge of
      //       computing the following quantities: potGrad, phaseFlux
      // It is easier to hard-code the if/else because it is difficult to address elem and face variables in a uniform way


      if( f >= 0 ) // the element is upstream
      {

        // Step 3.1.a: compute the derivative of phase flux wrt temperature
        real64 const dPhaseFlux_dT = m_phaseMob[er][esr][ei][ip] * dF_dT + m_dPhaseMob[er][esr][ei][ip][Deriv::dT] * f;

        // Step 3.2.a: compute the derivative of component flux wrt temperature

        // slice some constitutive arrays to avoid too much indexing in component loop
        arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 3 > phaseCompFracSub =
          m_phaseCompFrac[er][esr][ei][0][ip];
        arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 3 > dPhaseCompFracSub =
          m_dPhaseCompFrac[er][esr][ei][0][ip];

        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = phaseCompFracSub[ic];
          stack.dCompFlux_dT[ic] += dPhaseFlux_dT * ycp + phaseFlux * dPhaseCompFracSub[ic][Deriv::dT];
        }

        // Step 3.3.a: compute the enthalpy flux

        real64 const enthalpy = m_phaseEnthalpy[er][esr][ei][0][ip];
        stack.energyFlux += phaseFlux * enthalpy;
        stack.dEnergyFlux_dP += dPhaseFlux_dP * enthalpy + phaseFlux * m_dPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dP];
        stack.dEnergyFlux_dT += dPhaseFlux_dT * enthalpy + phaseFlux * m_dPhaseEnthalpy[er][esr][ei][0][ip][Deriv::dT];

        real64 dProp_dC[numComp]{};
        applyChainRule( numComp,
                        m_dCompFrac_dCompDens[er][esr][ei],
                        m_dPhaseEnthalpy[er][esr][ei][0][ip],
                        dProp_dC,
                        Deriv::dC );
        for( integer jc = 0; jc < numComp; ++jc )
        {
          stack.dEnergyFlux_dC[jc] += dPhaseFlux_dC[jc] * enthalpy + phaseFlux * dProp_dC[jc];
        }

      }
      else // the face is upstream
      {

        // Step 3.1.b: compute the derivative of phase flux wrt temperature
        real64 const dPhaseFlux_dT = facePhaseMob * dF_dT;

        // Step 3.2.b: compute the derivative of component flux wrt temperature

        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = facePhaseCompFrac[ip][ic];
          stack.dCompFlux_dT[ic] += dPhaseFlux_dT * ycp;
        }

        // Step 3.3.b: compute the enthalpy flux

        real64 const enthalpy = facePhaseEnthalpy[ip];
        stack.energyFlux += phaseFlux * enthalpy;
        stack.dEnergyFlux_dP += dPhaseFlux_dP * enthalpy;
        stack.dEnergyFlux_dT += dPhaseFlux_dT * enthalpy;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          stack.dEnergyFlux_dC[jc] += dPhaseFlux_dC[jc] * enthalpy;
        }

      }

    } );

    // *****************************************************
    // Computation of the conduction term in the energy flux
    // Note that the phase enthalpy term in the energy was computed above
    // Note that this term is computed using an explicit treatment of conductivity for now

    // Step 1: compute the thermal transmissibilities at this face
    // Below, the thermal conductivity used to compute (explicitly) the thermal conducivity
    // To avoid modifying the signature of the "computeWeights" function for now, we pass m_thermalConductivity twice
    // TODO: modify computeWeights to accomodate explicit coefficients
    real64 thermalTrans = 0.0;
    real64 dThermalTrans_dPerm[3]{}; // not used
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     thermalTrans,
                                     dThermalTrans_dPerm );

    // Step 2: compute temperature difference at the interface
    stack.energyFlux += thermalTrans
                        * ( m_temp[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )] - m_faceTemp[m_sei( iconn, Order::FACE )] );
    stack.dEnergyFlux_dT += thermalTrans;


    // **********************************************************************************
    // At this point, we have computed the energyFlux and the compFlux for all components
    // We have to do two things here:
    // 1) Add dCompFlux_dTemp to the localFluxJacobian of the component mass balance equations
    // 2) Add energyFlux and its derivatives to the localFlux(Jacobian) of the energy balance equation

    // Step 1: add dCompFlux_dTemp to localFluxJacobian
    for( integer ic = 0; ic < numComp; ++ic )
    {
      stack.localFluxJacobian[ic][numDof-1] =  m_dt * stack.dCompFlux_dT[ic];
    }

    // Step 2: add energyFlux and its derivatives to localFlux and localFluxJacobian
    integer const localRowIndexEnergy = numEqn-1;
    stack.localFlux[localRowIndexEnergy] =  m_dt * stack.energyFlux;

    stack.localFluxJacobian[localRowIndexEnergy][0] =  m_dt * stack.dEnergyFlux_dP;
    stack.localFluxJacobian[localRowIndexEnergy][numDof-1] =  m_dt * stack.dEnergyFlux_dT;
    for( integer jc = 0; jc < numComp; ++jc )
    {
      stack.localFluxJacobian[localRowIndexEnergy][jc+1] =  m_dt * stack.dEnergyFlux_dC[jc];
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Call Case::complete to assemble the component mass balance equations (i = 0 to i = numDof-2)
    // In the lambda, add contribution to residual and jacobian into the energy balance equation
    Base::complete( iconn, stack, [&] ( localIndex const localRow )
    {
      // beware, there is  volume balance eqn in m_localRhs and m_localMatrix!
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn], stack.localFlux[numEqn-1] );
      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
        ( localRow + numEqn,
        stack.dofColIndices,
        stack.localFluxJacobian[numEqn-1],
        numDof );

    } );
  }

protected:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on phase enthalpies
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_phaseEnthalpy;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dPhaseEnthalpy;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > const m_thermalConductivity;
  // for now, we treat thermal conductivity explicitly

};

/**
 * @class DirichletFluxComputeKernelFactory
 */
class DirichletFluxComputeKernelFactory
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
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] faceManager reference to the face manager
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidBase the multifluid constitutive model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const & dofKey,
                   string const & solverName,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   constitutive::MultiFluidBase & fluidBase,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    constitutive::constitutiveComponentUpdatePassThru< true >( fluidBase, numComps, [&]( auto & fluid, auto NC )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper const fluidWrapper = fluid.createKernelWrapper();

      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 2;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      // for now, we neglect capillary pressure in the kernel
      BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::KernelFlags::TotalMassEquation );

      using KernelType = DirichletFluxComputeKernel< NUM_COMP, NUM_DOF, typename FluidType::KernelWrapper >;
      typename KernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KernelType::ThermalCompFlowAccessors thermalCompFlowAccessors( elemManager, solverName );
      typename KernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KernelType::ThermalMultiFluidAccessors thermalMultiFluidAccessors( elemManager, solverName );
      typename KernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

      KernelType kernel( numPhases, rankOffset, faceManager, stencilWrapper, fluidWrapper,
                         dofNumberAccessor, compFlowAccessors, thermalCompFlowAccessors, multiFluidAccessors, thermalMultiFluidAccessors,
                         capPressureAccessors, permeabilityAccessors, thermalConductivityAccessors,
                         dt, localMatrix, localRhs, kernelFlags );
      KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace thermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP
