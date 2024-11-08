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
 * @file ThermalFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/FluxComputeKernel.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/thermalConductivity/MultiPhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseFVMKernels
{

/******************************** FluxComputeKernel ********************************/

/**
 * @class FluxComputeKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FluxComputeKernel : public isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
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
  using AbstractBase::m_phaseVolFrac;
  using AbstractBase::m_dPhaseVolFrac;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::numFluxSupportPoints;
  using Base::m_phaseMob;
  using Base::m_dPhaseMob;
  using Base::m_dPhaseMassDens;
  using Base::m_dPhaseCapPressure_dPhaseVolFrac;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

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
   * @param[in] stencilWrapper reference to the stencil wrapper
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
   * @param[in] kernelFlags flags packed all together
   */
  FluxComputeKernel( integer const numPhases,
                     globalIndex const rankOffset,
                     STENCILWRAPPER const & stencilWrapper,
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
                     BitFlags< isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelFlags > kernelFlags )
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
    m_temp( thermalCompFlowAccessors.get( fields::flow::temperature {} ) ),
    m_phaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::phaseEnthalpy {} ) ),
    m_dPhaseEnthalpy( thermalMultiFluidAccessors.get( fields::multifluid::dPhaseEnthalpy {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) )
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

    // Thermal transmissibility (for now, no derivatives)

    real64 thermalTransmissibility[maxNumConns][2]{};
  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives (including derivatives wrt temperature),
    //  2) enthalpy part of convectiveEnergyFlux  and its derivatives (including derivatives wrt temperature)
    //
    // Computing dCompFlux_dT and the enthalpy flux requires quantities already computed in the base computeFlux,
    // such as potGrad, phaseFlux, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
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
      // We are in the loop over phases, ip provides the current phase index.

      // Step 1: compute the derivatives of the mean density at the interface wrt temperature

      real64 dDensMean_dT[numFluxSupportPoints]{};

      real64 const trans[numFluxSupportPoints] = { stack.transmissibility[connectionIndex][0],
                                                   stack.transmissibility[connectionIndex][1] };

      real64 convectiveEnergyFlux = 0.0;
      real64 dConvectiveEnergyFlux_dP[numFluxSupportPoints]{};
      real64 dConvectiveEnergyFlux_dT[numFluxSupportPoints]{};
      real64 dConvectiveEnergyFlux_dC[numFluxSupportPoints][numComp]{};
      real64 dCompFlux_dT[numFluxSupportPoints][numComp]{};

      integer denom = 0;
      for( integer i = 0; i < numFluxSupportPoints; ++i )
      {
        localIndex const er  = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei  = sei[i];

        bool const phaseExists = (m_phaseVolFrac[er_up][esr_up][ei_up][ip] > 0);
        if( !phaseExists )
        {
          continue;
        }

        dDensMean_dT[i] = m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dT];
        denom++;
      }
      if( denom > 1 )
      {
        for( integer i = 0; i < numFluxSupportPoints; ++i )
        {
          dDensMean_dT[i] /= denom;
        }
      }

      // Step 2: compute the derivatives of the phase potential difference wrt temperature
      //***** calculation of flux *****

      real64 dPresGrad_dT[numFluxSupportPoints]{};
      real64 dGravHead_dT[numFluxSupportPoints]{};

      // compute potential difference MPFA-style
      for( integer i = 0; i < numFluxSupportPoints; ++i )
      {
        localIndex const er  = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei  = sei[i];

        // Step 2.1: compute derivative of capillary pressure wrt temperature
        real64 dCapPressure_dT = 0.0;
        if( AbstractBase::m_kernelFlags.isSet( isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelFlags::CapPressure ) )
        {
          for( integer jp = 0; jp < m_numPhases; ++jp )
          {
            real64 const dCapPressure_dS = m_dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
            dCapPressure_dT += dCapPressure_dS * m_dPhaseVolFrac[er][esr][ei][jp][Deriv::dT];
          }
        }

        // Step 2.2: compute derivative of phase pressure difference wrt temperature
        dPresGrad_dT[i] -= trans[i] * dCapPressure_dT;
        real64 const gravD = trans[i] * m_gravCoef[er][esr][ei];

        // Step 2.3: compute derivative of gravity potential difference wrt temperature
        for( integer j = 0; j < numFluxSupportPoints; ++j )
        {
          dGravHead_dT[j] += dDensMean_dT[j] * gravD;
        }
      }

      // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
      // *** upwinding ***

      // note: the upwinding is done in the base class, which is in charge of
      //       computing the following quantities: potGrad, phaseFlux, k_up, er_up, esr_up, ei_up

      real64 dPhaseFlux_dT[numFluxSupportPoints]{};

      // Step 3.1: compute the derivative of phase flux wrt temperature
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dT[ke] += dPresGrad_dT[ke];
      }
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dT[ke] -= dGravHead_dT[ke];
      }
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dPhaseFlux_dT[ke] *= m_phaseMob[er_up][esr_up][ei_up][ip];
      }
      dPhaseFlux_dT[k_up] += m_dPhaseMob[er_up][esr_up][ei_up][ip][Deriv::dT] * potGrad;

      // Step 3.2: compute the derivative of component flux wrt temperature

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 3 > phaseCompFracSub =
        m_phaseCompFrac[er_up][esr_up][ei_up][0][ip];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 3 > dPhaseCompFracSub =
        m_dPhaseCompFrac[er_up][esr_up][ei_up][0][ip];

      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const ycp = phaseCompFracSub[ic];
        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          dCompFlux_dT[ke][ic] += dPhaseFlux_dT[ke] * ycp;
        }
        dCompFlux_dT[k_up][ic] += phaseFlux * dPhaseCompFracSub[ic][Deriv::dT];
      }

      // Step 4: add dCompFlux_dTemp to localFluxJacobian
      for( integer ic = 0; ic < numComp; ++ic )
      {
        integer const eqIndex0 = k[0]* numEqn + ic;
        integer const eqIndex1 = k[1]* numEqn + ic;
        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          integer const localDofIndexTemp = k[ke] * numDof + numDof - 1;
          stack.localFluxJacobian[eqIndex0][localDofIndexTemp] += m_dt * dCompFlux_dT[ke][ic];
          stack.localFluxJacobian[eqIndex1][localDofIndexTemp] -= m_dt * dCompFlux_dT[ke][ic];
        }
      }

      // Step 5: compute the enthalpy flux
      real64 const enthalpy = m_phaseEnthalpy[er_up][esr_up][ei_up][0][ip];
      convectiveEnergyFlux += phaseFlux * enthalpy;

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        dConvectiveEnergyFlux_dP[ke] += dPhaseFlux_dP[ke] * enthalpy;
        dConvectiveEnergyFlux_dT[ke] += dPhaseFlux_dT[ke] * enthalpy;

        for( integer jc = 0; jc < numComp; ++jc )
        {
          dConvectiveEnergyFlux_dC[ke][jc] += dPhaseFlux_dC[ke][jc] * enthalpy;
        }
      }

      dConvectiveEnergyFlux_dP[k_up] += phaseFlux * m_dPhaseEnthalpy[er_up][esr_up][ei_up][0][ip][Deriv::dP];
      dConvectiveEnergyFlux_dT[k_up] += phaseFlux * m_dPhaseEnthalpy[er_up][esr_up][ei_up][0][ip][Deriv::dT];

      real64 dProp_dC[numComp]{};
      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er_up][esr_up][ei_up],
                      m_dPhaseEnthalpy[er_up][esr_up][ei_up][0][ip],
                      dProp_dC,
                      Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dConvectiveEnergyFlux_dC[k_up][jc] += phaseFlux * dProp_dC[jc];
      }

      // Step 6: add convectiveFlux and its derivatives to localFlux and localFluxJacobian
      integer const localRowIndexEnergy0 = k[0] * numEqn + numEqn - 1;
      integer const localRowIndexEnergy1 = k[1] * numEqn + numEqn - 1;
      stack.localFlux[localRowIndexEnergy0] +=  m_dt * convectiveEnergyFlux;
      stack.localFlux[localRowIndexEnergy1] -=  m_dt * convectiveEnergyFlux;

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        integer const localDofIndexPres = k[ke] * numDof;
        stack.localFluxJacobian[localRowIndexEnergy0][localDofIndexPres] += m_dt * dConvectiveEnergyFlux_dP[ke];
        stack.localFluxJacobian[localRowIndexEnergy1][localDofIndexPres] -= m_dt * dConvectiveEnergyFlux_dP[ke];
        integer const localDofIndexTemp = localDofIndexPres + numDof - 1;
        stack.localFluxJacobian[localRowIndexEnergy0][localDofIndexTemp] +=  m_dt * dConvectiveEnergyFlux_dT[ke];
        stack.localFluxJacobian[localRowIndexEnergy1][localDofIndexTemp] -=  m_dt * dConvectiveEnergyFlux_dT[ke];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          integer const localDofIndexComp = localDofIndexPres + jc + 1;
          stack.localFluxJacobian[localRowIndexEnergy0][localDofIndexComp] += m_dt * dConvectiveEnergyFlux_dC[ke][jc];
          stack.localFluxJacobian[localRowIndexEnergy1][localDofIndexComp] -= m_dt * dConvectiveEnergyFlux_dC[ke][jc];
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
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     m_thermalConductivity, // we have to pass something here, so we just use thermal conductivity
                                     stack.thermalTransmissibility,
                                     stack.dTrans_dPres ); // again, we have to pass something here, but this is unused for now



    localIndex k[2]{};
    localIndex connectionIndex = 0;

    for( k[0] = 0; k[0] < stack.numConnectedElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numConnectedElems; ++k[1] )
      {
        real64 const thermalTrans[2] = { stack.thermalTransmissibility[connectionIndex][0], stack.thermalTransmissibility[connectionIndex][1] };
        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        real64 conductiveEnergyFlux = 0.0;
        real64 dConductiveEnergyFlux_dT[numFluxSupportPoints]{};

        // Step 2: compute temperature difference at the interface
        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          localIndex const er  = seri[ke];
          localIndex const esr = sesri[ke];
          localIndex const ei  = sei[ke];

          conductiveEnergyFlux += thermalTrans[ke] * m_temp[er][esr][ei];
          dConductiveEnergyFlux_dT[ke] += thermalTrans[ke];
        }

        // Step 3: add conductiveFlux and its derivatives to localFlux and localFluxJacobian
        integer const localRowIndexEnergy0 = k[0] * numEqn + numEqn - 1;
        integer const localRowIndexEnergy1 = k[1] * numEqn + numEqn - 1;
        stack.localFlux[localRowIndexEnergy0] +=  m_dt * conductiveEnergyFlux;
        stack.localFlux[localRowIndexEnergy1] -=  m_dt * conductiveEnergyFlux;

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          integer const localDofIndexTemp = k[ke] * numDof + numDof - 1;
          stack.localFluxJacobian[localRowIndexEnergy0][localDofIndexTemp] +=  m_dt * dConductiveEnergyFlux_dT[ke];
          stack.localFluxJacobian[localRowIndexEnergy1][localDofIndexTemp] -=  m_dt * dConductiveEnergyFlux_dT[ke];
        }
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
  void complete( localIndex const iconn,
                 StackVariables & stack ) const
  {
    // Call Case::complete to assemble the component mass balance equations (i = 0 to i = numDof-2)
    // In the lambda, add contribution to residual and jacobian into the energy balance equation
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      // beware, there is  volume balance eqn in m_localRhs and m_localMatrix!
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn], stack.localFlux[i * numEqn + numEqn-1] );
      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
        ( localRow + numEqn,
        stack.dofColIndices.data(),
        stack.localFluxJacobian[i * numEqn + numEqn-1].dataIfContiguous(),
        stack.stencilSize * numDof );

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
 * @class FluxComputeKernelFactory
 */
class FluxComputeKernelFactory
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
      integer constexpr NUM_DOF = NC() + 2;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      BitFlags< isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelFlags > kernelFlags;
      if( hasCapPressure )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelFlags::CapPressure );
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelFlags::TotalMassEquation );

      using KernelType = FluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KernelType::ThermalCompFlowAccessors thermalCompFlowAccessors( elemManager, solverName );
      typename KernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KernelType::ThermalMultiFluidAccessors thermalMultiFluidAccessors( elemManager, solverName );
      typename KernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

      KernelType kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor,
                         compFlowAccessors, thermalCompFlowAccessors, multiFluidAccessors, thermalMultiFluidAccessors,
                         capPressureAccessors, permeabilityAccessors, thermalConductivityAccessors,
                         dt, localMatrix, localRhs, kernelFlags );
      KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace thermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALFLUXCOMPUTEKERNEL_HPP
