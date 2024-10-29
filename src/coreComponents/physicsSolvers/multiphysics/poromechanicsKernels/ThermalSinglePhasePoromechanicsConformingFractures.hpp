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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP

#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsConformingFractures.hpp"

namespace geos
{

namespace thermalSinglePhasePoromechanicsConformingFracturesKernels
{

using namespace fluxKernelsHelper;
using namespace constitutive;

template< integer NUM_EQN, integer NUM_DOF >
class ConnectorBasedAssemblyKernel : public singlePhasePoromechanicsConformingFracturesKernels::ConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF >
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

  using SinglePhaseFVMAbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = SinglePhaseFVMAbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = SinglePhaseFVMAbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = SinglePhaseFVMAbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = SinglePhaseFVMAbstractBase::PermeabilityAccessors;
  using FracturePermeabilityAccessors = StencilMaterialAccessors< PermeabilityBase,
                                                                  fields::permeability::dPerm_dDispJump >;
  using SinglePhaseFVMAbstractBase::m_dt;
  using SinglePhaseFVMAbstractBase::m_rankOffset;
  using SinglePhaseFVMAbstractBase::m_dofNumber;
  using SinglePhaseFVMAbstractBase::m_gravCoef;
  using SinglePhaseFVMAbstractBase::m_mob;
  using SinglePhaseFVMAbstractBase::m_dens;

  using SinglePhaseFVMBase = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >;
  using SinglePhaseFVMBase::numDof;
  using SinglePhaseFVMBase::numEqn;
  using SinglePhaseFVMBase::maxNumElems;
  using SinglePhaseFVMBase::maxNumConns;
  using SinglePhaseFVMBase::maxStencilSize;
  using SinglePhaseFVMBase::m_stencilWrapper;
  using SinglePhaseFVMBase::m_seri;
  using SinglePhaseFVMBase::m_sesri;
  using SinglePhaseFVMBase::m_sei;
  using SinglePhaseFVMBase::m_ghostRank;

  using Base = singlePhasePoromechanicsConformingFracturesKernels::ConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF >;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature,
                      fields::flow::dMobility_dTemperature >;

  using ThermalSinglePhaseFluidAccessors =
    StencilMaterialAccessors< SingleFluidBase,
                              fields::singlefluid::dDensity_dTemperature,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy_dPressure,
                              fields::singlefluid::dEnthalpy_dTemperature >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity >;



  ConnectorBasedAssemblyKernel( globalIndex const rankOffset,
                                SurfaceElementStencilWrapper const & stencilWrapper,
                                DofNumberAccessor const & flowDofNumberAccessor,
                                SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                                SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                                PermeabilityAccessors const & permeabilityAccessors,
                                FracturePermeabilityAccessors const & edfmPermeabilityAccessors,
                                ThermalConductivityAccessors const & thermalConductivityAccessors,
                                real64 const & dt,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs,
                                CRSMatrixView< real64, localIndex const > const & dR_dAper )
    : Base( rankOffset,
            stencilWrapper,
            flowDofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            edfmPermeabilityAccessors,
            dt,
            localMatrix,
            localRhs,
            dR_dAper ),
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_dMob_dTemp( thermalSinglePhaseFlowAccessors.get( fields::flow::dMobility_dTemperature {} ) ),
    m_dDens_dTemp( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dDensity_dTemperature {} ) ),
    m_enthalpy( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy_dPres( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy_dPressure {} ) ),
    m_dEnthalpy_dTemp( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy_dTemperature {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) )
  {}


  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
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
      : Base::StackVariables( size, numElems ),
      energyFlux( 0.0 ),
      dEnergyFlux_dTrans( 0.0 ),
      dEnergyFlux_dP( size ),
      dEnergyFlux_dT( size ),
      dEnergyFlux_dDispJump( size, 3 )
    {}
    using SinglePhaseFVMBase::StackVariables::stencilSize;
    using SinglePhaseFVMBase::StackVariables::numFluxElems;
    using SinglePhaseFVMBase::StackVariables::transmissibility;
    using SinglePhaseFVMBase::StackVariables::dTrans_dPres;
    using SinglePhaseFVMBase::StackVariables::dofColIndices;
    using SinglePhaseFVMBase::StackVariables::localFlux;
    using SinglePhaseFVMBase::StackVariables::localFluxJacobian;

    // Thermal transmissibility (for now, no derivatives)

    real64 thermalTransmissibility[maxNumConns][2]{};

    // Energy fluxes and derivatives

    /// Energy fluxes
    real64 energyFlux;
    /// Derivative of the Energy fluxes wrt transmissibility
    real64 dEnergyFlux_dTrans;
    /// Derivatives of energy fluxes wrt pressure
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dP;
    /// Derivatives of energy fluxes wrt temperature
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dT;
    /// Derivatives of energy fluxes wrt dispJump
    stackArray2d< real64, maxStencilSize *3 > dEnergyFlux_dDispJump{};

  };

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const

  {
    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives (including derivatives wrt temperature),
    //  2) enthalpy part of energyFlux  and its derivatives (including derivatives wrt temperature)
    //
    // Computing dFlux_dT and the enthalpy flux requires quantities already computed in the base computeFlux,
    // such as potGrad, fluxVal, and the indices of the upwind cell
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to access these variables
    Base::computeFlux( iconn, stack, [&] ( localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const,
                                           real64 const & alpha,
                                           real64 const & mobility,
                                           real64 const & potGrad,
                                           real64 const & massFlux,
                                           real64 const & dMassFlux_dTrans,
                                           real64 const (&dMassFlux_dP)[2] )
    {
      real64 trans[2] = {stack.transmissibility[0][0], stack.transmissibility[0][1]};
      real64 dMassFlux_dT[2]{};

      computeEnthalpyFlux( seri, sesri, sei,
                           trans,
                           m_enthalpy,
                           m_dEnthalpy_dPres,
                           m_dEnthalpy_dTemp,
                           m_gravCoef,
                           m_dDens_dTemp,
                           m_dMob_dTemp,
                           alpha,
                           mobility,
                           potGrad,
                           massFlux,
                           dMassFlux_dTrans,
                           dMassFlux_dP,
                           dMassFlux_dT,
                           stack.energyFlux,
                           stack.dEnergyFlux_dTrans,
                           stack.dEnergyFlux_dP,
                           stack.dEnergyFlux_dT );

      // add dMassFlux_dT to localFluxJacobian
      for( integer ke = 0; ke < 2; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numDof - 1;
        stack.localFluxJacobian[k[0]*numEqn][localDofIndexTemp] += m_dt * dMassFlux_dT[ke];
        stack.localFluxJacobian[k[1]*numEqn][localDofIndexTemp] -= m_dt * dMassFlux_dT[ke];
      }
    } );

    // *****************************************************
    // Computation of the conduction term in the energy flux
    // Note that the enthalpy term in the energy was computed above
    // Note that this term is computed using an explicit treatment of conductivity for now

    // Step 1: compute the thermal transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     m_thermalConductivity, // we have to pass something here, so we just use thermal conductivity
                                     stack.thermalTransmissibility,
                                     stack.dTrans_dPres ); // again, we have to pass something here, but this is unused for now

    localIndex k[2];
    localIndex connectionIndex = 0;

    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        real64 const thermalTrans[2] = { stack.thermalTransmissibility[connectionIndex][0], stack.thermalTransmissibility[connectionIndex][1] };

        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // Step 2: compute temperature difference at the interface
        computeConductiveFlux( seri, sesri, sei, m_temp, thermalTrans, stack.energyFlux, stack.dEnergyFlux_dT );

        // add energyFlux and its derivatives to localFlux and localFluxJacobian
        stack.localFlux[k[0]*numEqn + numEqn - 1] += m_dt * stack.energyFlux;
        stack.localFlux[k[1]*numEqn + numEqn - 1] -= m_dt * stack.energyFlux;

        for( integer ke = 0; ke < 2; ++ke )
        {
          integer const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[k[0]*numEqn + numEqn - 1][localDofIndexPres] =  m_dt * stack.dEnergyFlux_dP[ke];
          stack.localFluxJacobian[k[1]*numEqn + numEqn - 1][localDofIndexPres] = -m_dt * stack.dEnergyFlux_dP[ke];
          integer const localDofIndexTemp = localDofIndexPres + 1;
          stack.localFluxJacobian[k[0]*numEqn + numEqn - 1][localDofIndexTemp] =  m_dt * stack.dEnergyFlux_dT[ke];
          stack.localFluxJacobian[k[1]*numEqn + numEqn - 1][localDofIndexTemp] = -m_dt * stack.dEnergyFlux_dT[ke];
        }

        connectionIndex++;
      }
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
    // Call Base::complete to assemble the mass balance equations
    // In the lambda, add contribution to residual and jacobian into the energy balance equation
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      // The no. of fluxes is equal to the no. of equations in m_localRhs and m_localMatrix
      RAJA::atomicAdd( parallelDeviceAtomic{}, &SinglePhaseFVMAbstractBase::m_localRhs[localRow + numEqn-1], stack.localFlux[i * numEqn + numEqn-1] );

      SinglePhaseFVMAbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + numEqn-1,
                                                                                                      stack.dofColIndices.data(),
                                                                                                      stack.localFluxJacobian[i * numEqn + numEqn-1].dataIfContiguous(),
                                                                                                      stack.stencilSize * numDof );

    } );
  }


private:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on derivatives of fluid mobilities
  ElementViewConst< arrayView1d< real64 const > > const m_dMob_dTemp;

  /// Views on derivatives of fluid densities
  ElementViewConst< arrayView2d< real64 const > > const m_dDens_dTemp;

  /// Views on enthalpies
  ElementViewConst< arrayView2d< real64 const > > const m_enthalpy;
  ElementViewConst< arrayView2d< real64 const > > const m_dEnthalpy_dPres;
  ElementViewConst< arrayView2d< real64 const > > const m_dEnthalpy_dTemp;

  /// View on thermal conductivity
  ElementViewConst< arrayView3d< real64 const > > m_thermalConductivity;

};



/**
 * @class FaceBasedAssemblyKernelFactory
 */
class ConnectorBasedAssemblyKernelFactory
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
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & flowDofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   SurfaceElementStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, localIndex const > const & dR_dAper )
  {
    integer constexpr NUM_DOF = 2;   // pressure + temperature
    integer constexpr NUM_EQN = 2;   // mass balance + energy balance


    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > flowDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( flowDofKey );
    flowDofNumberAccessor.setName( solverName + "/accessors/" + flowDofKey );

    using kernelType = ConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors( elemManager, solverName );

    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::ThermalSinglePhaseFluidAccessors thermalFluidAccessors( elemManager, solverName );

    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
    typename kernelType::FracturePermeabilityAccessors edfmPermAccessors( elemManager, solverName );
    typename kernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper,
                       flowDofNumberAccessor,
                       flowAccessors, thermalFlowAccessors, fluidAccessors, thermalFluidAccessors,
                       permAccessors, edfmPermAccessors, thermalConductivityAccessors,
                       dt, localMatrix, localRhs, dR_dAper );

    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};



} // namespace SinglePhaseProppantFluxKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
