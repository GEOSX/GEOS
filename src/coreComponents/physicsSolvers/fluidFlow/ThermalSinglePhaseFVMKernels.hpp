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
 * @file ThermalSinglePhaseFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEFVMKERNELS_HPP

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"
#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"

namespace geos
{

namespace thermalSinglePhaseFVMKernels
{
/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_mob;
  using AbstractBase::m_dens;

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature,
                      fields::flow::dMobility_dTemperature >;

  using ThermalSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::dDensity_dTemperature,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy_dPressure,
                              fields::singlefluid::dEnthalpy_dTemperature >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< constitutive::SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity,
                              fields::thermalconductivity::dEffectiveConductivity_dT >;


  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] singlePhaseFlowAccessors accessor for wrappers registered by the solver
   * @param[in] thermalSinglePhaseFlowAccessors accessor for *thermal* wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors accessor for wrappers registered by the single fluid model
   * @param[in] thermalSinglePhaseFluidAccessors accessor for *thermal* wrappers registered by the single fluid model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] thermalConductivityAccessors accessor for wrappers registered by the thermal conductivity model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( globalIndex const rankOffset,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                           ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                           SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                           ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           ThermalConductivityAccessors const & thermalConductivityAccessors,
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
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_dMob_dTemp( thermalSinglePhaseFlowAccessors.get( fields::flow::dMobility_dTemperature {} ) ),
    m_dDens_dTemp( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dDensity_dTemperature {} ) ),
    m_enthalpy( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy_dPres( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy_dPressure {} ) ),
    m_dEnthalpy_dTemp( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy_dTemperature {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_dThermalCond_dT( thermalConductivityAccessors.get( fields::thermalconductivity::dEffectiveConductivity_dT {} ) )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems ),
      energyFlux( 0.0 ),
      dEnergyFlux_dP( size ),
      dEnergyFlux_dT( size )
    {}

    using Base::StackVariables::stencilSize;
    using Base::StackVariables::numFluxElems;
    using Base::StackVariables::transmissibility;
    using Base::StackVariables::dTrans_dPres;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;

    // Thermal transmissibility
    real64 thermalTransmissibility[maxNumConns][2]{};

    /// Derivatives of thermal transmissibility with respect to temperature
    real64 dThermalTrans_dT[maxNumConns][2]{};

    // Energy fluxes and derivatives

    /// Energy fluxes
    real64 energyFlux;
    /// Derivatives of energy fluxes wrt pressure
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dP;
    /// Derivatives of energy fluxes wrt temperature
    stackArray1d< real64, maxStencilSize > dEnergyFlux_dT;

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
                                           localIndex const connectionIndex,
                                           real64 const alpha,
                                           real64 const mobility,
                                           real64 const & potGrad,
                                           real64 const & fluxVal,
                                           real64 const (&dFlux_dP)[2] )
    {
      // Step 1: compute the derivatives of the mean density at the interface wrt temperature

      real64 dDensMean_dT[2]{0.0, 0.0};

      real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };

      for( integer ke = 0; ke < 2; ++ke )
      {
        real64 const dDens_dT = m_dDens_dTemp[seri[ke]][sesri[ke]][sei[ke]][0];
        dDensMean_dT[ke] = 0.5 * dDens_dT;
      }

      // Step 2: compute the derivatives of the potential difference wrt temperature
      //***** calculation of flux *****

      real64 dGravHead_dT[2]{0.0, 0.0};

      // compute potential difference
      for( integer ke = 0; ke < 2; ++ke )
      {
        localIndex const er  = seri[ke];
        localIndex const esr = sesri[ke];
        localIndex const ei  = sei[ke];

        // compute derivative of gravity potential difference wrt temperature
        real64 const gravD = trans[ke] * m_gravCoef[er][esr][ei];

        for( integer i = 0; i < 2; ++i )
        {
          dGravHead_dT[i] += dDensMean_dT[i] * gravD;
        }
      }

      // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
      // *** upwinding ***

      real64 dFlux_dT[2]{0.0, 0.0};

      // Step 3.1: compute the derivative of flux wrt temperature
      for( integer ke = 0; ke < 2; ++ke )
      {
        dFlux_dT[ke] -= dGravHead_dT[ke];
      }

      for( integer ke = 0; ke < 2; ++ke )
      {
        dFlux_dT[ke] *= mobility;
      }

      real64 dMob_dT[2]{};

      if( alpha <= 0.0 || alpha >= 1.0 )
      {
        localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );

        dMob_dT[k_up] = m_dMob_dTemp[seri[k_up]][sesri[k_up]][sei[k_up]];
      }
      else
      {
        real64 const mobWeights[2] = { alpha, 1.0 - alpha };
        for( integer ke = 0; ke < 2; ++ke )
        {
          dMob_dT[ke] = mobWeights[ke] * m_dMob_dTemp[seri[ke]][sesri[ke]][sei[ke]];
        }
      }

      // add contribution from upstream cell mobility derivatives
      for( integer ke = 0; ke < 2; ++ke )
      {
        dFlux_dT[ke] += dMob_dT[ke] * potGrad;
      }

      // add dFlux_dTemp to localFluxJacobian
      for( integer ke = 0; ke < 2; ++ke )
      {
        localIndex const localDofIndexTemp = k[ke] * numDof + numDof - 1;
        stack.localFluxJacobian[k[0]*numEqn][localDofIndexTemp] += m_dt * dFlux_dT[ke];
        stack.localFluxJacobian[k[1]*numEqn][localDofIndexTemp] -= m_dt * dFlux_dT[ke];
      }

      // Step 4: compute the enthalpy flux
      real64 enthalpy = 0.0;
      real64 dEnthalpy_dP[2]{0.0, 0.0};
      real64 dEnthalpy_dT[2]{0.0, 0.0};

      if( alpha <= 0.0 || alpha >= 1.0 )
      {
        localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );

        enthalpy = m_enthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0];
        dEnthalpy_dP[k_up] = m_dEnthalpy_dPres[seri[k_up]][sesri[k_up]][sei[k_up]][0];
        dEnthalpy_dT[k_up] = m_dEnthalpy_dTemp[seri[k_up]][sesri[k_up]][sei[k_up]][0];
      }
      else
      {
        real64 const mobWeights[2] = { alpha, 1.0 - alpha };
        for( integer ke = 0; ke < 2; ++ke )
        {
          enthalpy += mobWeights[ke] * m_enthalpy[seri[ke]][sesri[ke]][sei[ke]][0];
          dEnthalpy_dP[ke] = mobWeights[ke] * m_dEnthalpy_dPres[seri[ke]][sesri[ke]][sei[ke]][0];
          dEnthalpy_dT[ke] = mobWeights[ke] * m_dEnthalpy_dTemp[seri[ke]][sesri[ke]][sei[ke]][0];
        }
      }

      stack.energyFlux += fluxVal * enthalpy;

      for( integer ke = 0; ke < 2; ++ke )
      {
        stack.dEnergyFlux_dP[ke] += dFlux_dP[ke] * enthalpy;
        stack.dEnergyFlux_dT[ke] += dFlux_dT[ke] * enthalpy;
      }

      for( integer ke = 0; ke < 2; ++ke )
      {
        stack.dEnergyFlux_dP[ke] += fluxVal * dEnthalpy_dP[ke];
        stack.dEnergyFlux_dT[ke] += fluxVal * dEnthalpy_dT[ke];
      }

    } );

    // *****************************************************
    // Computation of the conduction term in the energy flux
    // Note that the enthalpy term in the energy was computed above
    // Note that this term is computed using an explicit treatment of conductivity for now

    // Step 1: compute the thermal transmissibilities at this face
    // We follow how the thermal compositional multi-phase solver does to update the thermal transmissibility
    m_stencilWrapper.computeWeights( iconn,
                                     m_thermalConductivity,
                                     m_dThermalCond_dT,
                                     stack.thermalTransmissibility,
                                     stack.dThermalTrans_dT );

    localIndex k[2];
    localIndex connectionIndex = 0;

    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        real64 const thermalTrans[2] = { stack.thermalTransmissibility[connectionIndex][0], stack.thermalTransmissibility[connectionIndex][1] };
        real64 const dThermalTrans_dT[2] = { stack.dThermalTrans_dT[connectionIndex][0], stack.dThermalTrans_dT[connectionIndex][1] };

        localIndex const seri[2]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // Step 2: compute temperature difference at the interface
        for( integer ke = 0; ke < 2; ++ke )
        {
          localIndex const er  = seri[ke];
          localIndex const esr = sesri[ke];
          localIndex const ei  = sei[ke];

          stack.energyFlux += thermalTrans[ke] * m_temp[er][esr][ei];
          stack.dEnergyFlux_dT[ke] += thermalTrans[ke] + dThermalTrans_dT[ke] * m_temp[er][esr][ei];
        }

        // add energyFlux and its derivatives to localFlux and localFluxJacobian
        stack.localFlux[k[0]*numEqn + numEqn - 1] += m_dt * stack.energyFlux;
        stack.localFlux[k[1]*numEqn + numEqn - 1] -= m_dt * stack.energyFlux;

        for( integer ke = 0; ke < 2; ++ke )
        {
          integer const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[k[0]*numEqn + numEqn - 1][localDofIndexPres] =  m_dt * stack.dEnergyFlux_dP[ke];
          stack.localFluxJacobian[k[1]*numEqn + numEqn - 1][localDofIndexPres] = -m_dt * stack.dEnergyFlux_dP[ke];
          integer const localDofIndexTemp = localDofIndexPres + numDof - 1;
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
    // Call Case::complete to assemble the mass balance equations
    // In the lambda, add contribution to residual and jacobian into the energy balance equation
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {
      // The no. of fluxes is equal to the no. of equations in m_localRhs and m_localMatrix
      // Different from the one in compositional multi-phase flow, which has a volume balance eqn.
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn-1], stack.localFlux[i * numEqn + numEqn-1] );

      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + numEqn-1,
                                                                                        stack.dofColIndices.data(),
                                                                                        stack.localFluxJacobian[i * numEqn + numEqn-1].dataIfContiguous(),
                                                                                        stack.stencilSize * numDof );

    } );
  }

protected:

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

  /// View on derivatives of thermal conductivity w.r.t. temperature
  ElementViewConst< arrayView3d< real64 const > > m_dThermalCond_dT;

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
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 2;
    integer constexpr NUM_EQN = 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using KernelType = FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename KernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename KernelType::ThermalSinglePhaseFlowAccessors thermalFlowAccessors( elemManager, solverName );
    typename KernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename KernelType::ThermalSinglePhaseFluidAccessors thermalFluidAccessors( elemManager, solverName );
    typename KernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
    typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

    KernelType kernel( rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, thermalFlowAccessors, fluidAccessors, thermalFluidAccessors,
                       permAccessors, thermalConductivityAccessors,
                       dt, localMatrix, localRhs );
    KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

/******************************** DirichletFaceBasedAssemblyKernel ********************************/

/**
 * @class DirichFaceBasedAssemblyKernel
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFaceBasedAssemblyKernel : public singlePhaseFVMKernels::DirichletFaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, FLUIDWRAPPER >
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

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_mob;
  using AbstractBase::m_pres;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;

  using AbstractBase::m_localMatrix;
  using AbstractBase::m_localRhs;

  using Base = singlePhaseFVMKernels::DirichletFaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, FLUIDWRAPPER >;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::m_facePres;
  using Base::m_faceGravCoef;

  using ThermalSinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::temperature,
                      fields::flow::dMobility_dTemperature >;

  using ThermalSinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::dDensity_dTemperature,
                              fields::singlefluid::enthalpy,
                              fields::singlefluid::dEnthalpy_dPressure,
                              fields::singlefluid::dEnthalpy_dTemperature >;

  using ThermalConductivityAccessors =
    StencilMaterialAccessors< constitutive::SinglePhaseThermalConductivityBase,
                              fields::thermalconductivity::effectiveConductivity,
                              fields::thermalconductivity::dEffectiveConductivity_dT >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of the MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor the degree of freedom number accessor
   * @param[in] singlePhaseFlowAccessors the single phase flow accessor
   * @param[in] thermalSinglePhaseFlowAccessors the thermal single phase flow accessor
   * @param[in] singlePhaseFluidAccessors the single phase fluid accessor
   * @param[in] thermalSinglePhaseFluidAccessors the thermal single phase fluid accessor
   * @param[in] permeabilityAccessors the permeability accessor
   * @param[in] thermalConductivityAccessors the thermal conductivity accessor
   * @param[in] dt the time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  DirichletFaceBasedAssemblyKernel( globalIndex const rankOffset,
                                    FaceManager const & faceManager,
                                    BoundaryStencilWrapper const & stencilWrapper,
                                    FLUIDWRAPPER const & fluidWrapper,
                                    DofNumberAccessor const & dofNumberAccessor,
                                    SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                    ThermalSinglePhaseFlowAccessors const & thermalSinglePhaseFlowAccessors,
                                    SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                    ThermalSinglePhaseFluidAccessors const & thermalSinglePhaseFluidAccessors,
                                    PermeabilityAccessors const & permeabilityAccessors,
                                    ThermalConductivityAccessors const & thermalConductivityAccessors,
                                    real64 const & dt,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )

    : Base( rankOffset,
            faceManager,
            stencilWrapper,
            fluidWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_temp( thermalSinglePhaseFlowAccessors.get( fields::flow::temperature {} ) ),
    m_faceTemp( faceManager.getField< fields::flow::faceTemperature >() ),
    m_dMob_dTemp( thermalSinglePhaseFlowAccessors.get( fields::flow::dMobility_dTemperature {} ) ),
    m_dDens_dTemp( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dDensity_dTemperature {} ) ),
    m_enthalpy( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::enthalpy {} ) ),
    m_dEnthalpy_dPres( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy_dPressure {} ) ),
    m_dEnthalpy_dTemp( thermalSinglePhaseFluidAccessors.get( fields::singlefluid::dEnthalpy_dTemperature {} ) ),
    m_thermalConductivity( thermalConductivityAccessors.get( fields::thermalconductivity::effectiveConductivity {} ) ),
    m_dThermalCond_dT( thermalConductivityAccessors.get( fields::thermalconductivity::dEffectiveConductivity_dT {} ) )
  {}


  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : Base::StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size,
                    localIndex numElems ):
      Base::StackVariables( size,
                            numElems )
    {}

    using Base::StackVariables::localFlux;
    using Base::StackVariables::localFluxJacobian;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::transmissibility;

    /// Energy fluxes and derivatives wrt pressure and temperature
    real64 energyFlux = 0.0;
    real64 dEnergyFlux_dP = 0.0;
    real64 dEnergyFlux_dT = 0.0;
  };

  /**
   * @brief Compute the local Dirichlet face flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack ) const
  {
    Base::computeFlux( iconn, stack, [&] ( localIndex const er,
                                           localIndex const esr,
                                           localIndex const ei,
                                           localIndex const kf,
                                           real64 const & f,
                                           real64 const & dF_dP,
                                           real64 const & mobility_up,
                                           real64 const & dMobility_dP_up )
    {

      // Compute the derivatives of the density wrt temperature

      real64 const dDens_dT = 0.5 * m_dDens_dTemp[er][esr][ei][0];

      // Compute the derivatives of the phase potential difference wrt temperature

      real64 const dF_dT = -stack.transmissibility * dDens_dT * ( m_gravCoef[er][esr][ei] - m_faceGravCoef[kf] );

      // Compute the (upwinded) energy flux

      real64 const flux = mobility_up * f;
      real64 const enthalpy = m_enthalpy[er][esr][ei][0];
      stack.energyFlux += flux * enthalpy;

      // Compute the derivatives of the (upwinded) energy flux wrt pressure and temperature

      if( f >= 0 ) // the element is upstream
      {
        real64 const dFlux_dP = mobility_up * dF_dP + dMobility_dP_up * f;
        real64 const dFlux_dT = mobility_up * dF_dT + m_dMob_dTemp[er][esr][ei] * f;

        stack.dEnergyFlux_dP += dFlux_dP * enthalpy + flux * m_dEnthalpy_dPres[er][esr][ei][0];
        stack.dEnergyFlux_dT += dFlux_dT * enthalpy + flux * m_dEnthalpy_dTemp[er][esr][ei][0];
      }
      else
      {
        real64 const dFlux_dP = mobility_up * dF_dP;
        real64 const dFlux_dT = mobility_up * dF_dT;

        stack.dEnergyFlux_dP += dFlux_dP * enthalpy;
        stack.dEnergyFlux_dT += dFlux_dT * enthalpy;
      }

      // Contribution of energy conduction through the solid phase
      real64 thermalTrans = 0.0;
      real64 dThermalTrans_dThermalCond[3]{};
      m_stencilWrapper.computeWeights( iconn,
                                       m_thermalConductivity,
                                       thermalTrans,
                                       dThermalTrans_dThermalCond );

      real64 const dThermalTrans_dT = LvArray::tensorOps::AiBi< 3 >( dThermalTrans_dThermalCond, m_dThermalCond_dT[er][esr][ei][0] );

      real64 const deltaT = m_temp[er][esr][ei] - m_faceTemp[kf];
      stack.energyFlux += thermalTrans * deltaT;
      stack.dEnergyFlux_dT += thermalTrans + dThermalTrans_dT * deltaT;

      // Add energyFlux and its derivatives to localFlux and localFluxJacobian
      integer const localRowIndexEnergy = numEqn - 1;
      stack.localFlux[localRowIndexEnergy] =  m_dt * stack.energyFlux;

      stack.localFluxJacobian[localRowIndexEnergy][0] =  m_dt * stack.dEnergyFlux_dP;
      stack.localFluxJacobian[localRowIndexEnergy][numDof-1] =  m_dt * stack.dEnergyFlux_dT;
    } );

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
    Base::complete( iconn, stack, [&] ( localIndex const localRow )
    {
      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + numEqn - 1], stack.localFlux[numEqn-1] );

      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
        ( localRow + numEqn - 1,
        stack.dofColIndices,
        stack.localFluxJacobian[numEqn-1],
        numDof );
    } );
  }

protected:

  /// Views on temperature
  ElementViewConst< arrayView1d< real64 const > > const m_temp;

  /// Views on face temperature
  arrayView1d< real64 const > const m_faceTemp;

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

  /// View on derivatives of thermal conductivity w.r.t. temperature
  ElementViewConst< arrayView3d< real64 const > > m_dThermalCond_dT;

};


/**
 * @class DirichletFaceBasedAssemblyKernelFactory
 */
class DirichletFaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] faceManager reference to the face manager
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the boundary stencil wrapper
   * @param[in] fluidBase the single phase fluid constitutive model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   BoundaryStencilWrapper const & stencilWrapper,
                   constitutive::SingleFluidBase & fluidBase,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

      integer constexpr NUM_DOF = 2;
      integer constexpr NUM_EQN = 2;

      using kernelType = DirichletFaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, typename FluidType::KernelWrapper >;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );

      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      typename kernelType::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, solverName );
      typename kernelType::ThermalSinglePhaseFlowAccessors thermalSinglePhaseFlowAccessors( elemManager, solverName );
      typename kernelType::SinglePhaseFluidAccessors singlePhaseFluidAccessors( elemManager, solverName );
      typename kernelType::ThermalSinglePhaseFluidAccessors thermalSinglePhaseFluidAccessors( elemManager, solverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename kernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

      kernelType kernel( rankOffset,
                         faceManager,
                         stencilWrapper,
                         fluidWrapper,
                         dofNumberAccessor,
                         singlePhaseFlowAccessors,
                         thermalSinglePhaseFlowAccessors,
                         singlePhaseFluidAccessors,
                         thermalSinglePhaseFluidAccessors,
                         permeabilityAccessors,
                         thermalConductivityAccessors,
                         dt,
                         localMatrix,
                         localRhs );

      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }

};

} // namespace thermalSinglePhaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALSINGLEPHASEFVMKERNELS_HPP
