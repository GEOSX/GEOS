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

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernel.hpp"

#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "constitutive/thermalConductivity/ThermalConductivityFields.hpp"

namespace geos
{

namespace thermalSinglePhaseFVMKernels
{

/******************************** DirichletFluxComputeKernel ********************************/

/**
 * @class DirichletFluxComputeKernel
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFluxComputeKernel : public singlePhaseFVMKernels::DirichletFluxComputeKernel< NUM_EQN, NUM_DOF, FLUIDWRAPPER >
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

  using AbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
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

  using Base = singlePhaseFVMKernels::DirichletFluxComputeKernel< NUM_EQN, NUM_DOF, FLUIDWRAPPER >;
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
  DirichletFluxComputeKernel( globalIndex const rankOffset,
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
 * @class DirichletFluxComputeKernelFactory
 */
class DirichletFluxComputeKernelFactory
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

      using kernelType = DirichletFluxComputeKernel< NUM_EQN, NUM_DOF, typename FluidType::KernelWrapper >;

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

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALDIRICHLETFLUXCOMPUTEKERNEL_HPP
