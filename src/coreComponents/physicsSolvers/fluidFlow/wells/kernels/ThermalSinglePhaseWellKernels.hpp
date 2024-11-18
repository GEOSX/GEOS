/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALSINGLEPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALSINGLEPHASEWELLKERNELS_HPP

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"

namespace geos
{

namespace thermalSinglePhaseWellKernels
{



/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation and volume balance
 */
template< integer NUM_DOF >
class ElementBasedAssemblyKernel : public singlePhaseWellKernels::ElementBasedAssemblyKernel< NUM_DOF >
{
public:
  using Base = singlePhaseWellKernels::ElementBasedAssemblyKernel< NUM_DOF >;
  using Base::m_rankOffset;
  using Base::m_wellElemDofNumber;
  using Base::m_elemGhostRank;
  using Base::m_wellElemVolume;
  using Base::m_wellElemDensity;
  using Base::m_wellElemDensity_n;
  using Base::m_dWellElemDensity_dPressure;
  using Base::m_localMatrix;
  using Base::m_localRhs;
  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_DOF;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    :    Base( rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs ),
    m_dWellElemDensity_dTemperature( fluid.dDensity_dTemperature()  ),
    m_internalEnergy( fluid.internalEnergy() ),
    m_internalEnergy_n( fluid.internalEnergy_n() ),
    m_dInternalEnergy_dPres( fluid.dInternalEnergy_dPressure() ),
    m_dInternalEnergy_dTemp( fluid.dInternalEnergy_dTemperature() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    {}
    using Base::StackVariables::eqnRowIndices;
    using Base::StackVariables::dofColIndices;
    using Base::StackVariables::localJacobian;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localRow;
    using Base::StackVariables::volume;
    using Base::StackVariables::density;
    using Base::StackVariables::density_n;
    using Base::StackVariables::dDensity_dPres;

  };
  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }



  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] phaseAmountKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const iwelem,
                            StackVariables & stack ) const
  {
    Base::computeAccumulation( iwelem, stack, [&]( )
    {

      // Step 1: assemble the derivatives of the mass balance equation w.r.t temperature
      stack.localJacobian[0][numDof-1] =  stack.volume * m_dWellElemDensity_dTemperature[iwelem][0];

      // Step 2: assemble the fluid part of the accumulation term of the energy equation
      real64 const fluidEnergy = stack.volume   * stack.density  * m_internalEnergy[iwelem][0];
      real64 const fluidEnergy_n = stack.volume   * stack.density_n  * m_internalEnergy_n[iwelem][0];

      real64 const dFluidEnergy_dP =  stack.volume   * stack.dDensity_dPres  * m_internalEnergy[iwelem][0]
                                     + stack.volume   * stack.density  * m_dInternalEnergy_dPres[iwelem][0];


      real64 const dFluidEnergy_dT = stack.volume   * m_dWellElemDensity_dTemperature[iwelem][0] * m_internalEnergy[iwelem][0]
                                     + stack.volume  * stack.density  * m_dInternalEnergy_dTemp[iwelem][0];

      // local accumulation
      stack.localResidual[numEqn-1] = fluidEnergy - fluidEnergy_n;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn-1][0]        = dFluidEnergy_dP;
      stack.localJacobian[numEqn-1][numDof-1] = dFluidEnergy_dT;
    } );
  }



  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
    {
      if( kernelComponent.elemGhostRank( iwelem ) >= 0 )
      {
        return;
      }
      typename KERNEL_TYPE::StackVariables stack;
      kernelComponent.setup( iwelem, stack );
      kernelComponent.computeAccumulation( iwelem, stack );
      kernelComponent.complete( iwelem, stack );

    } );
  }

protected:

  /// View on derivative of fluid density w.r.t temperature
  arrayView2d< real64 const > const m_dWellElemDensity_dTemperature;

  /// Views on fluid internal energy
  arrayView2d< real64 const > const m_internalEnergy;
  arrayView2d< real64 const > const m_internalEnergy_n;
  arrayView2d< real64 const > const m_dInternalEnergy_dPres;
  arrayView2d< real64 const > const m_dInternalEnergy_dTemp;

};


/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:
  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 2;
    ElementBasedAssemblyKernel< NUM_DOF >
    kernel( rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs );
    ElementBasedAssemblyKernel< NUM_DOF >::template
    launch< POLICY, ElementBasedAssemblyKernel< NUM_DOF > >( subRegion.size(), kernel );

  }
};
} // end namespace singlePhaseWellKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
