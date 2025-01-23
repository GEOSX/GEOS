/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ThermalAccumulationKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALACCUMULATIONKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALACCUMULATIONKERNELS_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/AccumulationKernels.hpp"

namespace geos
{

namespace thermalSinglePhaseBaseKernels
{

/******************************** AccumulationKernel ********************************/

/**
 * @class AccumulationKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< typename SUBREGION_TYPE, integer NUM_DOF >
class AccumulationKernel : public singlePhaseBaseKernels::AccumulationKernel< SUBREGION_TYPE, NUM_DOF >
{

public:

  using Base = singlePhaseBaseKernels::AccumulationKernel< SUBREGION_TYPE, NUM_DOF >;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_localMatrix;
  using Base::m_localRhs;
  using Base::m_dMass;

  /// Note: Derivative lineup only supports dP & dT, not component terms
  static constexpr integer isThermal = NUM_DOF-1;
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< isThermal >;
  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  AccumulationKernel( globalIndex const rankOffset,
                      string const dofKey,
                      SUBREGION_TYPE const & subRegion,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, localMatrix, localRhs ),
    m_energy( subRegion.template getField< fields::flow::energy >() ),
    m_energy_n( subRegion.template getField< fields::flow::energy_n >() ),
    m_dEnergy( subRegion.template getField< fields::flow::dEnergy >() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {};

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack );

    // assemble the derivatives of the mass balance equation w.r.t temperature
    stack.localJacobian[0][numDof-1] = m_dMass[ei][DerivOffset::dT];

    // assemble the accumulation term of the energy equation
    stack.localResidual[numEqn-1] = m_energy[ei] - m_energy_n[ei];
    stack.localJacobian[numEqn-1][0] += m_dEnergy[ei][DerivOffset::dP];
    stack.localJacobian[numEqn-1][numDof-1] += m_dEnergy[ei][DerivOffset::dT];
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the mass balance equation
    Base::complete( ei, stack );

    // Step 2: assemble the energy equation
    m_localRhs[stack.localRow + numEqn-1] += stack.localResidual[numEqn-1];
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow + numEqn-1,
                                                     stack.dofIndices,
                                                     stack.localJacobian[numEqn-1],
                                                     numDof );
  }

protected:

  /// View on energy
  arrayView1d< real64 const > const m_energy;
  arrayView1d< real64 const > const m_energy_n;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_dEnergy;

};

/**
 * @class SurfaceElementAccumulationKernel
 * @brief Define the interface for the assembly kernel in charge of accumulation in SurfaceElementSubRegion
 */
class SurfaceElementAccumulationKernel : public AccumulationKernel< SurfaceElementSubRegion, 2 >
{

public:

  using Base = AccumulationKernel< SurfaceElementSubRegion, 2 >;

  /**
   * @brief Constructor
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  SurfaceElementAccumulationKernel( globalIndex const rankOffset,
                                    string const dofKey,
                                    SurfaceElementSubRegion const & subRegion,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, localMatrix, localRhs ),
    m_creationMass( subRegion.getField< fields::flow::massCreated >() )
  {}

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            Base::StackVariables & stack ) const
  {
    Base::computeAccumulation( ei, stack );
    if( Base::m_mass_n[ei] > 1.1 * m_creationMass[ei] )
    {
      stack.localResidual[0] += m_creationMass[ei] * 0.25;
    }
  }

protected:

  arrayView1d< real64 const > const m_creationMass;

};

/**
 * @class AccumulationKernelFactory
 */
class AccumulationKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename SUBREGION_TYPE >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const dofKey,
                   SUBREGION_TYPE const & subRegion,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    if constexpr ( std::is_base_of_v< CellElementSubRegion, SUBREGION_TYPE > )
    {
      integer constexpr NUM_DOF = 2;
      AccumulationKernel< CellElementSubRegion, NUM_DOF > kernel( rankOffset, dofKey, subRegion, localMatrix, localRhs );
      AccumulationKernel< CellElementSubRegion, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
    }
    else if constexpr ( std::is_base_of_v< SurfaceElementSubRegion, SUBREGION_TYPE > )
    {
      SurfaceElementAccumulationKernel kernel( rankOffset, dofKey, subRegion, localMatrix, localRhs );
      SurfaceElementAccumulationKernel::launch< POLICY >( subRegion.size(), kernel );
    }
    else
    {
      GEOS_UNUSED_VAR( rankOffset, dofKey, subRegion, localMatrix, localRhs );
      GEOS_ERROR( "Unsupported subregion type: " << typeid(SUBREGION_TYPE).name() );
    }
  }

};

} // namespace thermalSinglePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALACCUMULATIONKERNELS_HPP
