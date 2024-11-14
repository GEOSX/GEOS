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
  using Base::m_volume;
  using Base::m_deltaVolume;
  using Base::m_porosity;
  using Base::m_dPoro_dPres;
  using Base::m_density;
  using Base::m_dDensity_dPres;
  using Base::m_localMatrix;
  using Base::m_localRhs;

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
                      constitutive::SingleFluidBase const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs ),
    m_dDensity_dTemp( fluid.dDensity_dTemperature() ),
    m_dPoro_dTemp( solid.getDporosity_dTemperature() ),
    m_internalEnergy( fluid.internalEnergy() ),
    m_dInternalEnergy_dPres( fluid.dInternalEnergy_dPressure() ),
    m_dInternalEnergy_dTemp( fluid.dInternalEnergy_dTemperature() ),
    m_rockInternalEnergy( solid.getInternalEnergy() ),
    m_dRockInternalEnergy_dTemp( solid.getDinternalEnergy_dTemperature() ),
    m_energy_n( subRegion.template getField< fields::flow::energy_n >() )
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

    using Base::StackVariables::poreVolume;
    using Base::StackVariables::dPoreVolume_dPres;
    using Base::StackVariables::localRow;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localJacobian;

    /// Derivative of pore volume with respect to temperature
    real64 dPoreVolume_dTemp = 0.0;

    // Solid energy

    /// Solid energy at time n+1
    real64 solidEnergy = 0.0;

    /// Derivative of solid internal energy with respect to pressure
    real64 dSolidEnergy_dPres = 0.0;

    /// Derivative of solid internal energy with respect to temperature
    real64 dSolidEnergy_dTemp = 0.0;
  };


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.dPoreVolume_dTemp = ( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dTemp[ei][0];

    // initialize the solid volume
    real64 const solidVolume = ( m_volume[ei] + m_deltaVolume[ei] ) * ( 1.0 - m_porosity[ei][0] );
    real64 const dSolidVolume_dPres = -( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dPres[ei][0];
    real64 const dSolidVolume_dTemp = -( m_volume[ei] + m_deltaVolume[ei] ) * m_dPoro_dTemp[ei][0];

    // initialize the solid internal energy
    stack.solidEnergy = solidVolume * m_rockInternalEnergy[ei][0];
    stack.dSolidEnergy_dPres = dSolidVolume_dPres * m_rockInternalEnergy[ei][0];
    stack.dSolidEnergy_dTemp = solidVolume * m_dRockInternalEnergy_dTemp[ei][0] + dSolidVolume_dTemp * m_rockInternalEnergy[ei][0];
  }

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
    stack.localResidual[numEqn-1] = -m_energy_n[ei];

    Base::computeAccumulation( ei, stack, [&] ()
    {
      // Step 1: assemble the derivatives of the mass balance equation w.r.t temperature
      stack.localJacobian[0][numDof-1] = stack.poreVolume * m_dDensity_dTemp[ei][0] + stack.dPoreVolume_dTemp * m_density[ei][0];

      // Step 2: assemble the fluid part of the accumulation term of the energy equation
      real64 const fluidEnergy = stack.poreVolume * m_density[ei][0] * m_internalEnergy[ei][0];

      real64 const dFluidEnergy_dP = stack.dPoreVolume_dPres * m_density[ei][0] * m_internalEnergy[ei][0]
                                     + stack.poreVolume * m_dDensity_dPres[ei][0] * m_internalEnergy[ei][0]
                                     + stack.poreVolume * m_density[ei][0] * m_dInternalEnergy_dPres[ei][0];

      real64 const dFluidEnergy_dT = stack.poreVolume * m_dDensity_dTemp[ei][0] * m_internalEnergy[ei][0]
                                     + stack.poreVolume * m_density[ei][0] * m_dInternalEnergy_dTemp[ei][0]
                                     + stack.dPoreVolume_dTemp * m_density[ei][0] * m_internalEnergy[ei][0];

      // local accumulation
      stack.localResidual[numEqn-1] += fluidEnergy;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn-1][0]        = dFluidEnergy_dP;
      stack.localJacobian[numEqn-1][numDof-1] = dFluidEnergy_dT;
    } );

    // Step 3: assemble the solid part of the accumulation term of the energy equation
    stack.localResidual[numEqn-1] += stack.solidEnergy;
    stack.localJacobian[numEqn-1][0] += stack.dSolidEnergy_dPres;
    stack.localJacobian[numEqn-1][numDof-1] += stack.dSolidEnergy_dTemp;
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

  /// View on derivative of fluid density w.r.t temperature
  arrayView2d< real64 const > const m_dDensity_dTemp;

  /// View on derivative of porosity w.r.t temperature
  arrayView2d< real64 const > const m_dPoro_dTemp;

  /// Views on fluid internal energy
  arrayView2d< real64 const > const m_internalEnergy;
  arrayView2d< real64 const > const m_dInternalEnergy_dPres;
  arrayView2d< real64 const > const m_dInternalEnergy_dTemp;

  /// Views on rock internal energy
  arrayView2d< real64 const > const m_rockInternalEnergy;
  arrayView2d< real64 const > const m_dRockInternalEnergy_dTemp;

  /// View on energy
  arrayView1d< real64 const > const m_energy_n;

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
                                    constitutive::SingleFluidBase const & fluid,
                                    constitutive::CoupledSolidBase const & solid,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
    : Base( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs ),
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
                   constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    if constexpr ( std::is_base_of_v< CellElementSubRegion, SUBREGION_TYPE > )
    {
      integer constexpr NUM_DOF = 2;
      AccumulationKernel< CellElementSubRegion, NUM_DOF > kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      AccumulationKernel< CellElementSubRegion, NUM_DOF >::template launch< POLICY >( subRegion.size(), kernel );
    }
    else if constexpr ( std::is_base_of_v< SurfaceElementSubRegion, SUBREGION_TYPE > )
    {
      SurfaceElementAccumulationKernel kernel( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      SurfaceElementAccumulationKernel::launch< POLICY >( subRegion.size(), kernel );
    }
    else
    {
      GEOS_UNUSED_VAR( rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs );
      GEOS_ERROR( "Unsupported subregion type: " << typeid(SUBREGION_TYPE).name() );
    }
  }

};

} // namespace thermalSinglePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALACCUMULATIONKERNELS_HPP
