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
 * @file ThermalCompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/IsothermalCompositionalMultiphaseBaseKernels.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{


/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @class PhaseVolumeFractionKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the phase volume fractions
 */
template< integer NUM_COMP, integer NUM_PHASE >
class PhaseVolumeFractionKernel : public isothermalCompositionalMultiphaseBaseKernels::PhaseVolumeFractionKernel< NUM_COMP, NUM_PHASE >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PhaseVolumeFractionKernel< NUM_COMP, NUM_PHASE >;
  using Base::m_dPhaseDens;
  using Base::m_dPhaseFrac;
  using Base::m_dPhaseVolFrac;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  PhaseVolumeFractionKernel( ObjectManagerBase & subRegion,
                             constitutive::MultiFluidBase const & fluid )
    : Base( subRegion, fluid )
  {}

  /**
   * @brief Compute the phase volume fractions in an element
   * @param[in] ei the element index
   */
  GEOS_HOST_DEVICE
  real64 compute( localIndex const ei ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDens = m_dPhaseDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseFrac = m_dPhaseFrac[ei][0];

    arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dPhaseVolFrac[ei];

    // Call the base compute the compute the phase volume fraction
    return Base::compute( ei, [&] ( localIndex const ip,
                                    real64 const & phaseVolFrac,
                                    real64 const & phaseDensInv,
                                    real64 const & totalDensity )
    {
      // when this lambda is called, we are in the phase loop
      // for each phase ip, compute the derivative of phase volume fraction wrt temperature
      dPhaseVolFrac[ip][Deriv::dT] = (dPhaseFrac[ip][Deriv::dT] - phaseVolFrac * dPhaseDens[ip][Deriv::dT]) * phaseDensInv;
      dPhaseVolFrac[ip][Deriv::dT] *= totalDensity;
    } );
  }

};

/**
 * @class PhaseVolumeFractionKernelFactory
 */
class PhaseVolumeFractionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static real64
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   constitutive::MultiFluidBase const & fluid )
  {
    real64 maxDeltaPhaseVolFrac = 0.0;
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::
        internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        PhaseVolumeFractionKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        maxDeltaPhaseVolFrac = PhaseVolumeFractionKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    return maxDeltaPhaseVolFrac;
  }
};


/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of thermal accumulation and volume balance
 */
template< localIndex NUM_COMP, localIndex NUM_DOF >
class ElementBasedAssemblyKernel : public isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_numPhases;
  using Base::m_rankOffset;
  using Base::m_dofNumber;
  using Base::m_elemGhostRank;
  using Base::m_volume;
  using Base::m_porosity;
  using Base::m_dPoro_dPres;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_phaseVolFrac;
  using Base::m_dPhaseVolFrac;
  using Base::m_phaseDens;
  using Base::m_dPhaseDens;
  using Base::m_phaseCompFrac;
  using Base::m_dPhaseCompFrac;
  using Base::m_localMatrix;
  using Base::m_localRhs;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              constitutive::MultiFluidBase const & fluid,
                              constitutive::CoupledSolidBase const & solid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > const kernelFlags )
    : Base( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs, kernelFlags ),
    m_dPoro_dTemp( solid.getDporosity_dTemperature() ),
    m_phaseInternalEnergy( fluid.phaseInternalEnergy() ),
    m_dPhaseInternalEnergy( fluid.dPhaseInternalEnergy() ),
    m_rockInternalEnergy( solid.getInternalEnergy() ),
    m_dRockInternalEnergy_dTemp( solid.getDinternalEnergy_dTemperature() ),
    m_energy_n( subRegion.getField< fields::flow::energy_n >() )
  {}

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables()
      : Base::StackVariables()
    {}

    using Base::StackVariables::localRow;
    using Base::StackVariables::dofIndices;
    using Base::StackVariables::localResidual;
    using Base::StackVariables::localJacobian;

    // derivative of pore volume wrt temperature
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

    // derivative of pore volume wrt temperature
    stack.dPoreVolume_dTemp = m_volume[ei] * m_dPoro_dTemp[ei][0];

    // initialize the solid volume
    real64 const solidVolume = m_volume[ei] * ( 1.0 - m_porosity[ei][0] );
    real64 const dSolidVolume_dPres = -m_volume[ei] * m_dPoro_dPres[ei][0];
    real64 const dSolidVolume_dTemp = -stack.dPoreVolume_dTemp;

    // initialize the solid internal energy
    stack.solidEnergy = solidVolume * m_rockInternalEnergy[ei][0];
    stack.dSolidEnergy_dPres = dSolidVolume_dPres * m_rockInternalEnergy[ei][0];
    stack.dSolidEnergy_dTemp = solidVolume * m_dRockInternalEnergy_dTemp[ei][0]
                               + dSolidVolume_dTemp * m_rockInternalEnergy[ei][0];
  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    // start with old time step value
    stack.localResidual[numEqn-1] = -m_energy_n[ei];

    Base::computeAccumulation( ei, stack, [&] ( integer const ip,
                                                real64 const & phaseAmount,
                                                real64 const & dPhaseAmount_dP,
                                                real64 const (&dPhaseAmount_dC)[numComp] )
    {
      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[numComp]{};

      // construct the slices
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
      arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];
      arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseInternalEnergy = m_phaseInternalEnergy[ei][0];
      arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > dPhaseInternalEnergy = m_dPhaseInternalEnergy[ei][0];

      // Step 1: assemble the derivatives of the component mass balance equations with respect to temperature

      real64 const dPhaseAmount_dT = stack.dPoreVolume_dTemp * phaseVolFrac[ip] * phaseDens[ip]
                                     + stack.poreVolume * (dPhaseVolFrac[ip][Deriv::dT] * phaseDens[ip] + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dT] );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        stack.localJacobian[ic][numDof-1] += dPhaseAmount_dT * phaseCompFrac[ip][ic]
                                             + phaseAmount * dPhaseCompFrac[ip][ic][Deriv::dT];
      }

      // Step 2: assemble the phase-dependent part of the accumulation term of the energy equation

      real64 const phaseEnergy = phaseAmount * phaseInternalEnergy[ip];
      real64 const dPhaseEnergy_dP = dPhaseAmount_dP * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dP];
      real64 const dPhaseEnergy_dT = dPhaseAmount_dT * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dT];

      // local accumulation
      stack.localResidual[numEqn-1] += phaseEnergy;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn-1][0]        += dPhaseEnergy_dP;
      stack.localJacobian[numEqn-1][numDof-1] += dPhaseEnergy_dT;

      // derivatives w.r.t. component densities
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numEqn-1][jc + 1] += phaseInternalEnergy[ip] * dPhaseAmount_dC[jc]
                                                 + dPhaseInternalEnergy_dC[jc] * phaseAmount;
      }
    } );

    // Step 3: assemble the solid part of the accumulation term

    // local accumulation and derivatives w.r.t. pressure and temperature
    stack.localResidual[numEqn-1] += stack.solidEnergy;
    stack.localJacobian[numEqn-1][0] += stack.dSolidEnergy_dPres;
    stack.localJacobian[numEqn-1][numDof-1] += stack.dSolidEnergy_dTemp;

  }

  /**
   * @brief Compute the local volume balance contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    Base::computeVolumeBalance( ei, stack, [&] ( real64 const & oneMinusPhaseVolFraction )
    {
      GEOS_UNUSED_VAR( oneMinusPhaseVolFraction );

      arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

      for( integer ip = 0; ip < m_numPhases; ++ip )
      {
        stack.localJacobian[numEqn-2][numDof-1] -= dPhaseVolFrac[ip][Deriv::dT];
      }
    } );
  }

  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Step 1: assemble the component mass balance equations and volume balance equations
    Base::complete( ei, stack );

    // Step 2: assemble the energy equation
    m_localRhs[stack.localRow + numEqn-1] += stack.localResidual[numEqn-1];
    m_localMatrix.template addToRow< serialAtomic >( stack.localRow + numEqn-1,
                                                     stack.dofIndices,
                                                     stack.localJacobian[numEqn-1],
                                                     numDof );
  }

protected:

  /// View on derivative of porosity w.r.t temperature
  arrayView2d< real64 const > const m_dPoro_dTemp;

  /// Views on phase internal energy
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseInternalEnergy;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseInternalEnergy;

  /// Views on rock internal energy
  arrayView2d< real64 const > m_rockInternalEnergy;
  arrayView2d< real64 const > m_dRockInternalEnergy_dTemp;

  /// Views on energy
  arrayView1d< real64 const > m_energy_n;

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
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( localIndex const numComps,
                   localIndex const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&] ( auto NC )
    {
      localIndex constexpr NUM_COMP = NC();
      localIndex constexpr NUM_DOF = NC()+2;

      BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation );

      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, solid, localMatrix, localRhs, kernelFlags );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
      launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel );
    } );
  }

};

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & temp,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp[k], compFrac[k] );
      }
    } );
  }
};

/******************************** SolidInternalEnergyUpdateKernel ********************************/

struct SolidInternalEnergyUpdateKernel
{

  template< typename POLICY, typename SOLID_INTERNAL_ENERGY_WRAPPER >
  static void
  launch( localIndex const size,
          SOLID_INTERNAL_ENERGY_WRAPPER const & solidInternalEnergyWrapper,
          arrayView1d< real64 const > const & temp )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      solidInternalEnergyWrapper.update( k, temp[k] );
    } );
  }
};

/******************************** ScalingForSystemSolutionKernel ********************************/

/**
 * @class ScalingForSystemSolutionKernel
 * @brief Define the kernel for scaling the Newton update
 */
class ScalingForSystemSolutionKernel : public isothermalCompositionalMultiphaseBaseKernels::ScalingForSystemSolutionKernel
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::ScalingForSystemSolutionKernel;
  using Base::m_numComp;
  using Base::m_localSolution;

  /**
   * @brief Create a new kernel instance
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxRelativeTempChange the max allowed relative temperature change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompDensChange the max allowed comp density change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] temperature the temperature vector
   * @param[in] compDens the component density vector
   * @param[in] pressureScalingFactor the pressure local scaling factor
   * @param[in] compDensScalingFactor the component density local scaling factor
   * @param[in] temperatureFactor the temperature local scaling factor
   */
  ScalingForSystemSolutionKernel( real64 const maxRelativePresChange,
                                  real64 const maxAbsolutePresChange,
                                  real64 const maxRelativeTempChange,
                                  real64 const maxCompFracChange,
                                  real64 const maxRelativeCompDensChange,
                                  globalIndex const rankOffset,
                                  integer const numComp,
                                  string const dofKey,
                                  ElementSubRegionBase const & subRegion,
                                  arrayView1d< real64 const > const localSolution,
                                  arrayView1d< real64 const > const pressure,
                                  arrayView1d< real64 const > const temperature,
                                  arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                                  arrayView1d< real64 > pressureScalingFactor,
                                  arrayView1d< real64 > compDensScalingFactor,
                                  arrayView1d< real64 > temperatureScalingFactor )
    : Base( maxRelativePresChange,
            maxAbsolutePresChange,
            maxCompFracChange,
            maxRelativeCompDensChange,
            rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor ),
    m_maxRelativeTempChange( maxRelativeTempChange ),
    m_temperature( temperature ),
    m_temperatureScalingFactor( temperatureScalingFactor )
  {}

  /**
   * @brief Compute the local value
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                StackVariables & stack ) const
  {
    computeScalingFactor( ei, stack );
  }

  /**
   * @brief Compute the local value of the scaling factor
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeScalingFactor( localIndex const ei,
                             StackVariables & stack ) const
  {
    real64 constexpr eps = isothermalCompositionalMultiphaseBaseKernels::minDensForDivision;
    Base::computeScalingFactor( ei, stack, [&] ()
    {
      // compute the change in temperature
      real64 const temp = m_temperature[ei];
      real64 const absTempChange = LvArray::math::abs( m_localSolution[stack.localRow + m_numComp + 1] );
      if( stack.localMaxDeltaTemp < absTempChange )
      {
        stack.localMaxDeltaTemp = absTempChange;
      }

      m_temperatureScalingFactor[ei] = 1.0;

      if( temp > eps )
      {
        real64 const relativeTempChange = absTempChange / temp;
        if( relativeTempChange > m_maxRelativeTempChange )
        {
          real64 const tempScalingFactor = m_maxRelativeTempChange / relativeTempChange;
          m_temperatureScalingFactor[ei] = tempScalingFactor;
          if( stack.localMinVal > tempScalingFactor )
          {
            stack.localMinVal = tempScalingFactor;
          }
          if( stack.localMinTempScalingFactor > tempScalingFactor )
          {
            stack.localMinTempScalingFactor = tempScalingFactor;
          }
        }
      }
    } );
  }

protected:

  /// Max allowed changes in primary variables
  real64 const m_maxRelativeTempChange;

  /// View on the primary variables
  arrayView1d< real64 const > const m_temperature;

  /// View on the scaling factor
  arrayView1d< real64 > const m_temperatureScalingFactor;

};

/**
 * @class ScalingForSystemSolutionKernelFactory
 */
class ScalingForSystemSolutionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxRelativeTempChange the max allowed relative temperature change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompdensChange the max allowed relative component density change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   */
  template< typename POLICY >
  static ScalingForSystemSolutionKernel::StackVariables
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxAbsolutePresChange,
                   real64 const maxRelativeTempChange,
                   real64 const maxCompFracChange,
                   real64 const maxRelativeCompDensChange,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const temperature = subRegion.getField< fields::flow::temperature >();
    arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();
    arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
    arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
    arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();
    ScalingForSystemSolutionKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxRelativeTempChange,
                                           maxCompFracChange, maxRelativeCompDensChange,
                                           rankOffset, numComp, dofKey, subRegion, localSolution,
                                           pressure, temperature, compDens, pressureScalingFactor,
                                           temperatureScalingFactor, compDensScalingFactor );
    return thermalCompositionalMultiphaseBaseKernels::
             ScalingForSystemSolutionKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernel
 * @brief Define the kernel for checking the updated solution
 */
class SolutionCheckKernel : public isothermalCompositionalMultiphaseBaseKernels::SolutionCheckKernel
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::SolutionCheckKernel;
  using Base::m_numComp;
  using Base::m_localSolution;
  using Base::m_scalingFactor;

  static real64 constexpr minTemperature = constants::zeroDegreesCelsiusInKelvin;

  /**
   * @brief Create a new kernel instance
   * @param[in] allowCompDensChopping flag to allow the component density chopping
   * @param[in] scalingFactor the scaling factor
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] temperature the temperature vector
   * @param[in] compDens the component density vector
   */
  SolutionCheckKernel( integer const allowCompDensChopping,
                       integer const allowNegativePressure,
                       CompositionalMultiphaseFVM::ScalingType const scalingType,
                       real64 const scalingFactor,
                       globalIndex const rankOffset,
                       integer const numComp,
                       string const dofKey,
                       ElementSubRegionBase const & subRegion,
                       arrayView1d< real64 const > const localSolution,
                       arrayView1d< real64 const > const pressure,
                       arrayView1d< real64 const > const temperature,
                       arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                       arrayView1d< real64 > pressureScalingFactor,
                       arrayView1d< real64 > compDensScalingFactor,
                       arrayView1d< real64 > temperatureScalingFactor )
    : Base( allowCompDensChopping,
            allowNegativePressure,
            scalingType,
            scalingFactor,
            rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor ),
    m_temperature( temperature ),
    m_temperatureScalingFactor( temperatureScalingFactor )
  {}

  /**
   * @brief Compute the local value of the solution check
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeSolutionCheck( localIndex const ei,
                             StackVariables & stack ) const
  {
    Base::computeSolutionCheck( ei, stack, [&] ()
    {
      bool const localScaling = m_scalingType == CompositionalMultiphaseFVM::ScalingType::Local;
      // compute the change in temperature
      real64 const newTemp = m_temperature[ei] + (localScaling ? m_temperatureScalingFactor[ei] : m_scalingFactor * m_localSolution[stack.localRow + m_numComp + 1]);
      if( newTemp < minTemperature )
      {
        stack.localMinVal = 0;
      }
    } );
  }

protected:

  /// View on the primary variables
  arrayView1d< real64 const > const m_temperature;

  /// View on the scaling factor
  arrayView1d< real64 const > const m_temperatureScalingFactor;

};

/**
 * @class SolutionCheckKernelFactory
 */
class SolutionCheckKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxRelativeTempChange the max allowed relative temperature change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   */
  template< typename POLICY >
  static SolutionCheckKernel::StackVariables
  createAndLaunch( integer const allowCompDensChopping,
                   integer const allowNegativePressure,
                   CompositionalMultiphaseFVM::ScalingType const scalingType,
                   real64 const scalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    arrayView1d< real64 const > const pressure =
      subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const temperature =
      subRegion.getField< fields::flow::temperature >();
    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getField< fields::flow::globalCompDensity >();
    arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
    arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
    arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();
    SolutionCheckKernel kernel( allowCompDensChopping, allowNegativePressure, scalingType, scalingFactor,
                                rankOffset, numComp, dofKey, subRegion, localSolution,
                                pressure, temperature, compDens, pressureScalingFactor, temperatureScalingFactor, compDensScalingFactor );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};


/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 3 >
{
public:

  using Base = ResidualNormKernelBase< 3 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numComponents,
                      integer const numPhases,
                      ElementSubRegionBase const & subRegion,
                      constitutive::MultiFluidBase const & fluid,
                      constitutive::CoupledSolidBase const & solid,
                      constitutive::SolidInternalEnergy const & solidInternalEnergy,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numComponents( numComponents ),
    m_numPhases( numPhases ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_phaseVolFrac_n( subRegion.getField< fields::flow::phaseVolumeFraction_n >() ),
    m_totalDens_n( fluid.totalDensity_n() ),
    m_phaseDens_n( fluid.phaseDensity_n() ),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n() ),
    m_solidInternalEnergy_n( solidInternalEnergy.getInternalEnergy_n() )
  {}

  GEOS_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const ei,
                                     real64 & massNormalizer,
                                     real64 & energyNormalizer ) const
  {
    massNormalizer = LvArray::math::max( m_minNormalizer, m_totalDens_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );
    real64 const poreVolume = m_porosity_n[ei][0] * m_volume[ei];
    energyNormalizer = m_solidInternalEnergy_n[ei][0] * ( 1.0 - m_porosity_n[ei][0] ) * m_volume[ei];
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      energyNormalizer += m_phaseInternalEnergy_n[ei][0][ip] * m_phaseDens_n[ei][0][ip] * m_phaseVolFrac_n[ei][ip] * poreVolume;
    }
    // warning: internal energy can be negative
    energyNormalizer = LvArray::math::max( m_minNormalizer, LvArray::math::abs( energyNormalizer ) );
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );
    real64 const volumeNormalizer = LvArray::math::max( m_minNormalizer, m_porosity_n[ei][0] * m_volume[ei] );

    // step 1: mass residual

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / massNormalizer;
      if( valMass > stack.localValue[0] )
      {
        stack.localValue[0] = valMass;
      }
    }

    // step 2: volume residual

    real64 const valVol = LvArray::math::abs( m_localResidual[stack.localRow + m_numComponents] ) / volumeNormalizer;
    if( valVol > stack.localValue[1] )
    {
      stack.localValue[1] = valVol;
    }

    // step 3: energy residual

    real64 const valEnergy = LvArray::math::abs( m_localResidual[stack.localRow + m_numComponents + 1] ) / energyNormalizer;
    if( valEnergy > stack.localValue[2] )
    {
      stack.localValue[2] = valEnergy;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    // note: for the L2 norm, we bundle the volume and mass residuals/normalizers
    real64 massNormalizer = 0.0, energyNormalizer = 0.0;
    computeMassEnergyNormalizers( ei, massNormalizer, energyNormalizer );

    // step 1: mass residual

    for( integer idof = 0; idof < m_numComponents; ++idof )
    {
      stack.localValue[0] += m_localResidual[stack.localRow + idof] * m_localResidual[stack.localRow + idof];
      stack.localNormalizer[0] += massNormalizer;
    }

    // step 2: volume residual

    real64 const valVol = m_localResidual[stack.localRow + m_numComponents] * m_totalDens_n[ei][0]; // we need a mass here, hence the
                                                                                                    // multiplication
    stack.localValue[1] += valVol * valVol;
    stack.localNormalizer[1] += massNormalizer;

    // step 3: energy residual

    stack.localValue[2] += m_localResidual[stack.localRow + m_numComponents + 1] * m_localResidual[stack.localRow + m_numComponents + 1];
    stack.localNormalizer[2] += energyNormalizer;
  }

protected:

  /// Number of fluid components
  integer const m_numComponents;

  /// Number of fluid phases
  integer const m_numPhases;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  arrayView2d< real64 const > const m_porosity_n;

  /// View on phase properties at the previous converged time step
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac_n;
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const m_totalDens_n;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens_n;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseInternalEnergy_n;

  /// View on solid properties at the previous converged time step
  arrayView2d< real64 const > const m_solidInternalEnergy_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[in] solidInternalEnergy the solid internal energy model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( solverBaseKernels::NormType const normType,
                   integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   ElementSubRegionBase const & subRegion,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   constitutive::SolidInternalEnergy const & solidInternalEnergy,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[3],
                   real64 (& residualNormalizer)[3] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               numComps, numPhases, subRegion, fluid, solid, solidInternalEnergy, minNormalizer );
    if( normType == solverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }

  }
};


} // namespace thermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_THERMALCOMPOSITIONALMULTIPHASEBASEKERNELS_HPP
