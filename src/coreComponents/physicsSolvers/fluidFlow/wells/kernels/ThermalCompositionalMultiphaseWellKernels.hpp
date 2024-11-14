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
 * @file ThermalCompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP

#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseWellKernels
{

using namespace constitutive;

/******************************** TotalMassDensityKernel ****************************/

/**
 * @class TotalMassDensityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the total mass density
 */
template< integer NUM_COMP, integer NUM_PHASE >
class TotalMassDensityKernel : public compositionalMultiphaseWellKernels::TotalMassDensityKernel< NUM_COMP, NUM_PHASE >
{
public:
  using Base = compositionalMultiphaseWellKernels::TotalMassDensityKernel< NUM_COMP, NUM_PHASE >;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_dPhaseMassDens;
  using Base::m_dPhaseVolFrac;
  using Base::m_dTotalMassDens;
  using Base::m_phaseMassDens;
  using Base::m_phaseVolFrac;
  using Base::m_totalMassDens;
  using Base::numComp;
  using Base::numPhase;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  TotalMassDensityKernel( ObjectManagerBase & subRegion,
                          MultiFluidBase const & fluid )
    : Base( subRegion, fluid )
  {}

  /**
   * @brief Compute the total mass density in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] totalMassDensityKernelOp the function used to customize the kernel
   */
  GEOS_HOST_DEVICE inline void compute( localIndex const ei ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseMassDens = m_phaseMassDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseMassDens = m_dPhaseMassDens[ei][0];

    real64 & dTotalMassDens_dT = m_dTotalMassDens[ei][Deriv::dT];

    // Call the base compute the compute the total mass density and derivatives
    return Base::compute( ei, [&]( localIndex const ip )
    {
      dTotalMassDens_dT += dPhaseVolFrac[ip][Deriv::dT] * phaseMassDens[ip] + phaseVolFrac[ip] * dPhaseMassDens[ip][Deriv::dT];
    } );
  }

protected:
  // outputs
  arrayView1d< real64 > m_dTotalMassDens_dTemp;
};

/**
 * @class TotalMassDensityKernelFactory
 */
class TotalMassDensityKernelFactory
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
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&]( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&]( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
template< localIndex NUM_COMP >
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 2 >
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;


  using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NUM_COMP, 1 >;

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 2 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const targetPhaseIndex,
                      WellElementSubRegion const & subRegion,
                      MultiFluidBase const & fluid,
                      WellControls const & wellControls,
                      real64 const timeAtEndOfStep,
                      real64 const dt,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numPhases( fluid.numFluidPhases()),
    m_targetPhaseIndex( targetPhaseIndex ),
    m_dt( dt ),
    m_isLocallyOwned( subRegion.isLocallyOwned() ),
    m_iwelemControl( subRegion.getTopWellElementIndex() ),
    m_isProducer( wellControls.isProducer() ),
    m_currentControl( wellControls.getControl() ),
    m_targetBHP( wellControls.getTargetBHP( timeAtEndOfStep ) ),
    m_targetTotalRate( wellControls.getTargetTotalRate( timeAtEndOfStep ) ),
    m_targetPhaseRate( wellControls.getTargetPhaseRate( timeAtEndOfStep ) ),
    m_targetMassRate( wellControls.getTargetMassRate( timeAtEndOfStep ) ),
    m_volume( subRegion.getElementVolume() ),
    m_phaseDens_n( fluid.phaseDensity_n() ),
    m_totalDens_n( fluid.totalDensity_n() ),
    m_phaseVolFraction_n( subRegion.getField< fields::well::phaseVolumeFraction_n >()),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n() )
  {}


  GEOS_HOST_DEVICE
  void computeMassEnergyNormalizers( localIndex const iwelem,
                                     real64 & massNormalizer,
                                     real64 & energyNormalizer ) const
  {
    massNormalizer = LvArray::math::max( m_minNormalizer, m_totalDens_n[iwelem][0] *  m_volume[iwelem] );

    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      energyNormalizer += m_phaseInternalEnergy_n[iwelem][0][ip] * m_phaseDens_n[iwelem][0][ip] * m_phaseVolFraction_n[iwelem][ip] *   m_volume[iwelem];
    }
    // warning: internal energy can be negative
    energyNormalizer = LvArray::math::max( m_minNormalizer, LvArray::math::abs( energyNormalizer ) );
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const iwelem,
                            LinfStackVariables & stack ) const override
  {
    real64 normalizer = 0.0;
    for( integer idof = 0; idof < WJ_ROFFSET::nEqn; ++idof )
    {

      // Step 1: compute a normalizer for the control or pressure equation

      // for the control equation, we distinguish two cases
      if( idof == WJ_ROFFSET::CONTROL )
      {

        // for the top well element, normalize using the current control
        if( m_isLocallyOwned && iwelem == m_iwelemControl )
        {
          if( m_currentControl == WellControls::Control::BHP )
          {
            // the residual entry is in pressure units
            normalizer = m_targetBHP;
          }
          else if( m_currentControl == WellControls::Control::TOTALVOLRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetTotalRate ), m_minNormalizer );
          }
          else if( m_currentControl == WellControls::Control::PHASEVOLRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetPhaseRate ), m_minNormalizer );
          }
          else if( m_currentControl == WellControls::Control::MASSRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetMassRate ), m_minNormalizer );
          }
        }
        // for the pressure difference equation, always normalize by the BHP
        else
        {
          normalizer = m_targetBHP;
        }
      }
      // Step 2: compute a normalizer for the mass balance equations
      else if( idof >= WJ_ROFFSET::MASSBAL && idof < WJ_ROFFSET::MASSBAL + numComp )
      {
        if( m_isProducer ) // only PHASEVOLRATE is supported for now
        {
          // the residual is in mass units
          normalizer = m_dt * LvArray::math::abs( m_targetPhaseRate ) * m_phaseDens_n[iwelem][0][m_targetPhaseIndex];
        }
        else // Type::INJECTOR, only TOTALVOLRATE is supported for now
        {
          if( m_currentControl == WellControls::Control::MASSRATE )
          {
            normalizer = m_dt * LvArray::math::abs( m_targetMassRate );
          }
          else
          {
            // the residual is in mass units
            normalizer = m_dt * LvArray::math::abs( m_targetTotalRate ) * m_totalDens_n[iwelem][0];
          }

        }

        // to make sure that everything still works well if the rate is zero, we add this check
        normalizer = LvArray::math::max( normalizer, m_volume[iwelem] * m_totalDens_n[iwelem][0] );
      }
      // Step 3: compute a normalizer for the volume balance equations
      else if( idof == WJ_ROFFSET::VOLBAL )
      {
        if( m_isProducer ) // only PHASEVOLRATE is supported for now
        {
          // the residual is in volume units
          normalizer = m_dt * LvArray::math::abs( m_targetPhaseRate );
        }
        else // Type::INJECTOR, only TOTALVOLRATE is supported for now
        {
          if( m_currentControl == WellControls::Control::MASSRATE )
          {
            normalizer = m_dt * LvArray::math::abs( m_targetMassRate/  m_totalDens_n[iwelem][0] );
          }
          else
          {
            normalizer = m_dt * LvArray::math::abs( m_targetTotalRate );
          }

        }
        // to make sure that everything still works well if the rate is zero, we add this check
        normalizer = LvArray::math::max( normalizer, m_volume[iwelem] );
      }
      // step 3: energy residual
      if( idof == WJ_ROFFSET::ENERGYBAL )
      {
        real64 massNormalizer = 0.0, energyNormalizer = 0.0;
        computeMassEnergyNormalizers( iwelem, massNormalizer, energyNormalizer );
        real64 const valEnergy = LvArray::math::abs( m_localResidual[stack.localRow + WJ_ROFFSET::ENERGYBAL] ) / energyNormalizer;
        if( valEnergy > stack.localValue[1] )
        {
          stack.localValue[1] = valEnergy;
        }

      }
      else
      {
        normalizer = LvArray::math::max( m_minNormalizer, normalizer );
        // Step 4: compute the contribution to the residual
        real64 const val = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / normalizer;
        if( val > stack.localValue[0] )
        {
          stack.localValue[0] = val;
        }
      }
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const iwelem,
                          L2StackVariables & stack ) const override
  {
    GEOS_UNUSED_VAR( iwelem, stack );
    GEOS_ERROR( "The L2 norm is not implemented for CompositionalMultiphaseWell" );
  }


protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Index of the target phase
  integer const m_targetPhaseIndex;

  /// Time step size
  real64 const m_dt;

  /// Flag indicating whether the well is locally owned or not
  bool const m_isLocallyOwned;

  /// Index of the element where the control is enforced
  localIndex const m_iwelemControl;

  /// Flag indicating whether the well is a producer or an injector
  bool const m_isProducer;

  /// Controls
  WellControls::Control const m_currentControl;
  real64 const m_targetBHP;
  real64 const m_targetTotalRate;
  real64 const m_targetPhaseRate;
  real64 const m_targetMassRate;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on phase/total density at the previous converged time step
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDens_n;
  arrayView2d< real64 const, multifluid::USD_FLUID > const m_totalDens_n;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFraction_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseInternalEnergy_n;

};

/*
 *@class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp number of fluid components
   * @param[in] numDof number of dofs per well element
   * @param[in] targetPhaseIndex the index of the target phase (for phase volume control)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the well element subregion
   * @param[in] fluid the fluid model
   * @param[in] wellControls the controls
   * @param[in] timeAtEndOfStep the time at the end of the step (time_n + dt)
   * @param[in] dt the time step size
   * @param[out] residualNorm the residual norm on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const targetPhaseIndex,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   WellElementSubRegion const & subRegion,
                   MultiFluidBase const & fluid,
                   WellControls const & wellControls,
                   real64 const timeAtEndOfStep,
                   real64 const dt,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[2] )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch ( numComp, [&]( auto NC )
    {

      integer constexpr NUM_COMP = NC();
      using kernelType = ResidualNormKernel< NUM_COMP >;
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      kernelType kernel( rankOffset, localResidual, dofNumber, ghostRank,
                         targetPhaseIndex, subRegion, fluid, wellControls, timeAtEndOfStep, dt, minNormalizer );
      kernelType::template launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    } );
  }

};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam IS_THERMAL thermal flag
 * @brief Define the interface for the assembly kernel in charge of thermal accumulation and volume balance
 */
template< localIndex NUM_COMP >
class ElementBasedAssemblyKernel : public compositionalMultiphaseWellKernels::ElementBasedAssemblyKernel< NUM_COMP, 1 >
{
public:
  using Base = compositionalMultiphaseWellKernels::ElementBasedAssemblyKernel< NUM_COMP, 1 >;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_dofNumber;
  using Base::m_dPhaseCompFrac;
  using Base::m_dPhaseDens;
  using Base::m_dPhaseVolFrac;
  using Base::m_dPoro_dPres;
  using Base::m_elemGhostRank;
  using Base::m_localMatrix;
  using Base::m_localRhs;
  using Base::m_numPhases;
  using Base::m_phaseCompFrac;
  using Base::m_phaseCompFrac_n;
  using Base::m_phaseDens;
  using Base::m_phaseDens_n;
  using Base::m_phaseVolFrac;
  using Base::m_phaseVolFrac_n;
  using Base::m_porosity;
  using Base::m_porosity_n;
  using Base::m_rankOffset;
  using Base::m_volume;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;

  using FLUID_PROP_COFFSET = multifluid::DerivativeOffsetC< NUM_COMP, 1 >;


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
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              integer const isProducer,
                              globalIndex const rankOffset,
                              string const dofKey,
                              WellElementSubRegion const & subRegion,
                              MultiFluidBase const & fluid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > const kernelFlags )
    : Base( numPhases, isProducer, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs, kernelFlags ),
    m_phaseInternalEnergy_n( fluid.phaseInternalEnergy_n()),
    m_phaseInternalEnergy( fluid.phaseInternalEnergy()),
    m_dPhaseInternalEnergy( fluid.dPhaseInternalEnergy())
  {}

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
    using Deriv = multifluid::DerivativeOffset;

    Base::computeAccumulation( ei, stack, [&]( integer const ip
                                               , real64 const & phaseAmount
                                               , real64 const & phaseAmount_n
                                               , real64 const (&dPhaseAmount)[FLUID_PROP_COFFSET::nDer] )
    {
      // We are in the loop over phases, ip provides the current phase index.
      // We have to do two things:
      //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
      //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

      real64 dPhaseInternalEnergy_dC[numComp]{};

      // construct the slices
      arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy_n = m_phaseInternalEnergy_n[ei][0];
      arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseInternalEnergy = m_phaseInternalEnergy[ei][0];
      arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseInternalEnergy = m_dPhaseInternalEnergy[ei][0];

      // Step 1: assemble the phase-dependent part of the accumulation term of the energy equation

      real64 const phaseEnergy = phaseAmount * phaseInternalEnergy[ip];
      real64 const phaseEnergy_n =  phaseAmount_n * phaseInternalEnergy_n[ip];
      real64 const dPhaseEnergy_dP = dPhaseAmount[FLUID_PROP_COFFSET::dP] * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dP];
      real64 const dPhaseEnergy_dT = dPhaseAmount[FLUID_PROP_COFFSET::dT] * phaseInternalEnergy[ip]
                                     + phaseAmount * dPhaseInternalEnergy[ip][Deriv::dT];
      // local accumulation
      stack.localResidual[numEqn-1] += phaseEnergy - phaseEnergy_n;

      // derivatives w.r.t. pressure and temperature
      stack.localJacobian[numEqn-1][0]        += dPhaseEnergy_dP;
      stack.localJacobian[numEqn-1][numDof-1] += dPhaseEnergy_dT;

      // derivatives w.r.t. component densities
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numEqn-1][jc + 1] += phaseInternalEnergy[ip] * dPhaseAmount[FLUID_PROP_COFFSET::dC+jc]
                                                 + dPhaseInternalEnergy_dC[jc] * phaseAmount;
      }
    } );


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
    Base::computeVolumeBalance( ei, stack );

  }

  GEOS_HOST_DEVICE
  void complete( localIndex const ei,
                 StackVariables & stack ) const
  {
    // Assemble the component mass balance equations and volume balance equations
    // Energy balance equation updates to solver matrices included in Base class
    Base::complete( ei, stack );

  }

protected:
  /// View on derivative of porosity w.r.t temperature
  arrayView2d< real64 const > const m_dPoro_dTemp;

  /// Views on phase internal energy
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseInternalEnergy;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseInternalEnergy;

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
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( localIndex const numComps,
                   localIndex const numPhases,
                   integer const isProducer,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   WellElementSubRegion const & subRegion,
                   MultiFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      localIndex constexpr NUM_COMP = NC();


      BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );

      ElementBasedAssemblyKernel< NUM_COMP >
      kernel( numPhases, isProducer, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs, kernelFlags );
      ElementBasedAssemblyKernel< NUM_COMP >::template
      launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP > >( subRegion.size(), kernel );
    } );
  }
};
/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NC >
class FaceBasedAssemblyKernel : public compositionalMultiphaseWellKernels::FaceBasedAssemblyKernel< NC, 1 >
{
public:
  static constexpr integer IS_THERMAL = 1;
  using Base  = compositionalMultiphaseWellKernels::FaceBasedAssemblyKernel< NC, IS_THERMAL >;

  // Well jacobian column and row indicies
  using WJ_COFFSET = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
  using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;

  using CP_Deriv = multifluid::DerivativeOffsetC< NC, IS_THERMAL >;

  using TAG = compositionalMultiphaseWellKernels::ElemTag;


  using Base::m_isProducer;
  using Base::m_dt;
  using Base::m_localRhs;
  using Base::m_localMatrix;
  using Base::m_rankOffset;
  using Base::maxNumElems;
  using Base::maxStencilSize;
  using Base::m_useTotalMassEquation;

  /// Compile time value for the number of components
  static constexpr integer numComp = NC;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = WJ_COFFSET::nDer;

/// Compile time value for the number of equations except volume and momentum
  static constexpr integer numEqn = WJ_ROFFSET::nEqn - 2;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  FaceBasedAssemblyKernel( real64 const dt,
                           globalIndex const rankOffset,
                           string const wellDofKey,
                           WellControls const & wellControls,
                           WellElementSubRegion const & subRegion,
                           MultiFluidBase const & fluid,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags )
    : Base( dt
            , rankOffset
            , wellDofKey
            , wellControls
            , subRegion
            , localMatrix
            , localRhs
            , kernelFlags ),
    m_numPhases ( fluid.numFluidPhases()),
    m_globalWellElementIndex( subRegion.getGlobalWellElementIndex() ),
    m_phaseFraction( fluid.phaseFraction()),
    m_dPhaseFraction( fluid.dPhaseFraction()),
    m_phaseEnthalpy( fluid.phaseEnthalpy()),
    m_dPhaseEnthalpy( fluid.dPhaseEnthalpy())
  { }

  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( localIndex const size )
      : Base::StackVariables( size )
    {}

    /// Storage for the face local residual vector (energy equation)
    stackArray1d< real64, maxNumElems > localEnergyFlux;
    /// Storage for the face local energy Jacobian matrix dC dP dT
    stackArray2d< real64, maxNumElems * maxStencilSize * CP_Deriv::nDer > localEnergyFluxJacobian;
    /// Storage for the face local Jacobian matrix dQ only
    stackArray2d< real64, maxNumElems * maxStencilSize > localEnergyFluxJacobian_dQ;
  };

  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const iwelem, StackVariables & stack ) const
  {
    Base::setup ( iwelem, stack );
    stack.localEnergyFlux.resize( stack.numConnectedElems );
    stack.localEnergyFluxJacobian.resize( stack.numConnectedElems, stack.stencilSize * numDof );
    stack.localEnergyFluxJacobian_dQ.resize( stack.numConnectedElems, 1 );
    for( integer i=0; i<stack.numConnectedElems; i++ )
    {
      stack.localEnergyFlux[i]=0.0;
      stack.localEnergyFluxJacobian_dQ[i][0]=0.0;
      for( integer j=0; j< stack.stencilSize * numDof; j++ )
        stack.localEnergyFluxJacobian[i][j] = 0.0;
    }
  }

  GEOS_HOST_DEVICE
  inline
  void complete( localIndex const iwelem, StackVariables & stack ) const
  {
    Base::complete ( iwelem, stack );

    using namespace compositionalMultiphaseUtilities;
    if( stack.numConnectedElems ==1 )
    {
      // Setup Jacobian global row indicies for energy equation
      globalIndex oneSidedEqnRowIndices  = stack.offsetUp + WJ_ROFFSET::ENERGYBAL - m_rankOffset;

      if( oneSidedEqnRowIndices  >= 0 && oneSidedEqnRowIndices < m_localMatrix.numRows() )
      {

        if( !m_isProducer  && m_globalWellElementIndex[iwelem] == 0 )
        {
          // For top segment energy balance eqn replaced with  T(n+1) - T = 0
          // No other energy balance derivatives
          // Assumption is global index == 0 is top segment with fixed temp BC
          for( integer i=0; i< CP_Deriv::nDer; i++ )
          {
            stack.localEnergyFluxJacobian[0][i] = 0.0;
          }
          stack.localEnergyFluxJacobian_dQ[0][0]=0;
          stack.localEnergyFlux[0]=0;
        }


        // Setup Jacobian global col indicies  ( Mapping from local jac order to well jac order)
        globalIndex oneSidedDofColIndices_dRate =   stack.offsetCurrent + WJ_COFFSET::dQ;
        globalIndex oneSidedDofColIndices_dPresCompTempUp[CP_Deriv::nDer]{};

        int ioff=0;
        oneSidedDofColIndices_dPresCompTempUp[ioff++] = stack.offsetUp + WJ_COFFSET::dP;
        oneSidedDofColIndices_dPresCompTempUp[ioff++] = stack.offsetUp + WJ_COFFSET::dT;
        for( integer jdof = 0; jdof < NC; ++jdof )
        {
          oneSidedDofColIndices_dPresCompTempUp[ioff++] = stack.offsetUp + WJ_COFFSET::dC+ jdof;
        }

        m_localMatrix.template addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndices,
                                                                 &oneSidedDofColIndices_dRate,
                                                                 stack.localEnergyFluxJacobian_dQ[0],
                                                                 1 );
        m_localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( oneSidedEqnRowIndices,
                                                                                     oneSidedDofColIndices_dPresCompTempUp,
                                                                                     stack.localEnergyFluxJacobian[0],
                                                                                     CP_Deriv::nDer );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[oneSidedEqnRowIndices], stack.localEnergyFlux[0] );
      }
    }
    else   // if ( stack.numConnectedElems == 2 )
    {
      globalIndex row_current = stack.offsetCurrent + WJ_ROFFSET::ENERGYBAL   - m_rankOffset;
      globalIndex row_next = stack.offsetNext + WJ_ROFFSET::ENERGYBAL   - m_rankOffset;

      if( !m_isProducer )
      {
        if( row_next >= 0 && row_next <  m_localMatrix.numRows() )
        {
          if( m_globalWellElementIndex[stack.iwelemNext] == 0 )
          {
            for( integer i=0; i<CP_Deriv::nDer; i++ )
              stack.localEnergyFluxJacobian[TAG::NEXT][i] = 0.0;
            stack.localEnergyFluxJacobian_dQ[TAG::NEXT][0] =0;
            stack.localEnergyFlux[TAG::NEXT] =0;
          }

        }
      }
      // Setup Jacobian global row indicies
      // equations for COMPONENT  + ENERGY balances
      globalIndex eqnRowIndices[2]{};

      // energy balance equations
      eqnRowIndices[TAG::CURRENT ] = row_current;
      eqnRowIndices[TAG::NEXT ]    = row_next;


      // Setup Jacobian global col indicies  ( Mapping from local jac order to well jac order)
      globalIndex dofColIndices[CP_Deriv::nDer]{};
      globalIndex dofColIndices_dRate = stack.offsetCurrent   + WJ_COFFSET::dQ;

      int ioff=0;
      // Indice storage order reflects local jac col storage order CP::Deriv order P T DENS
      // well jacobian order is P DENS Q T
      dofColIndices[ioff++] = stack.offsetUp + WJ_COFFSET::dP;

      if constexpr ( IS_THERMAL )
      {
        dofColIndices[ioff++] = stack.offsetUp + WJ_COFFSET::dT;
      }
      for( integer jdof = 0; jdof < NC; ++jdof )
      {
        dofColIndices[ioff++] = stack.offsetUp + WJ_COFFSET::dC+ jdof;
      }
      // Note this updates diag and offdiag
      for( integer i = 0; i < 2; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
        {
          m_localMatrix.template addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                   &dofColIndices_dRate,
                                                                   stack.localEnergyFluxJacobian_dQ[i],
                                                                   1 );
          m_localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                                       dofColIndices,
                                                                                       stack.localEnergyFluxJacobian[i],
                                                                                       CP_Deriv::nDer );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[eqnRowIndices[i]], stack.localEnergyFlux[i] );
        }
      }
    }

  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] ie the element index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */

  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iwelem, StackVariables & stack ) const
  {
    Base::computeFlux ( iwelem, stack, [&] ( localIndex const & iwelemNext
                                             , localIndex const & iwelemUp
                                             , real64 const & currentConnRate
                                             , real64 const (&dCompFrac_dCompDens)[NC][NC] )
    {

      if( iwelemNext < 0 &&  !m_isProducer )  // exit connection, injector
      {
        real64 eflux=0;
        real64 eflux_dq=0;
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          eflux    += m_phaseEnthalpy[iwelemUp][0][ip]* m_phaseFraction[iwelemUp][0][ip];
          eflux_dq += m_phaseEnthalpy[iwelemUp][0][ip] * m_phaseFraction[iwelemUp][0][ip];

          stack.localEnergyFluxJacobian[0] [CP_Deriv::dP] += m_phaseEnthalpy[iwelemUp][0][ip]*m_dPhaseFraction[iwelemUp][0][ip][CP_Deriv::dP]
                                                             +  m_dPhaseEnthalpy[iwelemUp][0][ip][CP_Deriv::dP]*m_phaseFraction[iwelemUp][0][ip];
          stack.localEnergyFluxJacobian[0] [CP_Deriv::dT] += m_phaseEnthalpy[iwelemUp][0][ip]*m_dPhaseFraction[iwelemUp][0][ip][CP_Deriv::dT]
                                                             +  m_dPhaseEnthalpy[iwelemUp][0][ip][CP_Deriv::dT]*m_phaseFraction[iwelemUp][0][ip];

          real64 dProp1_dC[numComp]{};
          applyChainRule( numComp, dCompFrac_dCompDens, m_dPhaseEnthalpy[iwelemUp][0][ip], dProp1_dC, CP_Deriv::dC );
          real64 dProp2_dC[numComp]{};
          applyChainRule( numComp, dCompFrac_dCompDens, m_dPhaseFraction[iwelemUp][0][ip], dProp2_dC, CP_Deriv::dC );
          for( integer dof=0; dof < numComp; dof++ )
          {
            stack.localEnergyFluxJacobian[0] [CP_Deriv::dC+dof] += m_phaseEnthalpy[iwelemUp][0][ip]*dProp2_dC[dof]
                                                                   + dProp1_dC[dof]*m_phaseFraction[iwelemUp][0][ip];
          }
        }
        for( integer dof=0; dof < CP_Deriv::nDer; dof++ )
        {
          stack.localEnergyFluxJacobian[0] [dof] *= -m_dt*currentConnRate;
        }
        // Energy equation
        stack.localEnergyFlux[0]   =  -m_dt * eflux * currentConnRate;
        stack.localEnergyFluxJacobian_dQ[0][0]  = -m_dt * eflux_dq;
      }
      else if( (  iwelemNext < 0 && m_isProducer ) || currentConnRate < 0 )    // exit connection, producer
      {
        real64 eflux=0;
        real64 eflux_dq=0;
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          eflux    += m_phaseEnthalpy[iwelemUp][0][ip]* m_phaseFraction[iwelemUp][0][ip];
          eflux_dq += m_phaseEnthalpy[iwelemUp][0][ip] * m_phaseFraction[iwelemUp][0][ip];
          stack.localEnergyFluxJacobian[0] [CP_Deriv::dP] += m_phaseEnthalpy[iwelemUp][0][ip]*m_dPhaseFraction[iwelemUp][0][ip][CP_Deriv::dP]
                                                             +  m_dPhaseEnthalpy[iwelemUp][0][ip][CP_Deriv::dP]*m_phaseFraction[iwelemUp][0][ip];
          stack.localEnergyFluxJacobian[0] [CP_Deriv::dT] += m_phaseEnthalpy[iwelemUp][0][ip]*m_dPhaseFraction[iwelemUp][0][ip][CP_Deriv::dT]
                                                             +  m_dPhaseEnthalpy[iwelemUp][0][ip][CP_Deriv::dT]*m_phaseFraction[iwelemUp][0][ip];

          real64 dProp1_dC[numComp]{};
          applyChainRule( numComp, dCompFrac_dCompDens, m_dPhaseEnthalpy[iwelemUp][0][ip], dProp1_dC, CP_Deriv::dC );
          real64 dProp2_dC[numComp]{};
          applyChainRule( numComp, dCompFrac_dCompDens, m_dPhaseFraction[iwelemUp][0][ip], dProp2_dC, CP_Deriv::dC );
          for( integer dof=0; dof < numComp; dof++ )
          {
            stack.localEnergyFluxJacobian[0] [CP_Deriv::dC+dof] += m_phaseEnthalpy[iwelemUp][0][ip]*dProp2_dC[dof]
                                                                   + dProp1_dC[dof]*m_phaseFraction[iwelemUp][0][ip];
          }

        }

        for( integer dof=0; dof < CP_Deriv::nDer; dof++ )
        {
          stack.localEnergyFluxJacobian[0][dof] *= -m_dt*currentConnRate;
        }
        stack.localEnergyFlux[0]   =  -m_dt * eflux * currentConnRate;
        stack.localEnergyFluxJacobian_dQ[0][0]  = -m_dt*eflux_dq;
      }
      else
      {
        real64 eflux=0;
        real64 eflux_dq=0;
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          eflux    += m_phaseEnthalpy[iwelemUp][0][ip]* m_phaseFraction[iwelemUp][0][ip];
          eflux_dq += m_phaseEnthalpy[iwelemUp][0][ip] * m_phaseFraction[iwelemUp][0][ip];

          real64 dprop_dp = m_phaseEnthalpy[iwelemUp][0][ip]*m_dPhaseFraction[iwelemUp][0][ip][CP_Deriv::dP]
                            +  m_dPhaseEnthalpy[iwelemUp][0][ip][CP_Deriv::dP]*m_phaseFraction[iwelemUp][0][ip];
          real64 dprop_dt = m_phaseEnthalpy[iwelemUp][0][ip]*m_dPhaseFraction[iwelemUp][0][ip][CP_Deriv::dT]
                            +  m_dPhaseEnthalpy[iwelemUp][0][ip][CP_Deriv::dT]*m_phaseFraction[iwelemUp][0][ip];

          stack.localEnergyFluxJacobian[TAG::NEXT ] [CP_Deriv::dP] += dprop_dp;
          stack.localEnergyFluxJacobian[TAG::NEXT] [CP_Deriv::dT] += dprop_dt;

          stack.localEnergyFluxJacobian[TAG::CURRENT ] [CP_Deriv::dP] += dprop_dp;
          stack.localEnergyFluxJacobian[TAG::CURRENT] [CP_Deriv::dT] += dprop_dt;

          real64 dPE_dC[numComp]{};
          applyChainRule( numComp, dCompFrac_dCompDens, m_dPhaseEnthalpy[iwelemUp][0][ip], dPE_dC, CP_Deriv::dC );
          real64 dPF_dC[numComp]{};
          applyChainRule( numComp, dCompFrac_dCompDens, m_dPhaseFraction[iwelemUp][0][ip], dPF_dC, CP_Deriv::dC );

          for( integer dof=0; dof < numComp; dof++ )
          {
            stack.localEnergyFluxJacobian[TAG::NEXT   ][CP_Deriv::dC+dof] +=  m_phaseEnthalpy[iwelemUp][0][ip]*dPF_dC[dof]
                                                                             +dPE_dC[dof]*m_phaseFraction[iwelemUp][0][ip];
            stack.localEnergyFluxJacobian[TAG::CURRENT ][CP_Deriv::dC+dof] +=  m_phaseEnthalpy[iwelemUp][0][ip]*dPF_dC[dof]
                                                                              +dPE_dC[dof]*m_phaseFraction[iwelemUp][0][ip];
          }
        }
        stack.localEnergyFlux[TAG::NEXT   ]   =  m_dt * eflux * currentConnRate;
        stack.localEnergyFlux[TAG::CURRENT  ] = -m_dt * eflux * currentConnRate;
        stack.localEnergyFluxJacobian_dQ [TAG::NEXT   ][0] =  m_dt * eflux_dq;
        stack.localEnergyFluxJacobian_dQ [TAG::CURRENT][0] =  -m_dt * eflux_dq;
        for( integer dof=0; dof < CP_Deriv::nDer; dof++ )
        {
          stack.localEnergyFluxJacobian[TAG::NEXT      ][dof]  *=  m_dt*currentConnRate;
          stack.localEnergyFluxJacobian[TAG::CURRENT   ][dof]  *=  -m_dt*currentConnRate;
        }
      }

    } );

  }


  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElements the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack
   * variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElements,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const ie )
    {
      typename KERNEL_TYPE::StackVariables stack( 1 );

      kernelComponent.setup( ie, stack );
      kernelComponent.computeFlux( ie, stack );
      kernelComponent.complete( ie, stack );
    } );
  }

protected:
  /// Number of phases
  integer const m_numPhases;

  /// Global index of local element
  arrayView1d< globalIndex const > m_globalWellElementIndex;

  /// Element phase fraction
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseFraction;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const m_dPhaseFraction;

  /// Views on phase enthalpy
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseEnthalpy;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseEnthalpy;


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
   * @param[in] numComps the number of fluid components
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] useTotalMassEquation flag specifying whether to replace one component bal eqn with total mass eqn
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] wellControls object holding well control/constraint information
   * @param[in] subregion well subregion
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComps,
                   real64 const dt,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   WellControls const & wellControls,
                   WellElementSubRegion const & subRegion,
                   MultiFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();


      BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );


      using kernelType = FaceBasedAssemblyKernel< NUM_COMP >;


      kernelType kernel( dt, rankOffset, dofKey, wellControls, subRegion, fluid, localMatrix, localRhs, kernelFlags );
      kernelType::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }
};

}   // end namespace thermalCompositionalMultiphaseWellKernels

} // end namespace geos

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_THERMALCOMPOSITIONALMULTIPHASEWELLKERNELS_HPP
