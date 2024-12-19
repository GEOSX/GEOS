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
 * @file SolutionScalingZFormulationKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGZFORMULATIONKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGZFORMULATIONKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/zFormulation/SolutionScalingAndCheckingZFormulationKernelBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** SolutionScalingZFormulationKernel ********************************/

/**
 * @class SolutionScalingZFormulationKernel
 * @brief Define the kernel for scaling the Newton update
 */
class SolutionScalingZFormulationKernel : public SolutionScalingAndCheckingZFormulationKernelBase< real64 >
{
public:

  using Base = SolutionScalingAndCheckingZFormulationKernelBase< real64 >;
  using Base::m_rankOffset;
  using Base::m_numComp;
  using Base::m_dofNumber;
  using Base::m_ghostRank;
  using Base::m_localSolution;
  using Base::m_pressure;
  using Base::m_compFrac;
  using Base::m_pressureScalingFactor;
  using Base::m_compFracScalingFactor;

  /**
   * @brief Create a new kernel instance
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compFrac the component density vector
   * @param[in] pressureScalingFactor the pressure local scaling factor
   * @param[in] compFracScalingFactor the component density local scaling factor
   */
  SolutionScalingZFormulationKernel( real64 const maxRelativePresChange,
                                     real64 const maxAbsolutePresChange,
                                     real64 const maxCompFracChange,
                                     globalIndex const rankOffset,
                                     integer const numComp,
                                     string const dofKey,
                                     ElementSubRegionBase const & subRegion,
                                     arrayView1d< real64 const > const localSolution,
                                     arrayView1d< real64 const > const pressure,
                                     arrayView2d< real64 const, compflow::USD_COMP > const compFrac,
                                     arrayView1d< real64 > pressureScalingFactor,
                                     arrayView1d< real64 > compFracScalingFactor )
    : Base( rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compFrac,
            pressureScalingFactor,
            compFracScalingFactor ),
    m_maxRelativePresChange( maxRelativePresChange ),
    m_maxAbsolutePresChange( maxAbsolutePresChange ),
    m_maxCompFracChange( maxCompFracChange )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
    GEOS_HOST_DEVICE
    StackVariables()
    { }

    StackVariables( real64 _localMinVal,
                    real64 _localMaxDeltaPres,
                    real64 _localMaxDeltaTemp,
                    real64 _localMaxDeltaCompFrac,
                    real64 _localMinPresScalingFactor,
                    real64 _localMinTempScalingFactor,
                    real64 _localMinCompFracScalingFactor )
      :
      Base::StackVariables( _localMinVal ),
      localMaxDeltaPres( _localMaxDeltaPres ),
      localMaxDeltaTemp( _localMaxDeltaTemp ),
      localMaxDeltaCompFrac( _localMaxDeltaCompFrac ),
      localMinPresScalingFactor( _localMinPresScalingFactor ),
      localMinTempScalingFactor( _localMinTempScalingFactor ),
      localMinCompFracScalingFactor( _localMinCompFracScalingFactor )
    { }

    real64 localMaxDeltaPres;
    real64 localMaxDeltaTemp;
    real64 localMaxDeltaCompFrac;

    real64 localMinPresScalingFactor;
    real64 localMinTempScalingFactor;
    real64 localMinCompFracScalingFactor;

  };

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static StackVariables
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > globalScalingFactor( 1.0 );

    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaPres( 0.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaTemp( 0.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaCompFrac( 0.0 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPresScalingFactor( 1.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minTempScalingFactor( 1.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minCompFracScalingFactor( 1.0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );

      globalScalingFactor.min( stack.localMinVal );

      maxDeltaPres.max( stack.localMaxDeltaPres );
      maxDeltaTemp.max( stack.localMaxDeltaTemp );
      maxDeltaCompFrac.max( stack.localMaxDeltaCompFrac );

      minPresScalingFactor.min( stack.localMinPresScalingFactor );
      minTempScalingFactor.min( stack.localMinTempScalingFactor );
      minCompFracScalingFactor.min( stack.localMinCompFracScalingFactor );
    } );

    return StackVariables( globalScalingFactor.get(),
                           maxDeltaPres.get(),
                           maxDeltaTemp.get(),
                           maxDeltaCompFrac.get(),
                           minPresScalingFactor.get(),
                           minTempScalingFactor.get(),
                           minCompFracScalingFactor.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.localMaxDeltaPres = 0.0;
    stack.localMaxDeltaTemp = 0.0;
    stack.localMaxDeltaCompFrac = 0.0;

    stack.localMinPresScalingFactor = 1.0;
    stack.localMinTempScalingFactor = 1.0;
    stack.localMinCompFracScalingFactor = 1.0;
  }

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
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeScalingFactor( localIndex const ei,
                             StackVariables & stack,
                             FUNC && kernelOp = NoOpFunc{} ) const
  {
    real64 constexpr eps = minDensForDivision;

    // compute the change in pressure
    real64 const pres = m_pressure[ei];
    real64 const absPresChange = LvArray::math::abs( m_localSolution[stack.localRow] );
    if( stack.localMaxDeltaPres < absPresChange )
    {
      stack.localMaxDeltaPres = absPresChange;
    }

    // compute pressure scaling factor
    real64 presScalingFactor = 1.0;
    // when enabled, absolute change scaling has a priority over relative change
    if( m_maxAbsolutePresChange > 0.0 ) // maxAbsolutePresChange <= 0.0 means that absolute scaling is disabled
    {
      if( absPresChange > m_maxAbsolutePresChange )
      {
        presScalingFactor = m_maxAbsolutePresChange / absPresChange;
      }
    }
    else if( pres > eps )
    {
      real64 const relativePresChange = absPresChange / pres;
      if( relativePresChange > m_maxRelativePresChange )
      {
        presScalingFactor = m_maxRelativePresChange / relativePresChange;
      }
    }
    m_pressureScalingFactor[ei] = presScalingFactor;
    if( stack.localMinVal > presScalingFactor )
    {
      stack.localMinVal = presScalingFactor;
    }
    if( stack.localMinPresScalingFactor > presScalingFactor )
    {
      stack.localMinPresScalingFactor = presScalingFactor;
    }

    // Component Fraction Change Scaling
    m_compFracScalingFactor[ei] = 1.0;

    for( integer ic = 0; ic < m_numComp; ++ic )
    {
      // compute scaling factor based on relative change in component densities
      real64 const absCompFracChange = LvArray::math::abs( m_localSolution[stack.localRow + ic + 1] );
      if( stack.localMaxDeltaCompFrac < absCompFracChange )
      {
        stack.localMaxDeltaCompFrac = absCompFracChange;
      }

      if( absCompFracChange > m_maxCompFracChange && absCompFracChange > eps )
      {
        real64 const compScalingFactor = m_maxCompFracChange / absCompFracChange;
        m_compFracScalingFactor[ei] = LvArray::math::min( m_compFracScalingFactor[ei], compScalingFactor );
        if( stack.localMinVal > compScalingFactor )
        {
          stack.localMinVal = compScalingFactor;
        }
        if( stack.localMinCompFracScalingFactor > compScalingFactor )
        {
          stack.localMinCompFracScalingFactor = compScalingFactor;
        }
      }
    }

    // compute the scaling factor for other vars, such as temperature
    kernelOp();
  }

protected:

  /// Max allowed changes in primary variables
  real64 const m_maxRelativePresChange;
  real64 const m_maxAbsolutePresChange;
  real64 const m_maxCompFracChange;

};

/**
 * @class SolutionScalingZFormulationKernelFactory
 */
class SolutionScalingZFormulationKernelFactory
{
public:

  /*
   * @brief Create and launch the kernel computing the scaling factor
   * @tparam POLICY the kernel policy
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @return the scaling factor
   */
  template< typename POLICY >
  static SolutionScalingZFormulationKernel::StackVariables
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxAbsolutePresChange,
                   real64 const maxCompFracChange,
                   arrayView1d< real64 const > const pressure,
                   arrayView2d< real64 const, compflow::USD_COMP > const compFrac,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > compFracScalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    SolutionScalingZFormulationKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxCompFracChange, rankOffset,
                                              numComp, dofKey, subRegion, localSolution, pressure, compFrac, pressureScalingFactor, compFracScalingFactor );
    return SolutionScalingZFormulationKernel::launch< POLICY >( subRegion.size(), kernel );
  }
};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGZFORMULATIONKERNEL_HPP
