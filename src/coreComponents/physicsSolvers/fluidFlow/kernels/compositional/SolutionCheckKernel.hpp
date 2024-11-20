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
 * @file SolutionCheckKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONCHECKKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONCHECKKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingAndCheckingKernelBase.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernel
 * @brief Define the kernel for checking the updated solution
 */
class SolutionCheckKernel : public SolutionScalingAndCheckingKernelBase< integer >
{
public:

  using Base = SolutionScalingAndCheckingKernelBase< integer >;
  using Base::m_rankOffset;
  using Base::m_numComp;
  using Base::m_dofNumber;
  using Base::m_ghostRank;
  using Base::m_localSolution;
  using Base::m_pressure;
  using Base::m_compDens;

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
   * @param[in] compDens the component density vector
   */
  SolutionCheckKernel( integer const allowCompDensChopping,
                       integer const allowNegativePressure,
                       CompositionalMultiphaseFVM::ScalingType const scalingType,
                       real64 const scalingFactor,
                       arrayView1d< real64 const > const pressure,
                       arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                       arrayView1d< real64 > pressureScalingFactor,
                       arrayView1d< real64 > compDensScalingFactor,
                       globalIndex const rankOffset,
                       integer const numComp,
                       string const dofKey,
                       ElementSubRegionBase const & subRegion,
                       arrayView1d< real64 const > const localSolution )
    : Base( rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor ),
    m_allowCompDensChopping( allowCompDensChopping ),
    m_allowNegativePressure( allowNegativePressure ),
    m_scalingFactor( scalingFactor ),
    m_scalingType( scalingType )
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
                    real64 _localMinPres,
                    real64 _localMinDens,
                    real64 _localMinTotalDens,
                    integer _localNumNegPressures,
                    integer _localNumNegDens,
                    integer _localNumNegTotalDens )
      :
      Base::StackVariables( _localMinVal ),
      localMinPres( _localMinPres ),
      localMinDens( _localMinDens ),
      localMinTotalDens( _localMinTotalDens ),
      localNumNegPressures( _localNumNegPressures ),
      localNumNegDens( _localNumNegDens ),
      localNumNegTotalDens( _localNumNegTotalDens )
    { }

    real64 localMinPres;
    real64 localMinDens;
    real64 localMinTotalDens;

    integer localNumNegPressures;
    integer localNumNegDens;
    integer localNumNegTotalDens;

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
    RAJA::ReduceMin< ReducePolicy< POLICY >, integer > globalMinVal( 1 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPres( 0.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minDens( 0.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minTotalDens( 0.0 );

    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegPressures( 0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegDens( 0 );
    RAJA::ReduceSum< ReducePolicy< POLICY >, integer > numNegTotalDens( 0 );

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.ghostRank( ei ) >= 0 )
      {
        return;
      }

      StackVariables stack;
      kernelComponent.setup( ei, stack );
      kernelComponent.compute( ei, stack );

      globalMinVal.min( stack.localMinVal );

      minPres.min( stack.localMinPres );
      minDens.min( stack.localMinDens );
      minTotalDens.min( stack.localMinTotalDens );

      numNegPressures += stack.localNumNegPressures;
      numNegDens += stack.localNumNegDens;
      numNegTotalDens += stack.localNumNegTotalDens;
    } );

    return StackVariables( globalMinVal.get(),
                           minPres.get(),
                           minDens.get(),
                           minTotalDens.get(),
                           numNegPressures.get(),
                           numNegDens.get(),
                           numNegTotalDens.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.localMinPres = 0.0;
    stack.localMinDens = 0.0;
    stack.localMinTotalDens = 0.0;

    stack.localNumNegPressures = 0;
    stack.localNumNegDens = 0;
    stack.localNumNegTotalDens = 0;
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
    computeSolutionCheck( ei, stack );
  }

  /**
   * @brief Compute the local value of the check
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] kernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeSolutionCheck( localIndex const ei,
                             StackVariables & stack,
                             FUNC && kernelOp = NoOpFunc{} ) const
  {
    bool const localScaling = m_scalingType == CompositionalMultiphaseFVM::ScalingType::Local;

    real64 const newPres = m_pressure[ei] + (localScaling ? m_pressureScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow];
    if( newPres < 0 )
    {
      if( !m_allowNegativePressure )
      {
        stack.localMinVal = 0;
      }
      stack.localNumNegPressures += 1;
      if( newPres < stack.localMinPres )
        stack.localMinPres = newPres;
    }

    // if component density chopping is not allowed, the time step fails if a component density is negative
    // otherwise, we just check that the total density is positive, and negative component densities
    // will be chopped (i.e., set to zero) in ApplySystemSolution)
    if( !m_allowCompDensChopping )
    {
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + (localScaling ? m_compDensScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + ic + 1];
        if( newDens < 0 )
        {
          stack.localMinVal = 0;
          stack.localNumNegDens += 1;
          if( newDens < stack.localMinDens )
            stack.localMinDens = newDens;
        }
      }
    }
    else
    {
      real64 totalDens = 0.0;
      for( integer ic = 0; ic < m_numComp; ++ic )
      {
        real64 const newDens = m_compDens[ei][ic] + (localScaling ? m_compDensScalingFactor[ei] : m_scalingFactor) * m_localSolution[stack.localRow + ic + 1];
        totalDens += ( newDens > 0.0 ) ? newDens : 0.0;
      }
      if( totalDens < 0 )
      {
        stack.localMinVal = 0;
        stack.localNumNegTotalDens += 1;
        if( totalDens < stack.localMinTotalDens )
          stack.localMinTotalDens = totalDens;
      }
    }

    kernelOp();
  }

protected:

  /// flag to allow the component density chopping
  integer const m_allowCompDensChopping;

  /// flag to allow negative pressure values
  integer const m_allowNegativePressure;

  /// scaling factor
  real64 const m_scalingFactor;

  /// scaling type (global or local)
  CompositionalMultiphaseFVM::ScalingType const m_scalingType;

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
   * @param[in] allowCompDensChopping flag to allow the component density chopping
   * @param[in] scalingFactor the scaling factor
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
                   arrayView1d< real64 const > const pressure,
                   arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > compDensScalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    SolutionCheckKernel kernel( allowCompDensChopping, allowNegativePressure, scalingType, scalingFactor,
                                pressure, compDens, pressureScalingFactor, compDensScalingFactor, rankOffset,
                                numComp, dofKey, subRegion, localSolution );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONCHECKKERNEL_HPP
