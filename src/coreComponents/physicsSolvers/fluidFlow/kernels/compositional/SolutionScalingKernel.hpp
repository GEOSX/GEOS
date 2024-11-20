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
 * @file SolutionScalingKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingAndCheckingKernelBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** SolutionScalingKernel ********************************/

/**
 * @class SolutionScalingKernel
 * @brief Define the kernel for scaling the Newton update
 */
class SolutionScalingKernel : public SolutionScalingAndCheckingKernelBase< real64 >
{
public:

  using Base = SolutionScalingAndCheckingKernelBase< real64 >;
  using Base::m_rankOffset;
  using Base::m_numComp;
  using Base::m_dofNumber;
  using Base::m_ghostRank;
  using Base::m_localSolution;
  using Base::m_pressure;
  using Base::m_compDens;
  using Base::m_pressureScalingFactor;
  using Base::m_compDensScalingFactor;

  /**
   * @brief Create a new kernel instance
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompDensChange the max allowed comp density change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @param[in] pressure the pressure vector
   * @param[in] compDens the component density vector
   * @param[in] pressureScalingFactor the pressure local scaling factor
   * @param[in] compDensScalingFactor the component density local scaling factor
   */
  SolutionScalingKernel( real64 const maxRelativePresChange,
                         real64 const maxAbsolutePresChange,
                         real64 const maxCompFracChange,
                         real64 const maxRelativeCompDensChange,
                         globalIndex const rankOffset,
                         integer const numComp,
                         string const dofKey,
                         ElementSubRegionBase const & subRegion,
                         arrayView1d< real64 const > const localSolution,
                         arrayView1d< real64 const > const pressure,
                         arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                         arrayView1d< real64 > pressureScalingFactor,
                         arrayView1d< real64 > compDensScalingFactor )
    : Base( rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor ),
    m_maxRelativePresChange( maxRelativePresChange ),
    m_maxAbsolutePresChange( maxAbsolutePresChange ),
    m_maxCompFracChange( maxCompFracChange ),
    m_maxRelativeCompDensChange( maxRelativeCompDensChange )
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
                    localIndex _localMaxDeltaPresLoc,
                    real64 _localMaxDeltaTemp,
                    localIndex _localMaxDeltaTempLoc,
                    real64 _localMaxDeltaCompDens,
                    localIndex _localMaxDeltaCompDensLoc,
                    real64 _localMinPresScalingFactor,
                    real64 _localMinTempScalingFactor,
                    real64 _localMinCompDensScalingFactor )
      :
      Base::StackVariables( _localMinVal ),
      localMaxDeltaPres( _localMaxDeltaPres ),
      localMaxDeltaPresLoc( _localMaxDeltaPresLoc ),
      localMaxDeltaTemp( _localMaxDeltaTemp ),
      localMaxDeltaTempLoc( _localMaxDeltaTempLoc ),
      localMaxDeltaCompDens( _localMaxDeltaCompDens ),
      localMaxDeltaCompDensLoc( _localMaxDeltaCompDensLoc ),
      localMinPresScalingFactor( _localMinPresScalingFactor ),
      localMinTempScalingFactor( _localMinTempScalingFactor ),
      localMinCompDensScalingFactor( _localMinCompDensScalingFactor )
    { }

    real64 localMaxDeltaPres;
    localIndex localMaxDeltaPresLoc;
    real64 localMaxDeltaTemp;
    localIndex localMaxDeltaTempLoc;
    real64 localMaxDeltaCompDens;
    localIndex localMaxDeltaCompDensLoc;

    real64 localMinPresScalingFactor;
    real64 localMinTempScalingFactor;
    real64 localMinCompDensScalingFactor;

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

    RAJA::ReduceMaxLoc< ReducePolicy< POLICY >, real64 > maxDeltaPres( std::numeric_limits< real64 >::min(), -1 );
    RAJA::ReduceMaxLoc< ReducePolicy< POLICY >, real64 > maxDeltaTemp( std::numeric_limits< real64 >::min(), -1 );
    RAJA::ReduceMaxLoc< ReducePolicy< POLICY >, real64 > maxDeltaCompDens( std::numeric_limits< real64 >::min(), -1 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minPresScalingFactor( 1.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minTempScalingFactor( 1.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minCompDensScalingFactor( 1.0 );

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

      maxDeltaPres.maxloc( stack.localMaxDeltaPres, ei );
      maxDeltaTemp.maxloc( stack.localMaxDeltaTemp, ei );
      maxDeltaCompDens.maxloc( stack.localMaxDeltaCompDens, ei );

      minPresScalingFactor.min( stack.localMinPresScalingFactor );
      minTempScalingFactor.min( stack.localMinTempScalingFactor );
      minCompDensScalingFactor.min( stack.localMinCompDensScalingFactor );
    } );

    return StackVariables( globalScalingFactor.get(),
                           maxDeltaPres.get(),
                           maxDeltaPres.getLoc(),
                           maxDeltaTemp.get(),
                           maxDeltaTemp.getLoc(),
                           maxDeltaCompDens.get(),
                           maxDeltaCompDens.getLoc(),
                           minPresScalingFactor.get(),
                           minTempScalingFactor.get(),
                           minCompDensScalingFactor.get() );
  }

  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    Base::setup( ei, stack );

    stack.localMaxDeltaPres = 0.0;
    stack.localMaxDeltaPresLoc = -1;
    stack.localMaxDeltaTemp = 0.0;
    stack.localMaxDeltaTempLoc = -1;
    stack.localMaxDeltaCompDens = 0.0;
    stack.localMaxDeltaCompDensLoc =-1;

    stack.localMinPresScalingFactor = 1.0;
    stack.localMinTempScalingFactor = 1.0;
    stack.localMinCompDensScalingFactor = 1.0;
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

    real64 prevTotalDens = 0;
    for( integer ic = 0; ic < m_numComp; ++ic )
    {
      prevTotalDens += m_compDens[ei][ic];
    }

    m_compDensScalingFactor[ei] = 1.0;

    // compute the change in component densities and component fractions
    for( integer ic = 0; ic < m_numComp; ++ic )
    {
      // compute scaling factor based on relative change in component densities
      real64 const absCompDensChange = LvArray::math::abs( m_localSolution[stack.localRow + ic + 1] );
      if( stack.localMaxDeltaCompDens < absCompDensChange )
      {
        stack.localMaxDeltaCompDens = absCompDensChange;
      }

      // This actually checks the change in component fraction, using a lagged total density
      // Indeed we can rewrite the following check as:
      //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
      // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
      // because I found it more robust than using directly newTotalDens (which can vary also
      // wildly when the compDens change is large)
      real64 const maxAbsCompDensChange = m_maxCompFracChange * prevTotalDens;
      if( absCompDensChange > maxAbsCompDensChange && absCompDensChange > eps )
      {
        real64 const compScalingFactor = maxAbsCompDensChange / absCompDensChange;
        m_compDensScalingFactor[ei] = LvArray::math::min( m_compDensScalingFactor[ei], compScalingFactor );
        if( stack.localMinVal > compScalingFactor )
        {
          stack.localMinVal = compScalingFactor;
        }
        if( stack.localMinCompDensScalingFactor > compScalingFactor )
        {
          stack.localMinCompDensScalingFactor = compScalingFactor;
        }
      }

      // switch from relative to absolute when value is < 1.0
      real64 const maxRelCompDensChange = m_maxRelativeCompDensChange * LvArray::math::max( m_compDens[ei][ic], 1.0 );
      if( absCompDensChange > maxRelCompDensChange && absCompDensChange > eps )
      {
        real64 const compScalingFactor = maxRelCompDensChange / absCompDensChange;
        m_compDensScalingFactor[ei] = LvArray::math::min( m_compDensScalingFactor[ei], compScalingFactor );
        if( stack.localMinVal > compScalingFactor )
        {
          stack.localMinVal = compScalingFactor;
        }
        if( stack.localMinCompDensScalingFactor > compScalingFactor )
        {
          stack.localMinCompDensScalingFactor = compScalingFactor;
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
  real64 const m_maxRelativeCompDensChange;

};

/**
 * @class SolutionScalingKernelFactory
 */
class SolutionScalingKernelFactory
{
public:

  /*
   * @brief Create and launch the kernel computing the scaling factor
   * @tparam POLICY the kernel policy
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompDensChange the max allowed comp density change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   * @return the scaling factor
   */
  template< typename POLICY >
  static SolutionScalingKernel::StackVariables
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxAbsolutePresChange,
                   real64 const maxCompFracChange,
                   real64 const maxRelativeCompDensChange,
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
    SolutionScalingKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxCompFracChange, maxRelativeCompDensChange, rankOffset,
                                  numComp, dofKey, subRegion, localSolution, pressure, compDens, pressureScalingFactor, compDensScalingFactor );
    return SolutionScalingKernel::launch< POLICY >( subRegion.size(), kernel );
  }
};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_SOLUTIONSCALINGKERNEL_HPP
