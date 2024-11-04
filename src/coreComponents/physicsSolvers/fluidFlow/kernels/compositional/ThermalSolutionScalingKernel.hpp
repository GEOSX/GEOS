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
 * @file ThermalSolutionScalingKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALSOLUTIONSCALINGKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALSOLUTIONSCALINGKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingKernel.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

/******************************** SolutionScalingKernel ********************************/

/**
 * @class SolutionScalingKernel
 * @brief Define the kernel for scaling the Newton update
 */
class SolutionScalingKernel : public isothermalCompositionalMultiphaseBaseKernels::SolutionScalingKernel
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::SolutionScalingKernel;
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
  SolutionScalingKernel( real64 const maxRelativePresChange,
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
                         arrayView1d< real64 > temperatureScalingFactor,
                         integer const temperatureOffset )
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
    m_temperatureScalingFactor( temperatureScalingFactor ),
    m_temperatureOffset( temperatureOffset )
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
      real64 const absTempChange = LvArray::math::abs( m_localSolution[stack.localRow + m_temperatureOffset] );
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

  /// Temperature offset in solution array
  integer const m_temperatureOffset;
};

/**
 * @class SolutionScalingKernelFactory
 */
class SolutionScalingKernelFactory
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
  static SolutionScalingKernel::StackVariables
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxAbsolutePresChange,
                   real64 const maxRelativeTempChange,
                   real64 const maxCompFracChange,
                   real64 const maxRelativeCompDensChange,
                   arrayView1d< real64 const > const pressure,
                   arrayView1d< real64 const > const temperature,
                   arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > compDensScalingFactor,
                   arrayView1d< real64 > temperatureScalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution,
                   integer const temperatureOffset )
  {
    SolutionScalingKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxRelativeTempChange,
                                  maxCompFracChange, maxRelativeCompDensChange,
                                  rankOffset, numComp, dofKey, subRegion, localSolution,
                                  pressure, temperature, compDens, pressureScalingFactor,
                                  compDensScalingFactor, temperatureScalingFactor, temperatureOffset );
    return thermalCompositionalMultiphaseBaseKernels::
             SolutionScalingKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

} // namespace thermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALSOLUTIONSCALINGKERNEL_HPP
