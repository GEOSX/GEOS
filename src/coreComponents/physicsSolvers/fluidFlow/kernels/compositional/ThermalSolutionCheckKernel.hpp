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
 * @file ThermalSolutionCheckKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALSOLUTIONCHECKKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALSOLUTIONCHECKKERNEL_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionCheckKernel.hpp"

namespace geos
{

namespace thermalCompositionalMultiphaseBaseKernels
{

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
                       compositionalMultiphaseUtilities::ScalingType const scalingType,
                       real64 const scalingFactor,
                       arrayView1d< real64 const > const pressure,
                       arrayView1d< real64 const > const temperature,
                       arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                       arrayView1d< real64 > pressureScalingFactor,
                       arrayView1d< real64 > compDensScalingFactor,
                       arrayView1d< real64 > temperatureScalingFactor,
                       globalIndex const rankOffset,
                       integer const numComp,
                       string const dofKey,
                       ElementSubRegionBase const & subRegion,
                       arrayView1d< real64 const > const localSolution,
                       integer const temperatureOffset )
    : Base( allowCompDensChopping,
            allowNegativePressure,
            scalingType,
            scalingFactor,
            pressure,
            compDens,
            pressureScalingFactor,
            compDensScalingFactor,
            rankOffset,
            numComp,
            dofKey,
            subRegion,
            localSolution ),
    m_temperature( temperature ),
    m_temperatureScalingFactor( temperatureScalingFactor ),
    m_temperatureOffset( temperatureOffset )
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
      bool const localScaling = m_scalingType == compositionalMultiphaseUtilities::ScalingType::Local;
      // compute the change in temperature
      real64 const newTemp = m_temperature[ei] + (localScaling ? m_temperatureScalingFactor[ei] : m_scalingFactor * m_localSolution[stack.localRow + m_temperatureOffset]);
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

  /// Offset to temperature variable
  integer m_temperatureOffset;

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
                   compositionalMultiphaseUtilities::ScalingType const scalingType,
                   real64 const scalingFactor,
                   arrayView1d< real64 const > const pressure,
                   arrayView1d< real64 const > const temperature,
                   arrayView2d< real64 const, compflow::USD_COMP > const compDens,
                   arrayView1d< real64 > pressureScalingFactor,
                   arrayView1d< real64 > temperatureScalingFactor,
                   arrayView1d< real64 > compDensScalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution,
                   integer temperatureOffset )
  {
    SolutionCheckKernel kernel( allowCompDensChopping, allowNegativePressure, scalingType, scalingFactor,
                                pressure, temperature, compDens, pressureScalingFactor, compDensScalingFactor, temperatureScalingFactor,
                                rankOffset, numComp, dofKey, subRegion, localSolution,
                                temperatureOffset );
    return SolutionCheckKernel::launch< POLICY >( subRegion.size(), kernel );
  }

};

} // namespace thermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_THERMALSOLUTIONCHECKKERNEL_HPP
