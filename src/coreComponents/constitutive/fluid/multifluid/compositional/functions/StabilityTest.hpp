/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StabilityTest.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_

#include "KValueInitialization.hpp"
#include "FugacityCalculator.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct StabilityTest
{
private:
  static constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
public:
  /**
   * @brief Perform a two-phase stability test
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] equationOfState The equation of state
   * @param[out] tangentPlaneDistance the minimum tangent plane distance (TPD)
   * @param[out] kValues the k-values estimated from the stationary points
   * @return a flag indicating that 2 stationary points have been found
   */
  template< integer USD1 >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, USD1 > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       EquationOfStateType const & equationOfState,
                       real64 & tangentPlaneDistance,
                       arraySlice1d< real64 > const & kValues )
  {
    constexpr integer numTrials = 2;    // Trial compositions
    stackArray1d< real64, maxNumComps > logFugacity( numComps );
    stackArray1d< real64, maxNumComps > normalizedComposition( numComps );
    stackArray2d< real64, numTrials *maxNumComps > trialComposition( numTrials, numComps );
    stackArray1d< real64, maxNumComps > logTrialComposition( numComps );
    stackArray1d< real64, maxNumComps > hyperplane( numComps );     // h-parameter
    stackArray1d< integer, maxNumComps > availableComponents( numComps );

    calculatePresentComponents( numComps, composition, availableComponents );
    auto const presentComponents = availableComponents.toSliceConst();

    // Calculate the hyperplane parameter
    // h_i = log( z_i ) + log( phi_i )
    for( integer ic = 0; ic < numComps; ++ic )
    {
      hyperplane[ic] = 0.0;
    }
    FugacityCalculator::computeLogFugacity( numComps,
                                            pressure,
                                            temperature,
                                            composition,
                                            componentProperties,
                                            equationOfState,
                                            logFugacity );
    for( integer const ic : presentComponents )
    {
      hyperplane[ic] = LvArray::math::log( composition[ic] ) + logFugacity[ic];
    }

    // Initialise the trial compositions using Wilson k-Values
    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        kValues );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      trialComposition( 0, ic ) = composition[ic] / kValues[ic];
      trialComposition( 1, ic ) = composition[ic] * kValues[ic];
    }

    integer numberOfStationaryPoints = 0;
    tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;
    for( integer trialIndex = 0; trialIndex < numTrials; ++trialIndex )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        normalizedComposition[ic] = trialComposition( trialIndex, ic );
      }
      normalizeComposition( numComps, normalizedComposition.toSlice() );

      FugacityCalculator::computeLogFugacity( numComps,
                                              pressure,
                                              temperature,
                                              normalizedComposition.toSliceConst(),
                                              componentProperties,
                                              equationOfState,
                                              logFugacity );
      for( integer const ic : presentComponents )
      {
        logTrialComposition[ic] = LvArray::math::log( trialComposition( trialIndex, ic ) );
      }
      for( localIndex iterationCount = 0; iterationCount < MultiFluidConstants::maxSSIIterations; ++iterationCount )
      {
        for( integer const ic : presentComponents )
        {
          logTrialComposition[ic] = hyperplane[ic] - logFugacity[ic];
          trialComposition( trialIndex, ic ) = LvArray::math::exp( logTrialComposition[ic] );
          normalizedComposition[ic] = trialComposition( trialIndex, ic );
        }
        normalizeComposition( numComps, normalizedComposition.toSlice() );

        FugacityCalculator::computeLogFugacity( numComps,
                                                pressure,
                                                temperature,
                                                normalizedComposition.toSliceConst(),
                                                componentProperties,
                                                equationOfState,
                                                logFugacity );

        real64 error = 0.0;
        for( integer const ic : presentComponents )
        {
          real64 const dG =  logTrialComposition[ic] + logFugacity[ic] - hyperplane[ic];
          error += (dG*dG);
        }
        error = LvArray::math::sqrt( error );

        if( error < MultiFluidConstants::fugacityTolerance )
        {
          // Calculate modified tangent plane distance (Michelsen, 1982b) of trial composition relative to input composition
          real64 tpd = 1.0;
          for( integer const ic : presentComponents )
          {
            tpd += trialComposition( trialIndex, ic ) * (logTrialComposition[ic] + logFugacity[ic] - hyperplane[ic] - 1.0);
          }
          if( tpd < tangentPlaneDistance )
          {
            tangentPlaneDistance = tpd;
          }
          numberOfStationaryPoints++;
          break;
        }
      }
    }
    if( numberOfStationaryPoints == numTrials )
    {
      for( integer const ic : presentComponents )
      {
        kValues[ic] = trialComposition( 1, ic ) / trialComposition( 0, ic );
      }
    }
    return numberOfStationaryPoints == numTrials;
  }

private:
  /**
   * @brief Calculate which components are present.
   * @details Creates a list of indices whose components have non-zero mole fraction.
   * @param[in] numComps number of components
   * @param[in] composition the composition of the fluid
   * @param[out] presentComponents the list of present components
   * @return the number of present components
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static integer calculatePresentComponents( integer const numComps,
                                             arraySlice1d< real64 const > const & composition,
                                             stackArray1d< integer, maxNumComps > & presentComponents )
  {
    // Check for machine-zero feed values
    integer presentCount = 0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      if( MultiFluidConstants::epsilon < composition[ic] )
      {
        presentComponents[presentCount++] = ic;
      }
    }
    presentComponents.resize( presentCount );
    return presentCount;
  }

  /**
   * @brief Normalise a composition in place to ensure that the components add up to unity
   * @param[in] numComps number of components
   * @param[in/out] composition composition to be normalized
   * @return the sum of the given values
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 normalizeComposition( integer const numComps,
                                      arraySlice1d< real64, USD > const & composition )
  {
    real64 totalMoles = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      totalMoles += composition[ic];
    }
    GEOS_ASSERT( MultiFluidConstants::epsilon < totalMoles );
    real64 const oneOverTotalMoles = 1.0 / totalMoles;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      composition[ic] *= oneOverTotalMoles;
    }
    return totalMoles;
  }
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
