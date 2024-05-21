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
 * @file StabilityFlash.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYFLASH_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYFLASH_HPP_

#include "KValueInitialization.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct StabilityFlash
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Perform negative two-phase EOS flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @return the minimum tangent plane distance (TPD)
   */
  template< typename EOS_TYPE, integer USD1 >
  GEOS_HOST_DEVICE
  static real64 compute( integer const numComps,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, USD1 > const & composition,
                         ComponentProperties::KernelWrapper const & componentProperties )
  {
    constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
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
    hyperplane.zero();
    EOS_TYPE::computeLogFugacityCoefficients( numComps,
                                              pressure,
                                              temperature,
                                              composition,
                                              componentProperties,
                                              logFugacity );
    for( integer const ic : presentComponents )
    {
      hyperplane[ic] = LvArray::math::log( composition[ic] ) + logFugacity[ic];
    }

    // Initialise the trial compositions using Wilson k-Values
    // Use fugacity space as temporary storage
    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        logFugacity.toSlice() );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      trialComposition( 0, ic ) = composition[ic] / logFugacity[ic];
      trialComposition( 1, ic ) = composition[ic] * logFugacity[ic];
    }

    real64 tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;
    for( integer trialIndex = 0; trialIndex < numTrials; ++trialIndex )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        normalizedComposition[ic] = trialComposition( trialIndex, ic );
      }
      normalizeComposition( numComps, normalizedComposition.toSlice() );
      EOS_TYPE::computeLogFugacityCoefficients( numComps,
                                                pressure,
                                                temperature,
                                                normalizedComposition.toSliceConst(),
                                                componentProperties,
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
        EOS_TYPE::computeLogFugacityCoefficients( numComps,
                                                  pressure,
                                                  temperature,
                                                  normalizedComposition.toSliceConst(),
                                                  componentProperties,
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
          double tpd = 1.0;

          for( integer const ic : presentComponents )
          {
            tpd += trialComposition( trialIndex, ic ) * (logTrialComposition[ic] + logFugacity[ic] - hyperplane[ic] - 1.0);
          }
          if( tpd < tangentPlaneDistance )
          {
            tangentPlaneDistance = tpd;
          }
          break;
        }
      }
    }
    return tangentPlaneDistance;
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
  template< typename ARRAY >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static integer calculatePresentComponents( integer const numComps,
                                             arraySlice1d< real64 const > const & composition,
                                             ARRAY & presentComponents )
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
    real64 const oneOverTotalMoles = 1.0 / (totalMoles + MultiFluidConstants::epsilon);
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

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYFLASH_HPP_
