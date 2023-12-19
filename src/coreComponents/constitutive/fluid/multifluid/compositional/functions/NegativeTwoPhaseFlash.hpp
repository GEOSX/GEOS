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
 * @file NegativeTwoPhaseFlash.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_

#include "common/DataTypes.hpp"
#include "RachfordRice.hpp"
#include "KValueInitialization.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct NegativeTwoPhaseFlash
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
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   * @return an indicator of success of the flash
   */
  template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       real64 & vapourPhaseMoleFraction,
                       arraySlice1d< real64 > const & liquidComposition,
                       arraySlice1d< real64 > const & vapourComposition )
  {
    constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
    stackArray1d< real64, maxNumComps > logLiquidFugacity( numComps );
    stackArray1d< real64, maxNumComps > logVapourFugacity( numComps );
    stackArray1d< real64, maxNumComps > kVapourLiquid( numComps );
    stackArray1d< real64, maxNumComps > fugacityRatios( numComps );
    stackArray1d< integer, maxNumComps > presentComponents( numComps );

    // Initialise compositions to feed composition
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic];
      vapourComposition[ic] = composition[ic];
    }

    calculatePresentComponents( numComps, composition, presentComponents );

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        kVapourLiquid );

    bool converged = false;
    for( localIndex iterationCount = 0; iterationCount < MultiFluidConstants::maxSSIIterations; ++iterationCount )
    {
      // Solve Rachford-Rice Equation
      vapourPhaseMoleFraction = RachfordRice::solve( kVapourLiquid, composition, presentComponents );

      // Assign phase compositions
      for( integer const ic : presentComponents )
      {
        liquidComposition[ic] = composition[ic] / ( 1.0 + vapourPhaseMoleFraction * ( kVapourLiquid[ic] - 1.0 ) );
        vapourComposition[ic] = kVapourLiquid[ic] * liquidComposition[ic];
      }

      normalizeComposition( numComps, liquidComposition );
      normalizeComposition( numComps, vapourComposition );

      // Compute the phase fugacities
      EOS_TYPE_LIQUID::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       logLiquidFugacity );
      EOS_TYPE_VAPOUR::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       logVapourFugacity );

      // Compute fugacity ratios and check convergence
      converged = true;

      for( integer const ic : presentComponents )
      {
        fugacityRatios[ic] = exp( logLiquidFugacity[ic] - logVapourFugacity[ic] ) * liquidComposition[ic] / vapourComposition[ic];
        if( MultiFluidConstants::fugacityTolerance < LvArray::math::abs( fugacityRatios[ic] - 1.0 ) )
        {
          converged = false;
        }
      }

      if( converged )
      {
        break;
      }

      // Update K-values
      for( integer const ic : presentComponents )
      {
        kVapourLiquid[ic] *= fugacityRatios[ic];
      }
    }

    // Retrieve physical bounds from negative flash values
    if( vapourPhaseMoleFraction < MultiFluidConstants::epsilon )
    {
      vapourPhaseMoleFraction = 0.0;
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidComposition[ic] = composition[ic];
        vapourComposition[ic] = composition[ic];
      }
    }
    else if( 1.0 - vapourPhaseMoleFraction < MultiFluidConstants::epsilon )
    {
      vapourPhaseMoleFraction = 1.0;
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidComposition[ic] = composition[ic];
        vapourComposition[ic] = composition[ic];
      }
    }

    return converged;
  }

  /**
   * @brief Calculate derivatives from the two-phase negative flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] vapourFraction the calculated vapour (gas) mole fraction
   * @param[in] liquidComposition the calculated liquid phase composition
   * @param[in] vapourComposition the calculated vapour phase composition
   * @param[out] vapourFractionDerivs derivatives of the calculated vapour (gas) mole fraction
   * @param[out] liquidCompositionDerivs derivatives of the calculated liquid phase composition
   * @param[out] vapourCompositionDerivs derivatives of the calculated vapour phase composition
   */
  template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
  GEOS_HOST_DEVICE
  static void computeDerivatives( integer const numComps,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  real64 const vapourFraction,
                                  arraySlice1d< real64 const > const & liquidComposition,
                                  arraySlice1d< real64 const > const & vapourComposition,
                                  arraySlice1d< real64 > const & vapourFractionDerivs,
                                  arraySlice2d< real64 > const & liquidCompositionDerivs,
                                  arraySlice2d< real64 > const & vapourCompositionDerivs )
  {
    constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
    constexpr integer maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;

    integer const numDofs = numComps + 2;

    auto const setZero = []( real64 & val ) { val = 0.0; };
    LvArray::forValuesInSlice( vapourFractionDerivs, setZero );
    LvArray::forValuesInSlice( liquidCompositionDerivs, setZero );
    LvArray::forValuesInSlice( vapourCompositionDerivs, setZero );

    // Check if we are single or 2-phase
    if( vapourFraction < MultiFluidConstants::epsilon )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
        vapourCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
      }
    }
    else if( 1.0 - vapourFraction < MultiFluidConstants::epsilon )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
        vapourCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
      }
    }
    else
    {
      // Calculate the liquid and vapour fugacities and derivatives
      stackArray1d< real64, maxNumComps > logLiquidFugacity( numComps );
      stackArray1d< real64, maxNumComps > logVapourFugacity( numComps );
      stackArray2d< real64, maxNumComps * maxNumDofs > logLiquidFugacityDerivs( numComps, numDofs );
      stackArray2d< real64, maxNumComps * maxNumDofs > logVapourFugacityDerivs( numComps, numDofs );

      EOS_TYPE_LIQUID::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       logLiquidFugacity );
      EOS_TYPE_LIQUID::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       logLiquidFugacity,
                                                       logLiquidFugacityDerivs );
      EOS_TYPE_VAPOUR::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       logVapourFugacity );
      EOS_TYPE_VAPOUR::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       logVapourFugacity,
                                                       logVapourFugacityDerivs );

      constexpr integer maxNumVals = 2*MultiFluidConstants::MAX_NUM_COMPONENTS+1;
      integer const numVals = 2*numComps;
      stackArray1d< real64, maxNumVals > b( numVals + 1 );
      stackArray1d< real64, maxNumVals > x( numVals + 1 );
      stackArray2d< real64, maxNumVals * maxNumVals > A( numVals + 1, numVals + 1 );

      auto printSystem = [&](){
        std::cout << "--------------------------------------------------------------------------------\n";
        for( integer row = 0; row <= numVals; row++ )
        {
          for( integer col = 0; col <= numVals; col++ )
          {
            std::cout << std::setw( 15 ) << A( row, col ) << " ";
          }
          std::cout << "| ";
          std::cout << std::setw( 15 ) << x( row ) << " ";
          std::cout << "| ";
          std::cout << std::setw( 15 ) << b( row ) << " ";
          std::cout << "\n";
        }
        std::cout << "--------------------------------------------------------------------------------\n";
      };
      std::cout << std::scientific << std::setprecision( 8 );

      b.zero();
      A.zero();
      for( integer ic = 0; ic < numComps; ++ic )
      {
        integer const xi = ic;
        integer const yi = ic + numComps;
        integer const vi = numVals;

        integer e = ic;
        A( e, xi ) = 1.0 - vapourFraction;
        A( e, yi ) = vapourFraction;
        A( e, vi ) = vapourComposition[ic] - liquidComposition[ic];

        e = ic + numComps;
        real64 const phiL = exp( logLiquidFugacity( ic ) );
        real64 const phiV = exp( logVapourFugacity( ic ) );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          integer const xj = jc;
          integer const yj = jc + numComps;
          real64 const dPhiLdx = logLiquidFugacityDerivs( ic, Deriv::dC+jc );
          real64 const dPhiVdy = logVapourFugacityDerivs( ic, Deriv::dC+jc );
          A( e, xj ) =  liquidComposition[ic] * phiL * dPhiLdx;
          A( e, yj ) = -vapourComposition[ic] * phiV * dPhiVdy;
        }
        A( e, xi ) += phiL;
        A( e, yi ) -= phiV;

        e = numVals;
        A( e, xi ) = -1.0;
        A( e, yi ) =  1.0;
      }
      // Pressure and temperature derivatives
      for( integer const pc : {Deriv::dP, Deriv::dT} )
      {
        for( integer ic = 0; ic < numComps; ++ic )
        {
          real64 const phiL = exp( logLiquidFugacity( ic ) );
          real64 const phiV = exp( logVapourFugacity( ic ) );
          b( ic ) = 0.0;
          b( ic + numComps ) = -liquidComposition[ic] * phiL * logLiquidFugacityDerivs( ic, pc )
                               + vapourComposition[ic] * phiV * logVapourFugacityDerivs( ic, pc );
        }
        b( numVals ) = 0.0;
        BlasLapackLA::solveLinearSystem( A, b, x );
        printSystem();
        for( integer ic = 0; ic < numComps; ++ic )
        {
          liquidCompositionDerivs( ic, pc ) = x( ic );
          vapourCompositionDerivs( ic, pc ) = x( ic + numComps );
        }
        vapourFractionDerivs( pc ) = x( numVals );
      }
      // Composition derivatives
      for( integer kc = 0; kc < numComps; ++kc )
      {
        integer const pc = Deriv::dC + kc;

        for( integer ic = 0; ic < numComps; ++ic )
        {
          b( ic ) = -composition[ic];
          b( ic + numComps ) = 0.0;
        }
        b( kc ) += 1.0;
        b( numVals ) = 0.0;
        BlasLapackLA::solveLinearSystem( A, b, x );
        printSystem();
        for( integer ic = 0; ic < numComps; ++ic )
        {
          liquidCompositionDerivs( ic, pc ) = x( ic );
          vapourCompositionDerivs( ic, pc ) = x( ic + numComps );
        }
        vapourFractionDerivs( pc ) = x( numVals );
      }
    }

    //real64 displacedVapourFraction = -1.0;
    //stackArray1d< real64, maxNumComps > displacedLiquidComposition( numComps );
    //stackArray1d< real64, maxNumComps > displacedVapourComposition( numComps );

/****
    // Pressure derivatives
    real64 const dp = 1.0e-4 * pressure;
    compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >( numComps,
                                                 pressure + dp,
                                                 temperature,
                                                 composition,
                                                 componentProperties,
                                                 displacedVapourFraction,
                                                 displacedLiquidComposition,
                                                 displacedVapourComposition );

    //vapourFractionDerivs[Deriv::dP] = (displacedVapourFraction - vapourFraction) / dp;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidCompositionDerivs( ic, Deriv::dP ) = (displacedLiquidComposition[ic] - liquidComposition[ic]) / dp;
      vapourCompositionDerivs( ic, Deriv::dP ) = (displacedVapourComposition[ic] - vapourComposition[ic]) / dp;
    }

    // Temperature derivatives
    real64 const dT = 1.0e-6 * temperature;
    compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >( numComps,
                                                 pressure,
                                                 temperature + dT,
                                                 composition,
                                                 componentProperties,
                                                 displacedVapourFraction,
                                                 displacedLiquidComposition,
                                                 displacedVapourComposition );

    //vapourFractionDerivs[Deriv::dT] = (displacedVapourFraction - vapourFraction) / dT;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidCompositionDerivs( ic, Deriv::dT ) = (displacedLiquidComposition[ic] - liquidComposition[ic]) / dT;
      vapourCompositionDerivs( ic, Deriv::dT ) = (displacedVapourComposition[ic] - vapourComposition[ic]) / dT;
    }

    // Composition derivatives
    real64 constexpr dz = 1.0e-7;
    stackArray1d< real64, maxNumComps > displacedComposition( numComps );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      displacedComposition[ic] = composition[ic];
    }

    for( integer jc = 0; jc < numComps; ++jc )
    {
      displacedComposition[jc] += dz;
      compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >( numComps,
                                                   pressure,
                                                   temperature,
                                                   displacedComposition,
                                                   componentProperties,
                                                   displacedVapourFraction,
                                                   displacedLiquidComposition,
                                                   displacedVapourComposition );
      displacedComposition[jc] = composition[jc];

      integer const kc = Deriv::dC + jc;
      vapourFractionDerivs[kc] = (displacedVapourFraction - vapourFraction) / dz;
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidCompositionDerivs( ic, kc ) = (displacedLiquidComposition[ic] - liquidComposition[ic]) / dz;
        vapourCompositionDerivs( ic, kc ) = (displacedVapourComposition[ic] - vapourComposition[ic]) / dz;
      }
    }
 */
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
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 normalizeComposition( integer const numComps,
                                      arraySlice1d< real64 > const & composition )
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

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
