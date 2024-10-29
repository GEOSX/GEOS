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
 * @file CubicEOSPhaseModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct PengRobinsonEOS
{
  static constexpr real64 omegaA = 0.457235529;
  static constexpr real64 omegaB = 0.077796074;
  static constexpr real64 delta1 = 2.4142135624; // 1 + sqrt( 2 )
  static constexpr real64 delta2 = -0.4142135624; // 1 - sqrt( 2 )

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64
  evaluate( real64 const & omega )
  {
    return ( omega < 0.49 )
      ? 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
      : 0.3796 + 1.485 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;
  }

  static constexpr char const * catalogName(){ return "PengRobinson"; }
};

struct SoaveRedlichKwongEOS
{
  static constexpr real64 omegaA = 0.42748;
  static constexpr real64 omegaB = 0.08664;
  static constexpr real64 delta1 = 0.0;
  static constexpr real64 delta2 = 1.0;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64
  evaluate( real64 const & omega )
  {
    return 0.480 + 1.574 * omega - 0.176 * omega * omega;
  }

  static constexpr char const * catalogName(){ return "SoaveRedlichKwong"; }
};

template< typename EOS_TYPE >
struct CubicEOSPhaseModel
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Generate a catalog name
   */
  static constexpr char const * catalogName(){ return EOS_TYPE::catalogName(); }

  /**
   * @brief Main entry point of the cubic EOS model
   * @details Computes the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  real64 const & pressure,
                                  real64 const & temperature,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  arraySlice1d< real64 > const & logFugacityCoefficients );

  /**
   * @brief Secondary entry point of the cubic EOS model
   * @details Computes the derivatives of the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] logFugacityCoefficients log of the fugacity coefficients
   * @param[out] logFugacityCoefficientDerivs derivatives of the log of the fugacity coefficients
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  real64 const & pressure,
                                  real64 const & temperature,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  arraySlice1d< real64 const > const & logFugacityCoefficients,
                                  arraySlice2d< real64 > const & logFugacityCoefficientDerivs );

  /**
   * @brief Compute compressibility factor for the cubic EOS model
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] compressibilityFactor the current compressibility factor
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 & compressibilityFactor,
                                arraySlice1d< real64 > const & compressibilityFactorDerivs );

  /**
   * @brief Calculate the dimensional volume shift
   * @details Computes the dimensional form of the volume shifts given the user defined non-dimensional form.
   * @param[in] numComps The number of components
   * @param[in] componentProperties The compositional model properties
   * @param[out] dimensionalVolumeShift The calculated dimensional volume shifts
   */
  GEOS_FORCE_INLINE
  static void calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                               arraySlice1d< real64 > const & dimensionalVolumeShift );

  /**
   * @brief Calculate the pure coefficients
   * @details Computes the pure coefficients
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           real64 & aCoefficient,
                           real64 & bCoefficient );

  /**
   * @brief Calculate the pure coefficients derivatives
   * @details Computes the pure coefficients derivatives
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   * @param[out] daCoefficient_dp pure coefficient (A) derivative w.r.t. pressure
   * @param[out] dbCoefficient_dp pure coefficient (B) derivative w.r.t. pressure
   * @param[out] daCoefficient_dt pure coefficient (A) derivative w.r.t. temperature
   * @param[out] dbCoefficient_dt pure coefficient (B) derivative w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           real64 & aCoefficient,
                           real64 & bCoefficient,
                           real64 & daCoefficient_dp,
                           real64 & dbCoefficient_dp,
                           real64 & daCoefficient_dt,
                           real64 & dbCoefficient_dt );

  /**
   * @brief Compute the mixture coefficients using pressure, temperature, composition and input
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[out] aPureCoefficient pure coefficient (A)
   * @param[out] bPureCoefficient pure coefficient (B)
   * @param[out] aMixtureCoefficient mixture coefficient (A)
   * @param[out] bMixtureCoefficient mixture coefficient (B)
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              arraySlice1d< real64 > const & aPureCoefficient,
                              arraySlice1d< real64 > const & bPureCoefficient,
                              real64 & aMixtureCoefficient,
                              real64 & bMixtureCoefficient );

  /**
   * @brief Compute the mixture coefficients derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] aMixtureCoefficientDerivs derivatives of mixture coefficient (A)
   * @param[out] bMixtureCoefficientDerivs derivatives of mixture coefficient (B)
   * @note Assumes that pressure and temperature are strictly positive
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              arraySlice1d< real64 const > const & aPureCoefficient,
                              arraySlice1d< real64 const > const & bPureCoefficient,
                              real64 const aMixtureCoefficient,
                              real64 const bMixtureCoefficient,
                              arraySlice1d< real64 > const & aMixtureCoefficientDerivs,
                              arraySlice1d< real64 > const & bMixtureCoefficientDerivs );

  /**
   * @brief Compute the compressibility factor using compositions, BICs, and mixture coefficients
   * @param[in] numComps number of components
   * @param[in] composition composition of the phase
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] compressibilityFactor compressibility factor
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeCompressibilityFactor( integer const numComps,
                                arraySlice1d< real64 const, USD > const & composition,
                                arraySlice2d< real64 const > const & binaryInteractionCoefficients,
                                arraySlice1d< real64 const > const & aPureCoefficient,
                                arraySlice1d< real64 const > const & bPureCoefficient,
                                real64 const & aMixtureCoefficient,
                                real64 const & bMixtureCoefficient,
                                real64 & compressibilityFactor );

  /**
   * @brief Compute the compressibility factor derivatives using mixture coefficients
   * @param[in] numComps number of components
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[in] compressibilityFactor the current compressibility factor
   * @param[in] aMixtureCoefficientDerivs derivatives of mixture coefficient (A)
   * @param[in] bMixtureCoefficientDerivs derivatives of mixture coefficient (B)
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   * @note Assumes that pressure and temperature are strictly positive
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & aMixtureCoefficient,
                                real64 const & bMixtureCoefficient,
                                real64 const & compressibilityFactor,
                                arraySlice1d< real64 const > const & aMixtureCoefficientDerivs,
                                arraySlice1d< real64 const > const & bMixtureCoefficientDerivs,
                                arraySlice1d< real64 > const & compressibilityFactorDerivs );

  /**
   * @brief Compute the log of the fugacity coefficients using compositions, BICs, compressibility factor and mixture coefficients
   * @param[in] numComps number of components
   * @param[in] composition composition of the phase
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[in] compressibilityFactor compressibility factor
   * @param[in] aPureCoefficient pure coefficient (A)
   * @param[in] bPureCoefficient pure coefficient (B)
   * @param[in] aMixtureCoefficient mixture coefficient (A)
   * @param[in] bMixtureCoefficient mixture coefficient (B)
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  arraySlice2d< real64 const > const & binaryInteractionCoefficients,
                                  real64 const & compressibilityFactor,
                                  arraySlice1d< real64 const > const & aPureCoefficient,
                                  arraySlice1d< real64 const > const & bPureCoefficient,
                                  real64 const & aMixtureCoefficient,
                                  real64 const & bMixtureCoefficient,
                                  arraySlice1d< real64 > const & logFugacityCoefficients );

  /**
   * @brief Helper functions solving a cubic equation using trigonometry
   *        m3 * x^3 + m2 * x^2 + m1 *x + m0  = 0
   * @param[in] m3 first coefficient (in front of x^3)
   * @param[in] m2 second coefficient (in front of x^2)
   * @param[in] m1 third coefficient (in front of x)
   * @param[in] m0 fourth coefficient
   * @param[out] roots the roots of the polynomial
   * @param[out] numRoots the number of roots
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void
    solveCubicPolynomial( real64 const & m3,
                          real64 const & m2,
                          real64 const & m1,
                          real64 const & m0,
                          real64 ( &roots )[3],
                          integer & numRoots );

};

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                arraySlice1d< real64 > const & logFugacityCoefficients )
{
  // step 0: allocate the stack memory needed for the update
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > aPureCoefficient( numComps );
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > bPureCoefficient( numComps );
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  real64 compressibilityFactor = 0.0;

  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  // step 1: compute the mixture coefficients aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // step 2: use mixture coefficients to update the compressibility factor
  computeCompressibilityFactor( numComps, // number of components
                                composition, // cell input
                                binaryInteractionCoefficients, // user input
                                aPureCoefficient, // computed by computeMixtureCoefficients
                                bPureCoefficient,
                                aMixtureCoefficient,
                                bMixtureCoefficient,
                                compressibilityFactor ); // output

  // step 3: use mixture coefficients and compressibility factor to update fugacity coefficients
  computeLogFugacityCoefficients( numComps, // number of components
                                  composition, // cell input
                                  binaryInteractionCoefficients, // user input
                                  compressibilityFactor, // computed by computeCompressibilityFactor
                                  aPureCoefficient, // computed by computeMixtureCoefficients
                                  bPureCoefficient,
                                  aMixtureCoefficient,
                                  bMixtureCoefficient,
                                  logFugacityCoefficients ); // output
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                arraySlice1d< real64 const > const & logFugacityCoefficients,
                                arraySlice2d< real64 > const & logFugacityCoefficientDerivs )
{
  integer constexpr numMaxComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  integer constexpr numMaxDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;
  integer const numDofs = 2 + numComps;

  GEOS_UNUSED_VAR( logFugacityCoefficients );

  stackArray1d< real64, numMaxComps > aPureCoefficient( numComps );
  stackArray1d< real64, numMaxComps > bPureCoefficient( numComps );
  stackArray2d< real64, 2*numMaxComps > aPureCoefficientDerivs( numComps, 2 );
  stackArray2d< real64, 2*numMaxComps > bPureCoefficientDerivs( numComps, 2 );
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  real64 compressibilityFactor = 0.0;
  stackArray1d< real64, numMaxDofs > aMixtureCoefficientDerivs( numDofs );
  stackArray1d< real64, numMaxDofs > bMixtureCoefficientDerivs( numDofs );
  stackArray1d< real64, numMaxDofs > compressibilityFactorDerivs( numDofs );

  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  // 1.1: Compute the pure and mixture coefficients
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // 1.2: Compute pure coefficient derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    computePureCoefficients( ic,
                             pressure,
                             temperature,
                             componentProperties,
                             aPureCoefficient[ic],
                             bPureCoefficient[ic],
                             aPureCoefficientDerivs( ic, Deriv::dP ),
                             bPureCoefficientDerivs( ic, Deriv::dP ),
                             aPureCoefficientDerivs( ic, Deriv::dT ),
                             bPureCoefficientDerivs( ic, Deriv::dT ));
  }

  // 1.3: Compute mixture coefficient derivatives
  computeMixtureCoefficients( numComps,
                              pressure,
                              temperature,
                              composition,
                              componentProperties,
                              aPureCoefficient,
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient,
                              aMixtureCoefficientDerivs,
                              bMixtureCoefficientDerivs );

  // 2.1: Update the compressibility factor
  computeCompressibilityFactor( numComps, // number of components
                                composition, // cell input
                                binaryInteractionCoefficients, // user input
                                aPureCoefficient, // computed by computeMixtureCoefficients
                                bPureCoefficient,
                                aMixtureCoefficient,
                                bMixtureCoefficient,
                                compressibilityFactor ); // output
  // 2.2: Update the compressibility factor derivatives
  computeCompressibilityFactor( numComps,
                                aMixtureCoefficient,
                                bMixtureCoefficient,
                                compressibilityFactor,
                                aMixtureCoefficientDerivs,
                                bMixtureCoefficientDerivs,
                                compressibilityFactorDerivs );

  // 3. Calculate derivatives of the logarithm of the fugacity coefficients
  stackArray1d< real64, numMaxComps > ki( numComps );
  stackArray2d< real64, numMaxComps * numMaxDofs > dki( numComps, numDofs );

  // ki
  for( integer ic = 0; ic < numComps; ++ic )
  {
    ki[ic] = 0.0;
    dki( ic, Deriv::dP ) = 0.0;
    dki( ic, Deriv::dT ) = 0.0;
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const aCoeffI = sqrt( aPureCoefficient[ic] );
      real64 const aCoeffJ = sqrt( aPureCoefficient[jc] );
      real64 const kij = ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) * aCoeffI * aCoeffJ;
      ki[ic] += composition[jc] * kij;
      dki( ic, Deriv::dC + jc ) = kij;
      dki( ic, Deriv::dP ) += 0.5 * composition[jc] * kij * ( aPureCoefficientDerivs( ic, Deriv::dP )/aPureCoefficient[ic] + aPureCoefficientDerivs( jc, Deriv::dP )/aPureCoefficient[jc] );
      dki( ic, Deriv::dT ) += 0.5 * composition[jc] * kij * ( aPureCoefficientDerivs( ic, Deriv::dT )/aPureCoefficient[ic] + aPureCoefficientDerivs( jc, Deriv::dT )/aPureCoefficient[jc] );
    }
  }

  auto const calculateDerivatives = [&]( integer const kc ){
    real64 const E = log( compressibilityFactor + EOS_TYPE::delta1 * bMixtureCoefficient )
                     - log( compressibilityFactor + EOS_TYPE::delta2 * bMixtureCoefficient );

    real64 const dE_dX = (compressibilityFactorDerivs[kc] + EOS_TYPE::delta1*bMixtureCoefficientDerivs[kc])/( compressibilityFactor + EOS_TYPE::delta1 * bMixtureCoefficient )
                         -(compressibilityFactorDerivs[kc] + EOS_TYPE::delta2*bMixtureCoefficientDerivs[kc])/( compressibilityFactor + EOS_TYPE::delta2 * bMixtureCoefficient );

    //real64 const F = log( compressibilityFactor - bMixtureCoefficient );
    real64 const dF_dX = (compressibilityFactorDerivs[kc] - bMixtureCoefficientDerivs[kc])/(compressibilityFactor - bMixtureCoefficient);

    real64 const G = 1.0 / ( ( EOS_TYPE::delta1 - EOS_TYPE::delta2 ) * bMixtureCoefficient );
    real64 const dG_dX = -G * bMixtureCoefficientDerivs[kc] / bMixtureCoefficient;

    real64 const A = aMixtureCoefficient;
    real64 const dA_dX = aMixtureCoefficientDerivs[kc];

    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const B = bPureCoefficient[ic] / bMixtureCoefficient;
      real64 dB_dX = -B*bMixtureCoefficientDerivs[kc] / bMixtureCoefficient;
      if( kc < Deriv::dC )
      {
        dB_dX += bPureCoefficientDerivs( ic, kc ) / bMixtureCoefficient;
      }

      // lnPhi = ( compressibilityFactor - 1 ) * B - F - G * ( 2 * ki[ic] - A * B ) * E;
      logFugacityCoefficientDerivs( ic, kc ) =
        compressibilityFactorDerivs[kc]*B + ( compressibilityFactor - 1 ) * dB_dX
        - dF_dX
        - dG_dX * ( 2 * ki[ic] - A * B ) * E
        - G * ( 2 * dki( ic, kc ) - dA_dX * B - A * dB_dX ) * E
        - G * ( 2 * ki[ic] - A * B ) * dE_dX;
    }
  };

  calculateDerivatives( Deriv::dP );
  calculateDerivatives( Deriv::dT );

  for( integer jc = 0; jc < numComps; ++jc )
  {
    calculateDerivatives( Deriv::dC+jc );
  }
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 & compressibilityFactor,
                              arraySlice1d< real64 > const & compressibilityFactorDerivs )
{
  // step 0: allocate the stack memory needed for the update
  integer constexpr numMaxComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  integer constexpr numMaxDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;
  integer const numDofs = 2 + numComps;

  stackArray1d< real64, numMaxComps > aPureCoefficient( numComps );
  stackArray1d< real64, numMaxComps > bPureCoefficient( numComps );
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  stackArray1d< real64, numMaxDofs > aMixtureCoefficientDerivs( numDofs );
  stackArray1d< real64, numMaxDofs > bMixtureCoefficientDerivs( numDofs );

  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  // step 1: compute the mixture coefficients aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient
  // 1.1: Compute the pure and mixture coefficients
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // 1.2: Compute mixture coefficient derivatives
  computeMixtureCoefficients( numComps,
                              pressure,
                              temperature,
                              composition,
                              componentProperties,
                              aPureCoefficient,
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient,
                              aMixtureCoefficientDerivs,
                              bMixtureCoefficientDerivs );

  // 2.1: Update the compressibility factor
  computeCompressibilityFactor( numComps, // number of components
                                composition, // cell input
                                binaryInteractionCoefficients, // user input
                                aPureCoefficient, // computed by computeMixtureCoefficients
                                bPureCoefficient,
                                aMixtureCoefficient,
                                bMixtureCoefficient,
                                compressibilityFactor ); // output

  // 2.2: Update the compressibility factor derivatives
  computeCompressibilityFactor( numComps,
                                aMixtureCoefficient,
                                bMixtureCoefficient,
                                compressibilityFactor,
                                aMixtureCoefficientDerivs,
                                bMixtureCoefficientDerivs,
                                compressibilityFactorDerivs );
}

template< typename EOS_TYPE >
void
CubicEOSPhaseModel< EOS_TYPE >::
calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                 arraySlice1d< real64 > const & dimensionalVolumeShift )
{
  integer const numComps = componentProperties.getNumberOfComponents();
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const Vs = componentProperties.getComponentVolumeShift()[ic];
    real64 const Pc = componentProperties.getComponentCriticalPressure()[ic];
    real64 const Tc = componentProperties.getComponentCriticalTemperature()[ic];
    real64 constexpr omegaB = EOS_TYPE::omegaB;
    dimensionalVolumeShift[ic] = constants::gasConstant * Vs * omegaB * Tc / Pc;
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 & aCoefficient,
                         real64 & bCoefficient )
{
  real64 daCoefficient_dp = 0.0;
  real64 dbCoefficient_dp = 0.0;
  real64 daCoefficient_dt = 0.0;
  real64 dbCoefficient_dt = 0.0;
  computePureCoefficients( ic,
                           pressure,
                           temperature,
                           componentProperties,
                           aCoefficient,
                           bCoefficient,
                           daCoefficient_dp,
                           dbCoefficient_dp,
                           daCoefficient_dt,
                           dbCoefficient_dt );
}


template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 & aCoefficient,
                         real64 & bCoefficient,
                         real64 & daCoefficient_dp,
                         real64 & dbCoefficient_dp,
                         real64 & daCoefficient_dt,
                         real64 & dbCoefficient_dt )
{
  arraySlice1d< real64 const > const & criticalPressure = componentProperties.m_componentCriticalPressure;
  arraySlice1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
  arraySlice1d< real64 const > const & acentricFactor = componentProperties.m_componentAcentricFactor;

  real64 const m = EOS_TYPE::evaluate( acentricFactor[ic] );
  real64 const pr = pressure / criticalPressure[ic];
  real64 const tr = temperature / criticalTemperature[ic];

  real64 const sqrtTr = sqrt( tr );
  real64 const mt = 1.0 + m * (1.0 - sqrtTr);

  aCoefficient = EOS_TYPE::omegaA * pr / (tr*tr) * mt * mt;
  bCoefficient = EOS_TYPE::omegaB * pr / tr;

  daCoefficient_dp = aCoefficient / pressure;
  dbCoefficient_dp = bCoefficient / pressure;

  daCoefficient_dt = -aCoefficient * (2.0/temperature + m/(mt * sqrtTr * criticalTemperature[ic]));
  dbCoefficient_dt = -bCoefficient / temperature;
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d< real64 const, USD > const & composition,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            arraySlice1d< real64 > const & aPureCoefficient,
                            arraySlice1d< real64 > const & bPureCoefficient,
                            real64 & aMixtureCoefficient,
                            real64 & bMixtureCoefficient )
{
  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  // mixture coefficients
  for( integer ic = 0; ic < numComps; ++ic )
  {
    computePureCoefficients( ic, pressure, temperature, componentProperties, aPureCoefficient[ic], bPureCoefficient[ic] );
  }

  aMixtureCoefficient = 0.0;
  bMixtureCoefficient = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      aMixtureCoefficient += composition[ic] * composition[jc] * ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
    }
    bMixtureCoefficient += composition[ic] * bPureCoefficient[ic];
  }
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d< real64 const, USD > const & composition,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            arraySlice1d< real64 const > const & aPureCoefficient,
                            arraySlice1d< real64 const > const & bPureCoefficient,
                            real64 const aMixtureCoefficient,
                            real64 const bMixtureCoefficient,
                            arraySlice1d< real64 > const & aMixtureCoefficientDerivs,
                            arraySlice1d< real64 > const & bMixtureCoefficientDerivs )
{
  arraySlice2d< real64 const > const & binaryInteractionCoefficients = componentProperties.m_componentBinaryCoeff;

  // Calculate pressure derivatives
  aMixtureCoefficientDerivs[Deriv::dP] = aMixtureCoefficient / pressure;
  bMixtureCoefficientDerivs[Deriv::dP] = bMixtureCoefficient / pressure;

  // Calculate temperature derivatives
  real64 aCoefficient = 0.0;
  real64 bCoefficient = 0.0;
  real64 dummy = 0.0;
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > daPureCoefficient_dt( numComps );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    computePureCoefficients( ic, pressure, temperature, componentProperties,
                             aCoefficient, bCoefficient, dummy, dummy, daPureCoefficient_dt[ic], dummy );
  }
  aMixtureCoefficientDerivs[Deriv::dT] = 0.0;
  bMixtureCoefficientDerivs[Deriv::dT] = -bMixtureCoefficient / temperature;

  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const coeff = composition[ic] * composition[jc] * ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) / sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      aMixtureCoefficientDerivs[Deriv::dT] += 0.5 * coeff * (daPureCoefficient_dt[ic]*aPureCoefficient[jc] + daPureCoefficient_dt[jc]*aPureCoefficient[ic]);
    }
  }

  // Calculate composition derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    aMixtureCoefficientDerivs[Deriv::dC+ic] = 0.0;
    bMixtureCoefficientDerivs[Deriv::dC+ic] = 0.0;
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const coeff = ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      aMixtureCoefficientDerivs[Deriv::dC+ic] += coeff * composition[jc];
      aMixtureCoefficientDerivs[Deriv::dC+jc] += coeff * composition[ic];
    }
    bMixtureCoefficientDerivs[Deriv::dC+ic] = bPureCoefficient[ic];
  }
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              arraySlice1d< real64 const, USD > const & composition,
                              arraySlice2d< real64 const > const & binaryInteractionCoefficients,
                              arraySlice1d< real64 const > const & aPureCoefficient,
                              arraySlice1d< real64 const > const & bPureCoefficient,
                              real64 const & aMixtureCoefficient,
                              real64 const & bMixtureCoefficient,
                              real64 & compressibilityFactor )
{
  // a Z3 + b Z2 + cZ + d = 0
  real64 const a = 1.0;
  real64 const b = ( EOS_TYPE::delta1 + EOS_TYPE::delta2 - 1.0 ) * bMixtureCoefficient - 1.0;
  real64 const c = aMixtureCoefficient + EOS_TYPE::delta1 * EOS_TYPE::delta2 * bMixtureCoefficient * bMixtureCoefficient
                   - ( EOS_TYPE::delta1 + EOS_TYPE::delta2 ) * bMixtureCoefficient * ( bMixtureCoefficient + 1.0 );
  real64 const d = -( aMixtureCoefficient * bMixtureCoefficient
                      + EOS_TYPE::delta1 * EOS_TYPE::delta2 * bMixtureCoefficient * bMixtureCoefficient * ( bMixtureCoefficient + 1.0 ) );

  real64 roots[3]{};
  integer numRoots = 0;
  solveCubicPolynomial( a, b, c, d, roots, numRoots );

  if( numRoots == 1 )
  {
    compressibilityFactor = roots[0];
  }
  else
  {
    real64 zMin =  LvArray::NumericLimits< real64 >::max;
    real64 zMax = -LvArray::NumericLimits< real64 >::max;

    for( integer i = 0; i < 3; ++i )
    {
      // skip unphysical roots
      if( roots[i] > bMixtureCoefficient )
      {
        // choose the root according to Gibbs' free energy minimization
        if( zMin > roots[i] )
        {
          zMin = roots[i];
        }
        if( zMax < roots[i] )
        {
          zMax = roots[i];
        }
      }
    }

    stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > logFugacityCoefficientsMax( numComps );
    stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > logFugacityCoefficientsMin( numComps );
    computeLogFugacityCoefficients( numComps, composition, binaryInteractionCoefficients, zMin,
                                    aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient,
                                    logFugacityCoefficientsMin.toSlice() );
    computeLogFugacityCoefficients( numComps, composition, binaryInteractionCoefficients, zMax,
                                    aPureCoefficient, bPureCoefficient, aMixtureCoefficient, bMixtureCoefficient,
                                    logFugacityCoefficientsMax.toSlice() );

    real64 dG = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dG += composition[ic] * ( logFugacityCoefficientsMin[ic] - logFugacityCoefficientsMax[ic] );
    }
    compressibilityFactor = ( dG < 0 ) ? zMin : zMax;
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              real64 const & aMixtureCoefficient,
                              real64 const & bMixtureCoefficient,
                              real64 const & compressibilityFactor,
                              arraySlice1d< real64 const > const & aMixtureCoefficientDerivs,
                              arraySlice1d< real64 const > const & bMixtureCoefficientDerivs,
                              arraySlice1d< real64 > const & compressibilityFactorDerivs )
{
  real64 constexpr d1pd2 = EOS_TYPE::delta1 + EOS_TYPE::delta2;
  real64 constexpr d1xd2 = EOS_TYPE::delta1 * EOS_TYPE::delta2;

  real64 constexpr a = 1.0;
  real64 const b = ( d1pd2 - 1.0 ) * bMixtureCoefficient - 1.0;
  real64 const c = aMixtureCoefficient + d1xd2 * bMixtureCoefficient * bMixtureCoefficient
                   - d1pd2 * bMixtureCoefficient * ( bMixtureCoefficient + 1.0 );

  // Implicit differentiation scale
  real64 const denominator = (3.0*a*compressibilityFactor + 2.0*b)*compressibilityFactor + c;
  real64 const scalingFactor = LvArray::math::abs( denominator ) < MultiFluidConstants::epsilon ? 0.0 : -1.0 / denominator;

  integer const numDofs = numComps + 2;

  for( integer kc = 0; kc < numDofs; ++kc )
  {
    // Given derivative of the mixture parameters a and b w.r.t. variable X, calculate the derivative of the
    // compressibility factor (z-factor) w.r.t. X
    real64 const da_dX = aMixtureCoefficientDerivs[kc];
    real64 const db_dX = bMixtureCoefficientDerivs[kc];
    // a Z3 + b Z2 + cZ + d = 0
    // Derivatives for a,b,c,d
    real64 const dbdx = ( d1pd2 - 1.0 ) * db_dX;
    real64 const dcdx = da_dX + ( 2.0*(d1xd2-d1pd2) * bMixtureCoefficient - d1pd2 )*db_dX;
    real64 const dddx = -(aMixtureCoefficient*db_dX + da_dX*bMixtureCoefficient
                          + d1xd2*((3.0*bMixtureCoefficient+2.0)*bMixtureCoefficient*db_dX));
    compressibilityFactorDerivs[kc] = (((dbdx*compressibilityFactor) + dcdx)*compressibilityFactor + dddx) * scalingFactor;
  }
}

template< typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                arraySlice1d< real64 const, USD > const & composition,
                                arraySlice2d< real64 const > const & binaryInteractionCoefficients,
                                real64 const & compressibilityFactor,
                                arraySlice1d< real64 const > const & aPureCoefficient,
                                arraySlice1d< real64 const > const & bPureCoefficient,
                                real64 const & aMixtureCoefficient,
                                real64 const & bMixtureCoefficient,
                                arraySlice1d< real64 > const & logFugacityCoefficients )
{
  stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > ki( numComps );

  // ki
  for( integer ic = 0; ic < numComps; ++ic )
  {
    ki[ic] = 0.0;
    for( integer jc = 0; jc < numComps; ++jc )
    {
      ki[ic] += composition[jc] * ( 1.0 - binaryInteractionCoefficients( ic, jc ) ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
    }
  }

  // E
  real64 const expE = ( compressibilityFactor + EOS_TYPE::delta1 * bMixtureCoefficient ) /
                      ( compressibilityFactor + EOS_TYPE::delta2 * bMixtureCoefficient );
  real64 const expF = compressibilityFactor - bMixtureCoefficient;
  GEOS_ERROR_IF( expE < MultiFluidConstants::epsilon || expF < MultiFluidConstants::epsilon,
                 GEOS_FMT( "Cubic EOS failed with exp(E)={} and exp(F)={}", expE, expF ));
  real64 const E = log( expE );
  real64 const F = log( expF );
  real64 const G = 1.0 / ( ( EOS_TYPE::delta1 - EOS_TYPE::delta2 ) * bMixtureCoefficient );
  real64 const A = aMixtureCoefficient;

  // ln phi
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const B = bPureCoefficient[ic] / bMixtureCoefficient;
    logFugacityCoefficients[ic] = ( compressibilityFactor - 1 ) * B - F - G * ( 2 * ki[ic] - A * B ) * E;
  }
}

template< typename EOS_TYPE >
GEOS_HOST_DEVICE
void
CubicEOSPhaseModel< EOS_TYPE >::
solveCubicPolynomial( real64 const & m3,
                      real64 const & m2,
                      real64 const & m1,
                      real64 const & m0,
                      real64 (& roots)[3],
                      integer & numRoots )
{
  // cubic equation : m3 * x^3 + m2 * x^2 + m1 *x + m0  = 0
  real64 const a1 = m2 / m3;
  real64 const a2 = m1 / m3;
  real64 const a3 = m0 / m3;

  real64 const q = ( a1 * a1 - 3 * a2 ) / 9;
  real64 const r = ( 2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3 ) / 54;
  real64 const qCubed = q * q * q;
  real64 const d = qCubed - r * r;

  // three real roots
  if( d >= 0 )
  {
    real64 const theta = acos( r / sqrt( qCubed ) );
    real64 const qSqrt = sqrt( q );
    roots[0] = -2 * qSqrt * cos( theta / 3 ) - a1 / 3;
    roots[1] = -2 * qSqrt * cos( ( theta + 2 * constants::pi ) / 3 ) - a1 / 3;
    roots[2] = -2 * qSqrt * cos( ( theta + 4 * constants::pi ) / 3 ) - a1 / 3;
    numRoots = 3;
  }
  // one real root
  else
  {
    real64 e = pow( sqrt( -d ) + LvArray::math::abs( r ), 1. / 3. );
    if( r > 0 )
    {
      e = -e;
    }
    roots[0] = ( e + q / e ) - a1 / 3.;
    numRoots = 1;
  }
}

using CubicEOSPR = CubicEOSPhaseModel< PengRobinsonEOS >;
using CubicEOSSRK = CubicEOSPhaseModel< SoaveRedlichKwongEOS >;

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_CUBICEOSPHASEMODEL_HPP_
