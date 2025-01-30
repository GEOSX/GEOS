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
 * @file SoreideWhitsonEOSModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSMODEL_HPP_

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"
#include "CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

enum class SoreideWhitsonPhaseType : integer
{
  Vapour,
  Aqueous
};

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
struct SoreideWhitsonEOSModel
{
  static constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  using Deriv = multifluid::DerivativeOffset;

  /**
   * @brief Main entry point of the cubic EOS model
   * @details Computes the logarithm of the fugacity coefficients
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] logFugacityCoefficients log of the fugacity coefficients
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  real64 const & pressure,
                                  real64 const & temperature,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  real64 const & salinity,
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
  static void
  computeLogFugacityCoefficients( integer const numComps,
                                  real64 const & pressure,
                                  real64 const & temperature,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  real64 const & salinity,
                                  arraySlice1d< real64 const > const & logFugacityCoefficients,
                                  arraySlice2d< real64 > const & logFugacityCoefficientDerivs );

  /**
   * @brief Calculate the pure coefficients
   * @details Computes the pure coefficients
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   */
  GEOS_HOST_DEVICE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           real64 const & salinity,
                           real64 & aCoefficient,
                           real64 & bCoefficient );

  /**
   * @brief Calculate the pure coefficients derivatives
   * @details Computes the pure coefficients derivatives
   * @param[in] ic Component index
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[out] aCoefficient pure coefficient (A)
   * @param[out] bCoefficient pure coefficient (B)
   * @param[out] daCoefficient_dp pure coefficient (A) derivative w.r.t. pressure
   * @param[out] dbCoefficient_dp pure coefficient (B) derivative w.r.t. pressure
   * @param[out] daCoefficient_dt pure coefficient (A) derivative w.r.t. temperature
   * @param[out] dbCoefficient_dt pure coefficient (B) derivative w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  static void
  computePureCoefficients( integer const ic,
                           real64 const & pressure,
                           real64 const & temperature,
                           ComponentProperties::KernelWrapper const & componentProperties,
                           real64 const & salinity,
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
   * @param[in] salinity salinity
   * @param[out] aPureCoefficient pure coefficient (A)
   * @param[out] bPureCoefficient pure coefficient (B)
   * @param[out] aMixtureCoefficient mixture coefficient (A)
   * @param[out] bMixtureCoefficient mixture coefficient (B)
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 const & salinity,
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
   * @param[in] salinity salinity
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
  static void
  computeMixtureCoefficients( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 const & salinity,
                              arraySlice1d< real64 const > const & aPureCoefficient,
                              arraySlice1d< real64 const > const & bPureCoefficient,
                              real64 const aMixtureCoefficient,
                              real64 const bMixtureCoefficient,
                              arraySlice1d< real64 > const & aMixtureCoefficientDerivs,
                              arraySlice1d< real64 > const & bMixtureCoefficientDerivs );

  /**
   * @brief Compute compressibility factor
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] compressibilityFactor the current compressibility factor
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 const & salinity,
                                real64 & compressibilityFactor );

  /**
   * @brief Compute compressibility factor derivatives
   * @details Computes the compressibility factor (z-factor) for the cubic EOS model including derivatives
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the phase
   * @param[in] componentProperties The compositional component properties
   * @param[in] compressibilityFactor the current compressibility factor
   * @param[out] compressibilityFactorDerivs derivatives of the compressibility factor
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  static void
  computeCompressibilityFactor( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 const & salinity,
                                real64 const & compressibilityFactor,
                                arraySlice1d< real64 > const & compressibilityFactorDerivs );

  /**
   * @brief Get the binary interaction coefficient between two components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] componentProperties The compositional component properties
   * @param[in] salinity salinity
   * @param[in] i index of first component
   * @param[in] j index of second component
   * @param[out] kij the binary interaction coefficient
   * @param[out] dkij_dT derivative of the binary interaction coefficient w.r.t. temperature
   */
  GEOS_HOST_DEVICE
  static void getBinaryInteractionCiefficient( real64 const & pressure,
                                               real64 const & temperature,
                                               ComponentProperties::KernelWrapper const & componentProperties,
                                               real64 const & salinity,
                                               integer const i,
                                               integer const j,
                                               real64 & kij,
                                               real64 & dkij_dT );

  /**
   * @brief Safe pow function
   * @param[in] a The value
   * @param[in] b The exponent
   * @return returns a**b if a is positive else 0
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 power( real64 const a, real64 const b )
  {
    return a < MultiFluidConstants::minForSpeciesPresence ? 0.0 : LvArray::math::exp( b * LvArray::math::log( a ) );
  }

  /**
   * @brief Check if a component has a specific type
   * @param[in] type The type of the component
   * @param[in] targetType The target type to check
   * @return returns true if the type matches the target type
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static bool isType( integer const type, ComponentProperties::ComponentType const targetType )
  {
    return (type == static_cast< int >(targetType));
  }
};

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 const & salinity,
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
                           salinity,
                           aCoefficient,
                           bCoefficient,
                           daCoefficient_dp,
                           dbCoefficient_dp,
                           daCoefficient_dt,
                           dbCoefficient_dt );
}

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computePureCoefficients( integer const ic,
                         real64 const & pressure,
                         real64 const & temperature,
                         ComponentProperties::KernelWrapper const & componentProperties,
                         real64 const & salinity,
                         real64 & aCoefficient,
                         real64 & bCoefficient,
                         real64 & daCoefficient_dp,
                         real64 & dbCoefficient_dp,
                         real64 & daCoefficient_dt,
                         real64 & dbCoefficient_dt )
{
  if constexpr (PHASE_TYPE == SoreideWhitsonPhaseType::Aqueous)
  {
    arraySlice1d< real64 const > const & criticalPressure = componentProperties.m_componentCriticalPressure;
    arraySlice1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;

    real64 const pr = pressure / criticalPressure[ic];
    real64 const tr = temperature / criticalTemperature[ic];

    real64 const trToMinus2 = 1.0/(tr*tr);
    real64 const trToMinus3 = trToMinus2/tr;

    real64 const csw = power( salinity, 1.1 );

    real64 const sqrtAlpha = 1.0 - 0.4530*(1.0 - tr*(1.0 - 0.0103*csw)) + 0.0034*(trToMinus3 - 1.0);
    real64 const dSqrtAlpha_dtr = 0.4530*(1.0 - 0.0103*csw) - 0.0102*trToMinus3/tr;
    real64 const alpha = sqrtAlpha * sqrtAlpha;
    real64 const dAlpha_dt = 2.0 * sqrtAlpha * dSqrtAlpha_dtr / criticalTemperature[ic];

    aCoefficient = EOS_TYPE::omegaA * pr * trToMinus2 * alpha;
    bCoefficient = EOS_TYPE::omegaB * pr / tr;

    daCoefficient_dp = aCoefficient / pressure;
    dbCoefficient_dp = bCoefficient / pressure;

    daCoefficient_dt = EOS_TYPE::omegaA * pr * (-2.0*trToMinus3 * alpha / criticalTemperature[ic] + trToMinus2 * dAlpha_dt);
    dbCoefficient_dt = -bCoefficient / temperature;
  }
  else
  {
    // Vapour phase usese normal cubic EOS
    CubicEOSPhaseModel< EOS_TYPE >::computePureCoefficients( ic,
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
}


template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
getBinaryInteractionCiefficient( real64 const & pressure,
                                 real64 const & temperature,
                                 ComponentProperties::KernelWrapper const & componentProperties,
                                 real64 const & salinity,
                                 integer const ic,
                                 integer const jc,
                                 real64 & kij,
                                 real64 & dkij_dT )
{
  GEOS_UNUSED_VAR( pressure );

  // Initialise default values
  kij = componentProperties.m_componentBinaryCoeff( ic, jc );
  dkij_dT = 0.0;

  if constexpr (PHASE_TYPE == SoreideWhitsonPhaseType::Aqueous)
  {
    integer i = ic;
    integer j = jc;
    integer type_i = componentProperties.m_componentType[i];
    integer type_j = componentProperties.m_componentType[j];

    if ( isType( type_i, ComponentProperties::ComponentType::Water ) && !isType( type_j, ComponentProperties::ComponentType::Water ) )
    {
      integer t = type_i;
      type_i = type_j;
      type_j = t;
      t = i;
      i = j;
      j = t;
    }

    if( isType( type_j, ComponentProperties::ComponentType::Water ) )
    {
      real64 const Tci = componentProperties.m_componentCriticalTemperature( i );
      real64 const dTr_dt = 1.0 / Tci;
      real64 const Tr = temperature * dTr_dt;

      if( isType( type_i, ComponentProperties::ComponentType::CarbonDioxide ))
      {
        // We have options here:
        // Original Soreide & Whitson (1992); Yan et al. (2011); Chabab et al. (2019)
        // Equation (14) Soreide & Whitson (1992)
        real64 constexpr A0 = -0.31092;
        real64 constexpr A1 = 0.23580;
        real64 constexpr A2 = -21.2566;
        real64 const a0 = 1.0 + 0.15587 * power( salinity, 0.7505 );
        real64 const a1 = 1.0 + 0.17837 * power( salinity, 0.979 );
        real64 const a2 = LvArray::math::exp( -6.7222*Tr - salinity );
        kij = A0*a0 + A1*a1 + A2*a2;
        dkij_dT = -6.7222*A2*a2 * dTr_dt;
      }
      else if( isType( type_i, ComponentProperties::ComponentType::HydrogenSulphide ))
      {
        // Equation (15) Soreide & Whitson (1992)
        real64 constexpr A0 = -0.20441;
        real64 constexpr A1 = 0.23426;
        kij = A0 + A1*Tr;
        dkij_dT = A1*dTr_dt;
      }
      else if( isType( type_i, ComponentProperties::ComponentType::Nitrogen ))
      {
        // Equation (13) Soreide & Whitson (1992)
        real64 constexpr A0 = -1.70235;
        real64 constexpr A1 = 0.44338;
        real64 const csw = power( salinity, 0.75 );
        real64 const a0 = 1.0 + 0.25587*csw;
        real64 const a1 = 1.0 + 0.08126*csw;
        kij = A0*a0 + A1*a1*Tr;
        dkij_dT = A1*a1*dTr_dt;
      }
      else
      {
        // Hydrocarbon-water interaction
        real64 const omega = componentProperties.m_componentAcentricFactor( i );

        // Table 2 from Soreide & Whitson (1992)
        // See also equations (11) and (12)
        real64 const A0 = 1.1120 - 1.7369 * power( omega, -0.1 );
        real64 const A1 = 1.1001 + 0.8360 * omega;
        real64 const A2 = -0.15742 - 1.0988 * omega;
        real64 const a0 = 1.0 + 0.017407*salinity;
        real64 const a1 = 1.0 + 0.033516*salinity;
        real64 const a2 = 1.0 + 0.011478*salinity;
        kij = A0*a0 + A1*a1*Tr + A2*a2*Tr*Tr;
        real64 const dkij_dTr = A1*a1 + 2.0*A2*a2*Tr;
        dkij_dT = dkij_dTr * dTr_dt;
      }
    }
  }
}

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d< real64 const, USD > const & composition,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            real64 const & salinity,
                            arraySlice1d< real64 > const & aPureCoefficient,
                            arraySlice1d< real64 > const & bPureCoefficient,
                            real64 & aMixtureCoefficient,
                            real64 & bMixtureCoefficient )
{
  if constexpr (PHASE_TYPE == SoreideWhitsonPhaseType::Aqueous)
  {
    // Pure component coefficients
    for( integer ic = 0; ic < numComps; ++ic )
    {
      computePureCoefficients( ic, pressure, temperature, componentProperties, salinity, aPureCoefficient[ic], bPureCoefficient[ic] );
    }

    // Mixture coefficients
    aMixtureCoefficient = 0.0;
    for( integer jc = 0; jc < numComps; ++jc )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        real64 kij = 0.0;
        real64 dkij_dT = 0.0;
        getBinaryInteractionCiefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij, dkij_dT );
        aMixtureCoefficient += composition[ic] * composition[jc] * ( 1.0 - kij ) * sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
      }
    }

    bMixtureCoefficient = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      bMixtureCoefficient += composition[ic] * bPureCoefficient[ic];
    }
  }
  else
  {
    // Vapour phase usese normal cubic EOS
    CubicEOSPhaseModel< EOS_TYPE >::computeMixtureCoefficients( numComps,
                                                                pressure,
                                                                temperature,
                                                                composition,
                                                                componentProperties,
                                                                aPureCoefficient,
                                                                bPureCoefficient,
                                                                aMixtureCoefficient,
                                                                bMixtureCoefficient );
  }
}

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computeMixtureCoefficients( integer const numComps,
                            real64 const & pressure,
                            real64 const & temperature,
                            arraySlice1d< real64 const, USD > const & composition,
                            ComponentProperties::KernelWrapper const & componentProperties,
                            real64 const & salinity,
                            arraySlice1d< real64 const > const & aPureCoefficient,
                            arraySlice1d< real64 const > const & bPureCoefficient,
                            real64 const aMixtureCoefficient,
                            real64 const bMixtureCoefficient,
                            arraySlice1d< real64 > const & aMixtureCoefficientDerivs,
                            arraySlice1d< real64 > const & bMixtureCoefficientDerivs )
{
  if constexpr (PHASE_TYPE == SoreideWhitsonPhaseType::Aqueous)
  {
    // a parameter derivatives
    aMixtureCoefficientDerivs[Deriv::dP] = aMixtureCoefficient / pressure;

    real64 aCoefficient = 0.0;
    real64 bCoefficient = 0.0;
    real64 dummy = 0.0;
    stackArray1d< real64, maxNumComps > daPureCoefficient_dt( numComps );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      computePureCoefficients( ic, pressure, temperature, componentProperties, salinity,
                               aCoefficient, bCoefficient, dummy, dummy, daPureCoefficient_dt[ic], dummy );
    }

    aMixtureCoefficientDerivs[Deriv::dT] = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      aMixtureCoefficientDerivs[Deriv::dC+ic] = 0.0;
    }
    for( integer ic = 0; ic < numComps; ++ic )
    {
      for( integer jc = 0; jc < numComps; ++jc )
      {
        real64 kij = 0.0;
        real64 dkij_dT = 0.0;
        getBinaryInteractionCiefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij, dkij_dT );

        real64 const aij = sqrt( aPureCoefficient[ic] * aPureCoefficient[jc] );
        real64 const coeff = composition[ic] * composition[jc] * ( 1.0 - kij ) / aij;
        aMixtureCoefficientDerivs[Deriv::dT] += 0.5 * coeff * (daPureCoefficient_dt[ic]*aPureCoefficient[jc] + daPureCoefficient_dt[jc]*aPureCoefficient[ic]);
        aMixtureCoefficientDerivs[Deriv::dT] -= composition[ic] * composition[jc] * dkij_dT * aij;

        aMixtureCoefficientDerivs[Deriv::dC+ic] += composition[jc] * ( 1.0 - kij ) * aij;
        aMixtureCoefficientDerivs[Deriv::dC+jc] += composition[ic] * ( 1.0 - kij ) * aij;
      }
    }

    // b parameter derivatives
    bMixtureCoefficientDerivs[Deriv::dP] = bMixtureCoefficient / pressure;
    bMixtureCoefficientDerivs[Deriv::dT] = -bMixtureCoefficient / temperature;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      bMixtureCoefficientDerivs[Deriv::dC+ic] = bPureCoefficient[ic];
    }
  }
  else
  {
    // Vapour phase usese normal cubic EOS
    CubicEOSPhaseModel< EOS_TYPE >::computeMixtureCoefficients( numComps,
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
  }
}

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 const & salinity,
                              real64 & compressibilityFactor )
{
  stackArray2d< real64, 2*maxNumComps > tempData( 2, numComps );
  arraySlice1d< real64 > aPureCoefficient = tempData[0];
  arraySlice1d< real64 > bPureCoefficient = tempData[1];
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;

  // step 1: compute the mixture coefficients
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              salinity,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  stackArray2d< real64, maxNumComps *maxNumComps > kij( numComps, numComps );
  real64 dkij_dT = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      getBinaryInteractionCiefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij( ic, jc ), dkij_dT );
    }
  }
  CubicEOSPhaseModel< EOS_TYPE >::computeCompressibilityFactor( numComps,
                                                                composition,
                                                                kij.toSliceConst(),
                                                                aPureCoefficient,
                                                                bPureCoefficient,
                                                                aMixtureCoefficient,
                                                                bMixtureCoefficient,
                                                                compressibilityFactor );
}

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computeCompressibilityFactor( integer const numComps,
                              real64 const & pressure,
                              real64 const & temperature,
                              arraySlice1d< real64 const, USD > const & composition,
                              ComponentProperties::KernelWrapper const & componentProperties,
                              real64 const & salinity,
                              real64 const & compressibilityFactor,
                              arraySlice1d< real64 > const & compressibilityFactorDerivs )
{
  integer constexpr numMaxDofs = maxNumComps + 2;

  stackArray2d< real64, 4*numMaxDofs > tempData( 4, numMaxDofs );
  arraySlice1d< real64 > aPureCoefficient = tempData[0];
  arraySlice1d< real64 > bPureCoefficient = tempData[1];
  real64 aMixtureCoefficient = 0.0;
  real64 bMixtureCoefficient = 0.0;
  arraySlice1d< real64 > aMixtureCoefficientDerivs = tempData[2];
  arraySlice1d< real64 > bMixtureCoefficientDerivs = tempData[3];

  // step 1: compute the mixture coefficients
  // 1.1: Compute the pure and mixture coefficients
  computeMixtureCoefficients( numComps, // number of components
                              pressure, // cell input
                              temperature,
                              composition,
                              componentProperties, // user input,
                              salinity,
                              aPureCoefficient, // output
                              bPureCoefficient,
                              aMixtureCoefficient,
                              bMixtureCoefficient );

  // 1.2: Compute the pure and mixture coefficient derivatives
  computeMixtureCoefficients( numComps,
                              pressure,
                              temperature,
                              composition,
                              componentProperties,
                              salinity,
                              aPureCoefficient.toSliceConst(),
                              bPureCoefficient.toSliceConst(),
                              aMixtureCoefficient,
                              bMixtureCoefficient,
                              aMixtureCoefficientDerivs,    // output
                              bMixtureCoefficientDerivs );

  // 2. Calculate derivatives
  CubicEOSPhaseModel< EOS_TYPE >::computeCompressibilityFactor( numComps,
                                                                aMixtureCoefficient,
                                                                bMixtureCoefficient,
                                                                compressibilityFactor,
                                                                aMixtureCoefficientDerivs.toSliceConst(),
                                                                bMixtureCoefficientDerivs.toSliceConst(),
                                                                compressibilityFactorDerivs );
}

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
computeLogFugacityCoefficients( integer const numComps,
                                real64 const & pressure,
                                real64 const & temperature,
                                arraySlice1d< real64 const, USD > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 const & salinity,
                                arraySlice1d< real64 > const & logFugacityCoefficients )
{
  if constexpr (PHASE_TYPE == SoreideWhitsonPhaseType::Aqueous)
  {
    // step 0: allocate the stack memory needed for the update
    stackArray2d< real64, 2*maxNumComps > pureCoefficients( 2, numComps );
    arraySlice1d< real64 > aPureCoefficient = pureCoefficients[0];
    arraySlice1d< real64 > bPureCoefficient = pureCoefficients[1];
    real64 aMixtureCoefficient = 0.0;
    real64 bMixtureCoefficient = 0.0;
    real64 compressibilityFactor = 0.0;

    // step 1: calculate the dynamic binary interaction coefficients
    stackArray2d< real64, maxNumComps *maxNumComps > kij( numComps, numComps );
    real64 dkij_dT = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      for( integer jc = 0; jc < numComps; ++jc )
      {
        getBinaryInteractionCiefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij( ic, jc ), dkij_dT );
      }
    }

    // step 2: compute the mixture coefficients
    computeMixtureCoefficients( numComps,
                                pressure,
                                temperature,
                                composition,
                                componentProperties,
                                salinity,
                                aPureCoefficient,
                                bPureCoefficient,
                                aMixtureCoefficient,
                                bMixtureCoefficient );

    // step 3: compute the compressibility factor
    CubicEOSPhaseModel< EOS_TYPE >::computeCompressibilityFactor( numComps,
                                                                  composition,
                                                                  kij.toSliceConst(),
                                                                  aPureCoefficient,
                                                                  bPureCoefficient,
                                                                  aMixtureCoefficient,
                                                                  bMixtureCoefficient,
                                                                  compressibilityFactor );

    // step 4: use mixture coefficients and compressibility factor to update fugacity coefficients
    CubicEOSPhaseModel< EOS_TYPE >::computeLogFugacityCoefficients( numComps,
                                                                    composition,
                                                                    kij.toSliceConst(),
                                                                    compressibilityFactor,
                                                                    aPureCoefficient,
                                                                    bPureCoefficient,
                                                                    aMixtureCoefficient,
                                                                    bMixtureCoefficient,
                                                                    logFugacityCoefficients );
  }
  else
  {
    CubicEOSPhaseModel< EOS_TYPE >::computeLogFugacityCoefficients( numComps,
                                                                    pressure,
                                                                    temperature,
                                                                    composition,
                                                                    componentProperties,
                                                                    logFugacityCoefficients );
  }
}

template< SoreideWhitsonPhaseType PHASE_TYPE, typename EOS_TYPE >
template< integer USD >
GEOS_HOST_DEVICE
void
SoreideWhitsonEOSModel< PHASE_TYPE, EOS_TYPE >::
  computeLogFugacityCoefficients( integer const numComps,
                                  real64 const & pressure,
                                  real64 const & temperature,
                                  arraySlice1d< real64 const, USD > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  real64 const & salinity,
                                  arraySlice1d< real64 const > const & logFugacityCoefficients,
                                  arraySlice2d< real64 > const & logFugacityCoefficientDerivs )
                                  {
                                    GEOS_UNUSED_VAR(numComps);
                                    GEOS_UNUSED_VAR(pressure);
                                    GEOS_UNUSED_VAR(temperature);
                                    GEOS_UNUSED_VAR(composition);
                                    GEOS_UNUSED_VAR(componentProperties);
                                    GEOS_UNUSED_VAR(salinity);
                                    GEOS_UNUSED_VAR(logFugacityCoefficients);
                                    GEOS_UNUSED_VAR(logFugacityCoefficientDerivs);
                                    /**
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

    // step 1: calculate the dynamic binary interaction coefficients
    stackArray2d< real64, maxNumComps *maxNumComps > kij( numComps, numComps );
    stackArray2d< real64, maxNumComps *maxNumComps > dkij_dT( numComps, numComps );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      for( integer jc = 0; jc < numComps; ++jc )
      {
        getBinaryInteractionCiefficient( pressure, temperature, componentProperties, salinity, ic, jc, kij( ic, jc ), dkij_dT(ic, jc) );
      }
    }

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
      real64 const bicValue = ( 1.0 - kij( ic, jc ) ) * aCoeffI * aCoeffJ;
      ki[ic] += composition[jc] * bicValue;
      dki( ic, Deriv::dC + jc ) = bicValue;
      dki( ic, Deriv::dP ) += 0.5 * composition[jc] * bicValue * ( aPureCoefficientDerivs( ic, Deriv::dP )/aPureCoefficient[ic] + aPureCoefficientDerivs( jc, Deriv::dP )/aPureCoefficient[jc] );
      dki( ic, Deriv::dT ) += 0.5 * composition[jc] * bicValue * ( aPureCoefficientDerivs( ic, Deriv::dT )/aPureCoefficient[ic] + aPureCoefficientDerivs( jc, Deriv::dT )/aPureCoefficient[jc] );
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
*/
                                  }

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_SOREIDEWHITSONEOSMODEL_HPP_
