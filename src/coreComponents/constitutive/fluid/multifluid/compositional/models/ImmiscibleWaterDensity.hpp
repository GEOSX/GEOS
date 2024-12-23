/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ImmiscibleWaterDensity.cpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERDENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERDENSITY_HPP_

#include "FunctionBase.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ImmiscibleWaterDensityUpdate final : public FunctionBaseUpdate
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  ImmiscibleWaterDensityUpdate( real64 const waterMolecularWeight,
                                real64 const referencePressure,
                                real64 const referenceTemperature,
                                real64 const density,
                                real64 const compressibility,
                                real64 const expansionCoefficient );

  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & molarDensity,
                arraySlice1d< real64, USD2 > const & dMolarDensity,
                real64 & massDensity,
                arraySlice1d< real64, USD2 > const & dMassDensity,
                bool useMass ) const;

private:
  real64 const m_waterMolecularWeight;
  real64 const m_referencePressure;
  real64 const m_referenceTemperature;
  real64 const m_density;
  real64 const m_compressibility;
  real64 const m_expansionCoefficient;
};

class ImmiscibleWaterDensity : public FunctionBase
{
public:
  ImmiscibleWaterDensity( string const & name,
                          ComponentProperties const & componentProperties,
                          integer const phaseIndex,
                          ModelParameters const & modelParameters );

  static string catalogName() { return "ImmiscibleWaterDensity"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ImmiscibleWaterDensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  ModelParameters const & m_parameters;
  real64 m_waterMolecularWeight{0.0};
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void ImmiscibleWaterDensityUpdate::compute(
  ComponentProperties::KernelWrapper const & componentProperties,
  real64 const & pressure,
  real64 const & temperature,
  arraySlice1d< real64 const, USD1 > const & phaseComposition,
  real64 & molarDensity,
  arraySlice1d< real64, USD2 > const & dMolarDensity,
  real64 & massDensity,
  arraySlice1d< real64, USD2 > const & dMassDensity,
  bool useMass ) const
{
  GEOS_UNUSED_VAR( componentProperties );
  GEOS_UNUSED_VAR( phaseComposition );
  GEOS_UNUSED_VAR( useMass );

  LvArray::forValuesInSlice( dMolarDensity, setZero );
  LvArray::forValuesInSlice( dMassDensity, setZero );

  real64 const density = m_density *
                         LvArray::math::exp( m_compressibility * (pressure - m_referencePressure) ) *
                         LvArray::math::exp( -m_expansionCoefficient * (temperature - m_referenceTemperature) );
  real64 const dDensity_dp = m_compressibility * density;
  real64 const dDensity_dT = -m_expansionCoefficient * density;

  massDensity = density;
  dMassDensity[Deriv::dP] = dDensity_dp;
  dMassDensity[Deriv::dT] = dDensity_dT;

  molarDensity = density / m_waterMolecularWeight;
  dMolarDensity[Deriv::dP] = dDensity_dp / m_waterMolecularWeight;
  dMolarDensity[Deriv::dT] = dDensity_dT / m_waterMolecularWeight;
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERDENSITY_HPP_
