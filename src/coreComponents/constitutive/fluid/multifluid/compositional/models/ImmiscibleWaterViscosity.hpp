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
 * @file ImmiscibleWaterViscosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERVISCOSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERVISCOSITY_HPP_

#include "FunctionBase.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ImmiscibleWaterViscosityUpdate final : public FunctionBaseUpdate
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  ImmiscibleWaterViscosityUpdate( real64 const referencePressure,
                                  real64 const referenceTemperature,
                                  real64 const viscosity,
                                  real64 const compressibility,
                                  real64 const expansionCoefficient );

  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 const & density,
                arraySlice1d< real64 const, USD2 > const & dDensity,
                real64 & viscosity,
                arraySlice1d< real64, USD2 > const & dViscosity,
                bool useMass ) const;

private:
  real64 const m_referencePressure;
  real64 const m_referenceTemperature;
  real64 const m_viscosity;
  real64 const m_compressibility;
  real64 const m_expansionCoefficient;
};

class ImmiscibleWaterViscosity : public FunctionBase
{
public:
  ImmiscibleWaterViscosity( string const & name,
                            ComponentProperties const & componentProperties,
                            integer const phaseIndex,
                            ModelParameters const & modelParameters );

  static string catalogName() { return ""; }

  FunctionType functionType() const override
  {
    return FunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ImmiscibleWaterViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  ModelParameters const & m_parameters;
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ImmiscibleWaterViscosityUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                              real64 const & pressure,
                                              real64 const & temperature,
                                              arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                              real64 const & density,
                                              arraySlice1d< real64 const, USD2 > const & dDensity,
                                              real64 & viscosity,
                                              arraySlice1d< real64, USD2 > const & dViscosity,
                                              bool useMass ) const
{
  GEOS_UNUSED_VAR( componentProperties );
  GEOS_UNUSED_VAR( phaseComposition );
  GEOS_UNUSED_VAR( density );
  GEOS_UNUSED_VAR( dDensity );
  GEOS_UNUSED_VAR( useMass );

  LvArray::forValuesInSlice( dViscosity, setZero );

  viscosity = m_viscosity *
              LvArray::math::exp( m_compressibility * (pressure - m_referencePressure) ) *
              LvArray::math::exp( -m_expansionCoefficient * (temperature - m_referenceTemperature) );
  dViscosity[Deriv::dP] = m_compressibility * viscosity;
  dViscosity[Deriv::dT] = -m_expansionCoefficient * viscosity;
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERVISCOSITY_HPP_
