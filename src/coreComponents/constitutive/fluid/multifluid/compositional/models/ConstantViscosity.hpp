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
 * @file ConstantViscosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CONSTANTVISCOSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CONSTANTVISCOSITY_HPP_

#include "FunctionBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ConstantViscosityUpdate final : public FunctionBaseUpdate
{
public:
  explicit ConstantViscosityUpdate( real64 const constantPhaseViscosity );

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
  real64 const m_constantPhaseViscosity;
};

class ConstantViscosity : public FunctionBase
{
public:
  static constexpr real64 defaultViscosity = 0.001;
public:
  ConstantViscosity( string const & name,
                     ComponentProperties const & componentProperties,
                     integer const phaseIndex,
                     ModelParameters const & modelParameters );

  static string catalogName() { return ""; }

  FunctionType functionType() const override
  {
    return FunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ConstantViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Parameters for constant viscosity model
  class Parameters : public ModelParameters
  {
public:
    Parameters( std::unique_ptr< ModelParameters > parameters );
    ~Parameters() override = default;

    array1d< real64 > m_constantPhaseViscosity;

private:
    void registerParametersImpl( MultiFluidBase * fluid ) override;
    void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

    struct viewKeyStruct
    {
      static constexpr char const * constantPhaseViscosityString() { return "constantPhaseViscosity"; }
    };
  };

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  real64 m_constantPhaseViscosity{defaultViscosity};
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ConstantViscosityUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                       real64 const & pressure,
                                       real64 const & temperature,
                                       arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                       real64 const & density,
                                       arraySlice1d< real64 const, USD2 > const & dDensity,
                                       real64 & viscosity,
                                       arraySlice1d< real64, USD2 > const & dViscosity,
                                       bool useMass ) const
{
  GEOS_UNUSED_VAR( componentProperties, pressure, temperature, useMass );
  GEOS_UNUSED_VAR( phaseComposition );
  GEOS_UNUSED_VAR( density, dDensity );

  viscosity = m_constantPhaseViscosity;

  LvArray::forValuesInSlice( dViscosity, setZero );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CONSTANTVISCOSITY_HPP_
