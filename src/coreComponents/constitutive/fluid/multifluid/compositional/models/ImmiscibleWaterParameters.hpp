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
 * @file ImmiscibleWaterParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERPARAMETERS_HPP_

#include "ModelParameters.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ImmiscibleWaterParameters : public ModelParameters
{
  static constexpr char const * waterComponentName = "h2o";

public:
  ImmiscibleWaterParameters( std::unique_ptr< ModelParameters > parameters );
  ~ImmiscibleWaterParameters() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters );

  static integer getWaterComponentIndex( ComponentProperties const & componentProperties );

  struct viewKeyStruct
  {
    static constexpr char const * waterReferencePressureString() { return "waterReferencePressure"; }
    static constexpr char const * waterReferenceTemperatureString() { return "waterReferenceTemperature"; }
    static constexpr char const * waterDensityString() { return "waterDensity"; }
    static constexpr char const * waterViscosityString() { return "waterViscosity"; }
    static constexpr char const * waterCompressibilityString() { return "waterCompressibility"; }
    static constexpr char const * waterViscosityCompressibilityString() { return "waterViscosityCompressibility"; }
    static constexpr char const * waterExpansionCoefficientString() { return "waterExpansionCoefficient"; }
    static constexpr char const * waterViscosityExpansionCoefficientString() { return "waterViscosityExpansionCoefficient"; }
  };

  real64 m_waterReferencePressure;
  real64 m_waterReferenceTemperature{293.15};
  real64 m_waterDensity;
  real64 m_waterViscosity;
  real64 m_waterCompressibility;
  real64 m_waterViscosityCompressibility{0.0};
  real64 m_waterExpansionCoefficient{0.0};
  real64 m_waterViscosityExpansionCoefficient{0.0};

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override;

  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERPARAMETERS_HPP_
