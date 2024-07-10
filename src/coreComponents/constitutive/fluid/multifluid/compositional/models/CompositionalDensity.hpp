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
 * @file CompositionalDensity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_

#include "FunctionBase.hpp"
#include "EquationOfState.hpp"

#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CompositionalProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class CompositionalDensityUpdate final : public FunctionBaseUpdate
{
public:
  CompositionalDensityUpdate( arrayView1d< real64 const > const & volumeShift,
                              EquationOfStateType const equationOfState )
    : m_componentDimensionalVolumeShift( volumeShift ),
    m_equationOfState( equationOfState )
  {}

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
  template< integer USD >
  GEOS_HOST_DEVICE
  void computeCompressibilityFactor( integer const numComps,
                                     real64 const & pressure,
                                     real64 const & temperature,
                                     arraySlice1d< real64 const, USD > const & composition,
                                     ComponentProperties::KernelWrapper const & componentProperties,
                                     EquationOfStateType const equationOfState,
                                     real64 & compressibilityFactor,
                                     arraySlice1d< real64 > const & compressibilityFactorDerivs ) const;

private:
  arrayView1d< real64 const > m_componentDimensionalVolumeShift;
  EquationOfStateType const m_equationOfState;
};

class CompositionalDensity : public FunctionBase
{
public:
  CompositionalDensity( string const & name,
                        ComponentProperties const & componentProperties,
                        integer const phaseIndex,
                        ModelParameters const & modelParameters );

  static string catalogName() { return "CompositionalDensity"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CompositionalDensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  static void calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                               EquationOfStateType const & equationOfState,
                                               arraySlice1d< real64 > componentDimensionalVolumeShift );

private:
  array1d< real64 > m_componentDimensionalVolumeShift;
  EquationOfStateType m_equationOfState;
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void CompositionalDensityUpdate::compute(
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
  GEOS_UNUSED_VAR( useMass );

  integer const numComps = componentProperties.m_componentMolarWeight.size();
  integer const numDofs = 2 + numComps;

  real64 compressibilityFactor = 0.0;
  stackArray1d< real64, 2+MultiFluidConstants::MAX_NUM_COMPONENTS > tempDerivs( numDofs );

  computeCompressibilityFactor( numComps,
                                pressure,
                                temperature,
                                phaseComposition,
                                componentProperties,
                                m_equationOfState,
                                compressibilityFactor,
                                tempDerivs.toSlice() );

  CompositionalProperties::computeMolarDensity( numComps,
                                                pressure,
                                                temperature,
                                                phaseComposition,
                                                m_componentDimensionalVolumeShift.toSliceConst(),
                                                compressibilityFactor,
                                                tempDerivs.toSlice(),
                                                molarDensity,
                                                dMolarDensity );

  CompositionalProperties::computeMassDensity( numComps,
                                               phaseComposition,
                                               componentProperties.m_componentMolarWeight.toSliceConst(),
                                               molarDensity,
                                               dMolarDensity.toSliceConst(),
                                               massDensity,
                                               dMassDensity );
}

template< integer USD >
GEOS_HOST_DEVICE
void CompositionalDensityUpdate::computeCompressibilityFactor( integer const numComps,
                                                               real64 const & pressure,
                                                               real64 const & temperature,
                                                               arraySlice1d< real64 const, USD > const & composition,
                                                               ComponentProperties::KernelWrapper const & componentProperties,
                                                               EquationOfStateType const equationOfState,
                                                               real64 & compressibilityFactor,
                                                               arraySlice1d< real64 > const & compressibilityFactorDerivs ) const
{
  if( equationOfState == EquationOfStateType::PengRobinson )
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::
    computeCompressibilityFactor( numComps,
                                  pressure,
                                  temperature,
                                  composition,
                                  componentProperties,
                                  compressibilityFactor,
                                  compressibilityFactorDerivs );
  }
  else if( equationOfState == EquationOfStateType::SoaveRedlichKwong )
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::
    computeCompressibilityFactor( numComps,
                                  pressure,
                                  temperature,
                                  composition,
                                  componentProperties,
                                  compressibilityFactor,
                                  compressibilityFactorDerivs );
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALDENSITY_HPP_
