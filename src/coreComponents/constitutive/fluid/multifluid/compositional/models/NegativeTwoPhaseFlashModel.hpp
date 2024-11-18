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
 * @file NegativeTwoPhaseFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_

#include "FunctionBase.hpp"
#include "EquationOfState.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/StabilityTest.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ModelParameters;

class NegativeTwoPhaseFlashModelUpdate final : public FunctionBaseUpdate
{
public:
  using PhaseProp = MultiFluidVar< real64, 3, constitutive::multifluid::LAYOUT_PHASE, constitutive::multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, constitutive::multifluid::LAYOUT_PHASE_COMP, constitutive::multifluid::LAYOUT_PHASE_COMP_DC >;
  using Deriv = constitutive::multifluid::DerivativeOffset;

  static constexpr real64 stabilityTolerance = MultiFluidConstants::fugacityTolerance;

  NegativeTwoPhaseFlashModelUpdate( integer const numComponents,
                                    integer const liquidIndex,
                                    integer const vapourIndex,
                                    EquationOfStateType const liquidEos,
                                    EquationOfStateType const vapourEos,
                                    arrayView1d< real64 const > const componentCriticalVolume );

  // Mark as a 2-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return 2; }

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice2d< real64, USD2 > const & kValues,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const
  {
    integer const numDofs = 2 + m_numComponents;

    // Perform stability test to check that we have 2 phases
    real64 tangentPlaneDistance = 0.0;
    bool const stabilityStatus = StabilityTest::compute( m_numComponents,
                                                         pressure,
                                                         temperature,
                                                         compFraction,
                                                         componentProperties,
                                                         m_liquidEos,
                                                         tangentPlaneDistance,
                                                         kValues[0] );
    GEOS_ERROR_IF( !stabilityStatus,
                   GEOS_FMT( "Stability test failed at pressure {:.5e} and temperature {:.3f}", pressure, temperature ));

    if( tangentPlaneDistance < -stabilityTolerance )
    {
      // Unstable mixture
      // Iterative solve to converge flash
      bool const flashStatus = NegativeTwoPhaseFlash::compute( m_numComponents,
                                                               pressure,
                                                               temperature,
                                                               compFraction,
                                                               componentProperties,
                                                               m_liquidEos,
                                                               m_vapourEos,
                                                               kValues,
                                                               phaseFraction.value[m_vapourIndex],
                                                               phaseCompFraction.value[m_liquidIndex],
                                                               phaseCompFraction.value[m_vapourIndex] );

      GEOS_ERROR_IF( !flashStatus,
                     GEOS_FMT( "Negative two phase flash failed to converge at pressure {:.5e} and temperature {:.3f}",
                               pressure, temperature ));

      // Calculate derivatives
      NegativeTwoPhaseFlash::computeDerivatives( m_numComponents,
                                                 pressure,
                                                 temperature,
                                                 compFraction,
                                                 componentProperties,
                                                 m_liquidEos,
                                                 m_vapourEos,
                                                 phaseFraction.value[m_vapourIndex],
                                                 phaseCompFraction.value[m_liquidIndex].toSliceConst(),
                                                 phaseCompFraction.value[m_vapourIndex].toSliceConst(),
                                                 phaseFraction.derivs[m_vapourIndex],
                                                 phaseCompFraction.derivs[m_liquidIndex],
                                                 phaseCompFraction.derivs[m_vapourIndex] );
    }
    else
    {
      // Stable mixture - simply label
      calculateLiCorrelation( componentProperties,
                              temperature,
                              compFraction,
                              phaseFraction.value[m_vapourIndex] );

      LvArray::forValuesInSlice( phaseFraction.derivs[m_vapourIndex], setZero );
      LvArray::forValuesInSlice( phaseCompFraction.derivs[m_liquidIndex], setZero );
      LvArray::forValuesInSlice( phaseCompFraction.derivs[m_vapourIndex], setZero );
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        phaseCompFraction.value( m_vapourIndex, ic ) = compFraction[ic];
        phaseCompFraction.value( m_liquidIndex, ic ) = compFraction[ic];
        phaseCompFraction.derivs( m_vapourIndex, ic, Deriv::dC + ic ) = 1.0;
        phaseCompFraction.derivs( m_liquidIndex, ic, Deriv::dC + ic ) = 1.0;
      }
    }

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_liquidIndex] = 1.0 - phaseFraction.value[m_vapourIndex];
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
    }
  }

  template< int USD >
  GEOS_HOST_DEVICE
  void calculateLiCorrelation( ComponentProperties::KernelWrapper const & componentProperties,
                               real64 const & temperature,
                               arraySlice1d< real64 const, USD > const & composition,
                               real64 & vapourFraction ) const
  {
    real64 sumVz = 0.0;
    real64 sumVzt = 0.0;
    arrayView1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      real64 const Vz = m_componentCriticalVolume[ic] * composition[ic];
      sumVz += Vz;
      sumVzt += Vz * criticalTemperature[ic];
    }
    real64 const pseudoCritTemperature = sumVzt / sumVz;
    if( pseudoCritTemperature < temperature )
    {
      vapourFraction = 1.0;
    }
    else
    {
      vapourFraction = 0.0;
    }
  }

private:
  integer const m_numComponents;
  integer const m_liquidIndex;
  integer const m_vapourIndex;
  EquationOfStateType const m_liquidEos;
  EquationOfStateType const m_vapourEos;
  arrayView1d< real64 const > const m_componentCriticalVolume;
};

class NegativeTwoPhaseFlashModel : public FunctionBase
{
public:
  NegativeTwoPhaseFlashModel( string const & name,
                              ComponentProperties const & componentProperties,
                              ModelParameters const & modelParameters );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NegativeTwoPhaseFlashModelUpdate;

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

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
