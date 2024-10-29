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
 * @file NegativeTwoPhaseFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_

#include "FunctionBase.hpp"
#include "EquationOfState.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class EquationOfState;

class NegativeTwoPhaseFlashModelUpdate final : public FunctionBaseUpdate
{
public:

  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;

  NegativeTwoPhaseFlashModelUpdate( integer const numComponents,
                                    integer const liquidIndex,
                                    integer const vapourIndex,
                                    EquationOfStateType const liquidEos,
                                    EquationOfStateType const vapourEos );

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

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_liquidIndex] = 1.0 - phaseFraction.value[m_vapourIndex];
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
    }
  }

private:
  integer const m_numComponents;
  integer const m_liquidIndex;
  integer const m_vapourIndex;
  EquationOfStateType const m_liquidEos;
  EquationOfStateType const m_vapourEos;
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
  EquationOfState const * m_parameters{};
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
