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
 * @file KValueFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/RachfordRice.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< integer NUM_PHASE >
class KValueFlashParameters;

template< integer NUM_PHASE >
class KValueFlashModelUpdate final : public FunctionBaseUpdate
{
  static constexpr integer numPhases = NUM_PHASE;
  // Must be 2-phase or 3-phase
  static_assert( NUM_PHASE == 2 || NUM_PHASE == 3, "KValue flash must be 2-phase or 3-phase" );

public:
  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;
  using Deriv = constitutive::multifluid::DerivativeOffset;

  KValueFlashModelUpdate( integer const numComponents,
                          TableFunction const & pressureTable,
                          TableFunction const & temperatureTable,
                          arrayView1d< integer const > const & presentComponents,
                          arrayView4d< real64 const > const & kValues );

  // Mark as a 2-phase or 3-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return numPhases; }

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice2d< real64, USD2 > const & kValues,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

private:
  static integer constexpr m_liquidIndex = 0;
  static integer constexpr m_vapourIndex = 1;
  static integer constexpr m_aqueousIndex = 2;

  integer const m_numComponents{0};
  integer const m_numPressurePoints{0};
  integer const m_numTemperaturePoints{0};

  /// Table with pressure indices
  TableFunction::KernelWrapper m_pressureTable;

  /// Table with temperature indices
  TableFunction::KernelWrapper m_temperatureTable;

  /// Index of available components
  arrayView1d< integer const > m_presentComponents;

  /// Hypercube of k-values
  arrayView4d< real64 const > m_kValues;
};

template< integer NUM_PHASE >
class KValueFlashModel : public FunctionBase
{
public:
  KValueFlashModel( string const & name,
                    ComponentProperties const & componentProperties,
                    ModelParameters const & modelParameters );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = KValueFlashModelUpdate< NUM_PHASE >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  /// Index of present componenyt
  array1d< integer > m_presentComponents;

  KValueFlashParameters< NUM_PHASE > const * m_parameters{};
};

template< integer NUM_PHASE >
template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void
KValueFlashModelUpdate< NUM_PHASE >::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                              real64 const & pressure,
                                              real64 const & temperature,
                                              arraySlice1d< real64 const, USD1 > const & compFraction,
                                              arraySlice2d< real64, USD2 > const & kValues,
                                              PhaseProp::SliceType const phaseFraction,
                                              PhaseComp::SliceType const phaseCompFraction ) const
{
  GEOS_UNUSED_VAR( componentProperties );

  integer const numDofs = 2 + m_numComponents;

  // Calculate k-values at p,t
  real64 pa = 0.0;
  real64 pa_dp = 0.0;
  pa = m_pressureTable.compute( &pressure, &pa_dp );
  integer const pn = LvArray::math::min( static_cast< integer >(pa), m_numPressurePoints - 2 );
  pa -= pn;

  real64 ta = 0.0;
  real64 ta_dt = 0.0;
  ta = m_temperatureTable.compute( &temperature, &ta_dt );
  integer const tn = LvArray::math::min( static_cast< integer >(ta), m_numTemperaturePoints - 2 );
  ta -= tn;

  stackArray3d< real64, 2*(NUM_PHASE-1)*MultiFluidConstants::MAX_NUM_COMPONENTS > kValueDerivatives( 2, NUM_PHASE-1, m_numComponents );
  arraySlice2d< real64 > dKValue_dP = kValueDerivatives[0];
  arraySlice2d< real64 > dKValue_dT = kValueDerivatives[1];
  for( integer ip = 0; ip < NUM_PHASE-1; ++ip )
  {
    for( integer ic = 0; ic < m_numComponents; ic++ )
    {
      real64 const k00 = m_kValues( ip, ic, pn, tn );
      real64 const k01 = m_kValues( ip, ic, pn, tn+1 );
      real64 const k10 = m_kValues( ip, ic, pn+1, tn );
      real64 const k11 = m_kValues( ip, ic, pn+1, tn+1 );
      kValues( ip, ic ) = (1.0-pa)*(1.0-ta)*k00 + (1.0-pa)*ta*k01 + pa*(1.0-ta)*k10 + pa*ta*k11;
      dKValue_dP( ip, ic ) = pa_dp*( -(1.0-ta)*k00 - ta*k01 + (1.0-ta)*k10 + ta*k11 );
      dKValue_dT( ip, ic ) = ta_dt*( -(1.0-pa)*k00 + (1.0-pa)*k01 - pa*k10 + pa*k11 );
    }
  }

  if constexpr (NUM_PHASE == 2)
  {
    real64 vapourFraction = RachfordRice::solve( kValues[0].toSliceConst(), compFraction, m_presentComponents.toSliceConst() );

    // Test for single phase
    if( vapourFraction < MultiFluidConstants::epsilon || 1.0-vapourFraction < MultiFluidConstants::epsilon )
    {
      vapourFraction = LvArray::math::min( LvArray::math::max( vapourFraction, 0.0 ), 1.0 );
      phaseFraction.value[m_vapourIndex] = vapourFraction;

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
    else
    {
      // Calculate derivatives implicitly from the Rachford-Rice equation
      real64 denominator = 0.0;
      real64 pressureNumerator = 0.0;
      real64 temperatureNumerator = 0.0;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const k = kValues( 0, ic ) - 1.0;
        real64 const r = 1.0 / ( 1.0 + vapourFraction*k );
        pressureNumerator += compFraction[ic] * dKValue_dP( 0, ic ) * r * r;
        temperatureNumerator += compFraction[ic] * dKValue_dT( 0, ic ) * r * r;
        denominator += compFraction[ic] * k * k * r * r;
        phaseFraction.derivs( m_vapourIndex, Deriv::dC + ic ) = k * r;
      }
      GEOS_ERROR_IF( denominator < MultiFluidConstants::epsilon,
                     "Failed to calculate derivatives for the Rachford-Rice equation." );
      real64 const invDenominator = 1.0 / denominator;
      phaseFraction.derivs( m_vapourIndex, Deriv::dP ) = pressureNumerator * invDenominator;
      phaseFraction.derivs( m_vapourIndex, Deriv::dT ) = temperatureNumerator * invDenominator;
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        phaseFraction.derivs( m_vapourIndex, Deriv::dC + ic ) *= invDenominator;
      }

      // Calculate phase compositions
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const k = kValues( 0, ic ) - 1.0;
        real64 const r = 1.0 / ( 1.0 + vapourFraction*k );
        real64 const xi = compFraction[ic] * r *;
        real64 const yi = xi * kValues( 0, ic );
        phaseCompFraction.value( m_liquidIndex, ic ) = xi;
        phaseCompFraction.value( m_vapourIndex, ic ) = yi;
      }
    }

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_liquidIndex] = 1.0 - phaseFraction.value[m_vapourIndex];
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
    }
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Rachford-Rice solve for {} phases not implemented.", NUM_PHASE ));
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_
