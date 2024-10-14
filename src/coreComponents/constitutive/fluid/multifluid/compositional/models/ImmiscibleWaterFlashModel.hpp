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
 * @file ImmiscibleWaterFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERFLASHMODEL_HPP_

#include "FunctionBase.hpp"
#include "EquationOfState.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "NegativeTwoPhaseFlashModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ModelParameters;

class ImmiscibleWaterFlashModelUpdate final : public FunctionBaseUpdate
{
private:
  static constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
public:

  using PhaseProp = NegativeTwoPhaseFlashModelUpdate::PhaseProp;
  using PhaseComp = NegativeTwoPhaseFlashModelUpdate::PhaseComp;
  using Deriv = multifluid::DerivativeOffset;

  ImmiscibleWaterFlashModelUpdate( integer const numComponents,
                                   integer const liquidIndex,
                                   integer const vapourIndex,
                                   integer const aqueousIndex,
                                   integer const waterComponentIndex,
                                   EquationOfStateType const liquidEos,
                                   EquationOfStateType const vapourEos,
                                   arrayView1d< real64 const > const componentCriticalVolume );

  // Mark as a 3-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return 3; }

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
  template< int USD >
  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  void convertCompositionDerivatives( real64 const hcMoleFraction,
                                      arraySlice1d< real64 const > const & composition,
                                      arraySlice1d< real64, USD > const & derivatives ) const
  {
    real64 dvdzi = 0.0;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      dvdzi += derivatives[Deriv::dC+ic] * composition[ic];
    }
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      derivatives[Deriv::dC+ic] /= hcMoleFraction;
    }
    derivatives[Deriv::dC+m_waterComponentIndex] = dvdzi / hcMoleFraction;
  }

private:
  NegativeTwoPhaseFlashModel::KernelWrapper const m_twoPhaseModel;
  integer const m_numComponents;
  integer const m_liquidIndex;
  integer const m_vapourIndex;
  integer const m_aquoesIndex;
  integer const m_waterComponentIndex;
};

class ImmiscibleWaterFlashModel : public FunctionBase
{
public:
  ImmiscibleWaterFlashModel( string const & name,
                             ComponentProperties const & componentProperties,
                             ModelParameters const & modelParameters );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ImmiscibleWaterFlashModelUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  ModelParameters const & m_parameters;
  integer m_waterComponentIndex{-1};
};

template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void ImmiscibleWaterFlashModelUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                               real64 const & pressure,
                                               real64 const & temperature,
                                               arraySlice1d< real64 const, USD1 > const & compFraction,
                                               arraySlice2d< real64, USD2 > const & kValues,
                                               PhaseProp::SliceType const phaseFraction,
                                               PhaseComp::SliceType const phaseCompFraction ) const
{
  LvArray::forValuesInSlice( phaseFraction.value, setZero );
  LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.value, setZero );
  LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );

  // Water phase
  phaseFraction.value[m_aquoesIndex] = compFraction[m_waterComponentIndex];
  phaseFraction.derivs( m_aquoesIndex, Deriv::dC + m_waterComponentIndex ) = 1.0;
  phaseCompFraction.value( m_aquoesIndex, m_waterComponentIndex ) = 1.0;

  // Total hydrocarbon mole fraction
  real64 const z_hc = 1.0 - compFraction[m_waterComponentIndex];

  if( z_hc < MultiFluidConstants::minForSpeciesPresence )
  {
    // Single phase water
    real64 const constantComposition = 1.0 / (m_numComponents - 1);
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      phaseCompFraction.value( m_liquidIndex, ic ) = constantComposition;
      phaseCompFraction.value( m_vapourIndex, ic ) = constantComposition;
    }
    phaseCompFraction.value( m_liquidIndex, m_waterComponentIndex ) = 0.0;
    phaseCompFraction.value( m_vapourIndex, m_waterComponentIndex ) = 0.0;
  }
  else
  {
    // Hydrocarbon phases

    // Calculate normalised hyrdocarbon composition
    stackArray1d< real64, maxNumComps > composition( m_numComponents );
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      composition[ic] = compFraction[ic] / z_hc;
    }
    composition[m_waterComponentIndex] = 0.0;

    // Perform negative two-phase flash
    m_twoPhaseModel.compute( componentProperties,
                             pressure,
                             temperature,
                             composition.toSliceConst(),
                             kValues,
                             phaseFraction,
                             phaseCompFraction );

    for( integer const phaseIndex : {m_liquidIndex, m_vapourIndex} )
    {
      real64 const v = phaseFraction.value[phaseIndex];
      phaseFraction.value[phaseIndex] *= z_hc;
      LvArray::forValuesInSlice( phaseFraction.derivs[phaseIndex], [&]( real64 & a ){ a *= z_hc; } );
      convertCompositionDerivatives( z_hc, composition.toSliceConst(), phaseFraction.derivs[phaseIndex] );
      phaseFraction.derivs( phaseIndex, Deriv::dC+m_waterComponentIndex ) = -v;

      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        convertCompositionDerivatives( z_hc, composition.toSliceConst(), phaseCompFraction.derivs[phaseIndex][ic] );
      }
    }
  }
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_IMMISCIBLEWATERFLASHMODEL_HPP_
