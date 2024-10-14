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
 * @file CompositionalMultiphaseFluidUpdates.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_

#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

namespace geos
{
namespace constitutive
{

/**
 * @brief Kernel wrapper class for CompositionalMultiphaseFluid.
 * @tparam FLASH Class describing the phase equilibrium model
 * @tparam PHASE1 Class describing the phase property models for the first phase.
 * @tparam PHASE2 Class describing the phase property models for the second phase.
 * @tparam PHASE3 Class describing the phase property models for the possible third phase.
 */
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
class CompositionalMultiphaseFluidUpdates final : public MultiFluidBase::KernelWrapper
{
public:
  CompositionalMultiphaseFluidUpdates( compositional::ComponentProperties const & componentProperties,
                                       FLASH const & flash,
                                       PHASE1 const & phase1,
                                       PHASE2 const & phase2,
                                       PHASE3 const & phase3,
                                       arrayView1d< real64 const > const & componentMolarWeight,
                                       bool const useMass,
                                       MultiFluidBase::PhaseProp::ViewType phaseFrac,
                                       MultiFluidBase::PhaseProp::ViewType phaseDens,
                                       MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                       MultiFluidBase::PhaseProp::ViewType phaseVisc,
                                       MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                       MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                       MultiFluidBase::PhaseComp::ViewType phaseCompFrac,
                                       MultiFluidBase::FluidProp::ViewType totalDensity,
                                       MultiFluidBase::PhaseComp::ViewValueType kValues );

  GEOS_HOST_DEVICE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                        MultiFluidBase::PhaseProp::SliceType const phaseFrac,
                        MultiFluidBase::PhaseProp::SliceType const phaseDens,
                        MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
                        MultiFluidBase::PhaseProp::SliceType const phaseVisc,
                        MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
                        MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
                        MultiFluidBase::PhaseComp::SliceType const phaseCompFrac,
                        MultiFluidBase::FluidProp::SliceType const totalDensity ) const override;

  GEOS_HOST_DEVICE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

protected:
  GEOS_HOST_DEVICE
  void compute( real64 const pressure,
                real64 const temperature,
                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                MultiFluidBase::PhaseProp::SliceType const phaseFrac,
                MultiFluidBase::PhaseProp::SliceType const phaseDens,
                MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
                MultiFluidBase::PhaseProp::SliceType const phaseVisc,
                MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
                MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
                MultiFluidBase::PhaseComp::SliceType const phaseCompFrac,
                MultiFluidBase::FluidProp::SliceType const totalDensity,
                MultiFluidBase::PhaseComp::SliceType::ValueType const & kValues ) const;

  /**
   * @brief Convert derivatives from phase mole fraction to total mole fraction
   * @details Given property derivatives @c dProperty where composition derivatives are with
   *          respect to a phase compositions, this will transform that properties so that
   *          they the composition derivatives are with respect to total composition. The derivatives
   *          of the phase composition should be provided in @c dPhaseComposition.
   * @param[in] numComps The number of components
   * @param[in] dPhaseComposition Derivatives of the phase composition
   * @param[in,out] dProperty The derivatives of the property
   * @param[in] workSpace Temporary workspace
   */
  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  static void convertDerivativesToTotalMoleFraction( integer const numComps,
                                                     arraySlice2d< real64 const, USD1 > const & dPhaseComposition,
                                                     arraySlice1d< real64, USD2 > const & dProperty,
                                                     arraySlice1d< real64 > const & workSpace );

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setZero( real64 & val ){ val = 0.0; }

private:
  // The component properties
  compositional::ComponentProperties::KernelWrapper m_componentProperties;

  // Flash kernel wrapper
  typename FLASH::KernelWrapper m_flash;

  // Phase model kernel wrappers
  typename PHASE1::KernelWrapper m_phase1;
  typename PHASE2::KernelWrapper m_phase2;
  typename PHASE3::KernelWrapper m_phase3;

  // Backup variables
  MultiFluidBase::PhaseComp::ViewValueType m_kValues;
};

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::
CompositionalMultiphaseFluidUpdates( compositional::ComponentProperties const & componentProperties,
                                     FLASH const & flash,
                                     PHASE1 const & phase1,
                                     PHASE2 const & phase2,
                                     PHASE3 const & phase3,
                                     arrayView1d< real64 const > const & componentMolarWeight,
                                     bool const useMass,
                                     MultiFluidBase::PhaseProp::ViewType phaseFrac,
                                     MultiFluidBase::PhaseProp::ViewType phaseDens,
                                     MultiFluidBase::PhaseProp::ViewType phaseMassDensity,
                                     MultiFluidBase::PhaseProp::ViewType phaseVisc,
                                     MultiFluidBase::PhaseProp::ViewType phaseEnthalpy,
                                     MultiFluidBase::PhaseProp::ViewType phaseInternalEnergy,
                                     MultiFluidBase::PhaseComp::ViewType phaseCompFrac,
                                     MultiFluidBase::FluidProp::ViewType totalDensity,
                                     MultiFluidBase::PhaseComp::ViewValueType kValues ):
  MultiFluidBase::KernelWrapper( componentMolarWeight,
                                 useMass,
                                 std::move( phaseFrac ),
                                 std::move( phaseDens ),
                                 std::move( phaseMassDensity ),
                                 std::move( phaseVisc ),
                                 std::move( phaseEnthalpy ),
                                 std::move( phaseInternalEnergy ),
                                 std::move( phaseCompFrac ),
                                 std::move( totalDensity ) ),
  m_componentProperties( componentProperties.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() ),
  m_phase1( phase1.createKernelWrapper() ),
  m_phase2( phase2.createKernelWrapper() ),
  m_phase3( phase3.createKernelWrapper() ),
  m_kValues( kValues )
{}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::compute(
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
  MultiFluidBase::PhaseProp::SliceType const phaseFrac,
  MultiFluidBase::PhaseProp::SliceType const phaseDens,
  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseVisc,
  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
  MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
  MultiFluidBase::PhaseComp::SliceType const phaseCompFrac,
  MultiFluidBase::FluidProp::SliceType const totalDensity ) const
{
  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES - 1;
  MultiFluidBase::PhaseComp::StackValueType< maxNumPhase *maxNumComp > kValues( 1, 1, numPhases() - 1, numComponents() );

  LvArray::forValuesInSlice( kValues[0][0], setZero );   // Force initialisation of k-Values

  compute( pressure,
           temperature,
           composition,
           phaseFrac,
           phaseDens,
           phaseMassDensity,
           phaseVisc,
           phaseEnthalpy,
           phaseInternalEnergy,
           phaseCompFrac,
           totalDensity,
           kValues[0][0] );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::compute(
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
  MultiFluidBase::PhaseProp::SliceType const phaseFrac,
  MultiFluidBase::PhaseProp::SliceType const phaseDens,
  MultiFluidBase::PhaseProp::SliceType const phaseMassDensity,
  MultiFluidBase::PhaseProp::SliceType const phaseVisc,
  MultiFluidBase::PhaseProp::SliceType const phaseEnthalpy,
  MultiFluidBase::PhaseProp::SliceType const phaseInternalEnergy,
  MultiFluidBase::PhaseComp::SliceType const phaseCompFrac,
  MultiFluidBase::FluidProp::SliceType const totalDensity,
  MultiFluidBase::PhaseComp::SliceType::ValueType const & kValues ) const
{
  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumDof = MultiFluidBase::MAX_NUM_COMPONENTS + 2;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();
  integer const numDof = numComp + 2;

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, maxNumComp > compMoleFrac( numComp );
  real64 dCompMoleFrac_dCompMassFrac[maxNumComp][maxNumComp]{};

  if( m_useMass )
  {
    // convert mass fractions to mole fractions
    convertToMoleFractions( composition,
                            compMoleFrac,
                            dCompMoleFrac_dCompMassFrac );
  }
  else
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      compMoleFrac[ic] = composition[ic];
    }
  }

  // 2. Compute phase fractions and phase component fractions
  m_flash.compute( m_componentProperties,
                   pressure,
                   temperature,
                   compMoleFrac.toSliceConst(),
                   kValues,
                   phaseFrac,
                   phaseCompFrac );

  // 3. Calculate the phase densities
  m_phase1.density.compute( m_componentProperties,
                            pressure,
                            temperature,
                            phaseCompFrac.value[0].toSliceConst(),
                            phaseDens.value[0],
                            phaseDens.derivs[0],
                            phaseMassDensity.value[0],
                            phaseMassDensity.derivs[0],
                            m_useMass );
  m_phase2.density.compute( m_componentProperties,
                            pressure,
                            temperature,
                            phaseCompFrac.value[1].toSliceConst(),
                            phaseDens.value[1],
                            phaseDens.derivs[1],
                            phaseMassDensity.value[1],
                            phaseMassDensity.derivs[1],
                            m_useMass );
  if constexpr (2 < FLASH::KernelWrapper::getNumberOfPhases())
  {
    m_phase3.density.compute( m_componentProperties,
                              pressure,
                              temperature,
                              phaseCompFrac.value[2].toSliceConst(),
                              phaseDens.value[2],
                              phaseDens.derivs[2],
                              phaseMassDensity.value[2],
                              phaseMassDensity.derivs[2],
                              m_useMass );
  }

  // 4. Calculate the phase viscosities
  m_phase1.viscosity.compute( m_componentProperties,
                              pressure,
                              temperature,
                              phaseCompFrac.value[0].toSliceConst(),
                              phaseMassDensity.value[0],
                              phaseMassDensity.derivs[0].toSliceConst(),
                              phaseVisc.value[0],
                              phaseVisc.derivs[0],
                              m_useMass );
  m_phase2.viscosity.compute( m_componentProperties,
                              pressure,
                              temperature,
                              phaseCompFrac.value[1].toSliceConst(),
                              phaseMassDensity.value[1],
                              phaseMassDensity.derivs[1].toSliceConst(),
                              phaseVisc.value[1],
                              phaseVisc.derivs[1],
                              m_useMass );
  if constexpr (2 < FLASH::KernelWrapper::getNumberOfPhases())
  {
    m_phase3.viscosity.compute( m_componentProperties,
                                pressure,
                                temperature,
                                phaseCompFrac.value[2].toSliceConst(),
                                phaseMassDensity.value[2],
                                phaseMassDensity.derivs[2].toSliceConst(),
                                phaseVisc.value[2],
                                phaseVisc.derivs[2],
                                m_useMass );
  }

  // 5. Convert derivatives from phase composition to total composition
  stackArray1d< real64, maxNumDof > workSpace( numDof );
  for( integer ip = 0; ip < FLASH::KernelWrapper::getNumberOfPhases(); ++ip )
  {
    convertDerivativesToTotalMoleFraction( numComp,
                                           phaseCompFrac.derivs[ip].toSliceConst(),
                                           phaseDens.derivs[ip],
                                           workSpace );
    convertDerivativesToTotalMoleFraction( numComp,
                                           phaseCompFrac.derivs[ip].toSliceConst(),
                                           phaseMassDensity.derivs[ip],
                                           workSpace );
    convertDerivativesToTotalMoleFraction( numComp,
                                           phaseCompFrac.derivs[ip].toSliceConst(),
                                           phaseVisc.derivs[ip],
                                           workSpace );
  }

  // 6. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {
    real64 phaseMolecularWeight[maxNumPhase]{};
    real64 dPhaseMolecularWeight[maxNumPhase][maxNumComp+2]{};

    arrayView1d< real64 const > const & componentMolarWeight = m_componentProperties.m_componentMolarWeight;

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = 0.0;
      for( integer kc = 0; kc < numDof; ++kc )
      {
        dPhaseMolecularWeight[ip][kc] = 0.0;
      }

      auto const & phaseComposition = phaseCompFrac.value[ip].toSliceConst();
      auto const & dPhaseComposition = phaseCompFrac.derivs[ip].toSliceConst();

      for( integer ic = 0; ic < numComp; ++ic )
      {
        phaseMolecularWeight[ip] += phaseComposition[ic] * componentMolarWeight[ic];
        for( integer kc = 0; kc < numDof; ++kc )
        {
          dPhaseMolecularWeight[ip][kc] += dPhaseComposition( ic, kc ) * componentMolarWeight[ic];
        }
      }
    }
    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight,
                            phaseFrac,
                            phaseCompFrac,
                            phaseDens.derivs,
                            phaseVisc.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );
  }

  // 7. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFrac,
                       phaseDens,
                       totalDensity );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::
update( localIndex const k,
        localIndex const q,
        real64 const pressure,
        real64 const temperature,
        arraySlice1d< geos::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  compute( pressure,
           temperature,
           composition,
           m_phaseFraction( k, q ),
           m_phaseDensity( k, q ),
           m_phaseMassDensity( k, q ),
           m_phaseViscosity( k, q ),
           m_phaseEnthalpy( k, q ),
           m_phaseInternalEnergy( k, q ),
           m_phaseCompFraction( k, q ),
           m_totalDensity( k, q ),
           m_kValues[k][q] );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void
CompositionalMultiphaseFluidUpdates< FLASH, PHASE1, PHASE2, PHASE3 >::
convertDerivativesToTotalMoleFraction( integer const numComps,
                                       arraySlice2d< real64 const, USD1 > const & dPhaseComposition,
                                       arraySlice1d< real64, USD2 > const & dProperty,
                                       arraySlice1d< real64 > const & workSpace )
{
  using Deriv = constitutive::multifluid::DerivativeOffset;
  integer const numDofs = numComps + 2;
  for( integer kc = 0; kc < numDofs; ++kc )
  {
    workSpace[kc] = dProperty[kc];
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    dProperty[Deriv::dC+ic] = 0.0;
  }
  for( integer kc = 0; kc < numDofs; ++kc )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      dProperty[kc] += (dPhaseComposition( ic, kc ) * workSpace[Deriv::dC+ic]);
    }
  }
}

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDUPDATES_HPP_
