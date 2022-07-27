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
 * @file ReactiveBrineFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_REACTIVEBRINEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_REACTIVEBRINEFLUID_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "constitutive/fluid/ReactiveMultiFluid.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/PhaseModel.hpp"
#include "constitutive/fluid/PVTFunctions/BrineEnthalpy.hpp"
#include "constitutive/fluid/PVTFunctions/NoOpPVTFunction.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineViscosity.hpp"

#include <memory>

namespace geosx
{

namespace constitutive
{

template< typename PHASE >
class ReactiveBrineFluid : public ReactiveMultiFluid
{
public:

  using exec_policy = parallelDevicePolicy<>;

  ReactiveBrineFluid( string const & name,
                      Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  virtual bool isThermal() const override;

  /**
   * @brief Kernel wrapper class for ReactiveBrineFluid.
   */
  class KernelWrapper final : public ReactiveMultiFluid::KernelWrapper
  {
public:

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
                          arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
                          arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
                          real64 & totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void compute( real64 const pressure,
                          real64 const temperature,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
                          PhaseProp::SliceType const phaseFraction,
                          PhaseProp::SliceType const phaseDensity,
                          PhaseProp::SliceType const phaseMassDensity,
                          PhaseProp::SliceType const phaseViscosity,
                          PhaseProp::SliceType const phaseEnthalpy,
                          PhaseProp::SliceType const phaseInternalEnergy,
                          PhaseComp::SliceType const phaseCompFraction,
                          FluidProp::SliceType const totalDensity ) const override;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class ReactiveBrineFluid;

    KernelWrapper( PHASE const & phase,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool const useMass,
                   bool const isThermal,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity,
                   integer const numPrimarySpecies,
                   chemicalReactions::EquilibriumReactions const & equilibriumReactions,
                   chemicalReactions::KineticReactions const & kineticReactions,
                   arrayView2d< real64 > const & primarySpeciesConcentration,
                   arrayView2d< real64 > const & secondarySpeciesConcentration,
                   arrayView2d< real64 > const & primarySpeciesTotalConcentration,
                   arrayView2d< real64 > const & kineticReactionRates );


    /// Flag to specify whether the model is thermal or not
    bool m_isThermal;

    /// Brine constitutive kernel wrappers
    typename PHASE::KernelWrapper m_phase;

  };

  virtual integer getWaterPhaseIndex() const override final;

  /**
   * @brief Names of the submodels for input
   */
  enum class SubModelInputNames : integer
  {
    DENSITY,         ///< the keyword for the density model
    VISCOSITY,       ///< the keyword for the viscosity model
    ENTHALPY         ///< the keyword for the enthalpy model
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : ReactiveMultiFluid::viewKeyStruct
  {
    static constexpr char const * phasePVTParaFilesString() { return "phasePVTParaFiles"; }
  };

protected:

  virtual void postProcessInput() override;

private:

  void createPVTModels();

  /// Brine constitutive models
  std::unique_ptr< PHASE > m_phase;

};

// these aliases are useful in constitutive dispatch
using BrinePhillipsFluid =
  ReactiveBrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction >,
                 PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                 PVTProps::CO2Solubility >;
using BrinePhillipsThermalFluid =
  ReactiveBrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy >,
                 PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                 PVTProps::CO2Solubility >;

template< typename PHASE >
GEOSX_HOST_DEVICE
inline void
ReactiveBrineFluid< PHASE >::KernelWrapper::
  compute( real64 pressure,
           real64 temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseFraction,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseMassDensity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseViscosity,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseEnthalpy,
           arraySlice1d< real64, multifluid::USD_PHASE - 2 > const & phaseInternalEnergy,
           arraySlice2d< real64, multifluid::USD_PHASE_COMP-2 > const & phaseCompFraction,
           real64 & totalDensity ) const
{
  integer constexpr numComp = 2;
  integer constexpr numPhase = 1;

  // 1. Convert input mass fractions to mole fractions and keep derivatives
  stackArray1d< real64, numComp > compMoleFrac( numComp );

  for( integer ic = 0; ic < numComp; ++ic )
  {
    compMoleFrac[ic] = composition[ic];
  }

  // 2. Compute phase fractions and phase component fractions
  real64 const temperatureInCelsius = temperature - 273.15;

  // 3. Compute phase densities and phase viscosities
  m_phase.density.compute( pressure,
                            temperatureInCelsius,
                            phaseCompFraction[0].toSliceConst(),
                            phaseDensity[0],
                            m_useMass );
  m_phase.viscosity.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction[0].toSliceConst(),
                              phaseViscosity[0],
                              m_useMass );

  // for now, we have to compute the phase mass density here
  m_phase.density.compute( pressure,
                           temperatureInCelsius,
                           phaseCompFraction[0].toSliceConst(),
                           phaseMassDensity[0],
                           true );

  // 4. Compute enthalpy and internal energy
  if( m_isThermal )
  {

    m_phase.enthalpy.compute( pressure,
                               temperatureInCelsius,
                               phaseCompFraction[0].toSliceConst(),
                               phaseEnthalpy[0],
                               m_useMass );

    computeInternalEnergy< numComp, numPhase >( pressure,
                                                phaseFraction,
                                                phaseMassDensity,
                                                phaseEnthalpy,
                                                phaseInternalEnergy );

  }

  // 6. Compute total fluid mass/molar density
  computeTotalDensity< numComp, numPhase >( phaseFraction,
                                            phaseDensity,
                                            totalDensity );
}

template< typename PHASE >
GEOSX_HOST_DEVICE
inline void
ReactiveBrineFluid< PHASE, PHASE2, FLASH >::KernelWrapper::
  compute( real64 const pressure,
           real64 const temperature,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition,
           PhaseProp::SliceType const phaseFraction,
           PhaseProp::SliceType const phaseDensity,
           PhaseProp::SliceType const phaseMassDensity,
           PhaseProp::SliceType const phaseViscosity,
           PhaseProp::SliceType const phaseEnthalpy,
           PhaseProp::SliceType const phaseInternalEnergy,
           PhaseComp::SliceType const phaseCompFraction,
           FluidProp::SliceType const totalDensity ) const
{
  integer constexpr numComp = 2;
  integer constexpr numPhase = 1;

  stackArray1d< real64, numComp > compMoleFrac( numComp );
  for( integer ic = 0; ic < numComp; ++ic )
  {
    compMoleFrac[ic] = composition[ic];
  }

  // 2. Compute phase fractions and phase component fractions
  real64 const temperatureInCelsius = temperature - 273.15;

  // 3. Compute phase densities and phase viscosities
  m_phase.density.compute( pressure,
                            temperatureInCelsius,
                            phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                            phaseDensity.value[0], phaseDensity.derivs[0],
                            m_useMass );

  m_phase.viscosity.compute( pressure,
                             temperatureInCelsius,
                             phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                             phaseViscosity.value[0], phaseViscosity.derivs[0],
                             m_useMass );


  // for now, we have to compute the phase mass density here
  m_phase.density.compute( pressure,
                           temperatureInCelsius,
                           phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                           phaseMassDensity.value[0], phaseMassDensity.derivs[0],
                           true );

  // 4. Compute enthalpy and internal energy
  if( m_isThermal )
  {
    m_phase.enthalpy.compute( pressure,
                               temperatureInCelsius,
                               phaseCompFraction.value[0].toSliceConst(), phaseCompFraction.derivs[0].toSliceConst(),
                               phaseEnthalpy.value[0], phaseEnthalpy.derivs[0],
                               m_useMass );

    computeInternalEnergy( pressure,
                           phaseFraction,
                           phaseMassDensity,
                           phaseEnthalpy,
                           phaseInternalEnergy );
  }

  // 6. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );
}

template< typename PHASE >
GEOSX_HOST_DEVICE inline void
ReactiveBrineFluid< PHASE >::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< geosx::real64 const, compflow::USD_COMP - 1 > const & composition ) const
{
  
  computeChemistry( pressure,
                    temperature,
                    composition,
                    m_primarySpeciesConcentration[k],
                    m_secondarySpeciesConcentration[k],
                    m_primarySpeciesTotalConcentration[k],
                    m_kineticReactionRates[k] );

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
           m_totalDensity( k, q ) );
}


} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_REACTIVEBRINEFLUID_HPP_
