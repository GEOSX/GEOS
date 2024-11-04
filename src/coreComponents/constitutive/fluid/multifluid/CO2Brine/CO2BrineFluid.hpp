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
 * @file CO2BrineFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_CO2BRINEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_CO2BRINEFLUID_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/PhaseModel.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/BrineEnthalpy.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Enthalpy.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Solubility.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/EzrokhiBrineDensity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/EzrokhiBrineViscosity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/FenghourCO2Viscosity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/NoOpPVTFunction.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PhillipsBrineViscosity.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp"
#include "common/Units.hpp"


#include <memory>

namespace geos
{

namespace constitutive
{

template< typename PHASE1, typename PHASE2, typename FLASH >
class CO2BrineFluid : public MultiFluidBase
{
public:

  using exec_policy = parallelDevicePolicy<>;

  CO2BrineFluid( string const & name,
                 Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName();

  virtual string getCatalogName() const override { return catalogName(); }

  static constexpr bool isThermalType()
  {
    return !( std::is_same_v< typename PHASE1::Enthalpy, PVTProps::NoOpPVTFunction > ||
              std::is_same_v< typename PHASE2::Enthalpy, PVTProps::NoOpPVTFunction > );
  }

  static constexpr integer min_n_components = 2;
  static constexpr integer max_n_components = 2;

  virtual bool isThermal() const override
  {
    return isThermalType();
  }

  /**
   * @brief Kernel wrapper class for CO2BrineFluid.
   */
  class KernelWrapper final : public MultiFluidBase::KernelWrapper
  {
public:
    GEOS_HOST_DEVICE
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

    GEOS_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         real64 const pressure,
                         real64 const temperature,
                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & composition ) const override;

private:

    friend class CO2BrineFluid;

    KernelWrapper( integer p1Index,
                   integer p2Index,
                   PHASE1 const & phase1,
                   PHASE2 const & phase2,
                   FLASH const & flash,
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
                   FluidProp::ViewType totalDensity );

    /// Index of the liquid phase
    integer m_p1Index;

    /// Index of the gas phase
    integer m_p2Index;

    /// Flag to specify whether the model is thermal or not
    bool m_isThermal;

    /// Brine constitutive kernel wrappers
    typename PHASE1::KernelWrapper m_phase1;

    // CO2 constitutive kernel wrapper
    typename PHASE2::KernelWrapper m_phase2;

    // Flash kernel wrapper
    typename FLASH::KernelWrapper m_flash;
  };

  virtual integer getWaterPhaseIndex() const override final;

  /**
   * @copydoc MultiFluidBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  void checkTablesParameters( real64 pressure, real64 temperature ) const override final;

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

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * flashModelParaFileString() { return "flashModelParaFile"; }
    static constexpr char const * solubilityTablesString() { return "solubilityTableNames"; }
    static constexpr char const * phasePVTParaFilesString() { return "phasePVTParaFiles"; }
    static constexpr char const * writeCSVFlagString() { return "writeCSV"; }
  };

protected:

  virtual void postInputInitialization() override;

  virtual void initializePreSubGroups() override;

private:

  /**
   * @brief Create a PVT Model and output them
   */
  void createPVTModels();

  /// Names of the files defining the viscosity and density models
  path_array m_phasePVTParaFiles;

  /// Name of the file defining the flash model
  Path m_flashModelParaFile;

  /// Names of solubility tables for each phase
  string_array m_solubilityTables;

  /// Index of the liquid phase
  integer m_p1Index;

  /// Index of the gas phase
  integer m_p2Index;

  /// Output csv file containing informations about PVT
  integer m_writeCSV;

  /// Brine constitutive models
  std::unique_ptr< PHASE1 > m_phase1;

  // CO2 constitutive models
  std::unique_ptr< PHASE2 > m_phase2;

  // Flash model
  std::unique_ptr< FLASH > m_flash;
};

// these aliases are useful in constitutive dispatch
using CO2BrinePhillipsFluid =
  CO2BrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction >,
                 PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                 PVTProps::CO2Solubility >;
using CO2BrinePhillipsThermalFluid =
  CO2BrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy >,
                 PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                 PVTProps::CO2Solubility >;

using CO2BrineEzrokhiFluid =
  CO2BrineFluid< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::NoOpPVTFunction >,
                 PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                 PVTProps::CO2Solubility >;
using CO2BrineEzrokhiThermalFluid =
  CO2BrineFluid< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::BrineEnthalpy >,
                 PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                 PVTProps::CO2Solubility >;

template< typename PHASE1, typename PHASE2, typename FLASH >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CO2BrineFluid< PHASE1, PHASE2, FLASH >::KernelWrapper::
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
  integer constexpr numPhase = 2;
  integer const ip1 = m_p1Index;
  integer const ip2 = m_p2Index;

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  stackArray1d< real64, numComp > compMoleFrac( numComp );
  real64 dCompMoleFrac_dCompMassFrac[numComp][numComp]{};

  if( m_useMass )
  {
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

  real64 const temperatureInCelsius = units::convertKToC( temperature );
  m_flash.compute( pressure,
                   temperatureInCelsius,
                   compMoleFrac.toSliceConst(),
                   phaseFraction,
                   phaseCompFraction );

  // 3. Compute phase densities and phase viscosities

  m_phase1.density.compute( pressure,
                            temperatureInCelsius,
                            phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.derivs[ip1].toSliceConst(),
                            phaseDensity.value[ip1], phaseDensity.derivs[ip1],
                            m_useMass );
  m_phase1.viscosity.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.derivs[ip1].toSliceConst(),
                              phaseViscosity.value[ip1], phaseViscosity.derivs[ip1],
                              m_useMass );
  m_phase2.density.compute( pressure,
                            temperatureInCelsius,
                            phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.derivs[ip2].toSliceConst(),
                            phaseDensity.value[ip2], phaseDensity.derivs[ip2],
                            m_useMass );
  m_phase2.viscosity.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.derivs[ip2].toSliceConst(),
                              phaseViscosity.value[ip2], phaseViscosity.derivs[ip2],
                              m_useMass );

  // 4. Depending on the m_useMass flag, convert to mass variables or simply compute mass density

  // TODO: for now the treatment of molar/mass density requires too many interpolations in the tables, it needs to be fixed
  //       we should modify the PVT functions so that they can return phaseMassDens, phaseDens, and phaseMW in one call

  if( m_useMass )
  {

    // 4.1 Compute the phase molecular weights (ultimately, get that from the PVT function)

    real64 phaseMolecularWeight[numPhase]{};
    real64 dPhaseMolecularWeight[numPhase][numComp+2]{};

    real64 phaseMolarDens{};
    stackArray1d< real64, numComp+2 > dPhaseMolarDens( numComp+2 );

    m_phase1.density.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.derivs[ip1].toSliceConst(),
                              phaseMolarDens, dPhaseMolarDens.toSlice(),
                              false );
    phaseMolecularWeight[ip1] = phaseDensity.value[ip1] / phaseMolarDens;
    for( integer idof = 0; idof < numComp+2; ++idof )
    {
      dPhaseMolecularWeight[ip1][idof] = phaseDensity.derivs[ip1][idof] / phaseMolarDens - phaseMolecularWeight[ip1] * dPhaseMolarDens[idof] / phaseMolarDens;
    }

    m_phase2.density.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.derivs[ip2].toSliceConst(),
                              phaseMolarDens, dPhaseMolarDens.toSlice(),
                              false );
    phaseMolecularWeight[ip2] = phaseDensity.value[ip2] / phaseMolarDens;
    for( integer idof = 0; idof < numComp+2; ++idof )
    {
      dPhaseMolecularWeight[ip2][idof] = phaseDensity.derivs[ip2][idof] / phaseMolarDens - phaseMolecularWeight[ip2] * dPhaseMolarDens[idof] / phaseMolarDens;
    }

    // 4.2 Convert the mole fractions to mass fractions
    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight,
                            phaseFraction,
                            phaseCompFraction,
                            phaseDensity.derivs,
                            phaseViscosity.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );


    // 4.3 Copy the phase densities into the phase mass densities
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMassDensity.value[ip] = phaseDensity.value[ip];
      for( integer idof = 0; idof < numComp+2; ++idof )
      {
        phaseMassDensity.derivs[ip][idof] = phaseDensity.derivs[ip][idof];
      }
    }
  }
  else
  {
    // for now, we have to compute the phase mass density here
    m_phase1.density.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.derivs[ip1].toSliceConst(),
                              phaseMassDensity.value[ip1], phaseMassDensity.derivs[ip1],
                              true );
    m_phase2.density.compute( pressure,
                              temperatureInCelsius,
                              phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.derivs[ip2].toSliceConst(),
                              phaseMassDensity.value[ip2], phaseMassDensity.derivs[ip2],
                              true );
  }

  // 5. Compute enthalpy and internal energy

  if( m_isThermal )
  {
    m_phase1.enthalpy.compute( pressure,
                               temperatureInCelsius,
                               phaseCompFraction.value[ip1].toSliceConst(), phaseCompFraction.derivs[ip1].toSliceConst(),
                               phaseEnthalpy.value[ip1], phaseEnthalpy.derivs[ip1],
                               m_useMass );
    m_phase2.enthalpy.compute( pressure,
                               temperatureInCelsius,
                               phaseCompFraction.value[ip2].toSliceConst(), phaseCompFraction.derivs[ip2].toSliceConst(),
                               phaseEnthalpy.value[ip2], phaseEnthalpy.derivs[ip2],
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

template< typename PHASE1, typename PHASE2, typename FLASH >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CO2BrineFluid< PHASE1, PHASE2, FLASH >::KernelWrapper::
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
           m_totalDensity( k, q ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_CO2BRINEFLUID_HPP_
