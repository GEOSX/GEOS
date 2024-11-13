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
 * @file CompositionalMultiphaseFluidPVTPackage.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDPVTPACKAGE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDPVTPACKAGE_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

#include "constitutive/PVTPackage/PVTPackage/source/pvt/pvt.hpp"

namespace geos
{
namespace constitutive
{

class CompositionalMultiphaseFluidPVTPackage : public MultiFluidBase
{
public:

  using exec_policy = serialPolicy;

  CompositionalMultiphaseFluidPVTPackage( string const & name, Group * const parent );

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  static string catalogName() { return "CompositionalMultiphaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  static constexpr bool isThermalType(){ return false; }

  // TODO: This method should be implemented if an incorrect extrapolation of the pressure and temperature is encountered in the kernel
  /**
   * @copydoc MultiFluidBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  virtual void checkTablesParameters( real64 pressure, real64 temperature ) const override final
  {
    GEOS_UNUSED_VAR( pressure, temperature );
  }

  virtual integer getWaterPhaseIndex() const override final;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * equationsOfStateString() { return "equationsOfState"; }
    static constexpr char const * componentCriticalPressureString() { return "componentCriticalPressure"; }
    static constexpr char const * componentCriticalTemperatureString() { return "componentCriticalTemperature"; }
    static constexpr char const * componentAcentricFactorString() { return "componentAcentricFactor"; }
    static constexpr char const * componentVolumeShiftString() { return "componentVolumeShift"; }
    static constexpr char const * componentBinaryCoeffString() { return "componentBinaryCoeff"; }
    static constexpr char const * constantPhaseViscosityString() { return "constantPhaseViscosity"; }
  };

  /**
   * @brief Kernel wrapper class for CompositionalMultiphaseFluidPVTPackage.
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

    friend class CompositionalMultiphaseFluidPVTPackage;

    KernelWrapper( pvt::MultiphaseSystem & fluid,
                   arrayView1d< pvt::PHASE_TYPE > const & phaseTypes,
                   arrayView1d< real64 const > const & constantPhaseViscosity,
                   arrayView1d< real64 const > const & componentMolarWeight,
                   bool const useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    pvt::MultiphaseSystem & m_fluid;

    arrayView1d< pvt::PHASE_TYPE > m_phaseTypes;
    arrayView1d< real64 const > m_constantPhaseViscosity;
  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostSubGroups() override;

private:

  void createFluid();

  /// PVTPackage fluid object
  std::unique_ptr< pvt::MultiphaseSystem > m_fluid{};

  /// PVTPackage phase labels
  array1d< pvt::PHASE_TYPE > m_phaseTypes;

  // names of equations of state to use for each phase
  string_array m_equationsOfState;

  // Phase viscosity
  array1d< real64 > m_constantPhaseViscosity;

  // standard EOS component input
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidPVTPackage::KernelWrapper::
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
  GEOS_UNUSED_VAR( phaseEnthalpy, phaseInternalEnergy );
#if defined(GEOS_DEVICE_COMPILE)
  GEOS_ERROR( "This function cannot be used on GPU" );
  GEOS_UNUSED_VAR( m_fluid );
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( composition );
  GEOS_UNUSED_VAR( phaseFraction );
  GEOS_UNUSED_VAR( phaseDensity );
  GEOS_UNUSED_VAR( phaseMassDensity );
  GEOS_UNUSED_VAR( phaseViscosity );
  GEOS_UNUSED_VAR( phaseCompFraction );
  GEOS_UNUSED_VAR( totalDensity );
#else

  using Deriv = constitutive::multifluid::DerivativeOffset;

  integer constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  integer constexpr maxNumPhase = MultiFluidBase::MAX_NUM_PHASES;
  integer const numComp = numComponents();
  integer const numPhase = numPhases();

  // 1. Convert input mass fractions to mole fractions and keep derivatives

  std::vector< double > compMoleFrac( numComp );
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

  // 2. Trigger PVTPackage compute and get back phase split

  m_fluid.Update( pressure, temperature, compMoleFrac );

  GEOS_WARNING_IF( !m_fluid.hasSucceeded(),
                   "Phase equilibrium calculations not converged" );

  pvt::MultiphaseSystemProperties const & props = m_fluid.getMultiphaseSystemProperties();

  // 3. Extract phase split, phase properties and derivatives from PVTPackage

  for( integer ip = 0; ip < numPhase; ++ip )
  {
    pvt::PHASE_TYPE const & phaseType = m_phaseTypes[ip];

    auto const & frac = props.getPhaseMoleFraction( phaseType );
    auto const & comp = props.getMoleComposition( phaseType );
    auto const & dens = m_useMass ? props.getMassDensity( phaseType ) : props.getMoleDensity( phaseType );
    auto const & massDens = props.getMassDensity( phaseType );

    phaseFraction.value[ip] = frac.value;
    phaseFraction.derivs[ip][Deriv::dP] = frac.dP;
    phaseFraction.derivs[ip][Deriv::dT] = frac.dT;

    phaseDensity.value[ip] = dens.value;
    phaseDensity.derivs[ip][Deriv::dP] = dens.dP;
    phaseDensity.derivs[ip][Deriv::dT] = dens.dT;

    phaseMassDensity.value[ip] = massDens.value;
    phaseMassDensity.derivs[ip][Deriv::dP] = massDens.dP;
    phaseMassDensity.derivs[ip][Deriv::dT] = massDens.dT;

    // TODO
    phaseViscosity.value[ip] = m_constantPhaseViscosity[ip];
    phaseViscosity.derivs[ip][Deriv::dP] = 0.0;
    phaseViscosity.derivs[ip][Deriv::dT] = 0.0;

    for( integer jc = 0; jc < numComp; ++jc )
    {
      phaseFraction.derivs[ip][Deriv::dC+jc] = frac.dz[jc];
      phaseDensity.derivs[ip][Deriv::dC+jc] = dens.dz[jc];
      phaseMassDensity.derivs[ip][Deriv::dC+jc] = massDens.dz[jc];
      phaseViscosity.derivs[ip][Deriv::dC+jc] = 0.0; // TODO

      phaseCompFraction.value[ip][jc] = comp.value[jc];
      phaseCompFraction.derivs[ip][jc][Deriv::dP] = comp.dP[jc];
      phaseCompFraction.derivs[ip][jc][Deriv::dT] = comp.dT[jc];

      for( integer ic = 0; ic < numComp; ++ic )
      {
        phaseCompFraction.derivs[ip][ic][Deriv::dC+jc] = comp.dz[ic][jc];
      }
    }
  }

  // 4. if mass variables used instead of molar, perform the conversion
  if( m_useMass )
  {

    // unfortunately here, we have to copy the molecular weight coming from PVT package...
    real64 phaseMolecularWeight[maxNumPhase]{};
    real64 dPhaseMolecularWeight[maxNumPhase][maxNumComp+2]{};

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      phaseMolecularWeight[ip] = props.getMolecularWeight( m_phaseTypes[ip] ).value;
      dPhaseMolecularWeight[ip][Deriv::dP] = props.getMolecularWeight( m_phaseTypes[ip] ).dP;
      dPhaseMolecularWeight[ip][Deriv::dT] = props.getMolecularWeight( m_phaseTypes[ip] ).dT;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dPhaseMolecularWeight[ip][Deriv::dC+ic] = props.getMolecularWeight( m_phaseTypes[ip] ).dz[ic];
      }
    }

    convertToMassFractions( dCompMoleFrac_dCompMassFrac,
                            phaseMolecularWeight,
                            dPhaseMolecularWeight,
                            phaseFraction,
                            phaseCompFraction,
                            phaseDensity.derivs,
                            phaseViscosity.derivs,
                            phaseEnthalpy.derivs,
                            phaseInternalEnergy.derivs );
  }

  // 5. Compute total fluid mass/molar density and derivatives

  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );

#endif
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
CompositionalMultiphaseFluidPVTPackage::KernelWrapper::
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

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUIDPVTPACKAGE_HPP_
