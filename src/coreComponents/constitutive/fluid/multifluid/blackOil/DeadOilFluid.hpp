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
 * @file DeadOilFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_DEADOILFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_BLACKOIL_DEADOILFLUID_HPP_

#include "constitutive/fluid/multifluid/blackOil/BlackOilFluidBase.hpp"

namespace geos
{

namespace constitutive
{

class DeadOilFluid : public BlackOilFluidBase
{
public:

  DeadOilFluid( string const & name, Group * const parent );

  static string catalogName() { return "DeadOilFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Kernel wrapper class for DeadOilFluid
   *        This kernel can be called on the GPU
   */
  class KernelWrapper final : public BlackOilFluidBase::KernelWrapper
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

    friend class DeadOilFluid;

    /**
     * @brief Constructor for the class doing in-kernel Dead-Oil fluid updates
     * @param[in] phaseTypes type of phases
     * @param[in] phaseOrder order of phases
     * @param[in] hydrocarbonPhaseOrder order of the hydrocarbon phases in the model
     * @param[in] surfacePhaseMassDensity surface phase mass densities provided by the user
     * @param[in] formationVolFractionTables hydrocarbon formation volume tables
     * @param[in] viscosityTables hydrocarbon viscosity tables
     * @param[in] waterParams parameters (Bo, visc) for the water phase
     * @param[in] componentMolarWeight component molecular weights
     * @param[in] useMass flag to decide whether we return mass or molar densities
     * @param[in] phaseFraction phase fractions (+ derivatives) in the cell
     * @param[in] phaseDensity phase mass/molar densities (+ derivatives) in the cell
     * @param[in] phaseMassDensity phase mass densities (+ derivatives) in the cell
     * @param[in] phaseViscosity phase viscosities (+ derivatives) in the cell
     * @param[in] phaseCompFraction phase component fractions (+ derivatives) in the cell
     * @param[in] totalDensity total density in the cell
     */
    KernelWrapper( arrayView1d< integer const > phaseTypes,
                   arrayView1d< integer const > phaseOrder,
                   arrayView1d< integer const > hydrocarbonPhaseOrder,
                   arrayView1d< real64 const > surfacePhaseMassDensity,
                   arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                   arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                   BlackOilFluidBase::WaterParams const waterParams,
                   arrayView1d< real64 const > componentMolarWeight,
                   bool useMass,
                   PhaseProp::ViewType phaseFraction,
                   PhaseProp::ViewType phaseDensity,
                   PhaseProp::ViewType phaseMassDensity,
                   PhaseProp::ViewType phaseViscosity,
                   PhaseProp::ViewType phaseEnthalpy,
                   PhaseProp::ViewType phaseInternalEnergy,
                   PhaseComp::ViewType phaseCompFraction,
                   FluidProp::ViewType totalDensity );

    /**
     * @brief Utility function to compute mass densities as a function of pressure (keeping derivatives)
     * @param[in] pressure pressure in the cell
     * @param[out] phaseMassDens the phase mass densities in the cell (+ derivatives)
     */
    GEOS_HOST_DEVICE
    void computeDensities( real64 const pressure,
                           PhaseProp::SliceType const & phaseMassDens ) const;

    /**
     * @brief Utility function to compute viscosities as a function of pressure (keeping derivatives)
     * @param[in] pressure pressure in the cell
     * @param[out] phaseVisc the phase viscosities in the cell (+ derivatives)
     */
    GEOS_HOST_DEVICE
    void computeViscosities( real64 const pressure,
                             PhaseProp::SliceType const & phaseVisc ) const;

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

private:

  /**
   * @brief Use the TableFunctions provided by the user to get the PVT data
   */
  virtual void readInputDataFromTableFunctions() override;

  /**
   * @brief Read all the PVT table provided by the user in Eclipse format
   */
  virtual void readInputDataFromPVTFiles() override;

};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
DeadOilFluid::KernelWrapper::
  computeDensities( real64 const pressure,
                    PhaseProp::SliceType const & phaseMassDens ) const
{
  using Deriv = multifluid::DerivativeOffset;

  LvArray::forValuesInSlice( phaseMassDens.derivs, []( real64 & val ){ val = 0.0; } );

  // 1. Hydrocarbon phases: look up in the formation vol factor tables

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    integer const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase formation vol factor and its derivative wrt pressure
    real64 derivative;
    real64 const fvf = m_formationVolFactorTables[iph].compute( &pressure, &derivative );

    // we are ready to update the densities
    real64 const fvfInv = 1.0 / fvf;
    phaseMassDens.value[ip] = m_surfacePhaseMassDensity[ip] * fvfInv;
    phaseMassDens.derivs[ip][Deriv::dP] = -derivative * phaseMassDens.value[ip] * fvfInv;
  }

  // 2. Water phase: use the constant formation volume factor and compressibility provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // note: double check, but I think std::exp is allowed in kernel
    real64 const expCompDeltaPres = std::exp( -m_waterParams.compressibility * ( pressure - m_waterParams.referencePressure ) );
    real64 const dExpCompDeltaPres_dPres = -m_waterParams.compressibility * expCompDeltaPres;
    real64 const denom = m_waterParams.formationVolFactor * expCompDeltaPres;
    real64 const dDenom_dPres = m_waterParams.formationVolFactor * dExpCompDeltaPres_dPres;
    real64 const denomInv = 1.0 / denom;
    phaseMassDens.value[ipWater] = m_surfacePhaseMassDensity[ipWater] * denomInv;
    phaseMassDens.derivs[ipWater][Deriv::dP] = -dDenom_dPres * phaseMassDens.value[ipWater] * denomInv;
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
DeadOilFluid::KernelWrapper::
  computeViscosities( real64 const pressure,
                      PhaseProp::SliceType const & phaseVisc ) const
{
  using Deriv = multifluid::DerivativeOffset;

  LvArray::forValuesInSlice( phaseVisc.derivs, []( real64 & val ){ val = 0.0; } );

  // 1. Hydrocarbon phases: look up in the viscosity tables

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    // get the phase index
    integer const ip = m_hydrocarbonPhaseOrder[iph];
    // interpolate in the table to get the phase viscosity and derivatives
    real64 dPhaseVisc_dPres = 0.0;
    phaseVisc.value[ip] = m_viscosityTables[iph].compute( &pressure, &dPhaseVisc_dPres );
    phaseVisc.derivs[ip][Deriv::dP] = dPhaseVisc_dPres;
  }

  // 2. Water phase: use the constant viscosity provided by the user

  using PT = DeadOilFluid::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];

  // if water is present
  if( ipWater >= 0 )
  {
    // just assign the viscosity value
    phaseVisc.value[ipWater] = m_waterParams.viscosity;
    phaseVisc.derivs[ipWater][Deriv::dP] = 0.0;
  }
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
DeadOilFluid::KernelWrapper::
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
  GEOS_UNUSED_VAR( temperature, phaseEnthalpy, phaseInternalEnergy );

  using Deriv = multifluid::DerivativeOffset;

  integer const nComps = numComponents();
  integer const nPhases = numPhases();

  // 1. Read viscosities and formation volume factors from tables, update mass densities
  computeViscosities( pressure,
                      phaseViscosity );
  computeDensities( pressure,
                    phaseMassDensity );

  // 2. Update phaseDens (mass density if useMass == 1, molar density otherwise)
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    real64 const mult = m_useMass ? 1.0 : 1.0 / m_componentMolarWeight[ip];
    phaseDensity.value[ip] = phaseMassDensity.value[ip] * mult;
    phaseDensity.derivs[ip][Deriv::dP] = phaseMassDensity.derivs[ip][Deriv::dP] * mult;
    for( integer ic = 0; ic < nComps; ++ic )
    {
      phaseDensity.derivs[ip][Deriv::dC+ic] = 0.0;
    }
  }

  // 3. Update remaining variables: phaseFrac, phaseCompFrac using Dead-Oil assumptions
  for( integer ip = 0; ip < nPhases; ++ip )
  {
    phaseFraction.value[ip] = composition[ip];
    phaseFraction.derivs[ip][Deriv::dP] = 0.0;
    for( integer ic = 0; ic < nComps; ++ic )
    {
      phaseFraction.derivs[ip][Deriv::dC+ic] = (ip == ic) ? 1.0 : 0.0;

      phaseCompFraction.value[ip][ic] = (ip == ic) ? 1.0 : 0.0;
      phaseCompFraction.derivs[ip][ic][Deriv::dP] = 0.0;
      for( integer jc = 0; jc < nComps; ++jc )
      {
        phaseCompFraction.derivs[ip][ic][Deriv::dC+jc] = 0.0;
      }
    }
  }

  // 4. Compute total fluid mass/molar density and derivatives
  computeTotalDensity( phaseFraction,
                       phaseDensity,
                       totalDensity );
}

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void
DeadOilFluid::KernelWrapper::
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

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_DEADOILFLUID_HPP_
