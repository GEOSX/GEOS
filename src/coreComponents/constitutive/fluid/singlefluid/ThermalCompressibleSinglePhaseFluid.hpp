/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ThermalCompressibleSinglePhaseFluid.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_THERMALCOMPRESSIBLESINGLEPHASEFLUID_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_THERMALCOMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/CompressibleSinglePhaseFluid.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Update class for the model suitable for lambda capture.
 * @tparam DENS_EAT type of density exponent approximation for the pressure part
 * @tparam VISC_EAT type of viscosity exponent approximation
 * @tparam INTENERGY_EAT type of internal energy exponent approximation
 */
template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT, ExponentApproximationType INTENERGY_EAT >
class ThermalCompressibleSinglePhaseUpdate : public SingleFluidBaseUpdate
{
public:

  using DensRelationType      = ExponentialRelation< real64, DENS_EAT, 3 >;
  using ViscRelationType      = ExponentialRelation< real64, VISC_EAT >;
  using IntEnergyRelationType = ExponentialRelation< real64, INTENERGY_EAT >;
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;

  ThermalCompressibleSinglePhaseUpdate( DensRelationType const & densRelation,
                                        ViscRelationType const & viscRelation,
                                        IntEnergyRelationType const & intEnergyRelation,
                                        arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & density,
                                        arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > const & dDensity,
                                        arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & viscosity,
                                        arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > const & dViscosity,
                                        arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & internalEnergy,
                                        arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > const & dInternalEnergy,
                                        arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & enthalpy,
                                        arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > const & dEnthalpy,
                                        real64 const & refIntEnergy )
    : SingleFluidBaseUpdate( density,
                             dDensity,
                             viscosity,
                             dViscosity ),
    m_internalEnergy( internalEnergy ),
    m_dInternalEnergy( dInternalEnergy ),
    m_enthalpy( enthalpy ),
    m_dEnthalpy( dEnthalpy ),
    m_densRelation( densRelation ),
    m_viscRelation( viscRelation ),
    m_intEnergyRelation( intEnergyRelation ),
    m_refIntEnergy( refIntEnergy )
  {}

  /// Default copy constructor
  ThermalCompressibleSinglePhaseUpdate( ThermalCompressibleSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  ThermalCompressibleSinglePhaseUpdate( ThermalCompressibleSinglePhaseUpdate && ) = default;

  /// Deleted copy assignment operator
  ThermalCompressibleSinglePhaseUpdate & operator=( ThermalCompressibleSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  ThermalCompressibleSinglePhaseUpdate & operator=( ThermalCompressibleSinglePhaseUpdate && ) = delete;

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure ) const override
  {
    m_densRelation.compute( pressure, density, dDensity_dPressure );
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        real64 & density,
                        real64 & dDensity_dPressure,
                        real64 & dDensity_dTemperature,
                        real64 & viscosity,
                        real64 & dViscosity_dPressure,
                        real64 & dViscosity_dTemperature,
                        real64 & internalEnergy,
                        real64 & dInternalEnergy_dPressure,
                        real64 & dInternalEnergy_dTemperature,
                        real64 & enthalpy,
                        real64 & dEnthalpy_dPressure,
                        real64 & dEnthalpy_dTemperature ) const override
  {
    m_viscRelation.compute( pressure, viscosity, dViscosity_dPressure );
    dViscosity_dTemperature = 0.0;

    m_densRelation.compute( pressure, temperature, density, dDensity_dPressure, dDensity_dTemperature );

    /// Compute the internal energy (only sensitive to temperature)
    m_intEnergyRelation.compute( temperature, internalEnergy, dInternalEnergy_dTemperature );
    dInternalEnergy_dPressure = 0.0;

    enthalpy = internalEnergy - m_refIntEnergy;
    dEnthalpy_dPressure = 0.0;
    dEnthalpy_dTemperature = dInternalEnergy_dTemperature;

  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure ) const override
  {
    compute( pressure,
             m_density[k][q],
             m_dDensity[k][q][DerivOffset::dP],
             m_viscosity[k][q],
             m_dViscosity[k][q][DerivOffset::dP] );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature ) const override
  {
    compute( pressure,
             temperature,
             m_density[k][q],
             m_dDensity[k][q][DerivOffset::dP],
             m_dDensity[k][q][DerivOffset::dT],
             m_viscosity[k][q],
             m_dViscosity[k][q][DerivOffset::dP],
             m_dViscosity[k][q][DerivOffset::dT],
             m_internalEnergy[k][q],
             m_dInternalEnergy[k][q][DerivOffset::dP],
             m_dInternalEnergy[k][q][DerivOffset::dT],
             m_enthalpy[k][q],
             m_dEnthalpy[k][q][DerivOffset::dP],
             m_dEnthalpy[k][q][DerivOffset::dT] );
  }

private:

  /// Fluid internal energy and derivatives
  arrayView2d< real64, constitutive::singlefluid::USD_FLUID > m_internalEnergy;
  arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > m_dInternalEnergy;

  /// Fluid enthalpy and derivatives
  arrayView2d< real64, constitutive::singlefluid::USD_FLUID > m_enthalpy;
  arrayView3d< real64, constitutive::singlefluid::USD_FLUID_DER > m_dEnthalpy;

  /// Relationship between the fluid density and pressure & temperature
  DensRelationType m_densRelation;

  /// Relationship between the fluid viscosity and pressure
  ViscRelationType m_viscRelation;

  /// Relationship between the fluid internal energy and temperature
  IntEnergyRelationType m_intEnergyRelation;

  /// Reference internal energy of the fluid
  real64 const m_refIntEnergy;

};

class ThermalCompressibleSinglePhaseFluid : public CompressibleSinglePhaseFluid
{
public:

  ThermalCompressibleSinglePhaseFluid( string const & name, Group * const parent );

  virtual ~ThermalCompressibleSinglePhaseFluid() override;

  static string catalogName() { return "ThermalCompressibleSinglePhaseFluid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  using CompressibleSinglePhaseFluid::m_densityModelType;

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using KernelWrapper = ThermalCompressibleSinglePhaseUpdate< ExponentApproximationType::Full, ExponentApproximationType::Linear, ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  struct viewKeyStruct : public CompressibleSinglePhaseFluid::viewKeyStruct
  {
    static constexpr char const * thermalExpansionCoeffString() { return "thermalExpansionCoeff"; }
    static constexpr char const * specificHeatCapacityString() { return "specificHeatCapacity"; }
    static constexpr char const * referenceTemperatureString() { return "referenceTemperature"; }
    static constexpr char const * referenceInternalEnergyString() { return "referenceInternalEnergy"; }
    static constexpr char const * internalEnergyModelTypeString() { return "internalEnergyModelType"; }
  };

protected:

  virtual void postInputInitialization() override;

private:

  /// scalar fluid thermal expansion coefficient
  real64 m_thermalExpansionCoeff;

  /// scalar fluid volumetric heat capacity coefficient
  real64 m_specificHeatCapacity;

  /// reference temperature parameter
  real64 m_referenceTemperature;

  /// reference internal energy parameter
  real64 m_referenceInternalEnergy;

  /// type of internal energy model
  ExponentApproximationType m_internalEnergyModelType;
};

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_THERMALCOMPRESSIBLESINGLEPHASEFLUID_HPP_ */
