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
 * @file ThermalCompressibleSinglePhaseFluid.cpp
 */

#include "ThermalCompressibleSinglePhaseFluid.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ThermalCompressibleSinglePhaseFluid::ThermalCompressibleSinglePhaseFluid( string const & name, Group * const parent ):
  CompressibleSinglePhaseFluid( name, parent ),
  m_internalEnergyModelType( ExponentApproximationType::Linear )
{
  m_densityModelType = ExponentApproximationType::Full;

  registerWrapper( viewKeyStruct::thermalExpansionCoeffString(), &m_thermalExpansionCoeff ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid thermal expansion coefficient. Unit: 1/K" );

  registerWrapper( viewKeyStruct::specificHeatCapacityString(), &m_specificHeatCapacity ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid heat capacity. Unit: J/kg/K" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference temperature" );

  registerWrapper( viewKeyStruct::referenceInternalEnergyString(), &m_referenceInternalEnergy ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid internal energy" );

  registerWrapper( viewKeyStruct::internalEnergyModelTypeString(), &m_internalEnergyModelType ).
    setApplyDefaultValue( m_internalEnergyModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of internal energy model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

}

ThermalCompressibleSinglePhaseFluid::~ThermalCompressibleSinglePhaseFluid() = default;

void ThermalCompressibleSinglePhaseFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  CompressibleSinglePhaseFluid::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_internalEnergy.setValues< serialPolicy >( m_referenceInternalEnergy );
}

void ThermalCompressibleSinglePhaseFluid::postInputInitialization()
{
  CompressibleSinglePhaseFluid::postInputInitialization();

  auto const checkNonnegative = [&]( real64 const value, auto const & attribute )
  {
    GEOS_THROW_IF_LT_MSG( value, 0.0,
                          GEOS_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                          InputError );
  };

  checkNonnegative( m_thermalExpansionCoeff, viewKeyStruct::thermalExpansionCoeffString() );
  checkNonnegative( m_specificHeatCapacity, viewKeyStruct::specificHeatCapacityString() );
  checkNonnegative( m_referenceInternalEnergy, viewKeyStruct::referenceInternalEnergyString() );

  // Due to the way update wrapper is currently implemented, we can only support one model type
  auto const checkModelType = [&]( ExponentApproximationType const value, auto const & attribute )
  {
    GEOS_THROW_IF( value != ExponentApproximationType::Linear && value != ExponentApproximationType::Full,
                   GEOS_FMT( "{}: invalid model type in attribute '{}' (only linear or fully exponential currently supported)", getFullName(), attribute ),
                   InputError );
  };
  checkModelType( m_internalEnergyModelType, viewKeyStruct::internalEnergyModelTypeString() );
}

ThermalCompressibleSinglePhaseFluid::KernelWrapper
ThermalCompressibleSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensRelationType( m_referencePressure, m_referenceTemperature, m_referenceDensity, m_compressibility, -m_thermalExpansionCoeff ),
                        KernelWrapper::ViscRelationType( m_referencePressure, m_referenceViscosity, m_viscosibility ),
                        KernelWrapper::IntEnergyRelationType( m_referenceTemperature, m_referenceInternalEnergy, m_specificHeatCapacity/m_referenceInternalEnergy ),
                        m_density,
                        m_dDensity_dPressure,
                        m_dDensity_dTemperature,
                        m_viscosity,
                        m_dViscosity_dPressure,
                        m_dViscosity_dTemperature,
                        m_internalEnergy,
                        m_dInternalEnergy_dPressure,
                        m_dInternalEnergy_dTemperature,
                        m_enthalpy,
                        m_dEnthalpy_dPressure,
                        m_dEnthalpy_dTemperature,
                        m_referenceInternalEnergy );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermalCompressibleSinglePhaseFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
