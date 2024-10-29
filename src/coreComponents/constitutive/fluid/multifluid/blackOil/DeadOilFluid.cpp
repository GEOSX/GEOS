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

#include "DeadOilFluid.hpp"

#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;
using namespace stringutilities;
namespace constitutive
{

DeadOilFluid::DeadOilFluid( string const & name,
                            Group * const parent )
  :
  BlackOilFluidBase( name, parent )
{}

void DeadOilFluid::readInputDataFromPVTFiles()
{
  GEOS_THROW_IF_NE_MSG( m_tableFiles.size(), numFluidPhases(),
                        GEOS_FMT( "{}: the number of table files must be equal to the number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF( m_formationVolFactorTableNames.size() > 0.0 || m_viscosityTableNames.size() > 0.0,
                 GEOS_FMT( "{}: input is redundant (both TableFunction names and pvt files)", getFullName() ),
                 InputError );

  array1d< array1d< real64 > > tableValues;
  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    tableValues.clear();
    PVTProps::BlackOilTables::readTable( m_tableFiles[ip], 3, tableValues );

    if( m_phaseTypes[ip] == PhaseType::WATER )
    {
      fillWaterData( tableValues );
    }
    else
    {
      fillHydrocarbonData( ip, tableValues );
    }
  }
}

void DeadOilFluid::readInputDataFromTableFunctions()
{
  GEOS_THROW_IF( !m_tableFiles.empty(),
                 GEOS_FMT( "{}: input is redundant (both TableFunction names and pvt files)", getFullName() ),
                 InputError );

  integer const ipWater = m_phaseOrder[PhaseType::WATER];
  integer const ipGas = m_phaseOrder[PhaseType::GAS];
  if( ipWater >= 0 ) // if water is present
  {
    validateWaterParams();
  }
  else
  {
    auto const errorIfPositiveValue = [&]( real64 const value, auto const & attribute )
    {
      GEOS_THROW_IF_GT_MSG( value, 0.0,
                            GEOS_FMT( "{}: if water is absent, attribute '{}' is not redundant", getFullName(), attribute ),
                            InputError );
    };
    errorIfPositiveValue( m_waterParams.referencePressure, viewKeyStruct::waterRefPressureString() );
    errorIfPositiveValue( m_waterParams.formationVolFactor, viewKeyStruct::waterFormationVolumeFactorString() );
    errorIfPositiveValue( m_waterParams.viscosity, viewKeyStruct::waterViscosityString() );
    errorIfPositiveValue( m_waterParams.compressibility, viewKeyStruct::waterCompressibilityString() );
  }

  integer const numExpectedTables = (ipGas >= 0) ? 2 : 1;
  GEOS_THROW_IF_NE_MSG( m_formationVolFactorTableNames.size(), numExpectedTables,
                        GEOS_FMT( "{}: one formation volume factor table must be provided for each hydrocarbon phase", getFullName() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( m_viscosityTableNames.size(), numExpectedTables,
                        GEOS_FMT( "{}: one viscosity table must be provided for each hydrocarbon phase", getFullName() ),
                        InputError );

  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    if( m_phaseTypes[ip] == PhaseType::OIL || m_phaseTypes[ip] == PhaseType::GAS )
    {
      m_hydrocarbonPhaseOrder.emplace_back( ip );
    }
  }

  FunctionManager const & functionManager = FunctionManager::getInstance();
  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    GEOS_THROW_IF( !functionManager.hasGroup( m_formationVolFactorTableNames[iph] ),
                   GEOS_FMT( "{}: formation volume factor table '{}' not found", getFullName(), m_formationVolFactorTableNames[iph] ),
                   InputError );
    GEOS_THROW_IF( !functionManager.hasGroup( m_viscosityTableNames[iph] ),
                   GEOS_FMT( "{}: viscosity table '{}' not found", getFullName(), m_viscosityTableNames[iph] ),
                   InputError );
  }
}

DeadOilFluid::KernelWrapper::
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
                 FluidProp::ViewType totalDensity )
  : BlackOilFluidBase::KernelWrapper( std::move( phaseTypes ),
                                      std::move( phaseOrder ),
                                      std::move( hydrocarbonPhaseOrder ),
                                      std::move( surfacePhaseMassDensity ),
                                      std::move( formationVolFactorTables ),
                                      std::move( viscosityTables ),
                                      waterParams,
                                      std::move( componentMolarWeight ),
                                      useMass,
                                      std::move( phaseFraction ),
                                      std::move( phaseDensity ),
                                      std::move( phaseMassDensity ),
                                      std::move( phaseViscosity ),
                                      std::move( phaseEnthalpy ),
                                      std::move( phaseInternalEnergy ),
                                      std::move( phaseCompFraction ),
                                      std::move( totalDensity ) )
{}

DeadOilFluid::KernelWrapper
DeadOilFluid::createKernelWrapper()
{
  return KernelWrapper( m_phaseTypes,
                        m_phaseOrder,
                        m_hydrocarbonPhaseOrder,
                        m_surfacePhaseMassDensity,
                        m_formationVolFactorTableKernels,
                        m_viscosityTableKernels,
                        m_waterParams,
                        m_componentMolarWeight,
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DeadOilFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geos
