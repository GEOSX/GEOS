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

#include "BlackOilFluidBase.hpp"

#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

namespace constitutive
{

BlackOilFluidBase::BlackOilFluidBase( string const & name,
                                      Group * const parent )
  :
  MultiFluidBase( name, parent )
{
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::surfacePhaseMassDensitiesString(), &m_surfacePhaseMassDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of surface mass densities for each phase" );

  // 1) First option: specify PVT tables from one file per phase, read the files line by line, and populate the internal TableFunctions

  registerWrapper( viewKeyStruct::tableFilesString(), &m_tableFiles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of filenames with input PVT tables (one per phase)" );

  // 2) Second option: specify TableFunction names for oil and gas,
  //                   and then specify water reference pressure, water Bw, water compressibility and water viscosity

  registerWrapper( viewKeyStruct::formationVolumeFactorTableNamesString(), &m_formationVolFactorTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of formation volume factor TableFunction names from the Functions block. \n"
                    "The user must provide one TableFunction per hydrocarbon phase, in the order provided in \"phaseNames\". \n"
                    "For instance, if \"oil\" is before \"gas\" in \"phaseNames\", the table order should be: oilTableName, gasTableName" );
  registerWrapper( viewKeyStruct::viscosityTableNamesString(), &m_viscosityTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of viscosity TableFunction names from the Functions block. \n"
                    "The user must provide one TableFunction per hydrocarbon phase, in the order provided in \"phaseNames\". \n"
                    "For instance, if \"oil\" is before \"gas\" in \"phaseNames\", the table order should be: oilTableName, gasTableName" );
  registerWrapper( viewKeyStruct::waterRefPressureString(), &m_waterParams.referencePressure ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water reference pressure" );
  registerWrapper( viewKeyStruct::waterFormationVolumeFactorString(), &m_waterParams.formationVolFactor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water formation volume factor" );
  registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterParams.compressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water compressibility" );
  registerWrapper( viewKeyStruct::waterViscosityString(), &m_waterParams.viscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Water viscosity" );

  // Register extra wrappers to enable auto-cloning
  registerWrapper( "phaseTypes", &m_phaseTypes )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "phaseOrder", &m_phaseOrder )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "hydrocarbonPhaseOrder", &m_hydrocarbonPhaseOrder )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "formationVolFactorTableWrappers", &m_formationVolFactorTableKernels )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( "viscosityTableWrappers", &m_viscosityTableKernels )
    .setSizedFromParent( 0 )
    .setRestartFlags( RestartFlags::NO_WRITE );
}

void BlackOilFluidBase::fillWaterData( array1d< array1d< real64 > > const & tableValues )
{
  GEOS_THROW_IF_NE_MSG( tableValues.size(), 1,
                        getFullName() << ": the water table must contain one line and only one",
                        InputError );
  GEOS_THROW_IF_NE_MSG( tableValues[0].size(), 4,
                        getFullName() << ": four columns (pressure, formation volume factor, compressibility, and viscosity) are expected for water",
                        InputError );

  GEOS_THROW_IF( m_waterParams.referencePressure > 0.0 || m_waterParams.formationVolFactor > 0.0 ||
                 m_waterParams.compressibility > 0.0 || m_waterParams.viscosity > 0.0,
                 getFullName() << ": input is redundant (user provided both water data and a water pvt file)",
                 InputError );

  m_waterParams.referencePressure = tableValues[0][0];
  m_waterParams.formationVolFactor = tableValues[0][1];
  m_waterParams.compressibility = tableValues[0][2];
  m_waterParams.viscosity = tableValues[0][3];

  validateWaterParams();
}

void BlackOilFluidBase::fillHydrocarbonData( integer const ip,
                                             array1d< array1d< real64 > > const & tableValues )
{
  array1d< array1d< real64 > > pressureCoords( 1 );
  pressureCoords[0].resize( tableValues.size() );
  array1d< real64 > formationVolFactor( tableValues.size() );
  array1d< real64 > viscosity( tableValues.size() );

  for( localIndex i = 0; i < tableValues.size(); ++i )
  {
    GEOS_THROW_IF_NE_MSG( tableValues[i].size(), 3,
                          getFullName() << ": three columns (pressure, formation volume factor, and viscosity) are expected for oil and gas",
                          InputError );

    pressureCoords[0][i] = tableValues[i][0];
    formationVolFactor[i] = tableValues[i][1];
    viscosity[i] = tableValues[i][2];
  }

  m_hydrocarbonPhaseOrder.emplace_back( ip );
  string const formationVolFactorTableName = getName() + (m_phaseTypes[ip] == PhaseType::OIL ? "_PVDO_Bo" : "_PVDG_Bg");
  string const viscosityTableName = getName() +  (m_phaseTypes[ip] == PhaseType::OIL ? "_PVDO_visco" : "_PVDG_viscg");
  m_formationVolFactorTableNames.emplace_back( formationVolFactorTableName );
  m_viscosityTableNames.emplace_back( viscosityTableName );

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & tablePVDX_B =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", formationVolFactorTableName ) );
  tablePVDX_B.setTableCoordinates( pressureCoords, { units::Pressure } );
  tablePVDX_B.setTableValues( formationVolFactor, units::Dimensionless );
  tablePVDX_B.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDX_visc =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", viscosityTableName ) );
  tablePVDX_visc.setTableCoordinates( pressureCoords, { units::Pressure } );
  tablePVDX_visc.setTableValues( viscosity, units::Viscosity );
  tablePVDX_visc.setInterpolationMethod( TableFunction::InterpolationType::Linear );
}

integer BlackOilFluidBase::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] = { "water" };
  return PVTProps::PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

void BlackOilFluidBase::postInputInitialization()
{
  m_componentNames = m_phaseNames;

  MultiFluidBase::postInputInitialization();

  m_phaseTypes.resize( numFluidPhases() );
  m_phaseOrder.resizeDefault( MAX_NUM_PHASES, -1 );

  auto const toPhaseType = [&]( string const & lookup )
  {
    static unordered_map< string, integer > const phaseDict =
    {
      { "gas", BlackOilFluidBase::PhaseType::GAS },
      { "oil", BlackOilFluidBase::PhaseType::OIL },
      { "water", BlackOilFluidBase::PhaseType::WATER }
    };
    return findOption( phaseDict, lookup, viewKeyStruct::phaseNamesString(), getFullName() );
  };

  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    m_phaseTypes[ip] = toPhaseType( m_phaseNames[ip] );
    m_phaseOrder[m_phaseTypes[ip]] = ip;
  }

  // Number of components should be at least equal to number of phases
  GEOS_THROW_IF_LT_MSG( numFluidComponents(), numFluidPhases(),
                        GEOS_FMT( "{}: {} number of components ({}) must be at least number of phases ({})",
                                  getFullName(), viewKeyStruct::componentNamesString(),
                                  numFluidComponents(), numFluidPhases() ),
                        InputError );

  auto const checkInputSize = [&]( auto const & array, auto const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), m_phaseNames.size(),
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );
  };
  checkInputSize( m_surfacePhaseMassDensity, viewKeyStruct::surfacePhaseMassDensitiesString() );
  checkInputSize( m_componentMolarWeight, viewKeyStruct::componentMolarWeightString() );

  // we make the distinction between the two input options
  if( m_tableFiles.empty() )
  {
    readInputDataFromTableFunctions();
  }
  else
  {
    readInputDataFromPVTFiles();
  }
}

void BlackOilFluidBase::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();

  {
    FunctionManager const & functionManager = FunctionManager::getInstance();
    for( integer IP_HDRINCL = 0; IP_HDRINCL < m_hydrocarbonPhaseOrder.size(); ++IP_HDRINCL )
    {
      // grab the tables by name from the function manager
      TableFunction const & fvfTable = functionManager.getGroup< TableFunction const >( m_formationVolFactorTableNames[IP_HDRINCL] );
      TableFunction const & viscosityTable = functionManager.getGroup< TableFunction const >( m_viscosityTableNames[IP_HDRINCL] );
      // validate them, then add them in a list to create their table wrappers when needed
      validateTable( fvfTable, false );
      validateTable( viscosityTable, true );
      m_formationVolFactorTables.emplace_back( &fvfTable );
      m_viscosityTables.emplace_back( &viscosityTable );
    }
  }

  createAllKernelWrappers();
}

void BlackOilFluidBase::checkTablesParameters( real64 const pressure,
                                               real64 const temperature ) const
{
  constexpr std::string_view errorMsg = "{} {}: {} table reading error for hydrocarbon phase {}.\n";

  GEOS_UNUSED_VAR( temperature );

  if( !m_checkPVTTablesRanges )
  {
    return;
  }

  for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
  {
    try
    {
      m_formationVolFactorTables[iph]->checkCoord( pressure, 0 );
    } catch( SimulationError const & ex )
    {
      throw SimulationError( ex, GEOS_FMT( errorMsg, getCatalogName(), getDataContext(),
                                           "formation volume factor", iph ) );
    }

    try
    {
      m_viscosityTables[iph]->checkCoord( pressure, 0 );
    } catch( SimulationError const & ex )
    {
      throw SimulationError( ex, GEOS_FMT( errorMsg, getCatalogName(), getDataContext(),
                                           "viscosity", iph ) );
    }
  }
}

void BlackOilFluidBase::createAllKernelWrappers()
{
  GEOS_THROW_IF( m_hydrocarbonPhaseOrder.size() != 1 && m_hydrocarbonPhaseOrder.size() != 2,
                 GEOS_FMT( "{}: the number of hydrocarbon phases must be 1 (oil) or 2 (oil+gas)", getFullName() ),
                 InputError );

  if( m_formationVolFactorTableKernels.empty() && m_viscosityTableKernels.empty() )
  {
    // loop over the hydrocarbon phases
    for( integer iph = 0; iph < m_hydrocarbonPhaseOrder.size(); ++iph )
    {
      // create the table wrapper for the oil and (if present) the gas phases
      m_formationVolFactorTableKernels.emplace_back( m_formationVolFactorTables[iph]->createKernelWrapper() );
      m_viscosityTableKernels.emplace_back( m_viscosityTables[iph]->createKernelWrapper() );
    }
  }
}

void BlackOilFluidBase::validateTable( TableFunction const & table,
                                       bool warningIfDecreasing ) const
{
  arrayView1d< real64 const > const property = table.getValues();
  GEOS_THROW_IF_NE_MSG( table.getInterpolationMethod(), TableFunction::InterpolationType::Linear,
                        GEOS_FMT( "{}: in table '{}' interpolation method must be linear", getFullName(), table.getName() ),
                        InputError );
  GEOS_THROW_IF_LT_MSG( property.size(), 2,
                        GEOS_FMT( "{}: table '{}' must contain at least two values", getFullName(), table.getName() ),
                        InputError );

  // we don't check the first interval, as the first value may be used to specify surface conditions
  // we only issue a warning here, as we still want to allow this configuration
  for( localIndex i = 3; i < property.size(); ++i )
  {
    GEOS_THROW_IF( (property[i] - property[i-1]) * (property[i-1] - property[i-2]) < 0,
                   GEOS_FMT( "{}: in table '{}', viscosity values must be monotone", getFullName(), table.getName() ),
                   InputError );
  }

  // we don't check the first value, as it may be used to specify surface conditions
  for( localIndex i = 2; i < property.size(); ++i )
  {
    GEOS_LOG_RANK_0_IF( ( property[i] - property[i-1] < 0 ) && warningIfDecreasing,
                        GEOS_FMT( "{}: Warning! in table '{}', values must be increasing as a function of pressure, please check your PVT tables",
                                  getFullName(), table.getName() ) );
    GEOS_LOG_RANK_0_IF( ( property[i] - property[i-1] > 0 ) && !warningIfDecreasing,
                        GEOS_FMT( "{}: Warning! In table '{}', values must be decreasing as a function of pressure, please check your PVT tables",
                                  getFullName(), table.getName() ) );
  }
}


void BlackOilFluidBase::validateWaterParams() const
{
  auto const checkPositiveValue = [&]( real64 const value, auto const & attribute )
  {
    GEOS_THROW_IF_LE_MSG( value, 0.0,
                          GEOS_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                          InputError );
  };
  checkPositiveValue( m_waterParams.referencePressure, viewKeyStruct::waterRefPressureString() );
  checkPositiveValue( m_waterParams.formationVolFactor, viewKeyStruct::waterFormationVolumeFactorString() );
  checkPositiveValue( m_waterParams.compressibility, viewKeyStruct::waterCompressibilityString() );
  checkPositiveValue( m_waterParams.viscosity, viewKeyStruct::waterViscosityString() );
}

BlackOilFluidBase::KernelWrapper::
  KernelWrapper( arrayView1d< integer const > phaseTypes,
                 arrayView1d< integer const > phaseOrder,
                 arrayView1d< integer const > hydrocarbonPhaseOrder,
                 arrayView1d< real64 const > surfacePhaseMassDensity,
                 arrayView1d< TableFunction::KernelWrapper const > formationVolFactorTables,
                 arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
                 BlackOilFluidBase::WaterParams const waterParams,
                 arrayView1d< real64 const > componentMolarWeight,
                 bool const useMass,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( std::move( componentMolarWeight ),
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_phaseTypes( std::move( phaseTypes ) ),
  m_phaseOrder( std::move( phaseOrder ) ),
  m_hydrocarbonPhaseOrder( std::move( hydrocarbonPhaseOrder ) ),
  m_surfacePhaseMassDensity( std::move( surfacePhaseMassDensity ) ),
  m_formationVolFactorTables( std::move( formationVolFactorTables ) ),
  m_viscosityTables( std::move( viscosityTables ) ),
  m_waterParams( waterParams )
{}

} // namespace constitutive

} // namespace geos
