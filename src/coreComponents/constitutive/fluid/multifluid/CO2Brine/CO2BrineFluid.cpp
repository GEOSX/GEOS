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
 * @file CO2BrineFluid.cpp
 */
#include "CO2BrineFluid.hpp"

#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "common/Units.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

using namespace PVTProps;

namespace
{
template< typename PHASE1, typename PHASE2, typename FLASH > class
  TwoPhaseCatalogNames {};

template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrinePhillipsThermalFluid"; }
};

template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::NoOpPVTFunction >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiFluid"; }
};
template<> class
  TwoPhaseCatalogNames< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::BrineEnthalpy >,
                        PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                        PVTProps::CO2Solubility >
{
public:
  static string name() { return "CO2BrineEzrokhiThermalFluid"; }
};

} // end namespace

// provide a definition for catalogName()
template< typename PHASE1, typename PHASE2, typename FLASH >
string CO2BrineFluid< PHASE1, PHASE2, FLASH >::catalogName()
{
  return TwoPhaseCatalogNames< PHASE1, PHASE2, FLASH >::name();
}

template< typename PHASE1, typename PHASE2, typename FLASH >
CO2BrineFluid< PHASE1, PHASE2, FLASH >::
CO2BrineFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::phasePVTParaFilesString(), &m_phasePVTParaFiles ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Names of the files defining the parameters of the viscosity and density models" );

  registerWrapper( viewKeyStruct::flashModelParaFileString(), &m_flashModelParaFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Name of the file defining the parameters of the flash model" );

  registerWrapper( viewKeyStruct::solubilityTablesString(), &m_solubilityTables ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Names of solubility tables for each phase" );

  // if this is a thermal model, we need to make sure that the arrays will be properly displayed and saved to restart
  if( isThermal() )
  {
    getField< fields::multifluid::phaseEnthalpy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );

    getField< fields::multifluid::phaseInternalEnergy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );
  }
}

template< typename PHASE1, typename PHASE2, typename FLASH >
std::unique_ptr< ConstitutiveBase >
CO2BrineFluid< PHASE1, PHASE2, FLASH >::
deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  CO2BrineFluid & newConstitutiveRelation = dynamicCast< CO2BrineFluid & >( *clone );
  newConstitutiveRelation.m_p1Index = m_p1Index;
  newConstitutiveRelation.m_p2Index = m_p2Index;

  newConstitutiveRelation.createPVTModels();

  return clone;
}

template< typename PHASE1, typename PHASE2, typename FLASH >
integer CO2BrineFluid< PHASE1, PHASE2, FLASH >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] =  { "Water", "water", "Liquid", "liquid" };
  return PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::checkTablesParameters( real64 const pressure,
                                                                    real64 const temperature ) const
{
  if( !m_checkPVTTablesRanges )
  {
    return;
  }

  real64 const temperatureInCelsius = units::convertKToC( temperature );
  try
  {
    m_phase1->density.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase1->viscosity.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase1->enthalpy.checkTablesParameters( pressure, temperatureInCelsius );
  } catch( SimulationError const & ex )
  {
    string const errorMsg = GEOS_FMT( "Table input error for {} phase (in table from \"{}\").\n",
                                      m_phaseNames[m_p1Index], m_phasePVTParaFiles[m_p1Index] );
    throw SimulationError( ex, errorMsg );
  }

  try
  {
    m_phase2->density.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase2->viscosity.checkTablesParameters( pressure, temperatureInCelsius );
    m_phase2->enthalpy.checkTablesParameters( pressure, temperatureInCelsius );
  } catch( SimulationError const & ex )
  {
    string const errorMsg = GEOS_FMT( "Table input error for {} phase (in table from \"{}\").\n",
                                      m_phaseNames[m_p2Index], m_phasePVTParaFiles[m_p2Index] );
    throw SimulationError( ex, errorMsg );
  }

  try
  {
    m_flash->checkTablesParameters( pressure, temperatureInCelsius );
  } catch( SimulationError const & ex )
  {
    string const errorMsg = GEOS_FMT( "Table input error for flash phase (in table from \"{}\").\n",
                                      m_flashModelParaFile );
    throw SimulationError( ex, errorMsg );
  }
}


template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::initializePreSubGroups()
{
#if defined(GEOS_DEVICE_COMPILE)
  GEOS_THROW_IF( this->getCatalogName() == CO2BrineEzrokhiThermalFluid::catalogName(),
                 GEOS_FMT( "The `{}` model is disabled for now. Please use the other thermal CO2-brine model instead: `{}`",
                           CO2BrineEzrokhiThermalFluid::catalogName(),
                           CO2BrinePhillipsThermalFluid::catalogName() ),
                 InputError );
#endif
}

template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  GEOS_THROW_IF_NE_MSG( numFluidPhases(), 2,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( numFluidComponents(), 2,
                        GEOS_FMT( "{}: invalid number of components", getFullName() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( m_phasePVTParaFiles.size(), 2,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName() ),
                        InputError );

  // Make sure one (and only one) of m_flashModelParaFile or m_solubilityTables is provided
  bool const hasParamFile = !m_flashModelParaFile.empty();
  bool const hasTables = !m_solubilityTables.empty();
  GEOS_THROW_IF( hasParamFile == hasTables,
                 GEOS_FMT( "{}: One and only one of {} or {} should be specified", getFullName(),
                           viewKeyStruct::flashModelParaFileString(),
                           viewKeyStruct::solubilityTablesString() ),
                 InputError );

  // NOTE: for now, the names of the phases are still hardcoded here
  // Later, we could read them from the XML file and we would then have a general class here

  string const expectedWaterPhaseNames[] = { "Water", "water", "Liquid", "liquid" };
  m_p1Index = PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );

  string const expectedGasPhaseNames[] = { "CO2", "co2", "gas", "Gas" };
  m_p2Index = PVTFunctionHelpers::findName( m_phaseNames, expectedGasPhaseNames, viewKeyStruct::phaseNamesString() );

  createPVTModels();
}

template< typename PHASE1, typename PHASE2, typename FLASH >
void CO2BrineFluid< PHASE1, PHASE2, FLASH >::createPVTModels()
{

  // TODO: get rid of these external files and move into XML, this is too error prone
  // For now, to support the legacy input, we read all the input parameters at once in the arrays below, and then we create the models
  array1d< array1d< string > > phase1InputParams;
  phase1InputParams.resize( 3 );
  array1d< array1d< string > > phase2InputParams;
  phase2InputParams.resize( 3 );

  // 1) Create the viscosity, density, enthalpy models
  for( string const & filename : m_phasePVTParaFiles )
  {
    std::ifstream is( filename );
    string str;
    while( std::getline( is, str ) )
    {
      array1d< string > const strs = stringutilities::tokenizeBySpaces< array1d >( str );

      if( !strs.empty() )
      {
        GEOS_THROW_IF( strs.size() < 2,
                       GEOS_FMT( "{}: missing PVT model in line '{}'", getFullName(), str ),
                       InputError );

        if( strs[0] == "DensityFun" )
        {
          if( strs[1] == PHASE1::Density::catalogName() )
          {
            phase1InputParams[PHASE1::InputParamOrder::DENSITY] = strs;
          }
          else if( strs[1] == PHASE2::Density::catalogName() )
          {
            phase2InputParams[PHASE2::InputParamOrder::DENSITY] = strs;
          }
        }
        else if( strs[0] == "ViscosityFun" )
        {
          if( strs[1] == PHASE1::Viscosity::catalogName() )
          {
            phase1InputParams[PHASE1::InputParamOrder::VISCOSITY] = strs;
          }
          else if( strs[1] == PHASE2::Viscosity::catalogName() )
          {
            phase2InputParams[PHASE2::InputParamOrder::VISCOSITY] = strs;
          }
        }
        else if( strs[0] == "EnthalpyFun" )
        {
          if( strs[1] == PHASE1::Enthalpy::catalogName() )
          {
            phase1InputParams[PHASE1::InputParamOrder::ENTHALPY] = strs;
          }
          else if( strs[1] == PHASE2::Enthalpy::catalogName() )
          {
            phase2InputParams[PHASE2::InputParamOrder::ENTHALPY] = strs;
          }
        }
        else
        {
          GEOS_THROW( GEOS_FMT( "{}: invalid PVT function type '{}'", getFullName(), strs[0] ), InputError );
        }
      }
    }
    is.close();
  }

  // at this point, we have read the file and we check the consistency of non-thermal models
  GEOS_THROW_IF( phase1InputParams[PHASE1::InputParamOrder::DENSITY].empty(),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE1::Density::catalogName() ),
                 InputError );
  GEOS_THROW_IF( phase2InputParams[PHASE2::InputParamOrder::DENSITY].empty(),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE2::Density::catalogName() ),
                 InputError );
  GEOS_THROW_IF( phase1InputParams[PHASE1::InputParamOrder::VISCOSITY].empty(),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE1::Viscosity::catalogName() ),
                 InputError );
  GEOS_THROW_IF( phase2InputParams[PHASE2::InputParamOrder::VISCOSITY].empty(),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE2::Viscosity::catalogName() ),
                 InputError );

  // we also detect any inconsistency arising in the enthalpy models
  GEOS_THROW_IF( phase1InputParams[PHASE1::InputParamOrder::ENTHALPY].empty() &&
                 ( PHASE1::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE1::Enthalpy::catalogName() ),
                 InputError );
  GEOS_THROW_IF( phase2InputParams[PHASE2::InputParamOrder::ENTHALPY].empty() &&
                 ( PHASE2::Enthalpy::catalogName() != PVTProps::NoOpPVTFunction::catalogName() ),
                 GEOS_FMT( "{}: PVT model {} not found in input files", getFullName(), PHASE2::Enthalpy::catalogName() ),
                 InputError );

  // then, we are ready to instantiate the phase models
  m_phase1 = std::make_unique< PHASE1 >( getName() + "_phaseModel1", phase1InputParams, m_componentNames, m_componentMolarWeight,
                                         getLogLevel() > 0 && logger::internal::rank==0 );
  m_phase2 = std::make_unique< PHASE2 >( getName() + "_phaseModel2", phase2InputParams, m_componentNames, m_componentMolarWeight,
                                         getLogLevel() > 0 && logger::internal::rank==0 );

  // 2) Create the flash model
  if( !m_flashModelParaFile.empty())
  {
    std::ifstream is( m_flashModelParaFile );
    string str;
    while( std::getline( is, str ) )
    {
      string_array const strs = stringutilities::tokenizeBySpaces< array1d >( str );

      if( !strs.empty() )
      {
        GEOS_THROW_IF( strs.size() < 2,
                       GEOS_FMT( "{}: missing flash model in line '{}'", getFullName(), str ),
                       InputError );

        if( strs[0] == "FlashModel" )
        {
          if( strs[1] == FLASH::catalogName() )
          {
            m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                                 strs,
                                                 m_phaseNames,
                                                 m_componentNames,
                                                 m_componentMolarWeight,
                                                 getLogLevel() > 0 && logger::internal::rank==0 );
          }
        }
        else
        {
          GEOS_THROW( GEOS_FMT( "{}: invalid flash model type '{}'", getFullName(), strs[0] ), InputError );
        }
      }
    }
    is.close();
  }
  else
  {
    // The user must provide 1 or 2 tables.
    GEOS_THROW_IF( m_solubilityTables.size() != 1 && m_solubilityTables.size() != 2,
                   GEOS_FMT( "{}: The number of table names in {} must be 1 or 2", getFullName(), viewKeyStruct::solubilityTablesString() ),
                   InputError );

    // If 1 table is provided, it is the CO2 solubility table and water vapourisation is zero
    // If 2 tables are provided, they are the CO2 solubility and water vapourisation tables depending
    // on how phaseNames is arranged
    string const solubilityModel = EnumStrings< CO2Solubility::SolubilityModel >::toString( CO2Solubility::SolubilityModel::Tables );
    string_array strs;
    strs.emplace_back( "FlashModel" );
    strs.emplace_back( solubilityModel );   // Marker to indicate that tables are provided
    strs.emplace_back( "" );   // 2 empty strings for the 2 phase tables gas first, then water
    strs.emplace_back( "" );
    if( m_solubilityTables.size() == 2 )
    {
      strs[2] = m_solubilityTables[m_p2Index];
      strs[3] = m_solubilityTables[m_p1Index];
    }
    else
    {
      strs[2] = m_solubilityTables[0];
    }
    m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                         strs,
                                         m_phaseNames,
                                         m_componentNames,
                                         m_componentMolarWeight,
                                         getLogLevel() > 0 && logger::internal::rank==0 );
  }

  GEOS_THROW_IF( m_flash == nullptr,
                 GEOS_FMT( "{}: flash model {} not found in input files", getFullName(), FLASH::catalogName() ),
                 InputError );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
typename CO2BrineFluid< PHASE1, PHASE2, FLASH >::KernelWrapper
CO2BrineFluid< PHASE1, PHASE2, FLASH >::createKernelWrapper()
{
  return KernelWrapper( m_p1Index,
                        m_p2Index,
                        *m_phase1,
                        *m_phase2,
                        *m_flash,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        isThermal(),
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

template< typename PHASE1, typename PHASE2, typename FLASH >
CO2BrineFluid< PHASE1, PHASE2, FLASH >::KernelWrapper::
  KernelWrapper( integer const p1Index,
                 integer const p2Index,
                 PHASE1 const & phase1,
                 PHASE2 const & phase2,
                 FLASH const & flash,
                 arrayView1d< geos::real64 const > componentMolarWeight,
                 bool const useMass,
                 bool const isThermal,
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
  m_p1Index( p1Index ),
  m_p2Index( p2Index ),
  m_isThermal( isThermal ),
  m_phase1( phase1.createKernelWrapper() ),
  m_phase2( phase2.createKernelWrapper() ),
  m_flash( flash.createKernelWrapper() )
{}

// explicit instantiation of the model template; unfortunately we can't use the aliases for this
template class CO2BrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::NoOpPVTFunction >,
                              PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                              PVTProps::CO2Solubility >;
template class CO2BrineFluid< PhaseModel< PVTProps::PhillipsBrineDensity, PVTProps::PhillipsBrineViscosity, PVTProps::BrineEnthalpy >,
                              PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                              PVTProps::CO2Solubility >;

template class CO2BrineFluid< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::NoOpPVTFunction >,
                              PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::NoOpPVTFunction >,
                              PVTProps::CO2Solubility >;
template class CO2BrineFluid< PhaseModel< PVTProps::EzrokhiBrineDensity, PVTProps::EzrokhiBrineViscosity, PVTProps::BrineEnthalpy >,
                              PhaseModel< PVTProps::SpanWagnerCO2Density, PVTProps::FenghourCO2Viscosity, PVTProps::CO2Enthalpy >,
                              PVTProps::CO2Solubility >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrinePhillipsFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrinePhillipsThermalFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineEzrokhiFluid, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, CO2BrineEzrokhiThermalFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geos
