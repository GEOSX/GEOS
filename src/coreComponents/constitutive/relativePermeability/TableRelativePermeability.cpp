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
 * @file TableRelativePermeability.cpp
 */

#include "TableRelativePermeability.hpp"
#include "constitutive/relativePermeability/TableRelativePermeabilityHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

TableRelativePermeability::TableRelativePermeability( std::string const & name,
                                                      Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::wettingNonWettingRelPermTableNamesString(), &m_wettingNonWettingRelPermTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupOfGroupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (wetting phase, non-wetting phase)\n"
                    "The expected format is \"{ wettingPhaseRelPermTableName, nonWettingPhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for two-phase flow.\n"
                    "If you want to do a three-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingIntermediateRelPermTableNamesString() ) +
                    " and " +
                    string( viewKeyStruct::nonWettingIntermediateRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::wettingIntermediateRelPermTableNamesString(), &m_wettingIntermediateRelPermTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupOfGroupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (wetting phase, intermediate phase)\n"
                    "The expected format is \"{ wettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::nonWettingIntermediateRelPermTableNamesString(), &m_nonWettingIntermediateRelPermTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupOfGroupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of relative permeability tables for the pair (non-wetting phase, intermediate phase)\n"
                    "The expected format is \"{ nonWettingPhaseRelPermTableName, intermediatePhaseRelPermTableName }\", in that order\n"
                    "Note that this input is only used for three-phase flow.\n"
                    "If you want to do a two-phase simulation, please use instead " +
                    string( viewKeyStruct::wettingNonWettingRelPermTableNamesString() ) +
                    " to specify the table names" );

  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::waterOilMaxRelPermString(), &m_waterOilMaxRelPerm ).
    setInputFlag( InputFlags::FALSE ). // will be deduced from tables
    setApplyDefaultValue( 0.0 ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::relPermKernelWrappersString(), &m_relPermKernelWrappers ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::threePhaseInterpolatorString(), &m_threePhaseInterpolator ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( ThreePhaseInterpolator::BAKER ).
    setDescription( "Type of Three phase interpolator."
                    "Valid options \n* " + EnumStrings< ThreePhaseInterpolator >::concat( "\n* " ) );
}

void TableRelativePermeability::postInputInitialization()
{
  RelativePermeabilityBase::postInputInitialization();

  integer const numPhases = m_phaseNames.size();
  integer const numDir = m_wettingNonWettingRelPermTableNames.size(0);
  //reshape Name containers


  GEOS_THROW_IF( numPhases != 2 && numPhases != 3,
                 GEOS_FMT( "{}: the expected number of fluid phases is either two, or three",
                           getFullName() ),
                 InputError );

  for( int dir=0; dir < numDir; ++dir )
  {
    if( numPhases == 2 )
    {
      m_wettingNonWettingRelPermTableNames.resize( numDir, numPhases );

      GEOS_THROW_IF( m_wettingNonWettingRelPermTableNames[dir][0].empty() || m_wettingNonWettingRelPermTableNames[dir][1].empty(),
                     GEOS_FMT(
                       "{}: for a two-phase flow simulation, we must use {} to specify the relative permeability tables for the pair (wetting phase, non-wetting phase)",
                       getFullName(),
                       viewKeyStruct::wettingNonWettingRelPermTableNamesString()),
                     InputError );

      GEOS_THROW_IF( m_wettingNonWettingRelPermTableNames[dir].size() != 2,
                     GEOS_FMT(
                       "{}: for a two-phase flow simulation, we must use {} to specify exactly two names: first the name of the wetting phase relperm table, second the name on the non-wetting phase relperm table",
                       getFullName(),
                       viewKeyStruct::wettingNonWettingRelPermTableNamesString()),
                     InputError );

    }
    else if( numPhases == 3 )
    {
      m_wettingIntermediateRelPermTableNames.resize( numDir, 2 );
      m_nonWettingIntermediateRelPermTableNames.resize( numDir, 2 );
      GEOS_THROW_IF( m_wettingIntermediateRelPermTableNames[dir][0].empty() || m_wettingIntermediateRelPermTableNames[dir][1].empty()
                     || m_nonWettingIntermediateRelPermTableNames[dir][0].empty() || m_nonWettingIntermediateRelPermTableNames[dir][1].empty(),
                     GEOS_FMT(
                       "{}: for a three-phase flow simulation, we must use {} to specify the relative permeability tables for the pair (wetting phase, intermediate phase), and {} to specify the relative permeability tables for the pair (non-wetting phase, intermediate phase)",
                       getFullName(),
                       viewKeyStruct::wettingIntermediateRelPermTableNamesString(),
                       viewKeyStruct::nonWettingIntermediateRelPermTableNamesString()),
                     InputError );

      GEOS_THROW_IF( m_wettingIntermediateRelPermTableNames[dir].size() != 2,
                     GEOS_FMT(
                       "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: first the name of the wetting phase relperm table, second the name on the intermediate phase relperm table",
                       getFullName(),
                       viewKeyStruct::wettingIntermediateRelPermTableNamesString()),
                     InputError );

      GEOS_THROW_IF( m_nonWettingIntermediateRelPermTableNames[dir].size() != 2,
                     GEOS_FMT(
                       "{}: for a three-phase flow simulation, we must use {} to specify exactly two names: first the name of the non-wetting phase relperm table, second the name on the intermediate phase relperm table",
                       getFullName(),
                       viewKeyStruct::nonWettingIntermediateRelPermTableNamesString()),
                     InputError );
    }
  }
}

void TableRelativePermeability::resizeFields( localIndex const size, localIndex const numPts )
{
  RelativePermeabilityBase::resizeFields( size, numPts );

  integer const numPhases = numFluidPhases();
  integer const numDir =  m_wettingNonWettingRelPermTableNames.size(0);


  m_phaseRelPerm.resize( size, numPts, numPhases, numDir );
  m_phaseRelPerm_n.resize( size, numPts, numPhases, numDir );
  m_dPhaseRelPerm_dPhaseVolFrac.resize( size, numPts, numPhases, numPhases, numDir );
  //phase trapped for stats
  m_phaseTrappedVolFrac.resize( size, numPts, numPhases );
  m_phaseTrappedVolFrac.zero();


}


void TableRelativePermeability::initializePreSubGroups()
{
  RelativePermeabilityBase::initializePreSubGroups();
  
  integer const numDir = m_wettingNonWettingRelPermTableNames.size(0);

  integer const numPhases = m_phaseNames.size();
  m_phaseMinVolumeFraction.resize( numDir, MAX_NUM_PHASES );


  string const fullName = getFullName();
  real64 phaseMinVolFrac = 0.0;
  real64 phaseMaxVolFrac = 0.0;
  real64 phaseRelPermEndPoint = 0.0;

  //initialize STONE-II only used var to avoid discrepancies in baselines
  m_waterOilMaxRelPerm = 1.0;

  FunctionManager const & functionManager = FunctionManager::getInstance();

  for( int dir=0; dir<numDir; ++dir )
  {
    if( numPhases == 2 )
    {
      for( integer ip = 0; ip < m_wettingNonWettingRelPermTableNames[dir].size(); ++ip )
      {
        GEOS_THROW_IF( !functionManager.hasGroup( m_wettingNonWettingRelPermTableNames[dir][ip] ),
                       GEOS_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_wettingNonWettingRelPermTableNames[dir][ip] ),
                       InputError );
        TableFunction const & relPermTable = functionManager.getGroup< TableFunction >(
          m_wettingNonWettingRelPermTableNames[dir][ip] );
        TableRelativePermeabilityHelpers::
          validateRelativePermeabilityTable( relPermTable,      // input
                                             fullName,
                                             phaseMinVolFrac,      // output
                                             phaseMaxVolFrac,
                                             phaseRelPermEndPoint );
        if( ip == 0 )        // wetting phase is either water, or oil (for two-phase oil-gas systems)
        {
          integer const ipWetting = (m_phaseOrder[PhaseType::WATER] >= 0) ? m_phaseOrder[PhaseType::WATER]
                                                                                    : m_phaseOrder[PhaseType::OIL];
          m_phaseMinVolumeFraction[dir][ipWetting] = phaseMinVolFrac;
        }
        else if( ip == 1 )          // non-wetting phase is either oil (for two-phase oil-water systems), or gas
        {
          integer const ipNonWetting = (m_phaseOrder[PhaseType::GAS] >= 0) ? m_phaseOrder[PhaseType::GAS]
                                                                                     : m_phaseOrder[PhaseType::OIL];
          m_phaseMinVolumeFraction[dir][ipNonWetting] = phaseMinVolFrac;
        }
      }
    }
    else if( numPhases == 3 )
    {
      for( integer ip = 0; ip < m_wettingIntermediateRelPermTableNames[dir].size(); ++ip )
      {
        GEOS_THROW_IF( !functionManager.hasGroup( m_wettingIntermediateRelPermTableNames[dir][ip] ),
                       GEOS_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_wettingIntermediateRelPermTableNames[dir][ip] ),
                       InputError );
        TableFunction const & relPermTable = functionManager.getGroup< TableFunction >(
          m_wettingIntermediateRelPermTableNames[dir][ip] );
        TableRelativePermeabilityHelpers::
          validateRelativePermeabilityTable( relPermTable,      // input
                                             fullName,
                                             phaseMinVolFrac,      // output
                                             phaseMaxVolFrac,
                                             phaseRelPermEndPoint );


        if( ip == 0 )        // wetting phase is water
        {
          m_phaseMinVolumeFraction[dir][m_phaseOrder[PhaseType::WATER]] = phaseMinVolFrac;
        }
        else if( ip == 1 )          // intermediate phase is oil
        {
          m_phaseMinVolumeFraction[dir][m_phaseOrder[PhaseType::OIL]] = phaseMinVolFrac;
          m_waterOilMaxRelPerm = phaseRelPermEndPoint;
        }
      }
      for( integer ip = 0; ip < m_nonWettingIntermediateRelPermTableNames[dir].size(); ++ip )
      {
        GEOS_THROW_IF( !functionManager.hasGroup( m_nonWettingIntermediateRelPermTableNames[dir][ip] ),
                       GEOS_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_nonWettingIntermediateRelPermTableNames[dir][ip] ),
                       InputError );
        TableFunction const & relPermTable = functionManager.getGroup< TableFunction >(
          m_nonWettingIntermediateRelPermTableNames[dir][ip] );
        TableRelativePermeabilityHelpers::
          validateRelativePermeabilityTable( relPermTable,      // input
                                             fullName,
                                             phaseMinVolFrac,      // output
                                             phaseMaxVolFrac,
                                             phaseRelPermEndPoint );

        if( ip == 0 )        // non-wetting phase is gas
        {
          m_phaseMinVolumeFraction[dir][m_phaseOrder[PhaseType::GAS]] = phaseMinVolFrac;
        }
        else if( ip == 1 )          // intermediate phase is oil
        {
          m_phaseMinVolumeFraction[dir][m_phaseOrder[PhaseType::OIL]] = phaseMinVolFrac;
        }
      }
    }
  }
}

void TableRelativePermeability::createAllTableKernelWrappers()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  integer const numPhases = m_phaseNames.size();
  integer const numDir =  m_wettingNonWettingRelPermTableNames.size(0);
  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_relPermKernelWrappers.clear();

  for( int dir=0; dir<numDir; ++dir )
  {
    if( numPhases == 2 )
    {
      m_relPermKernelWrappers.resize( numDir, numPhases );
      for( integer ip = 0; ip < m_wettingNonWettingRelPermTableNames[dir].size(); ++ip )
      {
        TableFunction const & relPermTable = functionManager.getGroup< TableFunction >(
          m_wettingNonWettingRelPermTableNames[dir][ip] );
        m_relPermKernelWrappers[dir][ip] = relPermTable.createKernelWrapper();
      }
    }
    else if( numPhases == 3 )
    {
      m_relPermKernelWrappers.resize( numDir, 4 );     //because of TPT indirection
      for( integer ip = 0; ip < m_wettingIntermediateRelPermTableNames[dir].size(); ++ip )
      {
        TableFunction const & relPermTable = functionManager.getGroup< TableFunction >(
          m_wettingIntermediateRelPermTableNames[dir][ip] );
        m_relPermKernelWrappers[dir][ip] = relPermTable.createKernelWrapper();
      }
      for( integer ip = 0; ip < m_nonWettingIntermediateRelPermTableNames[dir].size(); ++ip )
      {
        TableFunction const & relPermTable = functionManager.getGroup< TableFunction >(
          m_nonWettingIntermediateRelPermTableNames[dir][ip] );
        m_relPermKernelWrappers[dir][2 + ip] = relPermTable.createKernelWrapper();
      }
    }
  }
}

TableRelativePermeability::KernelWrapper::
  KernelWrapper( arrayView2d< TableFunction::KernelWrapper const > const & relPermKernelWrappers,
                 arrayView2d< real64 const > const & phaseMinVolumeFraction,
                 real64 const & waterPhaseMaxVolumeFraction,
                 arrayView1d< integer const > const & phaseTypes,
                 arrayView1d< integer const > const & phaseOrder,
                 ThreePhaseInterpolator const & threePhaseInterpolator,
                 arrayView4d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                 arrayView5d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
                 arrayView3d< real64, relperm::USD_PHASE > const & phaseTrappedVolFrac )
  : RelativePermeabilityBaseUpdate( phaseTypes,
                                    phaseOrder,
                                    phaseRelPerm,
                                    dPhaseRelPerm_dPhaseVolFrac,
                                    phaseTrappedVolFrac ),
  m_relPermKernelWrappers( relPermKernelWrappers ),
  m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
  m_waterOilRelPermMaxValue( waterPhaseMaxVolumeFraction ),
  m_threePhaseInterpolator( threePhaseInterpolator ) {}

TableRelativePermeability::KernelWrapper
TableRelativePermeability::createKernelWrapper()
{

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
  createAllTableKernelWrappers();

  // then we create the actual TableRelativePermeability::KernelWrapper
  return KernelWrapper( m_relPermKernelWrappers,
                        m_phaseMinVolumeFraction,
                        m_waterOilMaxRelPerm,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_threePhaseInterpolator,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac,
                        m_phaseTrappedVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableRelativePermeability, std::string const &, Group * const )

} // namespace constitutive

} // namespace geos
