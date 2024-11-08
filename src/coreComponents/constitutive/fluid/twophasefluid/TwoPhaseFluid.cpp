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
 * @file TwoPhaseFluid.cpp
 */

#include "TwoPhaseFluid.hpp"
#include "TwoPhaseFluidFields.hpp"

#include "functions/FunctionManager.hpp"


namespace geos
{

using namespace dataRepository;

namespace constitutive
{


TwoPhaseFluid::TwoPhaseFluid( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of fluid phases" );

  // 1) First option: specify PVT tables from one file per phase, read the files line by line, and populate the internal TableFunctions
  registerWrapper( viewKeyStruct::tableFilesString(), &m_tableFiles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of filenames with input PVT tables (one per phase)" );

  // 2) Second option: specify TableFunction names for each phase,
  registerWrapper( viewKeyStruct::densityTableNamesString(), &m_densityTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of density TableFuncion names from the Function block. \n"
                    "The user must provide one TableFunction per phase, respecting the order provided in \"phaseNames\"." );

  registerWrapper( viewKeyStruct::viscosityTableNamesString(), &m_viscosityTableNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of viscosity TableFuncion names from the Function block. \n"
                    "The user must provide one TableFunction per phase, respecting the order provided in \"phaseNames\"." );

  registerField( fields::twophasefluid::phaseDensity{}, &m_phaseDensity.value );
  registerField( fields::twophasefluid::dPhaseDensity{}, &m_phaseDensity.derivs );
  registerField( fields::twophasefluid::phaseDensity_n{}, &m_phaseDensity_n );

  registerField( fields::twophasefluid::phaseViscosity{}, &m_phaseViscosity.value );
  registerField( fields::twophasefluid::dPhaseViscosity{}, &m_phaseViscosity.derivs );
}


std::unique_ptr< ConstitutiveBase >
TwoPhaseFluid::deliverClone( string const & name, Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}


void TwoPhaseFluid::resizeFields( localIndex const size, localIndex const numPts )
{
  // Assume sole dependency on pressure, i.e. one derivative
  m_phaseDensity.value.resize( size, numPts, 2 );
  m_phaseDensity.derivs.resize( size, numPts, 2, 1 );

  m_phaseDensity_n.resize( size, numPts, 2 );

  m_phaseViscosity.value.resize( size, numPts, 2 );
  m_phaseViscosity.derivs.resize( size, numPts, 2, 1 );
}


void TwoPhaseFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}


void TwoPhaseFluid::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  // Input relationships can be provided either as text files or TableFunctions.
  m_tableFiles.empty() ? readInputDataFromTableFunctions() : readInputDataFromFileTableFunctions();

  checkTableConsistency();
  /* TODO
     // call to correctly set member array tertiary sizes on the 'main' material object
     resizeFields( 0, 0 );

     // set labels on array wrappers for plottable fields
     setLabels();
   */
}


void
TwoPhaseFluid::readTable( string const & fileName,
                          integer minRowLength,
                          array1d< array1d< real64 > > & data )
{
  std::ifstream is( fileName );
  GEOS_THROW_IF( !is.is_open(),
                 GEOS_FMT( "{}: could not open file: {}", getFullName(), fileName ),
                 InputError );

  // Read line-by-line until eof
  string str;
  while( std::getline( is, str ) )
  {
    // Remove whitespace and end-of-line characters, if any
    str = stringutilities::trim( str, " \r" );

    // Remove # and -- (Eclipse-style) comments
    str = stringutilities::removeStringAndFollowingContent( str, "#" );
    str = stringutilities::removeStringAndFollowingContent( str, "--" );

    // Skip empty or comment-only strings
    if( str.empty() )
      continue;

    // Add and read a new line entry
    array1d< real64 > newLine = stringutilities::fromStringToArray< real64 >( str );
    if( !newLine.empty() )
    {
      data.emplace_back( std::move( newLine ) );
    }
  }

  is.close();

  for( localIndex i = 0; i < data.size(); ++i )
  {
    GEOS_THROW_IF( data[i].size() < minRowLength,
                   GEOS_FMT( "{}: too few entries in row {} of table {}, minimum {} required", getFullName(), i, fileName, minRowLength ),
                   InputError );
  }
}


void TwoPhaseFluid::fillData( integer const ip,
                              array1d< array1d< real64 > > const & tableValues )
{
  array1d< array1d< real64 > > pressureCoords( 1 );
  pressureCoords[0].resize( tableValues.size() );
  array1d< real64 > density( tableValues.size() );
  array1d< real64 > viscosity( tableValues.size() );

  for( localIndex i = 0; i < tableValues.size(); ++i )
  {
    GEOS_THROW_IF_NE_MSG( tableValues[i].size(), 3,
                          GEOS_FMT( "{}: three columns (pressure, density, and viscosity) are expected", getFullName() ),
                          InputError );

    pressureCoords[0][i] = tableValues[i][0];
    density[i] = tableValues[i][1];
    viscosity[i] = tableValues[i][2];
  }

  string const densityTableName = getName() + "DensityPhase" + GEOS_FMT( "{}", ip );
  string const viscosityTableName = getName() +  "ViscosityPhase" + GEOS_FMT( "{}", ip );
  m_densityTableNames.emplace_back( densityTableName );
  m_viscosityTableNames.emplace_back( viscosityTableName );

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & tableDensity =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", densityTableName ) );
  tableDensity.setTableCoordinates( pressureCoords, { units::Pressure } );
  tableDensity.setTableValues( density, units::Density );
  tableDensity.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tableViscosity =
    dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", viscosityTableName ) );
  tableViscosity.setTableCoordinates( pressureCoords, { units::Pressure } );
  tableViscosity.setTableValues( viscosity, units::Viscosity );
  tableViscosity.setInterpolationMethod( TableFunction::InterpolationType::Linear );
}


void TwoPhaseFluid::readInputDataFromFileTableFunctions()
{
  // Check for ambiguous definition
  GEOS_THROW_IF( !(m_densityTableNames.empty() && m_viscosityTableNames.empty()),
                 GEOS_FMT( "{}: input is redundant (both TableFunction names and text files)", getFullName() ),
                 InputError );


  // Check that we have exactly two table files (one per phase)
  GEOS_THROW_IF_NE_MSG( m_tableFiles.size(), 2,
                        GEOS_FMT( "{}: expecting two table files (one per phase)", getFullName() ),
                        InputError );

  array1d< array1d< real64 > > tableValues;
  for( integer ip = 0; ip < 2; ++ip )
  {
    tableValues.clear();
    readTable( m_tableFiles[ip], 3, tableValues );
    fillData( ip, tableValues );
  }
}


void TwoPhaseFluid::readInputDataFromTableFunctions()
{
  // Check for ambiguous definition
  GEOS_THROW_IF( !m_tableFiles.empty(),
                 GEOS_FMT( "{}: input is redundant (both TableFunction names and text files)", getFullName() ),
                 InputError );

  // Since we are considering a two phase fluid, we should have exactly 2 tables per property
  GEOS_THROW_IF_NE_MSG( m_densityTableNames.size(), 2,
                        GEOS_FMT( "{}: one density table must be provided for each phase", getFullName() ),
                        InputError );

  GEOS_THROW_IF_NE_MSG( m_viscosityTableNames.size(), 2,
                        GEOS_FMT( "{}: one viscosity table must be provided for each phase", getFullName() ),
                        InputError );


  FunctionManager const & functionManager = FunctionManager::getInstance();

  for( integer iph = 0; iph < 2; ++iph )
  {
    GEOS_THROW_IF( !functionManager.hasGroup( m_densityTableNames[iph] ),
                   GEOS_FMT( "{}: density table '{}' not found", getFullName(), m_densityTableNames[iph] ),
                   InputError );

    GEOS_THROW_IF( !functionManager.hasGroup( m_viscosityTableNames[iph] ),
                   GEOS_FMT( "{}: viscosity table '{}' not found", getFullName(), m_viscosityTableNames[iph] ),
                   InputError );
  }
}


void TwoPhaseFluid::initializePostSubGroups()
{
  ConstitutiveBase::initializePostSubGroups();

  FunctionManager const & functionManager = FunctionManager::getInstance();
  for( integer iph = 0; iph < 2; ++iph )
  {
    // Grab the tables by name from the function manager,
    // then add them in a list to create their table wrappers when needed
    TableFunction const & densityTable = functionManager.getGroup< TableFunction const >( m_densityTableNames[iph] );
    m_densityTables.emplace_back( &densityTable );
    m_densityTableKernels.emplace_back( m_densityTables[iph]->createKernelWrapper() );

    TableFunction const & viscosityTable = functionManager.getGroup< TableFunction const >( m_viscosityTableNames[iph] );
    m_viscosityTables.emplace_back( &viscosityTable );
    m_viscosityTableKernels.emplace_back( m_viscosityTables[iph]->createKernelWrapper() );
  }
}


void TwoPhaseFluid::checkTableConsistency() const
{
  FunctionManager const & functionManager = FunctionManager::getInstance();
  for( integer iph = 0; iph < 2; ++iph )
  {
    TableFunction const & densityTable = functionManager.getGroup< TableFunction const >( m_densityTableNames[iph] );
    arrayView1d< real64 const > const density = densityTable.getValues();

    for( localIndex i = 1; i < density.size(); ++i )
    {
      GEOS_THROW_IF( density[i] - density[i-1] < 0,
                     GEOS_FMT( "{}: in table '{}' density values must be increasing", getFullName(), densityTable.getName() ),
                     InputError );
    }
  }
}


TwoPhaseFluid::KernelWrapper
TwoPhaseFluid::createKernelWrapper()
{
  return KernelWrapper( m_densityTableKernels,
                        m_viscosityTableKernels,
                        m_phaseDensity.toView(),
                        m_phaseViscosity.toView());
}


TwoPhaseFluid::KernelWrapper::KernelWrapper(
  arrayView1d< TableFunction::KernelWrapper const > densityTables,
  arrayView1d< TableFunction::KernelWrapper const > viscosityTables,
  PhaseProp::ViewType phaseDensity,
  PhaseProp::ViewType phaseViscosity )
  : m_densityTables( std::move( densityTables )),
  m_viscosityTables( std::move( viscosityTables )),
  m_phaseDensity( std::move( phaseDensity )),
  m_phaseViscosity( std::move( phaseViscosity )) {}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, TwoPhaseFluid, string const &, Group * const )

}  // namespace constitutive
}  // namespace geos
