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
 * @file PVTDriver.cpp
 */

#include "PVTDriver.hpp"

#include "common/MpiWrapper.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#include "codingUtilities/StringUtilities.hpp"

#include <fstream>

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

PVTDriver::PVTDriver( const string & name,
                      Group * const parent ):
  TaskBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::fluidNameString(), &m_fluidName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fluid to test" );

  registerWrapper( viewKeyStruct::feedString(), &m_feed ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Feed composition array [mol fraction]" );

  registerWrapper( viewKeyStruct::pressureFunctionString(), &m_pressureFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling pressure time history" );

  registerWrapper( viewKeyStruct::temperatureFunctionString(), &m_temperatureFunctionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Function controlling temperature time history" );

  registerWrapper( viewKeyStruct::outputMassDensityString(), &m_outputMassDensity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate that the mass density of each phase should be output" );

  registerWrapper( viewKeyStruct::outputCompressibilityString(), &m_outputCompressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate that the total compressibility should be output" );

  registerWrapper( viewKeyStruct::outputPhaseCompositionString(), &m_outputPhaseComposition ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate that phase compositions should be output" );

  //todo refactor in mother class
  registerWrapper( viewKeyStruct::numStepsString(), &m_numSteps ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of load steps to take" );

  registerWrapper( viewKeyStruct::outputString(), &m_outputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Output file" );

  registerWrapper( viewKeyStruct::baselineString(), &m_baselineFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "none" ).
    setDescription( "Baseline file" );
}

void PVTDriver::postInputInitialization()
{
  // Validate some inputs
  GEOS_ERROR_IF( m_outputMassDensity != 0 && m_outputMassDensity != 1,
                 getWrapperDataContext( viewKeyStruct::outputMassDensityString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );

  GEOS_ERROR_IF( m_outputCompressibility != 0 && m_outputCompressibility != 1,
                 getWrapperDataContext( viewKeyStruct::outputCompressibilityString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );

  GEOS_ERROR_IF( m_outputPhaseComposition != 0 && m_outputPhaseComposition != 1,
                 getWrapperDataContext( viewKeyStruct::outputPhaseCompositionString() ) <<
                 ": option can be either 0 (false) or 1 (true)" );

  // get number of phases and components

  constitutive::MultiFluidBase & baseFluid = getFluid();

  m_numPhases = baseFluid.numFluidPhases();
  m_numComponents = baseFluid.numFluidComponents();

  // Number of rows in numSteps+1
  integer const numRows = m_numSteps+1;

  // Number of columns depends on options
  // Default column order = time, pressure, temp, totalDensity, phaseFraction_{1:NP}, phaseDensity_{1:NP}, phaseViscosity_{1:NP}
  integer numCols = 3*m_numPhases+4;

  // If the mass density is requested then add NP columns
  if( m_outputMassDensity != 0 )
  {
    numCols += m_numPhases;
  }

  // If the total compressibility is requested then add a column
  if( m_outputCompressibility != 0 )
  {
    numCols++;
  }

  // If phase compositions are required we add {1:NP*NC} phase compositions
  if( m_outputPhaseComposition != 0 )
  {
    numCols += m_numPhases * m_numComponents;
  }

  // resize data table to fit number of timesteps and fluid phases:
  m_table.resize( numRows, numCols );

  // initialize functions

  FunctionManager & functionManager = FunctionManager::getInstance();

  TableFunction & pressureFunction = functionManager.getGroup< TableFunction >( m_pressureFunctionName );
  TableFunction & temperatureFunction = functionManager.getGroup< TableFunction >( m_temperatureFunctionName );

  pressureFunction.initializeFunction();
  temperatureFunction.initializeFunction();

  // determine time increment

  ArrayOfArraysView< real64 > coordinates = pressureFunction.getCoordinates();
  real64 const minTime = coordinates[0][0];
  real64 const maxTime = coordinates[0][coordinates.sizeOfArray( 0 )-1];
  real64 const dt = (maxTime-minTime) / m_numSteps;

  // set input columns

  for( integer n=0; n<m_numSteps+1; ++n )
  {
    m_table( n, TIME ) = minTime + n*dt;
    m_table( n, PRES ) = pressureFunction.evaluate( &m_table( n, TIME ) );
    m_table( n, TEMP ) = temperatureFunction.evaluate( &m_table( n, TIME ) );
  }
}


bool PVTDriver::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                         real64 const GEOS_UNUSED_PARAM( dt ),
                         integer const GEOS_UNUSED_PARAM( cycleNumber ),
                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                         real64 const GEOS_UNUSED_PARAM( eventProgress ),
                         DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // this code only makes sense in serial

  GEOS_THROW_IF( MpiWrapper::commRank() > 0, "PVTDriver should only be run in serial", std::runtime_error );

  // get the fluid out of the constitutive manager.
  // for the moment it is of type MultiFluidBase.

  constitutive::MultiFluidBase & baseFluid = getFluid();

  // depending on logLevel, print some useful info

  if( getLogLevel() > 0 )
  {
    GEOS_LOG_RANK_0( "Launching PVT Driver" );
    GEOS_LOG_RANK_0( "  Fluid .................. " << m_fluidName );
    GEOS_LOG_RANK_0( "  Type ................... " << baseFluid.getCatalogName() );
    GEOS_LOG_RANK_0( "  No. of Phases .......... " << m_numPhases );
    GEOS_LOG_RANK_0( "  No. of Components ...... " << m_numComponents );
    GEOS_LOG_RANK_0( "  Pressure Control ....... " << m_pressureFunctionName );
    GEOS_LOG_RANK_0( "  Temperature Control .... " << m_temperatureFunctionName );
    GEOS_LOG_RANK_0( "  Steps .................. " << m_numSteps );
    GEOS_LOG_RANK_0( "  Output ................. " << m_outputFile );
    GEOS_LOG_RANK_0( "  Baseline ............... " << m_baselineFile );
    GEOS_LOG_RANK_0( "  Output Mass Density .... " << m_outputMassDensity );
    GEOS_LOG_RANK_0( "  Output Compressibility . " << m_outputCompressibility );
    GEOS_LOG_RANK_0( "  Output Phase Comp. ..... " << m_outputPhaseComposition );
  }

  // create a dummy discretization with one quadrature point for
  // storing constitutive data

  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  dataRepository::Group discretization( "discretization", &rootGroup );

  discretization.resize( 1 );   // one element
  baseFluid.allocateConstitutiveData( discretization, 1 );   // one quadrature point

  // pass the fluid through the ConstitutivePassThru to downcast from the
  // base type to a known model type.  the lambda here then executes the
  // appropriate test driver. note that these calls will move data to device if available.

  constitutiveUpdatePassThru( baseFluid, [&] ( auto & selectedFluid )
  {
    using FLUID_TYPE = TYPEOFREF( selectedFluid );
    runTest< FLUID_TYPE >( selectedFluid, m_table );
  } );

  // move table back to host for output
  m_table.move( hostMemorySpace );

  if( m_outputFile != "none" )
  {
    outputResults();
  }

  if( m_baselineFile != "none" )
  {
    compareWithBaseline();
  }

  return false;
}

void PVTDriver::outputResults()
{
  // TODO: improve file path output to grab command line -o directory
  //       for the moment, we just use the specified m_outputFile directly

  FILE * fp = fopen( m_outputFile.c_str(), "w" );

  integer columnIndex = 0;
  fprintf( fp, "# column %d = time\n", ++columnIndex );
  fprintf( fp, "# column %d = pressure\n", ++columnIndex );
  fprintf( fp, "# column %d = temperature\n", ++columnIndex );
  fprintf( fp, "# column %d = density\n", ++columnIndex );
  if( m_outputCompressibility != 0 )
  {
    fprintf( fp, "# column %d = total compressibility\n", ++columnIndex );
  }

  auto const phaseNames = getFluid().phaseNames();

  fprintf( fp, "# columns %d-%d = phase fractions\n", columnIndex+1, columnIndex + m_numPhases );
  columnIndex += m_numPhases;
  fprintf( fp, "# columns %d-%d = phase densities\n", columnIndex+1, columnIndex + m_numPhases );
  columnIndex += m_numPhases;
  if( m_outputMassDensity != 0 )
  {
    fprintf( fp, "# columns %d-%d = phase mass densities\n", columnIndex+1, columnIndex + m_numPhases );
    columnIndex += m_numPhases;
  }
  fprintf( fp, "# columns %d-%d = phase viscosities\n", columnIndex+1, columnIndex + m_numPhases );
  columnIndex += m_numPhases;

  if( m_outputPhaseComposition != 0 )
  {
    string const componentNames = stringutilities::join( getFluid().componentNames(), ", " );
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      fprintf( fp, "# columns %d-%d = %s phase fractions [%s]\n", columnIndex+1, columnIndex + m_numComponents,
               phaseNames[ip].c_str(), componentNames.c_str() );
      columnIndex += m_numComponents;
    }
  }

  for( integer n=0; n<m_table.size( 0 ); ++n )
  {
    for( integer col=0; col<m_table.size( 1 ); ++col )
    {
      fprintf( fp, "%.4e ", m_table( n, col ) );
    }
    fprintf( fp, "\n" );
  }
  fclose( fp );
}


void PVTDriver::compareWithBaseline()
{
  // open baseline file

  std::ifstream file( m_baselineFile.c_str() );
  GEOS_THROW_IF( !file.is_open(), "Can't seem to open the baseline file " << m_baselineFile, InputError );

  // discard file header

  integer headerRows = 7;
  if( m_outputCompressibility )
  {
    headerRows++;
  }
  if( m_outputMassDensity )
  {
    headerRows++;
  }
  if( m_outputPhaseComposition )
  {
    headerRows += getFluid().numFluidPhases();
  }

  string line;
  for( integer row=0; row < headerRows; ++row )
  {
    getline( file, line );
  }

  // read data block.  we assume the file size is consistent with m_table,
  // but check for a premature end-of-file. we then compare results value by value.
  // we ignore the newton iteration and residual columns, as those may be platform
  // specific.

  real64 value;
  for( integer row=0; row < m_table.size( 0 ); ++row )
  {
    for( integer col=0; col < m_table.size( 1 ); ++col )
    {
      GEOS_THROW_IF( file.eof(), "Baseline file appears shorter than internal results", std::runtime_error );
      file >> value;

      real64 const error = fabs( m_table[row][col]-value ) / ( fabs( value )+1 );
      GEOS_THROW_IF( error > MultiFluidConstants::baselineTolerance,
                     GEOS_FMT( "Results do not match baseline ({} vs {}) at data row {} (row {} with header) and column {}",
                               m_table[row][col], value, row+1, row+headerRows, col+1 ),
                     std::runtime_error );
    }
  }

  // check we actually reached the end of the baseline file

  file >> value;
  GEOS_THROW_IF( !file.eof(), "Baseline file appears longer than internal results", std::runtime_error );

  // success

  if( getLogLevel() > 0 )
  {
    GEOS_LOG_RANK_0( "  Comparison ............. Internal results consistent with baseline." );
  }

  file.close();
}

constitutive::MultiFluidBase &
PVTDriver::getFluid()
{
  ConstitutiveManager & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  MultiFluidBase & baseFluid = constitutiveManager.getGroup< MultiFluidBase >( m_fluidName );
  return baseFluid;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        PVTDriver,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
