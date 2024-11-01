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
 * @file PVTFunctionHelpers.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "LvArray/src/sortedArrayManipulation.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

void
BlackOilTables::readTable( string const & fileName,
                           integer minRowLength,
                           array1d< array1d< real64 > > & data )
{
  std::ifstream is( fileName );
  GEOS_ERROR_IF( !is.is_open(),
                 "BlackOilTables: could not open file: " << fileName );

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
    {
      continue;
    }

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
    GEOS_ERROR_IF( data[i].size() < minRowLength,
                   "BlackOilTables: too few entries in row " << i << " of table " << fileName
                                                             << ", minimum " << std::to_string( minRowLength ) << " required" );
  }
}

array1d< array1d< array1d< real64 > > >
BlackOilTables::readAllTables( localIndex const ipWater,
                               array1d< integer > const & phaseTypes,
                               path_array const & tableFileNames ) const
{
  localIndex const numPhases = phaseTypes.size();
  array1d< array1d< array1d< real64 > > > result( numPhases );

  for( localIndex i = 0; i < numPhases; ++i )
  {
    localIndex const rowMinLength = ( phaseTypes[i] == ipWater ) ? 4 : 3;
    readTable( tableFileNames[i], rowMinLength, result[i] );
  }
  return result;
}

void
BlackOilTables::buildAllTables( localIndex const ipOil,
                                localIndex const ipWater,
                                localIndex const ipGas,
                                array1d< integer > const & phaseTypes,
                                path_array const & tableFiles )
{

  // check if both oil and gas are defined
  auto lower = LvArray::sortedArrayManipulation::find( phaseTypes.begin(), phaseTypes.size(), int(ipOil) );
  bool const containsOil = lower != phaseTypes.size();
  GEOS_ERROR_IF( !containsOil, "The oil phase must be defined for all PVT models" );

  // reading data from files
  array1d< array1d< array1d< real64 > > > const phaseTables = readAllTables( ipWater, phaseTypes, tableFiles );

  // finalize table construction
  for( localIndex i = 0; i < phaseTypes.size(); ++i )
  {
    array1d< array1d< real64 > > const & phaseTable = phaseTables[i];

    // oil
    if( phaseTypes[i] == ipOil )
    {
      m_oilTable.resize( phaseTable.size());
      for( localIndex ii = 0; ii < phaseTable.size(); ++ii )
      {
        m_oilTable[ii].resize( phaseTable[ii].size());
        for( localIndex jj = 0; jj < phaseTable[ii].size(); ++jj )
        {
          m_oilTable[ii][jj] = phaseTable[ii][jj];
        }
      }
    }
    // gas
    else if( phaseTypes[i] == ipGas )
    {
      m_gasTable.resize( phaseTable.size());
      for( localIndex ii = 0; ii < phaseTable.size(); ++ii )
      {
        m_gasTable[ii].resize( phaseTable[ii].size());
        for( localIndex jj = 0; jj < phaseTable[ii].size(); ++jj )
        {
          m_gasTable[ii][jj] = phaseTable[ii][jj];
        }
      }
    }
    // water
    else if( phaseTypes[i] == ipWater )
    {
      m_waterTable.resize( phaseTable[0].size());
      for( localIndex ii = 0; ii < phaseTable[0].size(); ++ii )
      {
        m_waterTable[ii] = phaseTable[0][ii];
      }
    }
    else
    {
      GEOS_ERROR( "Phase type not supported for Black Oil model" );
    }
  }
}

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
