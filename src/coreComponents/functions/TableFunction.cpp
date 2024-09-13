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
 * @file TableFunction.cpp
 */

#include "TableFunction.hpp"
#include "codingUtilities/Parsing.hpp"
#include "common/DataTypes.hpp"

#include <algorithm>

namespace geos
{

using namespace dataRepository;

TableFunction::TableFunction( const string & name,
                              Group * const parent ):
  FunctionBase( name, parent ),
  m_interpolationMethod( InterpolationType::Linear ),
  m_valueUnit( units::Unknown ),
  m_kernelWrapper( createKernelWrapper() )
{
  registerWrapper( viewKeyStruct::coordinatesString(), &m_tableCoordinates1D ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coordinates inputs for 1D tables" );

  registerWrapper( viewKeyStruct::valuesString(), &m_values ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Values for 1D tables" );

  registerWrapper( viewKeyStruct::coordinateFilesString(), &m_coordinateFiles ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "List of coordinate file names for ND Table" );

  registerWrapper( viewKeyStruct::voxelFileString(), &m_voxelFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Voxel file name for ND Table" );

  registerWrapper( viewKeyStruct::interpolationString(), &m_interpolationMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Interpolation method. Valid options:\n* " + EnumStrings< InterpolationType >::concat( "\n* " ) ).
    setApplyDefaultValue( m_interpolationMethod );
}

void TableFunction::readFile( string const & filename, array1d< real64 > & target )
{
  auto const skipped = []( char const c ){ return std::isspace( c ) || c == ','; };
  try
  {
    parseFile( filename, target, skipped );
  }
  catch( std::runtime_error const & e )
  {
    GEOS_THROW( GEOS_FMT( "{} {}: {}", catalogName(), getDataContext(), e.what() ), InputError );
  }
}

void TableFunction::setInterpolationMethod( InterpolationType const method )
{
  m_interpolationMethod = method;
  reInitializeFunction();
}

void TableFunction::setTableCoordinates( array1d< real64_array > const & coordinates,
                                         std::vector< units::Unit > const & dimUnits )
{
  m_dimUnits = dimUnits;
  m_coordinates.resize( 0 );
  for( localIndex i = 0; i < coordinates.size(); ++i )
  {
    for( localIndex j = 1; j < coordinates[i].size(); ++j )
    {
      GEOS_THROW_IF( coordinates[i][j] - coordinates[i][j-1] <= 0,
                     GEOS_FMT( "{} {}: coordinates must be strictly increasing, but axis {} is not",
                               catalogName(), getDataContext(), i ),
                     InputError );
    }
    m_coordinates.appendArray( coordinates[i].begin(), coordinates[i].end() );
  }
  reInitializeFunction();
}

void TableFunction::setTableValues( real64_array values, units::Unit unit )
{
  m_values = std::move( values );
  m_valueUnit = unit;
  reInitializeFunction();
}

void TableFunction::initializeFunction()
{
  // Read in data
  if( m_coordinates.size() > 0 )
  {
    // This function appears to be already initialized
    // Apparently, this can be called multiple times during unit tests?
  }
  else if( m_coordinateFiles.empty() )
  {
    // 1D Table
    m_coordinates.appendArray( m_tableCoordinates1D.begin(), m_tableCoordinates1D.end() );
    GEOS_THROW_IF_NE_MSG( m_tableCoordinates1D.size(), m_values.size(),
                          GEOS_FMT( "{} {}: 1D table function coordinates and values must have the same length",
                                    catalogName(), getDataContext() ),
                          InputError );
  }
  else
  {
    array1d< real64 > tmp;
    localIndex numValues = 1;
    for( localIndex ii = 0; ii < m_coordinateFiles.size(); ++ii )
    {
      tmp.clear();
      readFile( m_coordinateFiles[ii], tmp );
      m_coordinates.appendArray( tmp.begin(), tmp.end() );
      numValues *= tmp.size();
    }
    // ND Table
    m_values.reserve( numValues );
    readFile( m_voxelFile, m_values );
  }

  reInitializeFunction();
}

void TableFunction::reInitializeFunction()
{
  // Setup index increment (assume data is in Fortran array order)
  localIndex increment = 1;
  for( localIndex ii = 0; ii < m_coordinates.size(); ++ii )
  {
    increment *= m_coordinates.sizeOfArray( ii );
    for( localIndex j = 1; j < m_coordinates[ii].size(); ++j )
    {
      GEOS_THROW_IF( m_coordinates[ii][j] - m_coordinates[ii][j-1] <= 0,
                     GEOS_FMT( "{} {}: coordinates must be strictly increasing, but axis {} is not",
                               catalogName(), getDataContext(), ii ),
                     InputError );
    }
  }
  if( m_coordinates.size() > 0 && !m_values.empty() ) // coordinates and values have been set
  {
    GEOS_THROW_IF_NE_MSG( increment, m_values.size(),
                          GEOS_FMT( "{} {}: number of values does not match total number of table coordinates",
                                    catalogName(), getDataContext() ),
                          InputError );
  }

  // Create the kernel wrapper
  m_kernelWrapper = createKernelWrapper();
}

void TableFunction::checkCoord( real64 const coord, localIndex const dim ) const
{
  GEOS_THROW_IF( dim >= m_coordinates.size() || dim < 0,
                 GEOS_FMT( "{}: The {} dimension ( no. {} ) doesn't exist in the table.",
                           getDataContext(), units::getDescription( getDimUnit( dim ) ), dim ),
                 SimulationError );
  real64 const lowerBound = m_coordinates[dim][0];
  real64 const upperBound = m_coordinates[dim][m_coordinates.sizeOfArray( dim ) - 1];
  GEOS_THROW_IF( coord > upperBound || coord < lowerBound,
                 GEOS_FMT( "{}: Requested {} is out of the table bounds ( lower bound: {} -> upper bound: {} ).",
                           getDataContext(),
                           units::formatValue( coord, getDimUnit( dim ) ),
                           units::formatValue( lowerBound, getDimUnit( dim ) ),
                           units::formatValue( upperBound, getDimUnit( dim ) ) ),
                 SimulationError );
}

void TableFunction::print( std::string const & filename ) const
{
  std::ofstream os( joinPath( FunctionBase::getOutputDirectory(), filename + ".csv" ) );

  integer const numDimensions = LvArray::integerConversion< integer >( m_coordinates.size() );

  if( numDimensions != 2 )
  {
    // print header

    for( integer d = 0; d < numDimensions; d++ )
    {
      os << units::getDescription( getDimUnit( d )) << ",";
    }
    os << units::getDescription( m_valueUnit ) << "\n";

    // print values

    // prepare dividers
    std::vector< integer > div( numDimensions );
    div[0] = 1;
    for( integer d = 1; d < numDimensions; d++ )
    {
      div[d] = div[d-1] * m_coordinates[d-1].size();
    }
    // loop through all the values
    for( integer v = 0; v < m_values.size(); v++ )
    {
      // find coords indices
      std::vector< integer > idx( numDimensions );
      integer r = v;
      for( integer d = numDimensions-1; d >= 0; d-- )
      {
        idx[d] = r / div[d];
        r = r % div[d];
      }
      // finally print out in right order
      for( integer d = 0; d < numDimensions; d++ )
      {
        arraySlice1d< real64 const > const coords = m_coordinates[d];
        os << coords[idx[d]] << ",";
      }
      os << m_values[v] << "\n";
    }
  }
  else // numDimensions == 2
  {
    arraySlice1d< real64 const > const coordsX = m_coordinates[0];
    arraySlice1d< real64 const > const coordsY = m_coordinates[1];
    integer const nX = coordsX.size();
    integer const nY = coordsY.size();
    os<<units::getDescription( getDimUnit( 0 ));
    for( integer j = 0; j < nY; j++ )
    {
      os << "," << units::getDescription( getDimUnit( 1 )) << "=" << coordsY[j];
    }
    os << "\n";
    for( integer i = 0; i < nX; i++ )
    {
      os << coordsX[i];
      for( integer j = 0; j < nY; j++ )
      {
        os << "," << m_values[ j*nX + i ];
      }
      os << "\n";
    }
  }

  os.close();
}

TableFunction::KernelWrapper TableFunction::createKernelWrapper() const
{
  return { m_interpolationMethod,
           m_coordinates.toViewConst(),
           m_values.toViewConst() };
}

real64 TableFunction::evaluate( real64 const * const input ) const
{
  return m_kernelWrapper.compute( input );
}

TableFunction::KernelWrapper::KernelWrapper( InterpolationType const interpolationMethod,
                                             ArrayOfArraysView< real64 const > const & coordinates,
                                             arrayView1d< real64 const > const & values )
  :
  m_interpolationMethod( interpolationMethod ),
  m_coordinates( coordinates ),
  m_values( values )
{}

REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, string const &, Group * const )

} // end of namespace geos
