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
 * @file MultivariableTableFunction.cpp
 */

#include "MultivariableTableFunction.hpp"

#include "common/DataTypes.hpp"
#include <algorithm>

namespace geos
{

using namespace dataRepository;

MultivariableTableFunction::MultivariableTableFunction( const string & name,
                                                        Group * const parent ):
  FunctionBase( name, parent )
{}

void MultivariableTableFunction::initializeFunctionFromFile( string const & filename )
{
  std::ifstream file( filename.c_str() );
  GEOS_THROW_IF( !file, catalogName() << " " << getDataContext() << ": could not read input file " << filename, InputError );

  integer numDims, numOps;
  globalIndex numPointsTotal = 1;
  real64_array axisMinimums, axisMaximums;
  integer_array axisPoints;

  // 1. Read numDims and numOps


  file >> numDims;
  GEOS_THROW_IF( !file, "Can`t read number of table dimensions", InputError );
  file >> numOps;
  GEOS_THROW_IF( !file, "Can`t read number of interpolatored operators", InputError );

  // assume no more than 10 dimensions
  GEOS_THROW_IF_LT_MSG( numDims, 1, catalogName() << " " << getDataContext() << ": positive integer value expected", InputError );
  GEOS_THROW_IF_GT_MSG( numDims, 10, catalogName() << " " << getDataContext() << ": maximum 10 dimensions expected", InputError );

  // assume no more than 100 operators
  GEOS_THROW_IF_LT_MSG( numOps, 1, catalogName() << " " << getDataContext() << ": positive integer value expected", InputError );
  GEOS_THROW_IF_GT_MSG( numOps, 100, catalogName() << " " << getDataContext() << ": maximum 100 operators expected", InputError );

  axisMinimums.resize( numDims );
  axisMaximums.resize( numDims );
  axisPoints.resize( numDims );

  // 2. Read axis parameters

  for( integer i = 0; i < numDims; i++ )
  {
    file >> axisPoints[i];
    GEOS_THROW_IF( !file, catalogName() << " " << getDataContext() << ": can`t read the number of points for axis " + std::to_string( i ), InputError );
    GEOS_THROW_IF_LE_MSG( axisPoints[i], 1, catalogName() << " " << getDataContext() << ": minimum 2 discretization point per axis are expected", InputError );
    file >> axisMinimums[i];
    GEOS_THROW_IF( !file, catalogName() << " " << getDataContext() << ": can`t read minimum value for axis " + std::to_string( i ), InputError );
    file >> axisMaximums[i];
    GEOS_THROW_IF( !file, catalogName() << " " << getDataContext() << ": can`t read maximum value for axis " + std::to_string( i ), InputError );
    GEOS_THROW_IF_LT_MSG( axisMaximums[i], axisMinimums[i], catalogName() << " " << getDataContext() << ": maximum axis value is expected to be larger than minimum", InputError );

    numPointsTotal *= axisPoints[i];
  }

  // lets limit the point storage size with 1 Gb (taking into account that hypercube storage is 2^numDim larger)
  real64 pointStorageMemoryLimitGB = 1;

  GEOS_THROW_IF_GT_MSG( numPointsTotal * numOps, pointStorageMemoryLimitGB * 1024 * 1024 * 1024 / 8,
                        catalogName() << " " << getDataContext() <<
                        ": point storage size exceeds " + std::to_string( pointStorageMemoryLimitGB ) +
                        " Gb, please reduce number of points",
                        InputError );

  m_pointData.resize( numPointsTotal * numOps );

  // 3. Read table values
  for( auto i = 0; i < numPointsTotal; i++ )
  {
    for( auto j = 0; j < numOps; j++ )
    {
      file >> m_pointData[i * numOps + j];
      GEOS_THROW_IF( !file, catalogName() << " " << getDataContext() << ": table file is shorter than expected", InputError );
    }
  }
  real64 value;

  file >> value;
  GEOS_THROW_IF( file, catalogName() << " " << getDataContext() << ": table file is longer than expected", InputError );

  file.close();

  setTableCoordinates( numDims, numOps, axisMinimums, axisMaximums, axisPoints );
  initializeFunction();
}


void MultivariableTableFunction::setTableCoordinates( integer const numDims,
                                                      integer const numOps,
                                                      real64_array const & axisMinimums,
                                                      real64_array const & axisMaximums,
                                                      integer_array const & axisPoints )
{
  m_numDims = numDims;
  m_numOps = numOps;
  m_numVerts =  1 << numDims;

  m_axisMinimums = axisMinimums;
  m_axisMaximums = axisMaximums;
  m_axisPoints = axisPoints;
}

void MultivariableTableFunction::setTableValues( real64_array values )
{
  m_pointData = std::move( values );
}
void MultivariableTableFunction::getHypercubePoints( globalIndex const hypercubeIndex, globalIndex_array & hypercubePoints ) const
{
  auto remainder = hypercubeIndex;
  auto pwr = m_numVerts;

  for( auto j = 0; j < m_numVerts; ++j )
  {
    hypercubePoints[j] = 0;
  }

  for( auto i = 0; i < m_numDims; ++i )
  {

    integer const axis_idx = remainder / m_axisHypercubeMults[i];
    remainder = remainder % m_axisHypercubeMults[i];

    pwr /= 2;

    for( auto j = 0; j < m_numVerts; ++j )
    {
      auto zero_or_one = (j / pwr) % 2;
      hypercubePoints[j] += (axis_idx + zero_or_one) * m_axisPointMults[i];
    }
  }
}


void MultivariableTableFunction::initializeFunction()
{
  // check input


  GEOS_THROW_IF_NE_MSG( m_numDims, m_axisMinimums.size(), catalogName() << " " << getDataContext() <<
                        ": single minimum value is expected for each of " + std::to_string( m_numDims ) + "dimensions",
                        InputError );
  GEOS_THROW_IF_NE_MSG( m_numDims, m_axisMaximums.size(), catalogName() << " " << getDataContext() <<
                        ": single maxumum value is expected for each of " + std::to_string( m_numDims ) + "dimensions",
                        InputError );
  GEOS_THROW_IF_NE_MSG( m_numDims, m_axisPoints.size(), catalogName() << " " << getDataContext() <<
                        "single number is expected for each of " + std::to_string( m_numDims ) + "dimensions",
                        InputError );

  m_axisSteps.resize( m_numDims );
  m_axisStepInvs.resize( m_numDims );
  m_axisPointMults.resize( m_numDims );
  m_axisHypercubeMults.resize( m_numDims );

  // compute service data

  for( integer dim = 0; dim < m_numDims; dim++ )
  {
    m_axisSteps[dim] = (m_axisMaximums[dim] - m_axisMinimums[dim]) / (m_axisPoints[dim] - 1);
    m_axisStepInvs[dim] = 1 / m_axisSteps[dim];
  }

  m_axisPointMults[m_numDims - 1] = 1;
  m_axisHypercubeMults[m_numDims - 1] = 1;
  for( integer dim = m_numDims - 2; dim >= 0; --dim )
  {
    m_axisPointMults[dim] = m_axisPointMults[dim + 1] * m_axisPoints[dim + 1];
    m_axisHypercubeMults[dim] = m_axisHypercubeMults[dim + 1] * (m_axisPoints[dim + 1] - 1);
  }

  // check for point index overflow
  // fp type is intentional - to prevent overflow during computation and detect it later
  real64 numTablePoints = 1.0;
  globalIndex numTableHypercubes = 1;
  for( int dim = 0; dim < m_numDims; dim++ )
  {
    numTablePoints *= m_axisPoints[dim];
    numTableHypercubes *= m_axisPoints[dim] - 1;
  }


  // check is point data size is correct
  GEOS_THROW_IF_NE_MSG( globalIndex( numTablePoints ) * m_numOps, m_pointData.size(), catalogName() << " " << getDataContext() <<
                        ": table values array is expected to have length of " + std::to_string( globalIndex( numTablePoints ) * m_numOps ), InputError );

  // lets limit the hypercube storage size with 16 Gb
  real64 hypercubeStorageMemoryLimitGB = 16;

  GEOS_THROW_IF_GT_MSG( numTableHypercubes * m_numVerts * m_numOps, hypercubeStorageMemoryLimitGB * 1024 * 1024 * 1024 / 8,
                        catalogName() << " " << getDataContext() <<
                        ": hypercube storage size exceeds " + std::to_string( hypercubeStorageMemoryLimitGB ) +
                        " Gb, please reduce number of points",
                        InputError );

  // initialize hypercube data storage
  m_hypercubeData.resize( numTableHypercubes * m_numVerts * m_numOps );
  globalIndex_array points( m_numVerts );

  // fill each hypercube directly on the device with corresponding data from m_pointData
  for( auto i = 0; i < numTableHypercubes; i++ )
  {
    getHypercubePoints( i, points );

    for( auto j = 0; j < m_numVerts; ++j )
    {
      std::copy( m_pointData.begin() + points[j] * m_numOps,
                 m_pointData.begin() + (points[j] + 1) * m_numOps,
                 m_hypercubeData.begin() + m_numOps * (i * m_numVerts + j));
    }
  }

}

REGISTER_CATALOG_ENTRY( FunctionBase, MultivariableTableFunction, string const &, Group * const )

} // end of namespace geos
