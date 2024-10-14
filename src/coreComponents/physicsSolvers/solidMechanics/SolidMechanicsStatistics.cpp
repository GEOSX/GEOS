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
 * @file SolidMechanicsStatistics.cpp
 */

#include "SolidMechanicsStatistics.hpp"

#include "common/MpiWrapper.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsStatistics::SolidMechanicsStatistics( const string & name,
                                                    Group * const parent ):
  Base( name, parent )
{}

void SolidMechanicsStatistics::registerDataOnMesh( Group & meshBodies )
{
  // the fields have to be registered in "registerDataOnMesh" (and not later)
  // otherwise they cannot be targeted by TimeHistory

  // for now, this guard is needed to avoid breaking the xml schema generation
  if( m_solver == nullptr )
  {
    return;
  }

  m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    nodeManager.registerWrapper< NodeStatistics >( viewKeyStruct::nodeStatisticsString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
    nodeManager.excludeWrappersFromPacking( { viewKeyStruct::nodeStatisticsString() } );
    NodeStatistics & nodeStatistics = nodeManager.getReference< NodeStatistics >( viewKeyStruct::nodeStatisticsString() );

    nodeStatistics.minDisplacement.resizeDimension< 0 >( 3 );
    nodeStatistics.maxDisplacement.resizeDimension< 0 >( 3 );

    // write output header
    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + mesh.getName() + "_node_statistics" + ".csv" );
      outputFile << "Time [s],Min displacement X [m],Min displacement Y [m],Min displacement Z [m],"
                 << "Max displacement X [m],Max displacement Y [m],Max displacement Z [m]" << std::endl;
      outputFile.close();
    }
  } );
}

bool SolidMechanicsStatistics::execute( real64 const time_n,
                                        real64 const dt,
                                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                        DomainPartition & domain )
{
  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                          MeshLevel & mesh,
                                                                          arrayView1d< string const > const & )
  {
    // current time is time_n + dt
    computeNodeStatistics( mesh, time_n + dt );
  } );
  return false;
}

void SolidMechanicsStatistics::computeNodeStatistics( MeshLevel & mesh, real64 const time ) const
{
  GEOS_MARK_FUNCTION;

  // Step 1: increment the min/max quantities

  NodeManager & nodeManager = mesh.getNodeManager();
  arrayView1d< integer const > const ghostRank = nodeManager.ghostRank();
  solidMechanics::arrayViewConst2dLayoutTotalDisplacement const & u =
    nodeManager.getField< solidMechanics::totalDisplacement >();

  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxDispX( -LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxDispY( -LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxDispZ( -LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minDispX( LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minDispY( LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minDispZ( LvArray::NumericLimits< real64 >::max );

  forAll< parallelDevicePolicy<> >( nodeManager.size(), [u,
                                                         ghostRank,
                                                         maxDispX,
                                                         maxDispY,
                                                         maxDispZ,
                                                         minDispX,
                                                         minDispY,
                                                         minDispZ]
                                    GEOS_HOST_DEVICE ( localIndex const a )
  {
    if( ghostRank[a] < 0 )
    {
      maxDispX.max( u[a][0] );
      maxDispY.max( u[a][1] );
      maxDispZ.max( u[a][2] );
      minDispX.min( u[a][0] );
      minDispY.min( u[a][1] );
      minDispZ.min( u[a][2] );
    }
  } );

  // Step 2: synchronize the results over the MPI ranks

  NodeStatistics & nodeStatistics = nodeManager.getReference< NodeStatistics >( viewKeyStruct::nodeStatisticsString() );

  nodeStatistics.maxDisplacement[0] = maxDispX.get();
  nodeStatistics.maxDisplacement[1] = maxDispY.get();
  nodeStatistics.maxDisplacement[2] = maxDispZ.get();
  nodeStatistics.minDisplacement[0] = minDispX.get();
  nodeStatistics.minDisplacement[1] = minDispY.get();
  nodeStatistics.minDisplacement[2] = minDispZ.get();

  MpiWrapper::allReduce( nodeStatistics.maxDisplacement.data(),
                         nodeStatistics.maxDisplacement.data(),
                         3,
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Max ),
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( nodeStatistics.minDisplacement.data(),
                         nodeStatistics.minDisplacement.data(),
                         3,
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Min ),
                         MPI_COMM_GEOS );

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{} (time {} s): Min displacement (X, Y, Z): {}, {}, {} m",
                                      getName(), time, nodeStatistics.minDisplacement[0],
                                      nodeStatistics.minDisplacement[1], nodeStatistics.minDisplacement[2] ) );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{} (time {} s): Max displacement (X, Y, Z): {}, {}, {} m",
                                      getName(), time, nodeStatistics.maxDisplacement[0],
                                      nodeStatistics.maxDisplacement[1], nodeStatistics.maxDisplacement[2] ) );

  if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
  {
    std::ofstream outputFile( m_outputDir + "/" + mesh.getName() + "_node_statistics" + ".csv", std::ios_base::app );
    outputFile << time;
    for( integer i = 0; i < 3; ++i )
      outputFile << "," << nodeStatistics.minDisplacement[i];
    for( integer i = 0; i < 3; ++i )
      outputFile << "," << nodeStatistics.maxDisplacement[i];
    outputFile << std::endl;
    outputFile.close();
  }
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SolidMechanicsStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
