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
 * @file SinglePhaseStatistics.cpp
 */

#include "SinglePhaseStatistics.hpp"

#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/StatisticsKernel.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

SinglePhaseStatistics::SinglePhaseStatistics( const string & name,
                                              Group * const parent ):
  Base( name, parent )
{
  addLogLevel< logInfo::Statistics >();
}

void SinglePhaseStatistics::registerDataOnMesh( Group & meshBodies )
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
                                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    for( integer i = 0; i < regionNames.size(); ++i )
    {
      ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
      region.registerWrapper< RegionStatistics >( viewKeyStruct::regionStatisticsString() ).
        setRestartFlags( RestartFlags::NO_WRITE );
      region.excludeWrappersFromPacking( { viewKeyStruct::regionStatisticsString() } );

      // write output header
      if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
      {
        std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv" );
        outputFile <<
          "Time [s],Min pressure [Pa],Average pressure [Pa],Max pressure [Pa],Min delta pressure [Pa],Max delta pressure [Pa]," <<
          "Min temperature [Pa],Average temperature [Pa],Max temperature [Pa],Total dynamic pore volume [rm^3],Total fluid mass [kg]";
        outputFile << std::endl;
        outputFile.close();
      }
    }
  } );
}

bool SinglePhaseStatistics::execute( real64 const time_n,
                                     real64 const dt,
                                     integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                     integer const GEOS_UNUSED_PARAM( eventCounter ),
                                     real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                     DomainPartition & domain )
{
  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                          MeshLevel & mesh,
                                                                          arrayView1d< string const > const & regionNames )
  {
    // current time is time_n + dt
    computeRegionStatistics( time_n + dt, mesh, regionNames );
  } );
  return false;
}

void SinglePhaseStatistics::computeRegionStatistics( real64 const time,
                                                     MeshLevel & mesh,
                                                     arrayView1d< string const > const & regionNames ) const
{
  GEOS_MARK_FUNCTION;
  // Step 1: initialize the average/min/max quantities
  ElementRegionManager & elemManager = mesh.getElemManager();
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & stats = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.averagePressure = 0.0;
    stats.maxPressure = -LvArray::NumericLimits< real64 >::max;
    stats.minPressure = LvArray::NumericLimits< real64 >::max;

    stats.maxDeltaPressure = -LvArray::NumericLimits< real64 >::max;
    stats.minDeltaPressure = LvArray::NumericLimits< real64 >::max;

    stats.averageTemperature = 0.0;
    stats.maxTemperature = -LvArray::NumericLimits< real64 >::max;
    stats.minTemperature = LvArray::NumericLimits< real64 >::max;

    stats.totalPoreVolume = 0.0;
    stats.totalUncompactedPoreVolume = 0.0;
    stats.totalMass = 0.0;
  }

  // Step 2: increment the average/min/max quantities for all the subRegions
  elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                      ElementSubRegionBase & subRegion )
  {

    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();

    string const & solidName = subRegion.getReference< string >( SinglePhaseBase::viewKeyStruct::solidNamesString() );
    Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );
    CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
    arrayView1d< real64 const > const refPorosity = solid.getReferencePorosity();
    arrayView2d< real64 const > const porosity = solid.getPorosity();

    string const & fluidName = subRegion.template getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
    SingleFluidBase const & fluid = constitutiveModels.getGroup< SingleFluidBase >( fluidName );
    arrayView2d< real64 const > const densities = fluid.density();

    real64 subRegionAvgPresNumerator = 0.0;
    real64 subRegionMinPres = 0.0;
    real64 subRegionMaxPres = 0.0;
    real64 subRegionMinDeltaPres = 0.0;
    real64 subRegionMaxDeltaPres = 0.0;
    real64 subRegionAvgTempNumerator = 0.0;
    real64 subRegionMinTemp = 0.0;
    real64 subRegionMaxTemp = 0.0;
    real64 subRegionTotalUncompactedPoreVol = 0.0;
    real64 subRegionTotalPoreVol = 0.0;
    real64 subRegionTotalMass = 0.0;

    singlePhaseBaseKernels::StatisticsKernel::
      launch( subRegion.size(),
              elemGhostRank,
              volume,
              pres,
              deltaPres,
              temp,
              refPorosity,
              porosity,
              densities,
              subRegionMinPres,
              subRegionAvgPresNumerator,
              subRegionMaxPres,
              subRegionMinDeltaPres,
              subRegionMaxDeltaPres,
              subRegionMinTemp,
              subRegionAvgTempNumerator,
              subRegionMaxTemp,
              subRegionTotalUncompactedPoreVol,
              subRegionTotalPoreVol,
              subRegionTotalMass );

    ElementRegionBase & region = elemManager.getRegion( ElementRegionBase::getParentRegion( subRegion ).getName() );
    RegionStatistics & stats = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.averagePressure += subRegionAvgPresNumerator;
    if( subRegionMinPres < stats.minPressure )
    {
      stats.minPressure = subRegionMinPres;
    }
    if( subRegionMaxPres > stats.maxPressure )
    {
      stats.maxPressure = subRegionMaxPres;
    }

    if( subRegionMinDeltaPres < stats.minDeltaPressure )
    {
      stats.minDeltaPressure = subRegionMinDeltaPres;
    }
    if( subRegionMaxDeltaPres > stats.maxDeltaPressure )
    {
      stats.maxDeltaPressure = subRegionMaxDeltaPres;
    }

    stats.averageTemperature += subRegionAvgTempNumerator;
    if( subRegionMinTemp < stats.minTemperature )
    {
      stats.minTemperature = subRegionMinTemp;
    }
    if( subRegionMaxTemp > stats.maxTemperature )
    {
      stats.maxTemperature = subRegionMaxTemp;
    }

    stats.totalUncompactedPoreVolume += subRegionTotalUncompactedPoreVol;
    stats.totalPoreVolume += subRegionTotalPoreVol;
    stats.totalMass += subRegionTotalMass;
  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & stats = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.minPressure = MpiWrapper::min( stats.minPressure );
    stats.averagePressure = MpiWrapper::sum( stats.averagePressure );
    stats.maxPressure = MpiWrapper::max( stats.maxPressure );

    stats.minDeltaPressure = MpiWrapper::min( stats.minDeltaPressure );
    stats.maxDeltaPressure = MpiWrapper::max( stats.maxDeltaPressure );

    stats.minTemperature = MpiWrapper::min( stats.minTemperature );
    stats.averageTemperature = MpiWrapper::sum( stats.averageTemperature );
    stats.maxTemperature = MpiWrapper::max( stats.maxTemperature );

    stats.totalUncompactedPoreVolume = MpiWrapper::sum( stats.totalUncompactedPoreVolume );
    stats.totalPoreVolume = MpiWrapper::sum( stats.totalPoreVolume );
    stats.totalMass = MpiWrapper::sum( stats.totalMass );

    if( stats.totalUncompactedPoreVolume > 0 )
    {
      float invTotalUncompactedPoreVolume = 1.0 / stats.totalUncompactedPoreVolume;
      stats.averagePressure *= invTotalUncompactedPoreVolume;
      stats.averageTemperature *= invTotalUncompactedPoreVolume;
    }
    else
    {
      stats.averagePressure = 0.0;
      stats.averageTemperature = 0.0;
      GEOS_WARNING( GEOS_FMT( "{}, {}: Cannot compute average pressure & temperature because region pore volume is zero.", getName(), regionNames[i] ) );
    }

    string_view massUnit = units::getSymbol( m_solver->getMassUnit() );

    TableData singPhaseStatsData;
    singPhaseStatsData.addRow( "Pressure[Pa]", stats.minPressure, stats.averagePressure, stats.maxPressure );
    singPhaseStatsData.addRow( "Delta pressure [Pa]", stats.minDeltaPressure, "/", stats.maxDeltaPressure );
    singPhaseStatsData.addRow( "Temperature [K]", stats.minTemperature, stats.averageTemperature, stats.maxTemperature );
    singPhaseStatsData.addSeparator();
    singPhaseStatsData.addSeparator();
    singPhaseStatsData.addRow( "statistics", "phase", CellType::MergeNext, "value" );
    singPhaseStatsData.addSeparator();

    singPhaseStatsData.addRow( "Total dynamic pore volume [rm^3]", "all", CellType::MergeNext, stats.totalPoreVolume );
    singPhaseStatsData.addSeparator();
    singPhaseStatsData.addRow( GEOS_FMT("Total fluid mass [{}]", massUnit), "all", CellType::MergeNext, stats.totalMass );

    string const title = GEOS_FMT( "{}, {} (time {} s):", getName(), regionNames[i], time );
    TableLayout const singPhaseStatsLayout( title, { "statistics", "min", "average", "max" } );
    TableTextFormatter tableFormatter( singPhaseStatsLayout );
    GEOS_LOG_RANK_0( tableFormatter.toString( singPhaseStatsData ) );

    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      outputFile << time << "," << stats.minPressure << "," << stats.averagePressure << "," << stats.maxPressure << "," <<
        stats.minDeltaPressure << "," << stats.maxDeltaPressure << "," <<
        stats.minTemperature << "," << stats.averageTemperature << "," << stats.maxTemperature << "," <<
        stats.totalPoreVolume << "," << stats.totalMass << std::endl;
      outputFile.close();
    }
  }
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SinglePhaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
