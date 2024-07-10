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
 * @file SinglePhaseStatistics.cpp
 */

#include "SinglePhaseStatistics.hpp"

#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

SinglePhaseStatistics::SinglePhaseStatistics( const string & name,
                                              Group * const parent ):
  Base( name, parent )
{}

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
    RegionStatistics & regionStatistics = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    regionStatistics.averagePressure = 0.0;
    regionStatistics.maxPressure = -LvArray::NumericLimits< real64 >::max;
    regionStatistics.minPressure = LvArray::NumericLimits< real64 >::max;

    regionStatistics.maxDeltaPressure = -LvArray::NumericLimits< real64 >::max;
    regionStatistics.minDeltaPressure = LvArray::NumericLimits< real64 >::max;

    regionStatistics.averageTemperature = 0.0;
    regionStatistics.maxTemperature = -LvArray::NumericLimits< real64 >::max;
    regionStatistics.minTemperature = LvArray::NumericLimits< real64 >::max;

    regionStatistics.totalPoreVolume = 0.0;
    regionStatistics.totalUncompactedPoreVolume = 0.0;
    regionStatistics.totalMass = 0.0;
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
    RegionStatistics & regionStatistics = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    regionStatistics.averagePressure += subRegionAvgPresNumerator;
    if( subRegionMinPres < regionStatistics.minPressure )
    {
      regionStatistics.minPressure = subRegionMinPres;
    }
    if( subRegionMaxPres > regionStatistics.maxPressure )
    {
      regionStatistics.maxPressure = subRegionMaxPres;
    }

    if( subRegionMinDeltaPres < regionStatistics.minDeltaPressure )
    {
      regionStatistics.minDeltaPressure = subRegionMinDeltaPres;
    }
    if( subRegionMaxDeltaPres > regionStatistics.maxDeltaPressure )
    {
      regionStatistics.maxDeltaPressure = subRegionMaxDeltaPres;
    }

    regionStatistics.averageTemperature += subRegionAvgTempNumerator;
    if( subRegionMinTemp < regionStatistics.minTemperature )
    {
      regionStatistics.minTemperature = subRegionMinTemp;
    }
    if( subRegionMaxTemp > regionStatistics.maxTemperature )
    {
      regionStatistics.maxTemperature = subRegionMaxTemp;
    }

    regionStatistics.totalUncompactedPoreVolume += subRegionTotalUncompactedPoreVol;
    regionStatistics.totalPoreVolume += subRegionTotalPoreVol;
    regionStatistics.totalMass += subRegionTotalMass;
  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & regionStatistics = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    regionStatistics.minPressure = MpiWrapper::min( regionStatistics.minPressure );
    regionStatistics.averagePressure = MpiWrapper::sum( regionStatistics.averagePressure );
    regionStatistics.maxPressure = MpiWrapper::max( regionStatistics.maxPressure );

    regionStatistics.minDeltaPressure = MpiWrapper::min( regionStatistics.minDeltaPressure );
    regionStatistics.maxDeltaPressure = MpiWrapper::max( regionStatistics.maxDeltaPressure );

    regionStatistics.minTemperature = MpiWrapper::min( regionStatistics.minTemperature );
    regionStatistics.averageTemperature = MpiWrapper::sum( regionStatistics.averageTemperature );
    regionStatistics.maxTemperature = MpiWrapper::max( regionStatistics.maxTemperature );

    regionStatistics.totalUncompactedPoreVolume = MpiWrapper::sum( regionStatistics.totalUncompactedPoreVolume );
    regionStatistics.totalPoreVolume = MpiWrapper::sum( regionStatistics.totalPoreVolume );
    regionStatistics.totalMass = MpiWrapper::sum( regionStatistics.totalMass );

    if( regionStatistics.totalUncompactedPoreVolume > 0 )
    {
      float invTotalUncompactedPoreVolume = 1.0 / regionStatistics.totalUncompactedPoreVolume;
      regionStatistics.averagePressure *= invTotalUncompactedPoreVolume;
      regionStatistics.averageTemperature *= invTotalUncompactedPoreVolume;
    }
    else
    {
      regionStatistics.averagePressure = 0.0;
      regionStatistics.averageTemperature = 0.0;
      GEOS_WARNING( GEOS_FMT( "{}, {}: Cannot compute average pressure & temperature because region pore volume is zero.", getName(), regionNames[i] ) );
    }

    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Pressure (min, average, max): {}, {}, {} Pa",
                                        getName(), regionNames[i], time, regionStatistics.minPressure, regionStatistics.averagePressure, regionStatistics.maxPressure ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Delta pressure (min, max): {}, {} Pa",
                                        getName(), regionNames[i], time, regionStatistics.minDeltaPressure, regionStatistics.maxDeltaPressure ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Temperature (min, average, max): {}, {}, {} K",
                                        getName(), regionNames[i], time, regionStatistics.minTemperature, regionStatistics.averageTemperature, regionStatistics.maxTemperature ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Total dynamic pore volume: {} rm^3",
                                        getName(), regionNames[i], time, regionStatistics.totalPoreVolume ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Total fluid mass: {} kg",
                                        getName(), regionNames[i], time, regionStatistics.totalMass ) );

    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      outputFile << time << "," << regionStatistics.minPressure << "," << regionStatistics.averagePressure << "," << regionStatistics.maxPressure << "," <<
        regionStatistics.minDeltaPressure << "," << regionStatistics.maxDeltaPressure << "," <<
        regionStatistics.minTemperature << "," << regionStatistics.averageTemperature << "," << regionStatistics.maxTemperature << "," <<
        regionStatistics.totalPoreVolume << "," << regionStatistics.totalMass << std::endl;
      outputFile.close();
    }
  }
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SinglePhaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
