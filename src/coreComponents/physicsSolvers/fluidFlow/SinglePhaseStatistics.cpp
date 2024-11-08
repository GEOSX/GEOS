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

#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"

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

RegionSingStatsClass::RegionStatistics( string const & name,
                                        Group * const parent )
  : Group( name, parent )
{
  registerWrapper( viewKeyStruct::averagePressureString(), &m_averagePressure ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "average region pressure" );

  registerWrapper( viewKeyStruct::minPressureString(), &m_minPressure ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "minimum region pressure" );

  registerWrapper( viewKeyStruct::maxPressureString(), &m_maxPressure ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "maximum region pressure" );

  registerWrapper( viewKeyStruct::minDeltaPressureString(), &m_minDeltaPressure ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "minimum region delta pressure" );

  registerWrapper( viewKeyStruct::maxDeltaPressureString(), &m_maxDeltaPressure ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "maximum region delta pressure" );

  registerWrapper( viewKeyStruct::totalMassString(), &m_totalMass ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "fluid mass" );


  registerWrapper( viewKeyStruct::averageTemperatureString(), &m_averageTemperature ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "average region temperature" );

  registerWrapper( viewKeyStruct::minTemperatureString(), &m_minTemperature ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "minimum region temperature" );

  registerWrapper( viewKeyStruct::maxTemperatureString(), &m_maxTemperature ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "maximum region temperature" );


  registerWrapper( viewKeyStruct::totalPoreVolumeString(), &m_totalPoreVolume ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "total region pore volume" );

  registerWrapper( viewKeyStruct::totalUncompactedPoreVolumeString(), &m_totalUncompactedPoreVolume ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "total region uncompacted pore volume" );

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
      region.registerGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

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
                                                     arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;

  // Step 1: initialize the average/min/max quantities
  ElementRegionManager & elemManager = mesh.getElemManager();
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );

    RegionStatistics & stats = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.m_averagePressure =  0.0;
    stats.m_maxPressure =  0.0;
    constexpr real64 max = LvArray::NumericLimits< real64 >::max;
    stats.m_minPressure = max;

    stats.m_maxDeltaPressure = -max;
    stats.m_minDeltaPressure = max;

    stats.m_totalMass = 0.0;

    stats.m_averageTemperature = 0.0;
    stats.m_maxTemperature = 0.0;
    stats.m_minTemperature = max;

    stats.m_totalPoreVolume = 0.0;
    stats.m_totalUncompactedPoreVolume = 0.0;
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

    ElementRegionBase & region =
      elemManager.getRegion( ElementRegionBase::getParentRegion( subRegion ).getName() );
    RegionStatistics & stats =
      region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.m_averagePressure += subRegionAvgPresNumerator;
    if( subRegionMinPres < stats.m_minPressure )
    {
      stats.m_minPressure = subRegionMinPres;
    }
    if( subRegionMaxPres > stats.m_maxPressure )
    {
      stats.m_maxPressure = subRegionMaxPres;
    }
    if( subRegionMinDeltaPres < stats.m_minDeltaPressure )
    {
      stats.m_minDeltaPressure = subRegionMinDeltaPres;
    }
    if( subRegionMaxDeltaPres > stats.m_maxDeltaPressure )
    {
      stats.m_maxDeltaPressure = subRegionMaxDeltaPres;
    }
    stats.m_averageTemperature += subRegionAvgTempNumerator;
    if( subRegionMinTemp < stats.m_minTemperature )
    {
      stats.m_minTemperature = subRegionMinTemp;
    }
    if( subRegionMaxTemp > stats.m_maxTemperature )
    {
      stats.m_maxTemperature = subRegionMaxTemp;
    }

    stats.m_totalUncompactedPoreVolume += subRegionAvgPresNumerator;
    stats.m_totalPoreVolume += subRegionTotalPoreVol;
    stats.m_totalMass += subRegionTotalMass;
  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & stats = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );


    stats.m_minPressure = MpiWrapper::min( stats.m_minPressure );
    stats.m_maxPressure = MpiWrapper::max( stats.m_maxPressure );
    stats.m_averagePressure = MpiWrapper::sum( stats.m_averagePressure );
    stats.m_minDeltaPressure = MpiWrapper::min( stats.m_minDeltaPressure );
    stats.m_maxDeltaPressure = MpiWrapper::max( stats.m_maxDeltaPressure );

    stats.m_minTemperature= MpiWrapper::min( stats.m_minTemperature );
    stats.m_maxTemperature = MpiWrapper::max( stats.m_maxTemperature );
    stats.m_averageTemperature = MpiWrapper::sum( stats.m_averageTemperature );

    stats.m_totalUncompactedPoreVolume = MpiWrapper::sum( stats.m_totalUncompactedPoreVolume );
    stats.m_totalPoreVolume = MpiWrapper::sum( stats.m_totalPoreVolume );
    stats.m_totalMass = MpiWrapper::sum( stats.m_totalMass );

    if( stats.m_totalUncompactedPoreVolume > 0 )
    {
      float invTotalUncompactedPoreVolume = 1.0 / stats.m_totalUncompactedPoreVolume;
      stats.m_averagePressure *= invTotalUncompactedPoreVolume;
      stats.m_averageTemperature *= invTotalUncompactedPoreVolume;
    }
    else
    {
      stats.m_averagePressure = 0.0;
      stats.m_averageTemperature = 0.0;
      GEOS_WARNING( GEOS_FMT( "{}, {}: Cannot compute average pressure & temperature because region pore volume is zero.",
                              getName(), regionNames[i] ) );
    }

    string statPrefix = GEOS_FMT( "{}, {} (time {} s):", getName(), regionNames[i], time );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Pressure (min, average, max): {}, {}, {} Pa",
                                          statPrefix,
                                          stats.m_minPressure, stats.m_averagePressure, stats.m_maxPressure ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Delta pressure (min, max): {}, {} Pa",
                                          statPrefix,
                                          stats.m_minDeltaPressure, stats.m_maxDeltaPressure ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Temperature (min, average, max): {}, {}, {} K",
                                          statPrefix,
                                          stats.m_minTemperature, stats.m_averageTemperature, stats.m_maxTemperature ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Total dynamic pore volume: {} rm^3",
                                          statPrefix,
                                          stats.m_totalPoreVolume ) );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Total fluid mass: {} kg",
                                          statPrefix,
                                          stats.m_totalMass ) );
    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      outputFile << time << "," <<
        stats.m_minPressure << "," << stats.m_averagePressure << "," << stats.m_maxPressure << "," <<
        stats.m_minDeltaPressure << "," << stats.m_maxDeltaPressure << "," <<
        stats.m_minTemperature << "," << stats.m_averageTemperature << "," << stats.m_maxTemperature << "," <<
        stats.m_totalPoreVolume << "," << stats.m_totalMass << std::endl;
      outputFile.close();
    }
  }
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SinglePhaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
