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
 * @file CompositionalMultiphaseStatistics.cpp
 */

#include "CompositionalMultiphaseStatistics.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;

CompositionalMultiphaseStatistics::CompositionalMultiphaseStatistics( const string & name,
                                                                      Group * const parent ):
  Base( name, parent ),
  m_computeCFLNumbers( 0 ),
  m_computeRegionStatistics( 1 )
{
  registerWrapper( viewKeyStruct::computeCFLNumbersString(), &m_computeCFLNumbers ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether CFL numbers are computed or not" );

  registerWrapper( viewKeyStruct::computeRegionStatisticsString(), &m_computeRegionStatistics ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether region statistics are computed or not" );

  registerWrapper( viewKeyStruct::relpermThresholdString(), &m_relpermThreshold ).
    setApplyDefaultValue( 1e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether a phase is considered mobile (when the relperm is above the threshold) or immobile (when the relperm is below the threshold) in metric 2" );

  addLogLevel< logInfo::CFL >();
  addLogLevel< logInfo::Statistics >();
}

RegionCompStatsClass::RegionStatistics( const string & name,
                                        Group * const parent ) :
  Group( name, parent )
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


  registerWrapper( viewKeyStruct::phasePoreVolumeString(), &m_phasePoreVolume ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Phase region phase pore volume" );


  registerWrapper( viewKeyStruct::phaseMassString(), &m_phaseMass ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Region phase mass (trapped and non-trapped, immobile and mobile)" );

  registerWrapper( viewKeyStruct::trappedPhaseMassString(), &m_trappedPhaseMass ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Trapped region phase mass" );

  registerWrapper( viewKeyStruct::immobilePhaseMassString(), &m_immobilePhaseMass ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Immobile region phase mass" );

  registerWrapper( viewKeyStruct::dissolvedComponentMassString(), &m_dissolvedComponentMass ).
    setApplyDefaultValue( 0 ).
    //setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Dissolved region component mass" );
}

void RegionCompStatsClass::init( integer const numPhases, integer const numComps )
{
  m_phasePoreVolume.resizeDimension< 0 >( numPhases );
  m_phaseMass.resizeDimension< 0 >( numPhases );
  m_trappedPhaseMass.resizeDimension< 0 >( numPhases );
  m_immobilePhaseMass.resizeDimension< 0 >( numPhases );
  m_dissolvedComponentMass.resizeDimension< 0, 1 >( numPhases, numComps );
}


void CompositionalMultiphaseStatistics::postInputInitialization()
{
  Base::postInputInitialization();

  if( dynamicCast< CompositionalMultiphaseHybridFVM * >( m_solver ) && m_computeCFLNumbers != 0 )
  {
    GEOS_THROW( GEOS_FMT( "{} {}: the option to compute CFL numbers is incompatible with CompositionalMultiphaseHybridFVM",
                          catalogName(), getDataContext() ),
                InputError );
  }
}

void CompositionalMultiphaseStatistics::registerDataOnMesh( Group & meshBodies )
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

    integer const numPhases = m_solver->numFluidPhases();
    integer const numComps = m_solver->numFluidComponents();

    // if we have to report region statistics, we have to register them first here
    if( m_computeRegionStatistics )
    {

      for( integer i = 0; i < regionNames.size(); ++i )
      {
        ElementRegionBase & region = elemManager.getRegion( regionNames[i] );

        region.registerGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() ).
          setRestartFlags( RestartFlags::NO_WRITE );

        RegionStatistics & stats =
          region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );
        stats.init( numPhases, numComps );
        // write output header
        if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
        {
          std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv" );
          string_view massUnit = units::getSymbol( m_solver->getMassUnit() );
          outputFile <<
            "Time [s],Min pressure [Pa],Average pressure [Pa],Max pressure [Pa],Min delta pressure [Pa],Max delta pressure [Pa]," <<
            "Min temperature [Pa],Average temperature [Pa],Max temperature [Pa],Total dynamic pore volume [rm^3]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Phase " << ip << " dynamic pore volume [rm^3]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Phase " << ip << " mass [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Trapped phase " << ip << " mass (metric 1) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Non-trapped phase " << ip << " mass (metric 1) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Immobile phase " << ip << " mass (metric 2) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
            outputFile << ",Mobile phase " << ip << " mass (metric 2) [" << massUnit << "]";
          for( integer ip = 0; ip < numPhases; ++ip )
          {
            for( integer ic = 0; ic < numComps; ++ic )
              outputFile << ",Component " << ic << " (phase " << ip << ") mass [" << massUnit << "]";
          }
          outputFile << std::endl;
          outputFile.close();
        }
      }
    }

    // if we have to compute CFL numbers later, we need to register additional variables
    if( m_computeCFLNumbers )
    {
      elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                          ElementSubRegionBase & subRegion )
      {
        subRegion.registerField< fields::flow::phaseOutflux >( getName() ).
          reference().resizeDimension< 1 >( numPhases );
        subRegion.registerField< fields::flow::componentOutflux >( getName() ).
          reference().resizeDimension< 1 >( numComps );
        subRegion.registerField< fields::flow::phaseCFLNumber >( getName() );
        subRegion.registerField< fields::flow::componentCFLNumber >( getName() );
      } );
    }
  } );
}

bool CompositionalMultiphaseStatistics::execute( real64 const time_n,
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
    if( m_computeRegionStatistics )
    {
      // current time is time_n + dt
      computeRegionStatistics( time_n + dt, mesh, regionNames );
    }
  } );

  if( m_computeCFLNumbers )
  {
    // current time is time_n + dt
    computeCFLNumbers( time_n + dt, dt, domain );
  }

  return false;
}

void CompositionalMultiphaseStatistics::computeRegionStatistics( real64 const time,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;
  integer const numPhases = m_solver->numFluidPhases();
  integer const numComps = m_solver->numFluidComponents();

  // Step 1: initialize the average/min/max quantities
  ElementRegionManager & elemManager = mesh.getElemManager();
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );

    RegionStatistics & stats =
      region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.m_averagePressure = 0.0;
    stats.m_maxPressure = 0.0;
    constexpr real64 max = LvArray::NumericLimits< real64 >::max;
    stats.m_minPressure = max;

    stats.m_maxDeltaPressure = -max;
    stats.m_minDeltaPressure = max;

    stats.m_averageTemperature = 0.0;
    stats.m_maxTemperature = 0.0;
    stats.m_minTemperature = max;

    stats.m_totalPoreVolume = 0.0;
    stats.m_totalUncompactedPoreVolume = 0.0;

    stats.m_phasePoreVolume.setValues< serialPolicy >( 0.0 );

    stats.m_phaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_trappedPhaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_immobilePhaseMass.setValues< serialPolicy >( 0.0 );
    stats.m_dissolvedComponentMass.setValues< serialPolicy >( 0.0 );
  }

  // Step 2: increment the average/min/max quantities for all the subRegions
  elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                      ElementSubRegionBase & subRegion )
  {

    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      subRegion.getField< fields::flow::phaseVolumeFraction >();
    arrayView1d< real64 const > const deltaPres = subRegion.getField< fields::flow::deltaPressure >();

    Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

    string const & solidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::solidNamesString() );
    CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
    arrayView1d< real64 const > const refPorosity = solid.getReferencePorosity();
    arrayView2d< real64 const > const porosity = solid.getPorosity();

    string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
    MultiFluidBase const & fluid = constitutiveModels.getGroup< MultiFluidBase >( fluidName );
    arrayView3d< real64 const, multifluid::USD_PHASE > const phaseDensity = fluid.phaseDensity();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const phaseCompFraction = fluid.phaseCompFraction();


    //get min vol fraction for each phase to dispactche immobile/mobile mass
    string const & relpermName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase const & relperm = constitutiveModels.getGroup< RelativePermeabilityBase >( relpermName );
    arrayView3d< real64 const, relperm::USD_RELPERM > const phaseTrappedVolFrac = relperm.phaseTrappedVolFraction();
    arrayView3d< real64 const, relperm::USD_RELPERM > const phaseRelperm = relperm.phaseRelPerm();

    real64 subRegionAvgPresNumerator = 0.0;
    real64 subRegionMinPres = 0.0;
    real64 subRegionMaxPres = 0.0;
    real64 subRegionMinDeltaPres = 0.0;
    real64 subRegionMaxDeltaPres = 0.0;
    real64 subRegionAvgTempNumerator = 0.0;
    real64 subRegionMinTemp = 0.0;
    real64 subRegionMaxTemp = 0.0;
    real64 subRegionTotalUncompactedPoreVol = 0.0;
    array1d< real64 > subRegionPhaseDynamicPoreVol( numPhases );
    array1d< real64 > subRegionPhaseMass( numPhases );
    array1d< real64 > subRegionTrappedPhaseMass( numPhases );
    array1d< real64 > subRegionImmobilePhaseMass( numPhases );
    array1d< real64 > subRegionRelpermPhaseMass( numPhases );
    array2d< real64 > subRegiondissolvedComponentMass( numPhases, numComps );

    isothermalCompositionalMultiphaseBaseKernels::
      StatisticsKernel::
      launch< parallelDevicePolicy<> >( subRegion.size(),
                                        numComps,
                                        numPhases,
                                        m_relpermThreshold,
                                        elemGhostRank,
                                        volume,
                                        pres,
                                        deltaPres,
                                        temp,
                                        refPorosity,
                                        porosity,
                                        phaseDensity,
                                        phaseCompFraction,
                                        phaseVolFrac,
                                        phaseTrappedVolFrac,
                                        phaseRelperm,
                                        subRegionMinPres,
                                        subRegionAvgPresNumerator,
                                        subRegionMaxPres,
                                        subRegionMinDeltaPres,
                                        subRegionMaxDeltaPres,
                                        subRegionMinTemp,
                                        subRegionAvgTempNumerator,
                                        subRegionMaxTemp,
                                        subRegionTotalUncompactedPoreVol,
                                        subRegionPhaseDynamicPoreVol.toView(),
                                        subRegionPhaseMass.toView(),
                                        subRegionTrappedPhaseMass.toView(),
                                        subRegionImmobilePhaseMass.toView(),
                                        subRegiondissolvedComponentMass.toView() );

    ElementRegionBase & region = elemManager.getRegion( ElementRegionBase::getParentRegion( subRegion ).getName() );
    RegionStatistics & stats = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.m_averagePressure += subRegionAvgPresNumerator;
    if( subRegionMinPres < stats.m_minPressure )
    {
      stats.m_minPressure = subRegionMinPres;
    }
    if( subRegionMaxPres > stats.m_maxPressure )
    {
      stats.m_maxPressure = subRegionMaxPres;
    }
    if( subRegionMinDeltaPres <  stats.m_minDeltaPressure )
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

    stats.m_totalUncompactedPoreVolume += subRegionTotalUncompactedPoreVol;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      stats.m_phasePoreVolume[ip] += subRegionPhaseDynamicPoreVol[ip];
      stats.m_phaseMass[ip] += subRegionPhaseMass[ip];
      stats.m_trappedPhaseMass[ip] += subRegionTrappedPhaseMass[ip];
      stats.m_immobilePhaseMass[ip] += subRegionImmobilePhaseMass[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        stats.m_dissolvedComponentMass[ip][ic] += subRegiondissolvedComponentMass[ip][ic];
      }
    }

  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & stats = region.getGroup< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.m_minPressure = MpiWrapper::min( stats.m_minPressure );
    stats.m_maxPressure = MpiWrapper::max( stats.m_maxPressure );
    stats.m_minDeltaPressure = MpiWrapper::min( stats.m_minDeltaPressure );
    stats.m_maxDeltaPressure = MpiWrapper::max( stats.m_maxDeltaPressure );
    stats.m_minTemperature = MpiWrapper::min( stats.m_minTemperature );
    stats.m_maxTemperature = MpiWrapper::max( stats.m_maxTemperature );
    stats.m_totalUncompactedPoreVolume = MpiWrapper::sum( stats.m_totalUncompactedPoreVolume );

    array1d< real64 > phasePoreVolume = stats.m_phasePoreVolume;
    array1d< real64 > phaseMass = stats.m_phaseMass;
    array1d< real64 > trappedPhaseMass = stats.m_trappedPhaseMass;
    array1d< real64 > immobilePhaseMass = stats.m_immobilePhaseMass;
    array2d< real64 > dissolvedComponentMass = stats.m_dissolvedComponentMass;

    real64 totalPoreVolume = 0.0;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phasePoreVolume[ip] = MpiWrapper::sum( phasePoreVolume[ip] );
      phaseMass[ip] = MpiWrapper::sum( phaseMass[ip] );
      trappedPhaseMass[ip] = MpiWrapper::sum( trappedPhaseMass[ip] );
      immobilePhaseMass[ip] = MpiWrapper::sum( immobilePhaseMass[ip] );
      totalPoreVolume += phasePoreVolume[ip];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        dissolvedComponentMass[ip][ic] = MpiWrapper::sum( dissolvedComponentMass[ip][ic] );
      }
    }
    stats.m_totalPoreVolume = totalPoreVolume;

    stats.m_averagePressure = MpiWrapper::sum( stats.m_averagePressure );

    stats.m_averageTemperature = MpiWrapper::sum( stats.m_averageTemperature );

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
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                  GEOS_FMT( "{}, {}: Cannot compute average pressure because region pore volume is zero.",
                                            getName(), regionNames[i] ) );
    }


    // helpers to report statistics
    array1d< real64 > nonTrappedPhaseMass( numPhases );
    array1d< real64 > mobilePhaseMass( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      nonTrappedPhaseMass[ip] = phaseMass[ip] - trappedPhaseMass[ip];
      mobilePhaseMass[ip] = phaseMass[ip] - immobilePhaseMass[ip];
    }

    string_view massUnit = units::getSymbol( m_solver->getMassUnit() );
    string statPrefix = GEOS_FMT( "{}, {} (time {} s):", getName(), regionNames[i], time );
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Pressure (min, average, max): {}, {}, {} Pa",
                                          statPrefix,
                                          stats.m_minPressure, stats.m_averagePressure, stats.m_maxPressure ));
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Delta pressure (min, max): {}, {} Pa",
                                          statPrefix, stats.m_minDeltaPressure, stats.m_maxDeltaPressure ));
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Temperature (min, average, max): {}, {}, {} K",
                                          statPrefix,
                                          stats.m_minTemperature, stats.m_averageTemperature,
                                          stats.m_maxTemperature ));
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Total dynamic pore volume: {} rm^3",
                                          statPrefix, stats.m_totalPoreVolume ));
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Phase dynamic pore volume: {} rm^3",
                                          statPrefix, stats.m_phasePoreVolume ));
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Phase mass: {} {}",
                                          statPrefix, stats.m_phaseMass, massUnit ));

    // metric 1: trapping computed with the Land trapping coefficient (similar to Eclipse)
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Trapped phase mass (metric 1): {} {}",
                                          statPrefix, stats.m_trappedPhaseMass, massUnit ));
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Non-trapped phase mass (metric 1): {} {}",
                                          statPrefix, nonTrappedPhaseMass, massUnit ));

    // metric 2: immobile phase mass computed with a threshold on relative permeability
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Immobile phase mass (metric 2): {} {}",
                                          statPrefix, stats.m_immobilePhaseMass, " ", massUnit ));
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Mobile phase mass (metric 2): {} {}",
                                          statPrefix, mobilePhaseMass, massUnit ));

    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                GEOS_FMT( "{} Component mass: {} {}",
                                          statPrefix, stats.m_dissolvedComponentMass, massUnit ));

    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      outputFile << time << "," <<
        stats.m_minPressure << "," << stats.m_averagePressure << "," << stats.m_maxPressure << "," <<
        stats.m_minDeltaPressure << "," << stats.m_maxDeltaPressure << "," << stats.m_minTemperature << "," <<
        stats.m_averageTemperature << "," << stats.m_maxTemperature << "," << stats.m_totalPoreVolume;
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," <<phasePoreVolume[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << phaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << trappedPhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << nonTrappedPhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << immobilePhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << mobilePhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        for( integer ic = 0; ic < numComps; ++ic )
          outputFile << "," << dissolvedComponentMass[ip][ic];
      }
      outputFile << std::endl;
      outputFile.close();
    }
  }
}

void CompositionalMultiphaseStatistics::computeCFLNumbers( real64 const time,
                                                           real64 const dt,
                                                           DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;
  real64 maxPhaseCFL, maxCompCFL;
  m_solver->computeCFLNumbers( domain, dt, maxPhaseCFL, maxCompCFL );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::CFL, GEOS_FMT( "{} (time {} s): Max phase CFL number: {}", getName(), time, maxPhaseCFL ) );
  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::CFL, GEOS_FMT( "{} (time {} s): Max component CFL number: {}", getName(), time, maxCompCFL ) );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        CompositionalMultiphaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
