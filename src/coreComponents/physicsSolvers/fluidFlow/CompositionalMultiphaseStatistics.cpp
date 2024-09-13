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

        region.registerWrapper< RegionStatistics >( viewKeyStruct::regionStatisticsString() ).
          setRestartFlags( RestartFlags::NO_WRITE );
        region.excludeWrappersFromPacking( { viewKeyStruct::regionStatisticsString() } );
        RegionStatistics & regionStatistics = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

        regionStatistics.phasePoreVolume.resizeDimension< 0 >( numPhases );
        regionStatistics.phaseMass.resizeDimension< 0 >( numPhases );
        regionStatistics.trappedPhaseMass.resizeDimension< 0 >( numPhases );
        regionStatistics.immobilePhaseMass.resizeDimension< 0 >( numPhases );
        regionStatistics.componentMass.resizeDimension< 0, 1 >( numPhases, numComps );

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
                                                                 arrayView1d< string const > const & regionNames ) const
{
  GEOS_MARK_FUNCTION;

  integer const numPhases = m_solver->numFluidPhases();
  integer const numComps = m_solver->numFluidComponents();

  // Step 1: initialize the average/min/max quantities
  ElementRegionManager & elemManager = mesh.getElemManager();
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & regionStatistics = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    regionStatistics.averagePressure = 0.0;
    regionStatistics.maxPressure = 0.0;
    regionStatistics.minPressure = LvArray::NumericLimits< real64 >::max;

    regionStatistics.maxDeltaPressure = -LvArray::NumericLimits< real64 >::max;
    regionStatistics.minDeltaPressure = LvArray::NumericLimits< real64 >::max;

    regionStatistics.averageTemperature = 0.0;
    regionStatistics.maxTemperature = 0.0;
    regionStatistics.minTemperature = LvArray::NumericLimits< real64 >::max;

    regionStatistics.totalPoreVolume = 0.0;
    regionStatistics.totalUncompactedPoreVolume = 0.0;
    regionStatistics.phasePoreVolume.setValues< serialPolicy >( 0.0 );

    regionStatistics.phaseMass.setValues< serialPolicy >( 0.0 );
    regionStatistics.trappedPhaseMass.setValues< serialPolicy >( 0.0 );
    regionStatistics.immobilePhaseMass.setValues< serialPolicy >( 0.0 );
    regionStatistics.componentMass.setValues< serialPolicy >( 0.0 );
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
    arrayView3d< real64 const, relperm::USD_PHASE > const phaseTrappedVolFrac = relperm.phaseTrappedVolFraction();
    arrayView4d< real64 const, relperm::USD_RELPERM > const phaseRelperm = relperm.phaseRelPerm();

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
    array2d< real64 > subRegionComponentMass( numPhases, numComps );

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
                                        subRegionComponentMass.toView() );

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
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      regionStatistics.phasePoreVolume[ip] += subRegionPhaseDynamicPoreVol[ip];
      regionStatistics.phaseMass[ip] += subRegionPhaseMass[ip];
      regionStatistics.trappedPhaseMass[ip] += subRegionTrappedPhaseMass[ip];
      regionStatistics.immobilePhaseMass[ip] += subRegionImmobilePhaseMass[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        regionStatistics.componentMass[ip][ic] += subRegionComponentMass[ip][ic];
      }
    }

  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & regionStatistics = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    regionStatistics.minPressure = MpiWrapper::min( regionStatistics.minPressure );
    regionStatistics.maxPressure = MpiWrapper::max( regionStatistics.maxPressure );
    regionStatistics.minDeltaPressure = MpiWrapper::min( regionStatistics.minDeltaPressure );
    regionStatistics.maxDeltaPressure = MpiWrapper::max( regionStatistics.maxDeltaPressure );
    regionStatistics.minTemperature = MpiWrapper::min( regionStatistics.minTemperature );
    regionStatistics.maxTemperature = MpiWrapper::max( regionStatistics.maxTemperature );
    regionStatistics.totalUncompactedPoreVolume = MpiWrapper::sum( regionStatistics.totalUncompactedPoreVolume );
    regionStatistics.totalPoreVolume = 0.0;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      regionStatistics.phasePoreVolume[ip] = MpiWrapper::sum( regionStatistics.phasePoreVolume[ip] );
      regionStatistics.phaseMass[ip] = MpiWrapper::sum( regionStatistics.phaseMass[ip] );
      regionStatistics.trappedPhaseMass[ip] = MpiWrapper::sum( regionStatistics.trappedPhaseMass[ip] );
      regionStatistics.immobilePhaseMass[ip] = MpiWrapper::sum( regionStatistics.immobilePhaseMass[ip] );
      regionStatistics.totalPoreVolume += regionStatistics.phasePoreVolume[ip];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        regionStatistics.componentMass[ip][ic] = MpiWrapper::sum( regionStatistics.componentMass[ip][ic] );
      }
    }
    regionStatistics.averagePressure = MpiWrapper::sum( regionStatistics.averagePressure );
    regionStatistics.averageTemperature = MpiWrapper::sum( regionStatistics.averageTemperature );
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
      GEOS_LOG_LEVEL_RANK_0( 1, getName() << ", " << regionNames[i]
                                          << ": Cannot compute average pressure because region pore volume is zero." );
    }


    // helpers to report statistics
    array1d< real64 > nonTrappedPhaseMass( numPhases );
    array1d< real64 > mobilePhaseMass( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      nonTrappedPhaseMass[ip] = regionStatistics.phaseMass[ip] - regionStatistics.trappedPhaseMass[ip];
      mobilePhaseMass[ip] = regionStatistics.phaseMass[ip] - regionStatistics.immobilePhaseMass[ip];
    }

    string_view massUnit = units::getSymbol( m_solver->getMassUnit() );

    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Pressure (min, average, max): {}, {}, {} Pa",
                                        getName(), regionNames[i], time, regionStatistics.minPressure, regionStatistics.averagePressure, regionStatistics.maxPressure ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Delta pressure (min, max): {}, {} Pa",
                                        getName(), regionNames[i], time, regionStatistics.minDeltaPressure, regionStatistics.maxDeltaPressure ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Temperature (min, average, max): {}, {}, {} K",
                                        getName(), regionNames[i], time, regionStatistics.minTemperature, regionStatistics.averageTemperature, regionStatistics.maxTemperature ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Total dynamic pore volume: {} rm^3",
                                        getName(), regionNames[i], time, regionStatistics.totalPoreVolume ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Phase dynamic pore volume: {} rm^3",
                                        getName(), regionNames[i], time, regionStatistics.phasePoreVolume ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Phase mass: {} {}",
                                        getName(), regionNames[i], time, regionStatistics.phaseMass, massUnit ) );

    // metric 1: trapping computed with the Land trapping coefficient (similar to Eclipse)
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Trapped phase mass (metric 1): {} {}",
                                        getName(), regionNames[i], time, regionStatistics.trappedPhaseMass, massUnit ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Non-trapped phase mass (metric 1): {} {}",
                                        getName(), regionNames[i], time, nonTrappedPhaseMass, massUnit ) );

    // metric 2: immobile phase mass computed with a threshold on relative permeability
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Immobile phase mass (metric 2): {} {}",
                                        getName(), regionNames[i], time, regionStatistics.immobilePhaseMass, massUnit ) );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Mobile phase mass (metric 2): {} {}",
                                        getName(), regionNames[i], time, mobilePhaseMass, massUnit ) );

    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}, {} (time {} s): Component mass: {} {}",
                                        getName(), regionNames[i], time, regionStatistics.componentMass, massUnit ) );

    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      outputFile << time << "," << regionStatistics.minPressure << "," << regionStatistics.averagePressure << "," << regionStatistics.maxPressure << "," <<
        regionStatistics.minDeltaPressure << "," << regionStatistics.maxDeltaPressure << "," << regionStatistics.minTemperature << "," <<
        regionStatistics.averageTemperature << "," << regionStatistics.maxTemperature << "," << regionStatistics.totalPoreVolume;
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << regionStatistics.phasePoreVolume[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << regionStatistics.phaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << regionStatistics.trappedPhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << nonTrappedPhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << regionStatistics.immobilePhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
        outputFile << "," << mobilePhaseMass[ip];
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        for( integer ic = 0; ic < numComps; ++ic )
          outputFile << "," << regionStatistics.componentMass[ip][ic];
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

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{} (time {} s): Max phase CFL number: {}", getName(), time, maxPhaseCFL ) );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{} (time {} s): Max component CFL number: {}", getName(), time, maxCompCFL ) );
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        CompositionalMultiphaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
