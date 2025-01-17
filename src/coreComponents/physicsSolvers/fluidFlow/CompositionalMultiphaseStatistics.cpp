/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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

#include "mesh/DomainPartition.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/StatisticsKernel.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"


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
        RegionStatistics & stats = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

        stats.phasePoreVolume.resizeDimension< 0 >( numPhases );
        stats.phaseMass.resizeDimension< 0 >( numPhases );
        stats.trappedPhaseMass.resizeDimension< 0 >( numPhases );
        stats.immobilePhaseMass.resizeDimension< 0 >( numPhases );
        stats.componentMass.resizeDimension< 0, 1 >( numPhases, numComps );

        if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
        {
          auto addStatsValue = []( std::ostringstream & pstatsLayout, TableLayout & ptableLayout,
                                   string const & description, string_view pmassUnit,
                                   integer pnumPhases, integer pnumComps = 0 )
          {
            for( int ip = 0; ip < pnumPhases; ++ip )
            {
              if( pnumComps == 0 )
              {
                pstatsLayout << description << " (phase " << ip << ") [" << pmassUnit << "]";
              }
              else
              {
                for( int ic = 0; ic < pnumComps; ++ic )
                {
                  pstatsLayout << "Component " << ic << " (phase " << ip << ") mass [" << pmassUnit << "]";
                  if( ic == 0 )
                  {
                    pstatsLayout << ",";
                  }
                }
              }
              if( ip == 0 )
              {
                pstatsLayout << ",";
              }
            }

            ptableLayout.addToColumns( pstatsLayout.str());
            pstatsLayout.str( "" );
          };

          string_view massUnit = units::getSymbol( m_solver->getMassUnit() );

          TableLayout tableLayout( {
              TableLayout::Column().setName( "Time [s]" ),
              TableLayout::Column().setName( "Min pressure [Pa]" ),
              TableLayout::Column().setName( "Average pressure [Pa]" ),
              TableLayout::Column().setName( "Max pressure [Pa]" ),
              TableLayout::Column().setName( "Min delta pressure [Pa]" ),
              TableLayout::Column().setName( "Max delta pressure [Pa]" ),
              TableLayout::Column().setName( "Min temperature [Pa]" ),
              TableLayout::Column().setName( "Average temperature [Pa]" ),
              TableLayout::Column().setName( "Max temperature [Pa]" ),
              TableLayout::Column().setName( "Total dynamic pore volume [rm^3]" ),
              TableLayout::Column().setName( GEOS_FMT( "Phase mass [{}] dynamic pore volume [rm^3]", massUnit ) ),
            } );

          std::ostringstream statsLayout;
          addStatsValue( statsLayout, tableLayout, "Phase dynamic pore volume", massUnit, numPhases );
          addStatsValue( statsLayout, tableLayout, "Phase", massUnit, numPhases );
          addStatsValue( statsLayout, tableLayout, "Trapped phase mass (metric 1)", massUnit, numPhases );
          addStatsValue( statsLayout, tableLayout, "Non-trapped phase mass (metric 1)", massUnit, numPhases );
          addStatsValue( statsLayout, tableLayout, "Immobile phase mass (metric 2)", massUnit, numPhases );
          addStatsValue( statsLayout, tableLayout, "Mobile phase mass (metric 2)", massUnit, numPhases );
          addStatsValue( statsLayout, tableLayout, "Component", massUnit, numPhases, numComps );

          std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv" );
          TableCSVFormatter csvFormatter( tableLayout );
          outputFile << csvFormatter.headerToString();
        }
      }
    }
  } );

  // if we have to compute CFL numbers later, we need to register additional variables
  if( m_computeCFLNumbers )
  {
    m_solver->registerDataForCFL( meshBodies );
  }
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
    RegionStatistics & stats = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.averagePressure = 0.0;
    stats.maxPressure = 0.0;
    stats.minPressure = LvArray::NumericLimits< real64 >::max;

    stats.maxDeltaPressure = -LvArray::NumericLimits< real64 >::max;
    stats.minDeltaPressure = LvArray::NumericLimits< real64 >::max;

    stats.averageTemperature = 0.0;
    stats.maxTemperature = 0.0;
    stats.minTemperature = LvArray::NumericLimits< real64 >::max;

    stats.totalPoreVolume = 0.0;
    stats.totalUncompactedPoreVolume = 0.0;
    stats.phasePoreVolume.setValues< serialPolicy >( 0.0 );

    stats.phaseMass.setValues< serialPolicy >( 0.0 );
    stats.trappedPhaseMass.setValues< serialPolicy >( 0.0 );
    stats.immobilePhaseMass.setValues< serialPolicy >( 0.0 );
    stats.componentMass.setValues< serialPolicy >( 0.0 );
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
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      stats.phasePoreVolume[ip] += subRegionPhaseDynamicPoreVol[ip];
      stats.phaseMass[ip] += subRegionPhaseMass[ip];
      stats.trappedPhaseMass[ip] += subRegionTrappedPhaseMass[ip];
      stats.immobilePhaseMass[ip] += subRegionImmobilePhaseMass[ip];

      for( integer ic = 0; ic < numComps; ++ic )
      {
        stats.componentMass[ip][ic] += subRegionComponentMass[ip][ic];
      }
    }

  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer i = 0; i < regionNames.size(); ++i )
  {
    ElementRegionBase & region = elemManager.getRegion( regionNames[i] );
    RegionStatistics & stats = region.getReference< RegionStatistics >( viewKeyStruct::regionStatisticsString() );

    stats.minPressure = MpiWrapper::min( stats.minPressure );
    stats.maxPressure = MpiWrapper::max( stats.maxPressure );
    stats.minDeltaPressure = MpiWrapper::min( stats.minDeltaPressure );
    stats.maxDeltaPressure = MpiWrapper::max( stats.maxDeltaPressure );
    stats.minTemperature = MpiWrapper::min( stats.minTemperature );
    stats.maxTemperature = MpiWrapper::max( stats.maxTemperature );
    stats.totalUncompactedPoreVolume = MpiWrapper::sum( stats.totalUncompactedPoreVolume );
    stats.totalPoreVolume = 0.0;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      stats.phasePoreVolume[ip] = MpiWrapper::sum( stats.phasePoreVolume[ip] );
      stats.phaseMass[ip] = MpiWrapper::sum( stats.phaseMass[ip] );
      stats.trappedPhaseMass[ip] = MpiWrapper::sum( stats.trappedPhaseMass[ip] );
      stats.immobilePhaseMass[ip] = MpiWrapper::sum( stats.immobilePhaseMass[ip] );
      stats.totalPoreVolume += stats.phasePoreVolume[ip];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        stats.componentMass[ip][ic] = MpiWrapper::sum( stats.componentMass[ip][ic] );
      }
    }
    stats.averagePressure = MpiWrapper::sum( stats.averagePressure );
    stats.averageTemperature = MpiWrapper::sum( stats.averageTemperature );
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
      GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Statistics,
                                  GEOS_FMT( "{}, {}: Cannot compute average pressure because region pore volume is zero.", getName(), regionNames[i] ) );
    }


    // helpers to report statistics
    array1d< real64 > nonTrappedPhaseMass( numPhases );
    array1d< real64 > mobilePhaseMass( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      nonTrappedPhaseMass[ip] = stats.phaseMass[ip] - stats.trappedPhaseMass[ip];
      mobilePhaseMass[ip] = stats.phaseMass[ip] - stats.immobilePhaseMass[ip];
    }

    string_view massUnit = units::getSymbol( m_solver->getMassUnit() );

    std::vector< string > phaseCompName;
    phaseCompName.reserve( numPhases*numComps );
    std::vector< string > massValues;
    phaseCompName.reserve( numPhases*numComps );

    ConstitutiveManager const & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
    MultiFluidBase const & fluid = constitutiveManager.getGroup< MultiFluidBase >( m_solver->referenceFluidModelName() );
    auto const phaseNames = fluid.phaseNames();
    auto const componentNames = fluid.componentNames();
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        std::stringstream ss;
        ss << phaseNames[ip]<< ", " <<componentNames[ic];
        phaseCompName.push_back( ss.str() );
        massValues.push_back( GEOS_FMT( "{}", stats.componentMass[ip][ic] ) );
      }
    }

    TableData compPhaseStatsData;
    compPhaseStatsData.addRow( "Pressure[Pa]", stats.minPressure, stats.averagePressure, stats.maxPressure );
    compPhaseStatsData.addRow( "Delta pressure [Pa]", stats.minDeltaPressure, "/", stats.maxDeltaPressure );
    compPhaseStatsData.addRow( "Temperature [K]", stats.minTemperature, stats.averageTemperature, stats.maxTemperature );
    compPhaseStatsData.addSeparator();

    compPhaseStatsData.addSeparator();
    compPhaseStatsData.addRow( "statistics", "phase/component", CellType::MergeNext, "value" );
    compPhaseStatsData.addSeparator();

    compPhaseStatsData.addRow( "Total dynamic pore volume [rm^3]", "all", CellType::MergeNext, stats.totalPoreVolume );
    compPhaseStatsData.addSeparator();
    compPhaseStatsData.addRow( "Phase dynamic pore volume: [rm^3]",
                               stringutilities::joinLambda( phaseNames, "\n", []( auto data ) { return data[0]; } ),
                               CellType::MergeNext,
                               stringutilities::joinLambda( stats.phasePoreVolume, "\n", []( auto data ) { return data[0]; } ) );
    compPhaseStatsData.addSeparator();

    compPhaseStatsData.addRow( GEOS_FMT( "Phase mass [{}]", massUnit ),
                               stringutilities::joinLambda( phaseNames, "\n", []( auto data ) { return data[0]; } ),
                               CellType::MergeNext,
                               stringutilities::joinLambda( stats.phaseMass, "\n", []( auto data ) { return data[0]; } ) );
    compPhaseStatsData.addSeparator();

    compPhaseStatsData.addRow( GEOS_FMT( "Trapped phase mass (metric 1) [{}]", massUnit ),
                               stringutilities::joinLambda( phaseNames, "\n", []( auto value ) { return value[0]; } ),
                               CellType::MergeNext,
                               stringutilities::joinLambda( stats.trappedPhaseMass, "\n", []( auto value ) { return value[0]; } ) );
    compPhaseStatsData.addSeparator();
    compPhaseStatsData.addRow( GEOS_FMT( "nonTrappedPhaseMass [{}]", massUnit ),
                               stringutilities::joinLambda( phaseNames, "\n", []( auto value ) { return value[0]; } ),
                               CellType::MergeNext,
                               stringutilities::joinLambda( nonTrappedPhaseMass, "\n", []( auto value ) { return value[0]; } ) );
    compPhaseStatsData.addSeparator();

    compPhaseStatsData.addRow( GEOS_FMT( "Immobile phase mass (metric 2) [{}]", massUnit ),
                               stringutilities::joinLambda( phaseNames, "\n", []( auto value ) { return value[0]; } ),
                               CellType::MergeNext,
                               stringutilities::joinLambda( stats.immobilePhaseMass, "\n", []( auto value ) { return value[0]; } )  );
    compPhaseStatsData.addSeparator();
    compPhaseStatsData.addRow( GEOS_FMT( "Mobile phase mass (metric 2) [{}]", massUnit ),
                               stringutilities::joinLambda( phaseNames, "\n", []( auto value ) { return value[0]; } ),
                               CellType::MergeNext,
                               stringutilities::joinLambda( mobilePhaseMass, "\n", []( auto value ) { return value[0]; } ) );
    compPhaseStatsData.addSeparator();

    compPhaseStatsData.addRow( GEOS_FMT( "Component mass [{}]", massUnit ),
                               stringutilities::join( phaseCompName, '\n' ),
                               CellType::MergeNext,
                               stringutilities::join( massValues, '\n' ) );

    string const title = GEOS_FMT( "{}, {} (time {} s):", getName(), regionNames[i], time );
    TableLayout const compPhaseStatsLayout( title, { "statistics", "min", "average", "max" } );
    TableTextFormatter tableFormatter( compPhaseStatsLayout );
    GEOS_LOG_RANK_0( tableFormatter.toString( compPhaseStatsData ) );

    if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
    {
      TableData tableData;
      tableData.addRow( time, stats.minPressure, stats.averagePressure, stats.maxPressure, stats.minDeltaPressure, stats.maxDeltaPressure,
                        stats.minTemperature, stats.averageTemperature, stats.maxTemperature, stats.totalPoreVolume,
                        stringutilities::joinLambda( stats.phasePoreVolume, "\n", []( auto data ) { return data[0]; } ),
                        stringutilities::joinLambda( stats.phaseMass, "\n", []( auto data ) { return data[0]; } ),
                        stringutilities::joinLambda( stats.trappedPhaseMass, "\n", []( auto value ) { return value[0]; } ),
                        stringutilities::joinLambda( nonTrappedPhaseMass, "\n", []( auto value ) { return value[0]; } ),
                        stringutilities::joinLambda( stats.immobilePhaseMass, "\n", []( auto value ) { return value[0]; } ),
                        stringutilities::joinLambda( mobilePhaseMass, "\n", []( auto value ) { return value[0]; } ),
                        stringutilities::join( massValues, '\n' ) );

      std::ofstream outputFile( m_outputDir + "/" + regionNames[i] + ".csv", std::ios_base::app );
      TableCSVFormatter const csvOutput;
      outputFile << csvOutput.dataToString( tableData );
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
