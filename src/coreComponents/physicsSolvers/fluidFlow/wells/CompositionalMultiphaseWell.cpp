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
 * @file CompositionalMultiphaseWell.cpp
 */

#include "CompositionalMultiphaseWell.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/ThermalCompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/PerforationFluxKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/AccumulationKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalAccumulationKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalSolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalSolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/GlobalComponentFractionKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/PhaseVolumeFractionKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalPhaseVolumeFractionKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/FluidUpdateKernel.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace compositionalMultiphaseWellKernels;

CompositionalMultiphaseWell::CompositionalMultiphaseWell( const string & name,
                                                          Group * const parent )
  :
  WellSolverBase( name, parent ),
  m_useMass( false ),
  m_useTotalMassEquation( 1 ),
  m_maxCompFracChange( 1.0 ),
  m_maxRelativePresChange( 0.2 ),
  m_maxAbsolutePresChange( -1 ), // disabled by default
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 ),
  m_targetPhaseIndex( -1 )
{
  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::useTotalMassEquationString(), &m_useTotalMassEquation ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use total mass equation" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (absolute) change in a component fraction between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::maxRelativeCompDensChangeString(), &m_maxRelativeCompDensChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( LvArray::NumericLimits< real64 >::max/1.0e100 ). // disabled by default
    setDescription( "Maximum (relative) change in a component density between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString(), &m_maxRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (relative) change in pressure between two Newton iterations (recommended with rate control)" );

  this->registerWrapper( viewKeyStruct::maxAbsolutePresChangeString(), &m_maxAbsolutePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).       // disabled by default
    setDescription( "Maximum (absolute) pressure change in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::maxRelativeTempChangeString(), &m_maxRelativeTempChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (relative) change in temperature between two Newton iterations  " );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );
}

void CompositionalMultiphaseWell::postInputInitialization()
{
  WellSolverBase::postInputInitialization();

  GEOS_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) <<
                        ": The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOS_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) <<
                        ": The maximum absolute change in component fraction must larger or equal to 0.0" );

  GEOS_ERROR_IF_LE_MSG( m_maxRelativeCompDensChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxRelativeCompDensChangeString() ) <<
                        ": The maximum relative change in component density must be larger than 0.0" );
}

void CompositionalMultiphaseWell::registerDataOnMesh( Group & meshBodies )
{
  WellSolverBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_referenceFluidModelName.empty() )
      {
        m_referenceFluidModelName = getConstitutiveName< MultiFluidBase >( subRegion );
      }
    } );
  } );

  // 1. Set key dimensions of the problem
  // Empty check needed to avoid errors when running in schema generation mode.
  if( !m_referenceFluidModelName.empty() )
  {
    MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );
    m_numPhases = fluid0.numFluidPhases();
    m_numComponents = fluid0.numFluidComponents();
  }
  m_numDofPerWellElement =  isThermal() ?    m_numComponents + 3 : m_numComponents + 2; // 1 pressure + NC compositions + 1 connectionRate +
                                                                                        // temp if thermal
  m_numDofPerResElement =  isThermal() ? m_numComponents + 2 : m_numComponents + 1; // 1 pressure + NC compositions + temp if thermal

  // loop over the wells
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      string const & fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
      GEOS_ERROR_IF( fluidName.empty(), GEOS_FMT( "{}: Fluid model not found on subregion {}",
                                                  getDataContext(), subRegion.getName() ) );

      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.

      subRegion.registerField< fields::well::globalCompDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< fields::well::globalCompDensity_n >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      subRegion.registerField< fields::well::mixtureConnectionRate >( getName() );
      subRegion.registerField< fields::well::mixtureConnectionRate_n >( getName() );

      subRegion.registerField< fields::well::globalCompFraction >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< fields::well::dGlobalCompFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

      subRegion.registerField< fields::well::phaseVolumeFraction >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< fields::well::dPhaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents + 2 ); // dP, dT, dC

      subRegion.registerField< fields::well::totalMassDensity >( getName() );
      subRegion.registerField< fields::well::dTotalMassDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents +2 ); // dP, dT, dC

      subRegion.registerField< fields::well::phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< fields::well::pressureScalingFactor >( getName() );
      subRegion.registerField< fields::well::temperatureScalingFactor >( getName() );
      subRegion.registerField< fields::well::globalCompDensityScalingFactor >( getName() );

      PerforationData & perforationData = *subRegion.getPerforationData();
      perforationData.registerField< fields::well::compPerforationRate >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      perforationData.registerField< fields::well::dCompPerforationRate >( getName() ).reference().resizeDimension< 1, 2, 3 >( 2, m_numComponents, m_numComponents+ 2 );
      if( fluid.isThermal() )
      {
        perforationData.registerField< fields::well::energyPerforationFlux >( getName() );
        perforationData.registerField< fields::well::dEnergyPerforationFlux >( getName() ).
          reference().resizeDimension< 1, 2 >( 2, m_numComponents+2 );
      }

      WellControls & wellControls = getWellControls( subRegion );
      wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString() );

      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentBHPString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numComponents + 2 );   // dP, dT, dC


      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numPhases );

      wellControls.registerWrapper< array2d< real64 > >( viewKeyStruct::dCurrentPhaseVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0, 1 >( m_numPhases, m_numComponents + 3 );   // dP, dT, dC, dQ

      wellControls.registerWrapper< real64 >( viewKeyStruct::massDensityString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );
      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentTotalVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numComponents + 3 );   // dP, dT, dC dQ

      wellControls.registerWrapper< real64 >( viewKeyStruct::massDensityString() );

      // write rates output header
      // the rank that owns the reference well element is responsible
      if( m_writeCSV > 0 && subRegion.isLocallyOwned() )
      {
        string const wellControlsName = wellControls.getName();
        string const massUnit = m_useMass ? "kg" : "mol";
        integer const useSurfaceConditions = wellControls.useSurfaceConditions();
        string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
        string const unitKey = useSurfaceConditions ? "s" : "r";
        integer const numPhase = m_numPhases;
        // format: time,bhp,total_rate,total_vol_rate,phase0_vol_rate,phase1_vol_rate,...
        std::ofstream outputFile( m_ratesOutputDir + "/" + wellControlsName + ".csv" );
        outputFile << "Time [s],BHP [Pa],Total rate [" << massUnit << "/s],Total " << conditionKey << " Volumetric rate [" << unitKey << "m3/s]";
        for( integer ip = 0; ip < numPhase; ++ip )
          outputFile << ",Phase" << ip << " " << conditionKey << " volumetric rate [" << unitKey << "m3/s]";
        outputFile << std::endl;
        outputFile.close();
      }
    } );
  } );

}

void CompositionalMultiphaseWell::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{

  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
  GEOS_THROW_IF( fluidName.empty(),
                 GEOS_FMT( "{}: Fluid model not found on subregion {}",
                           getDataContext(), subRegion.getName() ),
                 InputError );

  string & relPermName = subRegion.registerWrapper< string >( viewKeyStruct::relPermNamesString() ).
                           setPlotLevel( PlotLevel::NOPLOT ).
                           setRestartFlags( RestartFlags::NO_WRITE ).
                           setSizedFromParent( 0 ).
                           setDescription( "Name of the relative permeability constitutive model to use" ).
                           reference();

  relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );

  GEOS_THROW_IF( relPermName.empty(),
                 GEOS_FMT( "{}: Relative permeability model not found on subregion {}",
                           getDataContext(), subRegion.getName() ),
                 InputError );

}
namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                        GEOS_FMT( "Mismatch in number of phases between constitutive models {} and {}",
                                  lhs.getDataContext(), rhs.getDataContext() ),
                        InputError );

  for( integer ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOS_THROW_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                          GEOS_FMT( "Mismatch in phase names between constitutive models {} and {}",
                                    lhs.getDataContext(), rhs.getDataContext() ),
                          InputError );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                        GEOS_FMT( "Mismatch in number of components between constitutive models {} and {}",
                                  lhs.getDataContext(), rhs.getDataContext() ),
                        InputError );

  for( integer ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOS_THROW_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                          GEOS_FMT( "Mismatch in component names between constitutive models {} and {}",
                                    lhs.getDataContext(), rhs.getDataContext() ),
                          InputError );
  }
}

}

/**
 * @brief Checks if the WellControls parameters are within the fluid tables ranges
 * @param fluid the fluid to check
 */
void CompositionalMultiphaseWell::validateWellControlsForFluid( WellControls const & wellControls,
                                                                MultiFluidBase const & fluid ) const
{
  if( wellControls.useSurfaceConditions() )
  {
    try
    {
      real64 const & surfaceTemp = wellControls.getSurfaceTemperature();
      real64 const & surfacePres = wellControls.getSurfacePressure();
      fluid.checkTablesParameters( surfacePres, surfaceTemp );
    } catch( SimulationError const & ex )
    {
      string const errorMsg = GEOS_FMT( "{}: wrong surface pressure / temperature.\n", getDataContext() );
      throw SimulationError( ex, errorMsg );
    }
  }
}

void CompositionalMultiphaseWell::validateConstitutiveModels( DomainPartition const & domain ) const
{
  GEOS_MARK_FUNCTION;

  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
  string const referenceFluidName = flowSolver.referenceFluidModelName();
  MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion const & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      compareMultiphaseModels( fluid, referenceFluid );
      compareMulticomponentModels( fluid, referenceFluid );

      string const & relpermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPerm = getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
      compareMultiphaseModels( relPerm, referenceFluid );

      WellControls const & wellControls = getWellControls( subRegion );
      validateWellControlsForFluid( wellControls, fluid );
    } );

  } );
}

void CompositionalMultiphaseWell::validateInjectionStreams( WellElementSubRegion const & subRegion ) const
{
  WellControls const & wellControls = getWellControls( subRegion );

  // check well injection stream for injectors
  if( wellControls.isInjector())
  {
    arrayView1d< real64 const > const & injectionStream = wellControls.getInjectionStream();

    integer const streamSize = injectionStream.size();
    GEOS_THROW_IF( ( streamSize == 0 ),
                   "WellControls '" << wellControls.getName() << "'" <<
                   ": Injection stream not specified for well " << subRegion.getName(),
                   InputError );
    GEOS_THROW_IF( ( streamSize != m_numComponents ),
                   "WellControls '" << wellControls.getName() << "'" <<
                   ": Injection stream for well " << subRegion.getName() << " should have " <<
                   m_numComponents << " components.",
                   InputError );

    real64 compFracSum = 0;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      real64 const compFrac = injectionStream[ic];
      GEOS_THROW_IF( ( compFrac < 0.0 ) || ( compFrac > 1.0 ),
                     "WellControls " << wellControls.getDataContext() <<
                     ": Invalid injection stream for well " << subRegion.getName(),
                     InputError );
      compFracSum += compFrac;
    }
    GEOS_THROW_IF( ( compFracSum < 1.0 - std::numeric_limits< real64 >::epsilon() ) ||
                   ( compFracSum > 1.0 + std::numeric_limits< real64 >::epsilon() ),
                   "WellControls " << wellControls.getDataContext() <<
                   ": Invalid injection stream for well " << subRegion.getName(),
                   InputError );
  }
}

void CompositionalMultiphaseWell::validateWellConstraints( real64 const & time_n,
                                                           real64 const & dt,
                                                           WellElementSubRegion const & subRegion )
{
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  // now that we know we are single-phase, we can check a few things in the constraints
  WellControls const & wellControls = getWellControls( subRegion );
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const & targetTotalRate = wellControls.getTargetTotalRate( time_n + dt );
  real64 const & targetPhaseRate = wellControls.getTargetPhaseRate( time_n + dt );

  GEOS_THROW_IF( wellControls.isInjector() && currentControl == WellControls::Control::PHASEVOLRATE,
                 "WellControls " << wellControls.getDataContext() <<
                 ": Phase rate control is not available for injectors",
                 InputError );
  GEOS_THROW_IF( wellControls.isProducer() && currentControl == WellControls::Control::TOTALVOLRATE,
                 "WellControls " << wellControls.getDataContext() <<
                 ": Total rate control is not available for producers",
                 InputError );

  GEOS_THROW_IF( wellControls.isInjector() && targetTotalRate < 0.0,
                 "WellControls " << wellControls.getDataContext() <<
                 ": Target total rate cannot be negative for injectors",
                 InputError );
  GEOS_THROW_IF( wellControls.isInjector() && !isZero( targetPhaseRate ),
                 "WellControls " << wellControls.getDataContext() <<
                 ": Target phase rate cannot be used for injectors",
                 InputError );

  // The user always provides positive rates, but these rates are later multiplied by -1 internally for producers
  GEOS_THROW_IF( wellControls.isProducer() && targetPhaseRate > 0.0,
                 "WellControls " << wellControls.getDataContext() <<
                 ": Target phase rate cannot be negative for producers",
                 InputError );
  GEOS_THROW_IF( wellControls.isProducer() && !isZero( targetTotalRate ),
                 "WellControls " << wellControls.getDataContext() <<
                 ": Target total rate cannot be used for producers",
                 InputError );

  // Find target phase index for phase rate constraint
  for( integer ip = 0; ip < fluid.numFluidPhases(); ++ip )
  {
    if( fluid.phaseNames()[ip] == wellControls.getTargetPhaseName() )
    {
      m_targetPhaseIndex = ip;
    }
  }
  GEOS_THROW_IF( wellControls.isProducer() && m_targetPhaseIndex == -1,
                 "WellControls " << wellControls.getDataContext() <<
                 ": Phase " << wellControls.getTargetPhaseName() << " not found for well control " << wellControls.getName(),
                 InputError );
}

void CompositionalMultiphaseWell::initializePostSubGroups()
{
  WellSolverBase::initializePostSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  validateConstitutiveModels( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      string & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );
      GEOS_THROW_IF( relPermName.empty(),
                     GEOS_FMT( "{}: Relative permeability not found on subregion {}",
                               getDataContext(), subRegion.getName() ),
                     InputError );

      validateInjectionStreams( subRegion );
    } );
  } );
}

void CompositionalMultiphaseWell::initializePostInitialConditionsPreSubGroups()
{
  WellSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    // loop over the wells
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      fluid.setMassFlag( m_useMass );
    } );
  } );
}

void CompositionalMultiphaseWell::updateGlobalComponentFraction( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  isothermalCompositionalMultiphaseBaseKernels::
    GlobalComponentFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               subRegion );

}

void CompositionalMultiphaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }
  using Deriv = multifluid::DerivativeOffset;

  integer const numComp = m_numComponents;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  integer isThermal = fluid.isThermal();
  // subRegion data

  arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();

  arrayView1d< real64 > const & totalMassDens = subRegion.getField< fields::well::totalMassDensity >();
  arrayView2d< real64, compflow::USD_FLUID_DC > const & dTotalMassDens = subRegion.getField< fields::well::dTotalMassDensity >();

  arrayView1d< real64 const > const wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

  // control data

  WellControls & wellControls = getWellControls( subRegion );
  string const wellControlsName = wellControls.getName();
  integer const logLevel = wellControls.getLogLevel();
  real64 const & refGravCoef = wellControls.getReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  arrayView1d< real64 > const & dCurrentBHP =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHPString() );

  geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
  {
    integer constexpr IS_THERMAL = ISTHERMAL();
    GEOS_UNUSED_VAR( NC );
    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [&numComp,
                                pres,
                                totalMassDens,
                                dTotalMassDens,
                                wellElemGravCoef,
                                &currentBHP,
                                &dCurrentBHP,
                                &iwelemRef,
                                &refGravCoef] ( localIndex const )
    {
      real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
      currentBHP = pres[iwelemRef] + totalMassDens[iwelemRef] * diffGravCoef;
      dCurrentBHP[Deriv::dP] =   1 + dTotalMassDens[iwelemRef][Deriv::dP] * diffGravCoef;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dCurrentBHP[Deriv::dC+ic] = dTotalMassDens[iwelemRef][Deriv::dC+ic] * diffGravCoef;
      }
      if constexpr ( IS_THERMAL )
      {
        dCurrentBHP[Deriv::dT] =  dTotalMassDens[iwelemRef][Deriv::dT] * diffGravCoef;
      }
    } );
  } );

  if( logLevel >= 2 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "{}: BHP (at the specified reference elevation) = {} Pa",
                             wellControlsName, currentBHP ) );
  }

}

void CompositionalMultiphaseWell::updateVolRatesForConstraint( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  integer constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const & temp = subRegion.getField< fields::well::temperature >();
  arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::mixtureConnectionRate >();

  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();
  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens = subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >();

  // fluid data

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  integer isThermal = fluid.isThermal();
  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac = fluid.phaseFraction();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac = fluid.dPhaseFraction();

  arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens = fluid.totalDensity();
  arrayView3d< real64 const, multifluid::USD_FLUID_DC > const & dTotalDens = fluid.dTotalDensity();

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens = fluid.phaseDensity();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens = fluid.dPhaseDensity();

  // control data

  WellControls & wellControls = getWellControls( subRegion );
  string const wellControlsName = wellControls.getName();
  integer const logLevel = wellControls.getLogLevel();
  string const massUnit = m_useMass ? "kg" : "mol";

  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 const & surfacePres = wellControls.getSurfacePressure();
  real64 const & surfaceTemp = wellControls.getSurfaceTemperature();

  arrayView1d< real64 > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  arrayView2d< real64 > const & dCurrentPhaseVolRate =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRateString() );

  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
  arrayView1d< real64 > const & dCurrentTotalVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRateString() );

  real64 & massDensity =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::massDensityString() );
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr IS_THERMAL = ISTHERMAL();
      using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NUM_COMP, IS_THERMAL >;
      // bring everything back to host, capture the scalars by reference
      forAll< serialPolicy >( 1, [&numComp,
                                  &numPhase,
                                  fluidWrapper,
                                  pres,
                                  temp,
                                  compFrac,
                                  dCompFrac_dCompDens,
                                  connRate,
                                  totalDens,
                                  dTotalDens,
                                  phaseDens,
                                  dPhaseDens,
                                  phaseFrac,
                                  dPhaseFrac,
                                  &useSurfaceConditions,
                                  &surfacePres,
                                  &surfaceTemp,
                                  &currentTotalVolRate,
                                  dCurrentTotalVolRate,
                                  currentPhaseVolRate,
                                  dCurrentPhaseVolRate,
                                  &iwelemRef,
                                  &logLevel,
                                  &wellControlsName,
                                  &massUnit,
                                  &massDensity] ( localIndex const )
      {
        GEOS_UNUSED_VAR( massUnit );
        using Deriv = multifluid::DerivativeOffset;
        stackArray1d< real64, maxNumComp > work( numComp );
        // Step 1: evaluate the phase and total density in the reference element

        //    We need to evaluate the density as follows:
        //      - Surface conditions: using the surface pressure provided by the user
        //      - Reservoir conditions: using the pressure in the top element
        if( useSurfaceConditions )
        {
          // we need to compute the surface density
          fluidWrapper.update( iwelemRef, 0, surfacePres, surfaceTemp, compFrac[iwelemRef] );
          if( logLevel >= 2 )
          {
            GEOS_LOG_RANK( GEOS_FMT( "{}: surface density computed with P_surface = {} Pa and T_surface = {} K",
                                     wellControlsName, surfacePres, surfaceTemp ) );
#ifdef GEOS_USE_HIP
            GEOS_UNUSED_VAR( wellControlsName );
#endif
          }
        }
        else
        {
          real64 const refPres = pres[iwelemRef];
          fluidWrapper.update( iwelemRef, 0, refPres, temp[iwelemRef], compFrac[iwelemRef] );
        }

        // Step 2: update the total volume rate

        real64 const currentTotalRate = connRate[iwelemRef];

        // Step 2.1: compute the inverse of the total density and derivatives
        massDensity = totalDens[iwelemRef][0];
        real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];

        stackArray1d< real64, maxNumComp > dTotalDensInv_dCompDens( numComp );
        for( integer ic = 0; ic < numComp; ++ic )
        {
          dTotalDensInv_dCompDens[ic] = -dTotalDens[iwelemRef][0][Deriv::dC+ic] * totalDensInv * totalDensInv;
        }
        applyChainRuleInPlace( numComp, dCompFrac_dCompDens[iwelemRef], dTotalDensInv_dCompDens, work.data() );

        // Step 2.2: divide the total mass/molar rate by the total density to get the total volumetric rate
        currentTotalVolRate = currentTotalRate * totalDensInv;
        // Compute derivatives  dP dT
        real64 const dTotalDensInv_dPres = -dTotalDens[iwelemRef][0][Deriv::dP] * totalDensInv * totalDensInv;
        dCurrentTotalVolRate[COFFSET_WJ::dP] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dTotalDensInv_dPres;
        if constexpr ( IS_THERMAL )
        {
          dCurrentTotalVolRate[COFFSET_WJ::dT] = ( useSurfaceConditions ==  0 ) * currentTotalRate * -dTotalDens[iwelemRef][0][Deriv::dT] * totalDensInv * totalDensInv;
        }

        dCurrentTotalVolRate[COFFSET_WJ::dQ] = totalDensInv;
        for( integer ic = 0; ic < numComp; ++ic )
        {
          dCurrentTotalVolRate[COFFSET_WJ::dC+ic] = currentTotalRate * dTotalDensInv_dCompDens[ic];
        }

        if( logLevel >= 2 && useSurfaceConditions )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: total fluid density at surface conditions = {} {}/sm3, total rate = {} {}/s, total surface volumetric rate = {} sm3/s",
                                   wellControlsName, totalDens[iwelemRef][0], massUnit, currentTotalRate, massUnit, currentTotalVolRate ) );
        }

        // Step 3: update the phase volume rate
        for( integer ip = 0; ip < numPhase; ++ip )
        {

          // Step 3.1: compute the inverse of the (phase density * phase fraction) and derivatives

          // skip the rest of this function if phase ip is absent
          bool const phaseExists = (phaseFrac[iwelemRef][0][ip] > 0);
          if( !phaseExists )
          {
            continue;
          }

          real64 const phaseDensInv =  1.0 / phaseDens[iwelemRef][0][ip];
          real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;
          real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac[iwelemRef][0][ip][Deriv::dP] * phaseDensInv
                                                           - dPhaseDens[iwelemRef][0][ip][Deriv::dP] * phaseFracTimesPhaseDensInv * phaseDensInv;


          // Step 3.2: divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
          currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
          dCurrentPhaseVolRate[ip][COFFSET_WJ::dP] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
          dCurrentPhaseVolRate[ip][COFFSET_WJ::dQ] = phaseFracTimesPhaseDensInv;
          if constexpr (IS_THERMAL )
          {
            real64 const dPhaseFracTimesPhaseDensInv_dTemp = dPhaseFrac[iwelemRef][0][ip][Deriv::dT] * phaseDensInv
                                                             - dPhaseDens[iwelemRef][0][ip][Deriv::dT] * phaseFracTimesPhaseDensInv * phaseDensInv;
            dCurrentPhaseVolRate[ip][COFFSET_WJ::dT] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dTemp;
          }

          for( integer ic = 0; ic < numComp; ++ic )
          {
            dCurrentPhaseVolRate[ip][COFFSET_WJ::dC+ic] = -phaseFracTimesPhaseDensInv * dPhaseDens[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
            dCurrentPhaseVolRate[ip][COFFSET_WJ::dC+ic] += dPhaseFrac[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
            dCurrentPhaseVolRate[ip][COFFSET_WJ::dC+ic] *= currentTotalRate;
          }
          applyChainRuleInPlace( numComp, dCompFrac_dCompDens[iwelemRef], &dCurrentPhaseVolRate[ip][COFFSET_WJ::dC], work.data() );

          if( logLevel >= 2 && useSurfaceConditions )
          {
            GEOS_LOG_RANK( GEOS_FMT( "{}: density of phase {} at surface conditions = {} {}/sm3, phase surface volumetric rate = {} sm3/s",
                                     wellControlsName, ip, phaseDens[iwelemRef][0][ip], massUnit, currentPhaseVolRate[ip] ) );
          }
        }
      } );
    } );
  } );
}



void CompositionalMultiphaseWell::updateFluidModel( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const & temp = subRegion.getField< fields::well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    thermalCompositionalMultiphaseBaseKernels::
      FluidUpdateKernel::
      launch< ExecPolicy >( subRegion.size(),
                            fluidWrapper,
                            pres,
                            temp,
                            compFrac );
  } );

}

real64 CompositionalMultiphaseWell::updatePhaseVolumeFraction( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  real64 maxDeltaPhaseVolFrac  =
    m_isThermal ?
    thermalCompositionalMultiphaseBaseKernels::
      PhaseVolumeFractionKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 subRegion,
                                                 fluid )
:    isothermalCompositionalMultiphaseBaseKernels::
      PhaseVolumeFractionKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 subRegion,
                                                 fluid );

  return maxDeltaPhaseVolFrac;
}

void CompositionalMultiphaseWell::updateTotalMassDensity( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  fluid.isThermal() ?
  thermalCompositionalMultiphaseWellKernels::
    TotalMassDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid )
  :
  TotalMassDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid );

}

void CompositionalMultiphaseWell::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  real64 maxPhaseVolFrac = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      real64 const maxRegionPhaseVolFrac = updateSubRegionState( subRegion );
      maxPhaseVolFrac = LvArray::math::max( maxRegionPhaseVolFrac, maxPhaseVolFrac );
    } );
  } );
  maxPhaseVolFrac = MpiWrapper::max( maxPhaseVolFrac );

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Max well phase volume fraction change = {}", getName(), fmt::format( "{:.{}f}", maxPhaseVolFrac, 4 ) ) );

}

real64 CompositionalMultiphaseWell::updateSubRegionState( WellElementSubRegion & subRegion )
{
  // update properties
  updateGlobalComponentFraction( subRegion );

  // update volumetric rates for the well constraints
  // note: this must be called before updateFluidModel
  updateVolRatesForConstraint( subRegion );

  // update densities, phase fractions, phase volume fractions

  updateFluidModel( subRegion );   //  Calculate fluid properties;
  real64 maxPhaseVolChange = updatePhaseVolumeFraction( subRegion );
  updateTotalMassDensity( subRegion );
  // update the current BHP pressure
  updateBHPForConstraint( subRegion );
  return maxPhaseVolChange;
}

void CompositionalMultiphaseWell::initializeWells( DomainPartition & domain, real64 const & time_n, real64 const & dt )
{
  GEOS_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // TODO: change the way we access the flowSolver here
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  // loop over the wells
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();
    PresTempCompFracInitializationKernel::CompFlowAccessors resCompFlowAccessors( mesh.getElemManager(), flowSolver.getName() );
    PresTempCompFracInitializationKernel::MultiFluidAccessors resMultiFluidAccessors( mesh.getElemManager(), flowSolver.getName() );

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls const & wellControls = getWellControls( subRegion );

      if( time_n <= 0.0 ||
          ( !wellControls.isWellOpen( time_n ) && wellControls.isWellOpen( time_n + dt ) )  )
      {

        PerforationData const & perforationData = *subRegion.getPerforationData();

        // get well primary variables on well elements
        arrayView1d< real64 > const & wellElemPressure = subRegion.getField< fields::well::pressure >();
        arrayView1d< real64 > const & wellElemTemp = subRegion.getField< fields::well::temperature >();
        arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens = subRegion.getField< fields::well::globalCompDensity >();
        arrayView1d< real64 > const & connRate = subRegion.getField< fields::well::mixtureConnectionRate >();

        // get the info stored on well elements
        arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac = subRegion.getField< fields::well::globalCompFraction >();
        arrayView1d< real64 const > const & wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

        // get the element region, subregion, index
        arrayView1d< localIndex const > const resElementRegion = perforationData.getField< fields::perforation::reservoirElementRegion >();
        arrayView1d< localIndex const > const resElementSubRegion = perforationData.getField< fields::perforation::reservoirElementSubRegion >();
        arrayView1d< localIndex const > const resElementIndex = perforationData.getField< fields::perforation::reservoirElementIndex >();

        arrayView1d< real64 const > const & perfGravCoef = perforationData.getField< fields::well::gravityCoefficient >();

        // 1) Loop over all perforations to compute an average mixture density and component fraction
        // 2) Initialize the reference pressure
        // 3) Estimate the pressures in the well elements using the average density
        PresTempCompFracInitializationKernel::
          launch( perforationData.size(),
                  subRegion.size(),
                  numComp,
                  numPhase,
                  perforationData.getNumPerforationsGlobal(),
                  wellControls,
                  0.0, // initialization done at t = 0
                  resCompFlowAccessors.get( fields::flow::pressure{} ),
                  resCompFlowAccessors.get( fields::flow::temperature{} ),
                  resCompFlowAccessors.get( fields::flow::globalCompDensity{} ),
                  resCompFlowAccessors.get( fields::flow::phaseVolumeFraction{} ),
                  resMultiFluidAccessors.get( fields::multifluid::phaseMassDensity{} ),
                  resElementRegion,
                  resElementSubRegion,
                  resElementIndex,
                  perfGravCoef,
                  wellElemGravCoef,
                  wellElemPressure,
                  wellElemTemp,
                  wellElemCompFrac );

        // get well secondary variables on well elements
        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
        MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
        arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens = fluid.phaseDensity();
        arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens = fluid.totalDensity();

        // 4) Back calculate component densities
        constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
        {
          typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

          thermalCompositionalMultiphaseBaseKernels::
            FluidUpdateKernel::
            launch< serialPolicy >( subRegion.size(),
                                    fluidWrapper,
                                    wellElemPressure,
                                    wellElemTemp,
                                    wellElemCompFrac );
        } );

        CompDensInitializationKernel::launch( subRegion.size(),
                                              numComp,
                                              wellElemCompFrac,
                                              wellElemTotalDens,
                                              wellElemCompDens );

        // 5) Recompute the pressure-dependent properties
        updateSubRegionState( subRegion );

        // 6) Estimate the well rates
        // TODO: initialize rates using perforation rates
        compositionalMultiphaseWellKernels::
          RateInitializationKernel::
          launch( subRegion.size(),
                  m_targetPhaseIndex,
                  wellControls,
                  0.0, // initialization done at t = 0
                  wellElemPhaseDens,
                  wellElemTotalDens,
                  connRate );
      }
    } );

  } );
}

void CompositionalMultiphaseWell::assembleFluxTerms( real64 const & time,
                                                     real64 const & dt,
                                                     DomainPartition & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;


  string const wellDofKey = dofManager.getKey( wellElementDofName());
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {
      WellControls const & well_controls = getWellControls( subRegion );
      if( well_controls.isWellOpen( time+ dt ) )
      {
        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
        MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
        int numComponents = fluid.numFluidComponents();

        if( isThermal() )
        {
          thermalCompositionalMultiphaseWellKernels::
            FaceBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                       dt,
                                                       dofManager.rankOffset(),
                                                       m_useTotalMassEquation,
                                                       wellDofKey,
                                                       well_controls,
                                                       subRegion,
                                                       fluid,
                                                       localMatrix,
                                                       localRhs );
        }
        else
        {
          compositionalMultiphaseWellKernels::
            FaceBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                       dt,
                                                       dofManager.rankOffset(),
                                                       m_useTotalMassEquation,
                                                       wellDofKey,
                                                       well_controls,
                                                       subRegion,
                                                       localMatrix,
                                                       localRhs );
        }
      }
    } );
  } );

}

void CompositionalMultiphaseWell::assembleAccumulationTerms( real64 const & time,
                                                             real64 const & dt,
                                                             DomainPartition & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time );
  GEOS_UNUSED_VAR( dt );
  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      int numPhases = fluid.numFluidPhases();
      int numComponents = fluid.numFluidComponents();
      WellControls const & wellControls = getWellControls( subRegion );
      if( wellControls.isWellOpen( time+ dt ) )
      {
        if( isThermal() )
        {

          thermalCompositionalMultiphaseWellKernels::
            ElementBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                       numPhases,
                                                       wellControls.isProducer(),
                                                       dofManager.rankOffset(),
                                                       m_useTotalMassEquation,
                                                       wellDofKey,
                                                       subRegion,
                                                       fluid,
                                                       localMatrix,
                                                       localRhs );
        }
        else
        {
          compositionalMultiphaseWellKernels::
            ElementBasedAssemblyKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                       numPhases,
                                                       wellControls.isProducer(),
                                                       dofManager.rankOffset(),
                                                       m_useTotalMassEquation,
                                                       wellDofKey,
                                                       subRegion,
                                                       fluid,
                                                       localMatrix,
                                                       localRhs );
        }
      }
      else
      {
        //wellControls.setWellOpen(false);
        // get the degrees of freedom and ghosting info
        arrayView1d< globalIndex const > const & wellElemDofNumber =
          subRegion.getReference< array1d< globalIndex > >( wellDofKey );
        arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();
        localIndex rank_offset = dofManager.rankOffset();
        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
        {
          if( wellElemGhostRank[ei] < 0 )
          {
            globalIndex const dofIndex = wellElemDofNumber[ei];
            localIndex const localRow = dofIndex - rank_offset;

            real64 unity = 1.0;
            for( integer i=0; i < m_numDofPerWellElement; i++ )
            {
              globalIndex const rindex =  localRow+i;
              globalIndex const cindex =dofIndex + i;
              localMatrix.template addToRow< serialAtomic >( rindex,
                                                             &cindex,
                                                             &unity,
                                                             1 );
              localRhs[rindex] = 0.0;
            }
          }
        } );
      }
    } );
  } );


}


real64
CompositionalMultiphaseWell::calculateResidualNorm( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer numNorm = 1; // mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;

  if( isThermal() )
  {
    numNorm = 2;  // mass balance and energy balance
  }
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );


  globalIndex const rankOffset = dofManager.rankOffset();
  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {


    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {


      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      WellControls const & wellControls = getWellControls( subRegion );

      // step 1: compute the norm in the subRegion

      if( isThermal() )
      {
        real64 subRegionResidualNorm[2]{};

        thermalCompositionalMultiphaseWellKernels::ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     m_targetPhaseIndex,
                                                     rankOffset,
                                                     wellDofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     wellControls,
                                                     time_n + dt,
                                                     dt,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionResidualNorm );
        // step 2: reduction across meshBodies/regions/subRegions

        for( integer i=0; i<numNorm; i++ )
        {
          if( subRegionResidualNorm[i] > localResidualNorm[i] )
          {
            localResidualNorm[i] = subRegionResidualNorm[i];
          }
        }

      }
      else
      {
        real64 subRegionResidualNorm[1]{};
        ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                     numDofPerWellElement(),
                                                     m_targetPhaseIndex,
                                                     rankOffset,
                                                     wellDofKey,
                                                     localRhs,
                                                     subRegion,
                                                     fluid,
                                                     wellControls,
                                                     time_n + dt,
                                                     dt,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionResidualNorm );



        // step 2: reduction across meshBodies/regions/subRegions

        if( subRegionResidualNorm[0] > localResidualNorm[0] )
        {
          localResidualNorm[0] = subRegionResidualNorm[0];
        }
      }
    } );
  } );

  // step 3: second reduction across MPI ranks
  real64 resNorm=localResidualNorm[0];
  if( isThermal() )
  {
    real64 globalResidualNorm[2]{};
    globalResidualNorm[0] = MpiWrapper::max( localResidualNorm[0] );
    globalResidualNorm[1] = MpiWrapper::max( localResidualNorm[1] );
    resNorm=sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );
    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                             coupledSolverAttributePrefix(), globalResidualNorm[0], globalResidualNorm[1] );
    }
  }
  else
  {
    resNorm= MpiWrapper::max( resNorm );
    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), resNorm );
    }
  }
  return resNorm;
}

real64
CompositionalMultiphaseWell::scalingForSystemSolution( DomainPartition & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  real64 scalingFactor = 1.0;
  real64 maxDeltaPres = 0.0, maxDeltaCompDens = 0.0, maxDeltaTemp = 0.0;
  real64 minPresScalingFactor = 1.0, minCompDensScalingFactor = 1.0, minTempScalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure = subRegion.getField< fields::well::pressure >();
      arrayView1d< real64 const > const temperature = subRegion.getField< fields::well::temperature >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::well::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::well::pressureScalingFactor >();
      arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::well::temperatureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::well::globalCompDensityScalingFactor >();
      const integer temperatureOffset = m_numComponents+2;
      auto const subRegionData =
        m_isThermal
  ? thermalCompositionalMultiphaseBaseKernels::
          SolutionScalingKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxAbsolutePresChange,
                                                     m_maxRelativeTempChange,
                                                     m_maxCompFracChange,
                                                     m_maxRelativeCompDensChange,
                                                     pressure,
                                                     temperature,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     temperatureScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     wellDofKey,
                                                     subRegion,
                                                     localSolution,
                                                     temperatureOffset )
  : isothermalCompositionalMultiphaseBaseKernels::
          SolutionScalingKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxAbsolutePresChange,
                                                     m_maxCompFracChange,
                                                     m_maxRelativeCompDensChange,
                                                     pressure,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     wellDofKey,
                                                     subRegion,
                                                     localSolution );


      scalingFactor = std::min( subRegionData.localMinVal, scalingFactor );

      maxDeltaPres  = std::max( maxDeltaPres, subRegionData.localMaxDeltaPres );
      maxDeltaCompDens = std::max( maxDeltaCompDens, subRegionData.localMaxDeltaCompDens );
      maxDeltaTemp = std::max( maxDeltaTemp, subRegionData.localMaxDeltaTemp );
      minPresScalingFactor = std::min( minPresScalingFactor, subRegionData.localMinPresScalingFactor );
      minCompDensScalingFactor = std::min( minCompDensScalingFactor, subRegionData.localMinCompDensScalingFactor );
      minTempScalingFactor = std::min( minTempScalingFactor, subRegionData.localMinTempScalingFactor );
    } );
  } );

  scalingFactor = MpiWrapper::min( scalingFactor );
  maxDeltaPres  = MpiWrapper::max( maxDeltaPres );
  maxDeltaCompDens = MpiWrapper::max( maxDeltaCompDens );
  minPresScalingFactor = MpiWrapper::min( minPresScalingFactor );
  minCompDensScalingFactor = MpiWrapper::min( minCompDensScalingFactor );

  string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Max well pressure change: {} Pa (before scaling)",
                                      getName(), GEOS_FMT( "{:.{}f}", maxDeltaPres, 3 ) ) );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Max well component density change: {} {} (before scaling)",
                                      getName(), GEOS_FMT( "{:.{}f}", maxDeltaCompDens, 3 ), massUnit ) );

  if( m_isThermal )
  {
    maxDeltaTemp = MpiWrapper::max( maxDeltaTemp );
    minTempScalingFactor = MpiWrapper::min( minTempScalingFactor );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Max well temperature change: {} K (before scaling)",
                                        getName(), GEOS_FMT( "{:.{}f}", maxDeltaTemp, 3 ) ) );
  }


  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Min well pressure scaling factor: {}", getName(), minPresScalingFactor ) );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Min well component density scaling factor: {}", getName(), minCompDensScalingFactor ) );
  if( m_isThermal )
  {
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Min well temperature scaling factor: {}", getName(), minTempScalingFactor ) );
  }


  return LvArray::math::max( scalingFactor, m_minScalingFactor );

}

bool
CompositionalMultiphaseWell::checkSystemSolution( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  integer localCheck = 1;
  if( 0 )
  {

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                          [&]( localIndex const,
                                                                               WellElementSubRegion & subRegion )
      {
        arrayView1d< real64 const > const pressure =
          subRegion.getField< fields::well::pressure >();
        arrayView2d< real64 const, compflow::USD_COMP > const compDens =
          subRegion.getField< fields::well::globalCompDensity >();
        arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::well::pressureScalingFactor >();
        arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::well::globalCompDensityScalingFactor >();

        auto const subRegionData =
          compositionalMultiphaseWellKernels::
            SolutionCheckKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                       CompositionalMultiphaseFVM::ScalingType::Global,
                                                       scalingFactor,
                                                       pressure,
                                                       compDens,
                                                       pressureScalingFactor,
                                                       compDensScalingFactor,
                                                       dofManager.rankOffset(),
                                                       m_numComponents,
                                                       wellDofKey,
                                                       subRegion,
                                                       localSolution );

        if( !subRegionData.localMinVal )
        {
          GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Solution,
                                      GEOS_FMT( "Solution is invalid in well {} (either a negative pressure or a negative component density was found)", subRegion.getName()) );
        }

        localCheck = std::min( localCheck, subRegionData.localMinVal );
      } );
    } );
  }
  else
  {

    real64 minPres = 0.0, minDens = 0.0, minTotalDens = 0.0;
    integer numNegPres = 0, numNegDens = 0, numNegTotalDens = 0;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions( regionNames,
                                                  [&]( localIndex const,
                                                       ElementSubRegionBase & subRegion )
      {
        //integer const m_allowCompDensChopping(true);
        integer const m_allowNegativePressure( false );
        CompositionalMultiphaseFVM::ScalingType const m_scalingType( CompositionalMultiphaseFVM::ScalingType::Global );
        arrayView1d< real64 const > const pressure =
          subRegion.getField< fields::well::pressure >();
        arrayView1d< real64 const > const temperature =
          subRegion.getField< fields::well::temperature >();
        arrayView2d< real64 const, compflow::USD_COMP > const compDens =
          subRegion.getField< fields::well::globalCompDensity >();
        arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::well::pressureScalingFactor >();
        arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::well::temperatureScalingFactor >();
        arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::well::globalCompDensityScalingFactor >();

        // check that pressure and component densities are non-negative
        // for thermal, check that temperature is above 273.15 K
        const integer temperatureOffset = m_numComponents+2;
        auto const subRegionData =
          m_isThermal
  ? thermalCompositionalMultiphaseBaseKernels::
            SolutionCheckKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                       m_allowNegativePressure,
                                                       m_scalingType,
                                                       scalingFactor,
                                                       pressure,
                                                       temperature,
                                                       compDens,
                                                       pressureScalingFactor,
                                                       temperatureScalingFactor,
                                                       compDensScalingFactor,
                                                       dofManager.rankOffset(),
                                                       m_numComponents,
                                                       wellDofKey,
                                                       subRegion,
                                                       localSolution,
                                                       temperatureOffset )
  : isothermalCompositionalMultiphaseBaseKernels::
            SolutionCheckKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                       m_allowNegativePressure,
                                                       m_scalingType,
                                                       scalingFactor,
                                                       pressure,
                                                       compDens,
                                                       pressureScalingFactor,
                                                       compDensScalingFactor,
                                                       dofManager.rankOffset(),
                                                       m_numComponents,
                                                       wellDofKey,
                                                       subRegion,
                                                       localSolution );

        localCheck = std::min( localCheck, subRegionData.localMinVal );

        minPres  = std::min( minPres, subRegionData.localMinPres );
        minDens = std::min( minDens, subRegionData.localMinDens );
        minTotalDens = std::min( minTotalDens, subRegionData.localMinTotalDens );
        numNegPres += subRegionData.localNumNegPressures;
        numNegDens += subRegionData.localNumNegDens;
        numNegTotalDens += subRegionData.localNumNegTotalDens;
      } );
    } );

    minPres  = MpiWrapper::min( minPres );
    minDens = MpiWrapper::min( minDens );
    minTotalDens = MpiWrapper::min( minTotalDens );
    numNegPres = MpiWrapper::sum( numNegPres );
    numNegDens = MpiWrapper::sum( numNegDens );
    numNegTotalDens = MpiWrapper::sum( numNegTotalDens );

    if( numNegPres > 0 )
      GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Number of negative well pressure values: {}, minimum value: {} Pa",
                                          getName(), numNegPres, fmt::format( "{:.{}f}", minPres, 3 ) ) );
    string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
    if( numNegDens > 0 )
      GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Number of negative well component density values: {}, minimum value: {} {} ",
                                          getName(), numNegDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );
    if( minTotalDens > 0 )
      GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Number of negative total well density values: {}, minimum value: {} {} ",
                                          getName(), minTotalDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );

  }

  return MpiWrapper::min( localCheck );
}

void CompositionalMultiphaseWell::computePerforationRates( real64 const & time_n, real64 const & dt, DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
  {

    // TODO: change the way we access the flowSolver here
    CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion & subRegion )
    {
      PerforationData * const perforationData = subRegion.getPerforationData();
      WellControls const & wellControls = getWellControls( subRegion );
      if( wellControls.isWellOpen( time_n+ dt ) )
      {

        bool const disableReservoirToWellFlow = wellControls.isInjector() and !wellControls.isCrossflowEnabled();

        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
        MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
        bool isThermal = fluid.isThermal();

        if( isThermal )
        {
          thermalPerforationFluxKernels::
            PerforationFluxKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       flowSolver.getName(),
                                                       perforationData,
                                                       subRegion,
                                                       fluid,
                                                       elemManager,
                                                       disableReservoirToWellFlow );
        }
        else
        {
          isothermalPerforationFluxKernels::
            PerforationFluxKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       m_numPhases,
                                                       flowSolver.getName(),
                                                       perforationData,
                                                       subRegion,
                                                       elemManager,
                                                       disableReservoirToWellFlow );
        }
      }
      else
      {
        // Zero completion flow rate
        arrayView2d< real64 > const compPerfRate = perforationData->getField< fields::well::compPerforationRate >();
        for( integer iperf=0; iperf<perforationData->size(); iperf++ )
        {
          for( integer ic = 0; ic < m_numComponents; ++ic )
          {
            compPerfRate[iperf][ic] = 0.0;
          }
        }
      }
    } );

  } );

}


void
CompositionalMultiphaseWell::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  real64 const dt,
                                                  DomainPartition & domain )
{


  DofManager::CompMask pressureMask( m_numDofPerWellElement, 0, 1 );
  DofManager::CompMask componentMask( m_numDofPerWellElement, 1, numFluidComponents()+1 );
  DofManager::CompMask connRateMask( m_numDofPerWellElement, numFluidComponents()+1, numFluidComponents()+2 );
  GEOS_UNUSED_VAR( dt );
  // update all the fields using the global damping coefficients
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::globalCompDensity::key(),
                               scalingFactor,
                               componentMask );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::mixtureConnectionRate::key(),
                               scalingFactor,
                               connRateMask );

  if( isThermal() )
  {
    DofManager::CompMask temperatureMask( m_numDofPerWellElement, numFluidComponents()+2, numFluidComponents()+3 );

    dofManager.addVectorToField( localSolution,
                                 wellElementDofName(),
                                 fields::well::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );

  }
  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    // synchronize
    FieldIdentifiers fieldsToBeSync;
    if( isThermal() )
    {
      fieldsToBeSync.addElementFields( { fields::well::pressure::key(),
                                         fields::well::globalCompDensity::key(),
                                         fields::well::mixtureConnectionRate::key(),
                                         fields::well::temperature::key() },
                                       regionNames );
    }
    else
    {
      fieldsToBeSync.addElementFields( { fields::well::pressure::key(),
                                         fields::well::globalCompDensity::key(),
                                         fields::well::mixtureConnectionRate::key() },
                                       regionNames );
    }
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );


}

void CompositionalMultiphaseWell::chopNegativeDensities( DomainPartition & domain )
{
  integer const numComp = m_numComponents;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getField< fields::well::globalCompDensity >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        if( wellElemGhostRank[iwelem] < 0 )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            // we allowed for some densities to be slightly negative in CheckSystemSolution
            // if the new density is negative, chop back to zero
            if( wellElemCompDens[iwelem][ic] < 0 )
            {
              wellElemCompDens[iwelem][ic] = 0;
            }
          }
        }
      } );
    } );

  } );
}


void CompositionalMultiphaseWell::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      // get a reference to the primary variables on well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getField< fields::well::pressure >();
      arrayView1d< real64 const > const & wellElemPressure_n =
        subRegion.getField< fields::well::pressure_n >();
      wellElemPressure.setValues< parallelDevicePolicy<> >( wellElemPressure_n );

      if( isThermal() )
      {
        // get a reference to the primary variables on well elements
        arrayView1d< real64 > const & wellElemTemperature =
          subRegion.getField< fields::well::temperature >();
        arrayView1d< real64 const > const & wellElemTemperature_n =
          subRegion.getField< fields::well::temperature_n >();
        wellElemTemperature.setValues< parallelDevicePolicy<> >( wellElemTemperature_n );
      }
      arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity =
        subRegion.getField< fields::well::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemGlobalCompDensity_n =
        subRegion.getField< fields::well::globalCompDensity_n >();
      wellElemGlobalCompDensity.setValues< parallelDevicePolicy<> >( wellElemGlobalCompDensity_n );

      arrayView1d< real64 > const & connRate =
        subRegion.getField< fields::well::mixtureConnectionRate >();
      arrayView1d< real64 const > const & connRate_n =
        subRegion.getField< fields::well::mixtureConnectionRate_n >();
      connRate.setValues< parallelDevicePolicy<> >( connRate_n );

      updateSubRegionState( subRegion );
    } );
  } );
}

void CompositionalMultiphaseWell::assemblePressureRelations( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel const & mesh,
                                                                 arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      WellControls & wellControls = getWellControls( subRegion );

      if( wellControls.isWellOpen( time_n+ dt ) )
      {
        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
        MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
        bool isThermal = fluid.isThermal();
        // get the degrees of freedom, depth info, next welem index
        string const wellDofKey = dofManager.getKey( wellElementDofName() );
        arrayView1d< globalIndex const > const & wellElemDofNumber =
          subRegion.getReference< array1d< globalIndex > >( wellDofKey );
        arrayView1d< real64 const > const & wellElemGravCoef =
          subRegion.getField< fields::well::gravityCoefficient >();
        arrayView1d< localIndex const > const & nextWellElemIndex =
          subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

        // get primary variables on well elements
        arrayView1d< real64 const > const & wellElemPres =
          subRegion.getField< fields::well::pressure >();

        // get total mass density on well elements (for potential calculations)
        arrayView1d< real64 const > const & wellElemTotalMassDens =
          subRegion.getField< fields::well::totalMassDensity >();
        arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens =
          subRegion.getField< fields::well::dTotalMassDensity >();

        bool controlHasSwitched = false;
        isothermalCompositionalMultiphaseBaseKernels::
          KernelLaunchSelectorCompTherm< PressureRelationKernel >( numFluidComponents(),
                                                                   isThermal,
                                                                   subRegion.size(),
                                                                   dofManager.rankOffset(),
                                                                   subRegion.isLocallyOwned(),
                                                                   subRegion.getTopWellElementIndex(),
                                                                   m_targetPhaseIndex,
                                                                   wellControls,
                                                                   time_n + dt, // controls evaluated with BHP/rate of the end of step
                                                                   wellElemDofNumber,
                                                                   wellElemGravCoef,
                                                                   nextWellElemIndex,
                                                                   wellElemPres,
                                                                   wellElemTotalMassDens,
                                                                   dWellElemTotalMassDens,
                                                                   controlHasSwitched,
                                                                   localMatrix,
                                                                   localRhs );

        if( controlHasSwitched )
        {
          // TODO: move the switch logic into wellControls
          // TODO: implement a more general switch when more then two constraints per well type are allowed

          real64 const timeAtEndOfStep = time_n + dt;

          if( wellControls.getControl() == WellControls::Control::BHP )
          {
            if( wellControls.isProducer() )
            {
              wellControls.switchToPhaseRateControl( wellControls.getTargetPhaseRate( timeAtEndOfStep ) );
              GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::WellControl,
                                          GEOS_FMT( "Control switch for well {} from BHP constraint to phase volumetric rate constraint", subRegion.getName() ) );
            }
            else
            {
              wellControls.switchToTotalRateControl( wellControls.getTargetTotalRate( timeAtEndOfStep ) );
              GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::WellControl,
                                          GEOS_FMT( "Control switch for well {} from BHP constraint to total volumetric rate constraint", subRegion.getName()) );
            }
          }
          else
          {
            wellControls.switchToBHPControl( wellControls.getTargetBHP( timeAtEndOfStep ) );
            GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::WellControl,
                                        GEOS_FMT( "Control switch for well {} from rate constraint to BHP constraint", subRegion.getName() ) );
          }
        }
      }

    } );
  } );
}

void CompositionalMultiphaseWell::shutDownWell( real64 const time_n,
                                                real64 const dt,
                                                DomainPartition const & domain,
                                                DofManager const & dofManager,
                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      // if the well is open, we don't have to do anything, so we just return
      WellControls const & wellControls = getWellControls( subRegion );
      if( wellControls.isWellOpen( time_n + dt ) )
      {
        return;
      }

      globalIndex const rankOffset = dofManager.rankOffset();

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      arrayView1d< real64 const > const pres =
        subRegion.getField< fields::well::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::well::globalCompDensity >();
      arrayView1d< real64 const > const connRate =
        subRegion.getField< fields::well::mixtureConnectionRate >();

      integer const numComp = m_numComponents;

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 4.1. Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    pres[ei],   // freeze the current pressure value
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 4.2. For each component, apply target global density value
        for( integer ic = 0; ic < numComp; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      compDens[ei][ic],   // freeze the current component density values
                                                      compDens[ei][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }

        // 4.3. Apply rate value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex + numComp + 1,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    connRate[ei],   // freeze the current pressure value
                                                    connRate[ei] );
        localRhs[localRow + numComp + 1] = rhsValue;

      } );
    } );
  } );
}

void CompositionalMultiphaseWell::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  WellSolverBase::implicitStepSetup( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getField< fields::well::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemGlobalCompDensity =
        subRegion.getField< fields::well::globalCompDensity >();
      arrayView1d< real64 const > const & wellElemTemperature =
        subRegion.getField< fields::well::temperature >();

      arrayView1d< real64 > const & wellElemPressure_n =
        subRegion.getField< fields::well::pressure_n >();
      wellElemPressure_n.setValues< parallelDevicePolicy<> >( wellElemPressure );

      if( isThermal() )
      {

        arrayView1d< real64 > const & wellElemTemperature_n =
          subRegion.getField< fields::well::temperature_n >();
        wellElemTemperature_n.setValues< parallelDevicePolicy<> >( wellElemTemperature );
      }

      arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity_n =
        subRegion.getField< fields::well::globalCompDensity_n >();
      wellElemGlobalCompDensity_n.setValues< parallelDevicePolicy<> >( wellElemGlobalCompDensity );

      arrayView1d< real64 const > const & connRate =
        subRegion.getField< fields::well::mixtureConnectionRate >();
      arrayView1d< real64 > const & connRate_n =
        subRegion.getField< fields::well::mixtureConnectionRate_n >();
      connRate_n.setValues< parallelDevicePolicy<> >( connRate );

      arrayView2d< real64 const, compflow::USD_PHASE > const wellElemPhaseVolFrac =
        subRegion.getField< fields::well::phaseVolumeFraction >();
      arrayView2d< real64, compflow::USD_PHASE > const wellElemPhaseVolFrac_n =
        subRegion.getField< fields::well::phaseVolumeFraction_n >();
      wellElemPhaseVolFrac_n.setValues< parallelDevicePolicy<> >( wellElemPhaseVolFrac );

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      fluid.saveConvergedState();

      validateWellConstraints( time_n, dt, subRegion );

      updateSubRegionState( subRegion );
    } );
  } );
}

void CompositionalMultiphaseWell::implicitStepComplete( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  WellSolverBase::implicitStepComplete( time_n, dt, domain );

  if( getLogLevel() > 0 )
  {
    printRates( time_n, dt, domain );
  }
}

void CompositionalMultiphaseWell::printRates( real64 const & time_n,
                                              real64 const & dt,
                                              DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      // the rank that owns the reference well element is responsible for the calculations below.
      if( !subRegion.isLocallyOwned() )
      {
        return;
      }

      integer const numPhase = m_numPhases;

      // control data

      WellControls const & wellControls = getWellControls( subRegion );
      string const wellControlsName = wellControls.getName();

      // format: time,total_rate,total_vol_rate,phase0_vol_rate,phase1_vol_rate,...
      std::ofstream outputFile;
      if( m_writeCSV > 0 )
      {
        outputFile.open( m_ratesOutputDir + "/" + wellControlsName + ".csv", std::ios_base::app );
        outputFile << time_n + dt;
      }

      if( !wellControls.isWellOpen( time_n + dt ) )
      {
        GEOS_LOG( GEOS_FMT( "{}: well is shut", wellControlsName ) );
        if( outputFile.is_open())
        {
          // print all zeros in the rates file
          outputFile << ",0.0,0.0,0.0";
          for( integer ip = 0; ip < numPhase; ++ip )
            outputFile << ",0.0";
          outputFile << std::endl;
          outputFile.close();
        }
        return;
      }

      localIndex const iwelemRef = subRegion.getTopWellElementIndex();
      string const massUnit = m_useMass ? "kg" : "mol";

      // subRegion data

      arrayView1d< real64 const > const & connRate =
        subRegion.getField< fields::well::mixtureConnectionRate >();

      integer const useSurfaceConditions = wellControls.useSurfaceConditions();

      real64 const & currentBHP =
        wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
      arrayView1d< real64 const > const & currentPhaseVolRate =
        wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
      real64 const & currentTotalVolRate =
        wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

      // bring everything back to host, capture the scalars by reference
      forAll< serialPolicy >( 1, [&numPhase,
                                  &useSurfaceConditions,
                                  &currentBHP,
                                  connRate,
                                  &currentTotalVolRate,
                                  currentPhaseVolRate,
                                  &iwelemRef,
                                  &wellControlsName,
                                  &massUnit,
                                  &outputFile] ( localIndex const )
      {
        string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
        string const unitKey = useSurfaceConditions ? "s" : "r";

        real64 const currentTotalRate = connRate[iwelemRef];
        GEOS_LOG( GEOS_FMT( "{}: BHP (at the specified reference elevation): {} Pa",
                            wellControlsName, currentBHP ) );
        GEOS_LOG( GEOS_FMT( "{}: Total rate: {} {}/s; total {} volumetric rate: {} {}m3/s",
                            wellControlsName, currentTotalRate, massUnit, conditionKey, currentTotalVolRate, unitKey ) );
        for( integer ip = 0; ip < numPhase; ++ip )
          GEOS_LOG( GEOS_FMT( "{}: Phase {} {} volumetric rate: {} {}m3/s",
                              wellControlsName, ip, conditionKey, currentPhaseVolRate[ip], unitKey ) );
        if( outputFile.is_open())
        {
          outputFile << "," << currentBHP;
          outputFile << "," << currentTotalRate << "," << currentTotalVolRate;
          for( integer ip = 0; ip < numPhase; ++ip )
            outputFile << "," << currentPhaseVolRate[ip];
          outputFile << std::endl;
          outputFile.close();
        }
      } );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseWell, string const &, Group * const )
}     // namespace geos
