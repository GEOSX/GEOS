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
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/solid/SolidInternalEnergy.hpp"
#include "constitutive/contact/HydraulicApertureBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseKernels.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

template< typename POROUSWRAPPER_TYPE >
void updatePorosityAndPermeabilityFromPressureAndTemperature( POROUSWRAPPER_TYPE porousWrapper,
                                                              CellElementSubRegion & subRegion,
                                                              arrayView1d< real64 const > const & pressure,
                                                              arrayView1d< real64 const > const & pressure_k,
                                                              arrayView1d< real64 const > const & pressure_n,
                                                              arrayView1d< real64 const > const & temperature,
                                                              arrayView1d< real64 const > const & temperature_k,
                                                              arrayView1d< real64 const > const & temperature_n )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromPressureAndTemperature( k, q,
                                                           pressure[k],
                                                           pressure_k[k],
                                                           pressure_n[k],
                                                           temperature[k],
                                                           temperature_k[k],
                                                           temperature_n[k] );
    }
  } );
}


template< typename POROUSWRAPPER_TYPE >
void updatePorosityAndPermeabilityFromPressureAndAperture( POROUSWRAPPER_TYPE porousWrapper,
                                                           SurfaceElementSubRegion & subRegion,
                                                           arrayView1d< real64 const > const & pressure,
                                                           arrayView1d< real64 const > const & oldHydraulicAperture,
                                                           arrayView1d< real64 const > const & newHydraulicAperture )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromPressureAndAperture( k, q,
                                                        pressure[k],
                                                        oldHydraulicAperture[k],
                                                        newHydraulicAperture[k] );
    }
  } );
}

FlowSolverBase::FlowSolverBase( string const & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_numDofPerCell( 0 ),
  m_isThermal( 0 ),
  m_keepFlowVariablesConstantDuringInitStep( 0 ),
  m_isFixedStressPoromechanicsUpdate( false ),
  m_isJumpStabilized( false ),
  m_isLaggingFractureStencilWeightsUpdate( 0 )
{
  this->registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not." );

  this->registerWrapper( viewKeyStruct::allowNegativePressureString(), &m_allowNegativePressure ).
    setApplyDefaultValue( 1 ). // negative pressure is allowed by default
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating if negative pressure is allowed" );

  this->registerWrapper( viewKeyStruct::maxAbsolutePresChangeString(), &m_maxAbsolutePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).       // disabled by default
    setDescription( "Maximum (absolute) pressure change in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::maxSequentialPresChangeString(), &m_maxSequentialPresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e5 ).           // 0.1 bar = 1e5 Pa
    setDescription( "Maximum (absolute) pressure change in a sequential iteration, used for outer loop convergence check" );

  this->registerWrapper( viewKeyStruct::maxSequentialTempChangeString(), &m_maxSequentialTempChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "Maximum (absolute) temperature change in a sequential iteration, used for outer loop convergence check" );

  // allow the user to select a norm
  getNonlinearSolverParameters().getWrapper< solverBaseKernels::NormType >( NonlinearSolverParameters::viewKeysStruct::normTypeString() ).setInputFlag( InputFlags::OPTIONAL );
}

void FlowSolverBase::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< fields::flow::gravityCoefficient >( getName() ).
        setApplyDefaultValue( 0.0 );
      subRegion.registerField< fields::flow::netToGross >( getName() );

      subRegion.registerField< fields::flow::pressure >( getName() );
      subRegion.registerField< fields::flow::pressure_n >( getName() );
      subRegion.registerField< fields::flow::initialPressure >( getName() );
      subRegion.registerField< fields::flow::deltaPressure >( getName() ); // for reporting/stats purposes
      subRegion.registerField< fields::flow::bcPressure >( getName() ); // needed for the application of boundary conditions
      if( m_isFixedStressPoromechanicsUpdate )
      {
        subRegion.registerField< fields::flow::pressure_k >( getName() ); // needed for the fixed-stress porosity update
      }

      subRegion.registerField< fields::flow::temperature >( getName() );
      subRegion.registerField< fields::flow::temperature_n >( getName() );
      subRegion.registerField< fields::flow::initialTemperature >( getName() );
      subRegion.registerField< fields::flow::bcTemperature >( getName() ); // needed for the application of boundary conditions
      if( m_isFixedStressPoromechanicsUpdate )
      {
        subRegion.registerField< fields::flow::temperature_k >( getName() ); // needed for the fixed-stress porosity update
      }
      if( m_isThermal )
      {
        subRegion.registerField< fields::flow::energy >( getName() );
        subRegion.registerField< fields::flow::energy_n >( getName() );
      }
    } );

    elemManager.forElementSubRegionsComplete< SurfaceElementSubRegion >( [&]( localIndex const,
                                                                              localIndex const,
                                                                              ElementRegionBase & region,
                                                                              SurfaceElementSubRegion & subRegion )
    {
      SurfaceElementRegion & faceRegion = dynamicCast< SurfaceElementRegion & >( region );

      subRegion.registerField< fields::flow::gravityCoefficient >( getName() );

      subRegion.registerField< fields::flow::aperture0 >( getName() ).
        setApplyDefaultValue( faceRegion.getDefaultAperture() );

      subRegion.registerField< fields::flow::hydraulicAperture >( getName() ).
        setApplyDefaultValue( faceRegion.getDefaultAperture() );

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::flow::gravityCoefficient >( getName() ).
      setApplyDefaultValue( 0.0 );
    faceManager.registerField< fields::flow::transMultiplier >( getName() );

  } );

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {

    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
    fluxApprox.addFieldName( fields::flow::pressure::key() );
    fluxApprox.setCoeffName( fields::permeability::permeability::key() );
    if( m_isThermal )
    {
      fluxApprox.addFieldName( fields::flow::temperature::key() );
    }
  }
}

void FlowSolverBase::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  arrayView1d< real64 const > const pres = subRegion.template getField< fields::flow::pressure >();
  arrayView1d< real64 > const pres_n = subRegion.template getField< fields::flow::pressure_n >();
  pres_n.setValues< parallelDevicePolicy<> >( pres );

  arrayView1d< real64 const > const temp = subRegion.template getField< fields::flow::temperature >();
  arrayView1d< real64 > const temp_n = subRegion.template getField< fields::flow::temperature_n >();
  temp_n.setValues< parallelDevicePolicy<> >( temp );

  if( m_isThermal )
  {
    arrayView1d< real64 const > const energy = subRegion.template getField< fields::flow::energy >();
    arrayView1d< real64 > const energy_n = subRegion.template getField< fields::flow::energy_n >();
    energy_n.setValues< parallelDevicePolicy<> >( energy );
  }

  if( m_isFixedStressPoromechanicsUpdate )
  {
    arrayView1d< real64 > const pres_k = subRegion.template getField< fields::flow::pressure_k >();
    arrayView1d< real64 > const temp_k = subRegion.template getField< fields::flow::temperature_k >();
    pres_k.setValues< parallelDevicePolicy<> >( pres );
    temp_k.setValues< parallelDevicePolicy<> >( temp );
  }
}

void FlowSolverBase::saveSequentialIterationState( DomainPartition & domain )
{
  GEOS_ASSERT( m_isFixedStressPoromechanicsUpdate );

  real64 maxPresChange = 0.0;
  real64 maxTempChange = 0.0;
  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&]( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions ( regionNames,
                                                 [&]( localIndex const,
                                                      ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 > const pres_k = subRegion.getField< fields::flow::pressure_k >();
      arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
      arrayView1d< real64 > const temp_k = subRegion.getField< fields::flow::temperature_k >();

      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPresChange( 0.0 );
      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxTempChange( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          subRegionMaxPresChange.max( LvArray::math::abs( pres[ei] - pres_k[ei] ) );
          pres_k[ei] = pres[ei];
          subRegionMaxTempChange.max( LvArray::math::abs( temp[ei] - temp_k[ei] ) );
          temp_k[ei] = temp[ei];
        }
      } );

      maxPresChange = LvArray::math::max( maxPresChange, subRegionMaxPresChange.get() );
      maxTempChange = LvArray::math::max( maxTempChange, subRegionMaxTempChange.get() );
    } );
  } );

  // store to be later used in convergence check
  m_sequentialPresChange = MpiWrapper::max( maxPresChange );
  m_sequentialTempChange = m_isThermal ? MpiWrapper::max( maxTempChange ) : 0.0;
}

void FlowSolverBase::enableFixedStressPoromechanicsUpdate()
{
  m_isFixedStressPoromechanicsUpdate = true;
}

void FlowSolverBase::enableJumpStabilization()
{
  m_isJumpStabilized = true;
}

void FlowSolverBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  SolverBase::setConstitutiveNamesCallSuper( subRegion );

  subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  subRegion.registerWrapper< string >( viewKeyStruct::solidNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  solidName = getConstitutiveName< CoupledSolidBase >( subRegion );
  GEOS_ERROR_IF( solidName.empty(), GEOS_FMT( "{}: Solid model not found on subregion {}",
                                              getDataContext(), subRegion.getName() ) );

  subRegion.registerWrapper< string >( viewKeyStruct::permeabilityNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
  permName = getConstitutiveName< PermeabilityBase >( subRegion );
  GEOS_ERROR_IF( permName.empty(), GEOS_FMT( "{}: Permeability model not found on subregion {}",
                                             getDataContext(), subRegion.getName() ) );

  if( m_isThermal )
  {
    string & solidInternalEnergyName = subRegion.registerWrapper< string >( viewKeyStruct::solidInternalEnergyNamesString() ).
                                         setPlotLevel( PlotLevel::NOPLOT ).
                                         setRestartFlags( RestartFlags::NO_WRITE ).
                                         setSizedFromParent( 0 ).
                                         setDescription( "Name of the solid internal energy constitutive model to use" ).
                                         reference();

    solidInternalEnergyName = getConstitutiveName< SolidInternalEnergy >( subRegion );
    GEOS_THROW_IF( solidInternalEnergyName.empty(),
                   GEOS_FMT( "{}: Solid internal energy model not found on subregion {}",
                             getDataContext(), subRegion.getName() ),
                   InputError );
  }
}

void FlowSolverBase::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  GEOS_UNUSED_VAR( subRegion );
}

void FlowSolverBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                  MeshLevel &,
                                                                  arrayView1d< string const > const & regionNames )
    {
      array1d< string > & stencilTargetRegions = fluxApprox.targetRegions( meshBodyName );
      std::set< string > stencilTargetRegionsSet( stencilTargetRegions.begin(), stencilTargetRegions.end() );
      stencilTargetRegionsSet.insert( regionNames.begin(), regionNames.end() );
      stencilTargetRegions.clear();
      for( auto const & targetRegion: stencilTargetRegionsSet )
      {
        stencilTargetRegions.emplace_back( targetRegion );
      }
    } );
  }
}

void FlowSolverBase::validatePoreVolumes( DomainPartition const & domain ) const
{
  real64 minPoreVolume = LvArray::NumericLimits< real64 >::max;
  real64 maxPorosity = -LvArray::NumericLimits< real64 >::max;
  globalIndex numElemsBelowPoreVolumeThreshold = 0;
  globalIndex numElemsAbovePorosityThreshold = 0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase const & subRegion )
    {

      string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      real64 minPoreVolumeInSubRegion = 0.0;
      real64 maxPorosityInSubRegion = 0.0;
      localIndex numElemsBelowPoreVolumeThresholdInSubRegion = 0;
      localIndex numElemsAbovePorosityThresholdInSubRegion = 0;

      flowSolverBaseKernels::MinPoreVolumeMaxPorosityKernel::
        computeMinPoreVolumeMaxPorosity( subRegion.size(),
                                         ghostRank,
                                         porosity,
                                         volume,
                                         minPoreVolumeInSubRegion,
                                         maxPorosityInSubRegion,
                                         numElemsBelowPoreVolumeThresholdInSubRegion,
                                         numElemsAbovePorosityThresholdInSubRegion );

      if( minPoreVolumeInSubRegion < minPoreVolume )
      {
        minPoreVolume = minPoreVolumeInSubRegion;
      }
      if( maxPorosityInSubRegion > maxPorosity )
      {
        maxPorosity = maxPorosityInSubRegion;
      }

      numElemsBelowPoreVolumeThreshold += numElemsBelowPoreVolumeThresholdInSubRegion;
      numElemsAbovePorosityThreshold += numElemsAbovePorosityThresholdInSubRegion;
    } );
  } );

  minPoreVolume = MpiWrapper::min( minPoreVolume );
  maxPorosity = MpiWrapper::max( maxPorosity );
  numElemsBelowPoreVolumeThreshold = MpiWrapper::sum( numElemsBelowPoreVolumeThreshold );
  numElemsAbovePorosityThreshold = MpiWrapper::sum( numElemsAbovePorosityThreshold );

  GEOS_LOG_RANK_0_IF( numElemsBelowPoreVolumeThreshold > 0,
                      GEOS_FMT( "\nWarning! The mesh contains {} elements with a pore volume below {} m^3."
                                "\nThe minimum pore volume is {} m^3."
                                "\nOur recommendation is to check the validity of mesh and/or increase the porosity in these elements.\n",
                                numElemsBelowPoreVolumeThreshold, flowSolverBaseKernels::poreVolumeThreshold, minPoreVolume ) );
  GEOS_LOG_RANK_0_IF( numElemsAbovePorosityThreshold > 0,
                      GEOS_FMT( "\nWarning! The mesh contains {} elements with a porosity above 1."
                                "\nThe maximum porosity is {}.\n",
                                numElemsAbovePorosityThreshold, maxPorosity ) );
}

void FlowSolverBase::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    precomputeData( mesh, regionNames );
  } );
}

void FlowSolverBase::precomputeData( MeshLevel & mesh,
                                     arrayView1d< string const > const & regionNames )
{
  FaceManager & faceManager = mesh.getFaceManager();
  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                        ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();

    arrayView1d< real64 > const gravityCoef =
      subRegion.getField< fields::flow::gravityCoefficient >();

    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      gravityCoef[ ei ] = LvArray::tensorOps::AiBi< 3 >( elemCenter[ ei ], gravVector );
    } );
  } );

  {
    arrayView2d< real64 const > const faceCenter = faceManager.faceCenter();

    arrayView1d< real64 > const gravityCoef =
      faceManager.getField< fields::flow::gravityCoefficient >();

    forAll< parallelHostPolicy >( faceManager.size(), [=] ( localIndex const kf )
    {
      gravityCoef[ kf ] = LvArray::tensorOps::AiBi< 3 >( faceCenter[ kf ], gravVector );
    } );
  }
}

void FlowSolverBase::updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const & pressure_n = subRegion.getField< fields::flow::pressure_n >();

  arrayView1d< real64 const > const & temperature = subRegion.getField< fields::flow::temperature >();
  arrayView1d< real64 const > const & temperature_n = subRegion.getField< fields::flow::temperature_n >();

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( solidName );

  constitutive::ConstitutivePassThru< CoupledSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
  {
    typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();
    if( m_isFixedStressPoromechanicsUpdate )
    {
      arrayView1d< real64 const > const & pressure_k = subRegion.getField< fields::flow::pressure_k >();
      arrayView1d< real64 const > const & temperature_k = subRegion.getField< fields::flow::temperature_k >();
      updatePorosityAndPermeabilityFromPressureAndTemperature( porousWrapper, subRegion, pressure, pressure_k, pressure_n, temperature, temperature_k, temperature_n );
    }
    else
    {
      updatePorosityAndPermeabilityFromPressureAndTemperature( porousWrapper, subRegion, pressure, pressure_n, pressure_n, temperature, temperature_n, temperature_n );
    }
  } );
}

void FlowSolverBase::updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const & pressure = subRegion.getField< fields::flow::pressure >();

  arrayView1d< real64 const > const newHydraulicAperture = subRegion.getField< fields::flow::hydraulicAperture >();
  arrayView1d< real64 const > const oldHydraulicAperture = subRegion.getField< fields::flow::aperture0 >();

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  CoupledSolidBase & porousSolid = subRegion.getConstitutiveModel< CoupledSolidBase >( solidName );

  constitutive::ConstitutivePassThru< CompressibleSolidBase >::execute( porousSolid, [=, &subRegion] ( auto & castedPorousSolid )
  {
    typename TYPEOFREF( castedPorousSolid ) ::KernelWrapper porousWrapper = castedPorousSolid.createKernelUpdates();

    updatePorosityAndPermeabilityFromPressureAndAperture( porousWrapper, subRegion, pressure, oldHydraulicAperture, newHydraulicAperture );

  } );
}


void FlowSolverBase::findMinMaxElevationInEquilibriumTarget( DomainPartition & domain, // cannot be const...
                                                             std::map< string, localIndex > const & equilNameToEquilId,
                                                             arrayView1d< real64 > const & maxElevation,
                                                             arrayView1d< real64 > const & minElevation ) const
{
  array1d< real64 > localMaxElevation( equilNameToEquilId.size() );
  array1d< real64 > localMinElevation( equilNameToEquilId.size() );
  localMaxElevation.setValues< parallelHostPolicy >( -1e99 );
  localMinElevation.setValues< parallelHostPolicy >( 1e99 );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply< ElementSubRegionBase,
                   EquilibriumInitialCondition >( 0.0,
                                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                                  EquilibriumInitialCondition::catalogName(),
                                                  [&] ( EquilibriumInitialCondition const & fs,
                                                        string const &,
                                                        SortedArrayView< localIndex const > const & targetSet,
                                                        ElementSubRegionBase & subRegion,
                                                        string const & )
  {
    RAJA::ReduceMax< parallelDeviceReduce, real64 > targetSetMaxElevation( -1e99 );
    RAJA::ReduceMin< parallelDeviceReduce, real64 > targetSetMinElevation( 1e99 );

    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();

    // TODO: move to FlowSolverBaseKernels to make this function "protected"
    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      targetSetMaxElevation.max( elemCenter[k][2] );
      targetSetMinElevation.min( elemCenter[k][2] );
    } );

    localIndex const equilIndex = equilNameToEquilId.at( fs.getName() );
    localMaxElevation[equilIndex] = LvArray::math::max( targetSetMaxElevation.get(), localMaxElevation[equilIndex] );
    localMinElevation[equilIndex] = LvArray::math::min( targetSetMinElevation.get(), localMinElevation[equilIndex] );

  } );

  MpiWrapper::allReduce( localMaxElevation.data(),
                         maxElevation.data(),
                         localMaxElevation.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Max ),
                         MPI_COMM_GEOS );
  MpiWrapper::allReduce( localMinElevation.data(),
                         minElevation.data(),
                         localMinElevation.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Min ),
                         MPI_COMM_GEOS );
}

void FlowSolverBase::computeSourceFluxSizeScalingFactor( real64 const & time,
                                                         real64 const & dt,
                                                         DomainPartition & domain, // cannot be const...
                                                         std::map< string, localIndex > const & bcNameToBcId,
                                                         arrayView1d< globalIndex > const & bcAllSetsSize ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&]( SourceFluxBoundaryCondition const & fs,
                                                         string const &,
                                                         SortedArrayView< localIndex const > const & targetSet,
                                                         ElementSubRegionBase & subRegion,
                                                         string const & )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // TODO: move to FlowSolverBaseKernels to make this function "protected"
      // loop over all the elements of this target set
      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > localSetSize( 0 );
      forAll< parallelDevicePolicy<> >( targetSet.size(),
                                        [targetSet, ghostRank, localSetSize] GEOS_HOST_DEVICE ( localIndex const k )
      {
        localIndex const ei = targetSet[k];
        if( ghostRank[ei] < 0 )
        {
          localSetSize += 1;
        }
      } );

      // increment the set size for this source flux boundary conditions
      bcAllSetsSize[bcNameToBcId.at( fs.getName())] += localSetSize.get();
    } );
  } );

  // synchronize the set size over all the MPI ranks
  MpiWrapper::allReduce( bcAllSetsSize.data(),
                         bcAllSetsSize.data(),
                         bcAllSetsSize.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOS );
}

void FlowSolverBase::saveAquiferConvergedState( real64 const & time,
                                                real64 const & dt,
                                                DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  // This step requires three passes:
  //    - First we count the number of individual aquifers
  //    - Second we loop over all the stencil entries to compute the sum of aquifer influxes
  //    - Third we loop over the aquifers to save the sums of each individual aquifer

  // Step 1: count individual aquifers

  std::map< string, localIndex > aquiferNameToAquiferId;
  localIndex aquiferCounter = 0;

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
  {
    aquiferNameToAquiferId[bc.getName()] = aquiferCounter;
    aquiferCounter++;
  } );

  // Step 2: sum the aquifer fluxes for each individual aquifer

  array1d< real64 > globalSumFluxes( aquiferNameToAquiferId.size() );
  array1d< real64 > localSumFluxes( aquiferNameToAquiferId.size() );

  fsManager.apply< FaceManager, AquiferBoundaryCondition >( time + dt,
                                                            mesh,
                                                            AquiferBoundaryCondition::catalogName(),
                                                            [&] ( AquiferBoundaryCondition const & bc,
                                                                  string const & setName,
                                                                  SortedArrayView< localIndex const > const &,
                                                                  FaceManager &,
                                                                  string const & )
  {
    BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
    if( stencil.size() == 0 )
    {
      return;
    }

    AquiferBoundaryCondition::KernelWrapper aquiferBCWrapper = bc.createKernelWrapper();

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > pressure =
      elemManager.constructFieldAccessor< fields::flow::pressure >();
    pressure.setName( getName() + "/accessors/" + fields::flow::pressure::key() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > pressure_n =
      elemManager.constructFieldAccessor< fields::flow::pressure_n >();
    pressure_n.setName( getName() + "/accessors/" + fields::flow::pressure_n::key() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > gravCoef =
      elemManager.constructFieldAccessor< fields::flow::gravityCoefficient >();
    gravCoef.setName( getName() + "/accessors/" + fields::flow::gravityCoefficient::key() );

    real64 const targetSetSumFluxes =
      fluxKernelsHelper::AquiferBCKernel::sumFluxes( stencil,
                                                     aquiferBCWrapper,
                                                     pressure.toNestedViewConst(),
                                                     pressure_n.toNestedViewConst(),
                                                     gravCoef.toNestedViewConst(),
                                                     time,
                                                     dt );

    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );
    localSumFluxes[aquiferIndex] += targetSetSumFluxes;
  } );

  MpiWrapper::allReduce( localSumFluxes.data(),
                         globalSumFluxes.data(),
                         localSumFluxes.size(),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOS );

  // Step 3: we are ready to save the summed fluxes for each individual aquifer

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition & bc )
  {
    localIndex const aquiferIndex = aquiferNameToAquiferId.at( bc.getName() );

    if( bc.getLogLevel() >= 1 )
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "{} {}: at time {} s, the boundary condition produces a volume of {} m3.",
                                 bc.getCatalogName(), bc.getName(),
                                 time + dt, dt * globalSumFluxes[aquiferIndex] ) );
    }
    bc.saveConvergedState( dt * globalSumFluxes[aquiferIndex] );
  } );
}

void FlowSolverBase::prepareStencilWeights( DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( getDiscretizationName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
      m_isLaggingFractureStencilWeightsUpdate ?
      mesh.getElemManager().constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::aperture0::key() ) :
      mesh.getElemManager().constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

    fluxApprox.forStencils< SurfaceElementStencil, FaceElementToCellStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
    {
      using STENCILWRAPPER_TYPE =  typename TYPEOFREF( stencil ) ::KernelWrapper;

      STENCILWRAPPER_TYPE stencilWrapper = stencil.createKernelWrapper();

      flowSolverBaseKernels::stencilWeightsUpdateKernel< STENCILWRAPPER_TYPE >::prepareStencilWeights( stencilWrapper, hydraulicAperture.toNestedViewConst() );
    } );
  } );
}

void FlowSolverBase::updateStencilWeights( DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( getDiscretizationName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > hydraulicAperture =
      mesh.getElemManager().constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( fields::flow::hydraulicAperture::key() );

    fluxApprox.forStencils< SurfaceElementStencil, FaceElementToCellStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
    {
      using STENCILWRAPPER_TYPE =  typename TYPEOFREF( stencil ) ::KernelWrapper;

      STENCILWRAPPER_TYPE stencilWrapper = stencil.createKernelWrapper();

      flowSolverBaseKernels::stencilWeightsUpdateKernel< STENCILWRAPPER_TYPE >::updateStencilWeights( stencilWrapper, hydraulicAperture.toNestedViewConst() );
    } );
  } );
}

bool FlowSolverBase::checkSequentialSolutionIncrements( DomainPartition & GEOS_UNUSED_PARAM( domain ) ) const
{

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Convergence, GEOS_FMT( "    {}: Max pressure change during outer iteration: {} Pa",
                                                              getName(), GEOS_FMT( "{:.{}f}", m_sequentialPresChange, 3 ) ) );

  if( m_isThermal )
  {
    GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Convergence, GEOS_FMT( "    {}: Max temperature change during outer iteration: {} K",
                                                                getName(), GEOS_FMT( "{:.{}f}", m_sequentialTempChange, 3 ) ) );
  }

  return (m_sequentialPresChange < m_maxSequentialPresChange) && (m_sequentialTempChange < m_maxSequentialTempChange);
}

string FlowSolverBase::BCMessage::generateMessage( string_view baseMessage,
                                                   string_view fieldName, string_view setName )
{
  return GEOS_FMT( "{} Check if you have added or applied the appropriate fields to "
                   "the FieldSpecification component with fieldName=\"{}\" "
                   "and setNames=\"{}\"\n", baseMessage, fieldName, setName );
}

string FlowSolverBase::BCMessage::pressureConflict( string_view regionName, string_view subRegionName,
                                                    string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Conflicting pressure boundary conditions on set {}/{}/{}\n",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::temperatureConflict( string_view regionName, string_view subRegionName,
                                                       string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Conflicting temperature boundary conditions on set {}/{}/{}\n",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::missingPressure( string_view regionName, string_view subRegionName,
                                                   string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Pressure boundary condition not prescribed on set {}/{}/{}\n",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::missingTemperature( string_view regionName, string_view subRegionName,
                                                      string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Temperature boundary condition not prescribed on set {}/{}/{}\n",
                                    regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::conflictingComposition( int comp,
                                                          string_view regionName, string_view subRegionName,
                                                          string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Conflicting composition[{}] boundary conditions on set {}/{}/{}\n",
                                    comp, regionName, subRegionName, setName ), fieldName, setName );
}

string FlowSolverBase::BCMessage::invalidComponentIndex( int comp,
                                                         string_view fsName,
                                                         string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Invalid component index [{}] in composition boundary condition {}\n",
                                    comp, fsName ), fieldName, fsName );
}

string FlowSolverBase::BCMessage::inconsistency( int componentIndex, string_view componentName,
                                                 string_view regionName, string_view subRegionName,
                                                 string_view setName, string_view fieldName )
{
  return generateMessage( GEOS_FMT( "Boundary condition not applied to composition[{}] ({}) "
                                    "on region {}/{}/{}\n",
                                    componentIndex, componentName, regionName, subRegionName, setName ),
                          fieldName, setName );
}

} // namespace geos
