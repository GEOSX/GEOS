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
 * @file ImmiscibleMultiphaseFlow.cpp
 */

#include "ImmiscibleMultiphaseFlow.hpp"

#include "FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseKernels.hpp"
#include "physicsSolvers/SolverBaseKernels.hpp"

#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"  // should be removed eventually
//#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"

#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"


#include "constitutive/ConstitutivePassThru.hpp"

#include "constitutive/fluid/twophasefluid/TwoPhaseFluid.hpp"

#include <cmath>

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields::immiscibleMultiphaseFlow;
using namespace immiscibleMultiphaseKernels;


ImmiscibleMultiphaseFlow::ImmiscibleMultiphaseFlow( const string & name,
                                                    Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 2 ),
  m_hasCapPressure( 0 ),
  m_useTotalMassEquation ( 0 )
{
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Temperature" );
}

void ImmiscibleMultiphaseFlow::postInputInitialization()
{
  FlowSolverBase::postInputInitialization();
}

void ImmiscibleMultiphaseFlow::registerDataOnMesh( Group & meshBodies )
{
  FlowSolverBase::registerDataOnMesh( meshBodies );

  // DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  // ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 0. Find a "reference" fluid model name (at this point, models are already attached to subregions)
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // If at least one region has a capillary pressure model, consider it enabled for all
      string const capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
      if( !capPresName.empty() )
      {
        m_hasCapPressure = true;
      }
    } );
  } );

  m_numDofPerCell = m_numPhases;

  // 2. Register and resize all fields as necessary
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_hasCapPressure )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::capPressureNamesString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 ).
          setDescription( "Name of the capillary pressure constitutive model to use" ).
          reference();

        string & capPresName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
        GEOS_THROW_IF( capPresName.empty(),
                       GEOS_FMT( "{}: Capillary pressure model not found on subregion {}",
                                 getDataContext(), subRegion.getDataContext() ),
                       InputError );
      }

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.
      subRegion.registerField< phaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< bcPhaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMass >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMass_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMobility >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< dPhaseMobility >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numPhases ); // dP, dS

    } );

    // FaceManager & faceManager = mesh.getFaceManager();
    // {
    //   faceManager.registerField< totalMassFlux >( getName() );
    // }
  } );
}

void ImmiscibleMultiphaseFlow::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{

  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< TwoPhaseFluid >( subRegion );
  GEOS_ERROR_IF( fluidName.empty(), GEOS_FMT( "{}: Fluid model not found on subregion {}",
                                              getDataContext(), subRegion.getName() ) );



  string & relPermName = subRegion.registerWrapper< string >( viewKeyStruct::relPermNamesString() ).
                           setPlotLevel( PlotLevel::NOPLOT ).
                           setRestartFlags( RestartFlags::NO_WRITE ).
                           setSizedFromParent( 0 ).
                           setDescription( "Name of the relative permeability constitutive model to use" ).
                           reference();

  relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );

  GEOS_THROW_IF( relPermName.empty(),
                 GEOS_FMT( "{}: Relative permeability model not found on subregion {}",
                           getDataContext(), subRegion.getDataContext() ),
                 InputError );
}

void ImmiscibleMultiphaseFlow::initializePreSubGroups()
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::immiscibleMultiphaseFVM;

  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  // ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // // 1. Validate various models against each other (must have same phases and components)
  // validateConstitutiveModels( domain );

  // 2. Set the value of temperature
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getField< fields::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );
}


void ImmiscibleMultiphaseFlow::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();

  TwoPhaseFluid & fluid = getConstitutiveModel< TwoPhaseFluid >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    FluidUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(), fluidWrapper, pres );
  } );
}



void ImmiscibleMultiphaseFlow::updateRelPermModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;


  GEOS_UNUSED_VAR( dataGroup );

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
    dataGroup.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

  string const & relPermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase & relPerm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relPermName );

  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

    isothermalCompositionalMultiphaseBaseKernels::
      RelativePermeabilityUpdateKernel::
      launch< parallelDevicePolicy<> >( dataGroup.size(),
                                        relPermWrapper,
                                        phaseVolFrac );
  } );
}

void ImmiscibleMultiphaseFlow::updateCapPressureModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  if( m_hasCapPressure )
  {
    arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
      dataGroup.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

    string const & cappresName = dataGroup.getReference< string >( viewKeyStruct::capPressureNamesString() );
    CapillaryPressureBase & capPressure = getConstitutiveModel< CapillaryPressureBase >( dataGroup, cappresName );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      isothermalCompositionalMultiphaseBaseKernels::
        CapillaryPressureUpdateKernel::
        launch< parallelDevicePolicy<> >( dataGroup.size(),
                                          capPresWrapper,
                                          phaseVolFrac );
    } );
  }
}


void ImmiscibleMultiphaseFlow::updateFluidState( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  updateFluidModel( subRegion );
  updatePhaseVolumeFraction( subRegion );
  updatePhaseMass( subRegion );
  updateRelPermModel( subRegion );
  updatePhaseMobility( subRegion );
  updateCapPressureModel( subRegion );
}


void ImmiscibleMultiphaseFlow::updatePhaseMass( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );

  TwoPhaseFluid const & fluid = getConstitutiveModel< TwoPhaseFluid >( subRegion, fluidName );
  CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView2d< real64 const > const porosity = solid.getPorosity();

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac= subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDens = fluid.phaseDensity();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseMass = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass >();

  // Might be needed for geomechanics????? if so, need to change the accumulation as well?
  //arrayView1d< real64 > const deltaVolume = subRegion.getField< fields::flow::deltaVolume >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const poreVolume = volume[ei] * porosity[ei][0];
    for( integer ip = 0; ip < 2; ++ip )
    {
      phaseMass[ei][ip] = poreVolume * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip];
    }
  } );
}


void ImmiscibleMultiphaseFlow::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  // note that the phase mobility computed here also includes phase density
  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  TwoPhaseFluid const & fluid = getConstitutiveModel< TwoPhaseFluid >( dataGroup, fluidName );

  string const & relpermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relpermName );

  immiscibleMultiphaseKernels::
    PhaseMobilityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                               dataGroup,
                                               fluid,
                                               relperm );
}

void ImmiscibleMultiphaseFlow::initializeFluidState( MeshLevel & mesh,
                                                     DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                     arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;

  mesh.getElemManager().forElementSubRegions( regionNames,
                                              [&]( localIndex const,
                                                   ElementSubRegionBase & subRegion )
  {
    // 2. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    updateFluidModel( subRegion );

  } );

  // for some reason CUDA does not want the host_device lambda to be defined inside the generic lambda
  // I need the exact type of the subRegion for updateSolidflowProperties to work well.
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                              SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                           auto & subRegion )
  {
    // 4. Initialize/update dependent state quantities

    // 4.1 Update the constitutive models that only depend on
    //      - the primary variables
    //      - the fluid constitutive quantities (as they have already been updated)
    // We postpone the other constitutive models for now
    // In addition, to avoid multiplying permeability/porosity bay netToGross in the assembly kernel, we do it once and for all here
    arrayView1d< real64 const > const netToGross = subRegion.template getField< fields::flow::netToGross >();
    CoupledSolidBase const & porousSolid =
      getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
    PermeabilityBase const & permeabilityModel =
      getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() ) );
    permeabilityModel.scaleHorizontalPermeability( netToGross );
    porousSolid.scaleReferencePorosity( netToGross );
    saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
    updatePorosityAndPermeability( subRegion );

    // Now, we initialize and update each constitutive model one by one

    // 4.2 Save the computed porosity into the old porosity
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity
    string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
    porousMaterial.initializeState();

    // 4.3 Initialize/update the relative permeability model using the initial phase volume fraction
    //     This is needed to handle relative permeability hysteresis
    //     Also, initialize the fluid model
    //
    // Note:
    // - This must be called after updatePhaseVolumeFraction
    // - This step depends on phaseVolFraction

    // initialized phase volume fraction
    arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
      subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

    string const & relpermName = subRegion.template getReference< string >( viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase & relPermMaterial =
      getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
    relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac ); // this needs to happen before calling updateRelPermModel
    updateRelPermModel( subRegion );
    relPermMaterial.saveConvergedState(); // this needs to happen after calling updateRelPermModel

    // 4.4 Then, we initialize/update the capillary pressure model
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity and permeability
    if( m_hasCapPressure )
    {
      // initialized porosity
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      string const & permName = subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeabilityMaterial = getConstitutiveModel< PermeabilityBase >( subRegion, permName );
      // initialized permeability
      arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

      string const & capPressureName = subRegion.template getReference< string >( viewKeyStruct::capPressureNamesString() );
      CapillaryPressureBase const & capPressureMaterial =
        getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressureName );
      capPressureMaterial.initializeRockState( porosity, permeability ); // this needs to happen before calling updateCapPressureModel
      updateCapPressureModel( subRegion );
    }

    // 4.5 Update the phase mobility
    //
    // Note:
    // - This must be called after updateRelPermModel
    // - This step depends phaseRelPerm
    updatePhaseMobility( subRegion );

  } );

  // 5. Save initial pressure
  mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 > const initPres = subRegion.getField< fields::flow::initialPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
    arrayView1d< real64 > const initTemp = subRegion.template getField< fields::flow::initialTemperature >();
    initPres.setValues< parallelDevicePolicy<> >( pres );
    initTemp.setValues< parallelDevicePolicy<> >( temp );

    // TODO: Missing updatePhaseMass?
  } );
}

void ImmiscibleMultiphaseFlow::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::pressure::key(),
                                       fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames,
                                                                                                 [&]( localIndex const,
                                                                                                      auto & subRegion )
    {
      saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
      updatePorosityAndPermeability( subRegion );

      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion,
                                                  subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      porousSolid.initializeState();
    } );

    // Initialize primary variables from applied initial conditions
    initializeFluidState( mesh, domain, regionNames );

    mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames,
                                                                     [&]( localIndex const,
                                                                          SurfaceElementRegion & region )
    {
      region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
      {
        subRegion.getWrapper< real64_array >( fields::flow::hydraulicAperture::key() ).
          setApplyDefaultValue( region.getDefaultAperture() );
      } );
    } );

  } );

  // report to the user if some pore volumes are very small
  // note: this function is here because: 1) porosity has been initialized and 2) NTG has been applied
  validatePoreVolumes( domain );
}

void
ImmiscibleMultiphaseFlow::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                             real64 const & GEOS_UNUSED_PARAM( dt ),
                                             DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      saveConvergedState( subRegion );

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updatePhaseVolumeFraction( subRegion );
      updateFluidState( subRegion );

      // after the update, save the new saturation
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
      phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );

      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const & phaseMass =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const & phaseMass_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();
      phaseMass_n.setValues< parallelDevicePolicy<> >( phaseMass );

    } );
  } );
}

void ImmiscibleMultiphaseFlow::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                               real64 const dt,
                                               DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationTerm( domain,
                            dofManager,
                            localMatrix,
                            localRhs );


  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}



void ImmiscibleMultiphaseFlow::assembleAccumulationTerm( DomainPartition & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );

      TwoPhaseFluid const & fluid = getConstitutiveModel< TwoPhaseFluid >( subRegion, fluidName );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      //arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      integer const numofPhases = 2;
      globalIndex const rankOffset = dofManager.rankOffset();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const elemGhostRank= subRegion.ghostRank();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView2d< real64 const > const porosity = solid.getPorosity();
      arrayView2d< real64 const > const dPoro_dPres = solid.getDporosity_dPressure();
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac= subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      arrayView3d< real64 const, multifluid::USD_PHASE > phaseDens = fluid.phaseDensity();
      arrayView4d< real64 const, multifluid::USD_PHASE_DC > dPhaseDens = fluid.dPhaseDensity();
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const PhaseMass_n = subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();

      // The line below is to be used if we want to use the total mass flux formulation
      //BitFlags< ElementBasedAssemblyKernelFlags > m_kernelFlags = kernelFlags;


      integer const numElems = subRegion.size();
      forAll< parallelDevicePolicy<> >(
        numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank( ei ) >= 0 )
        {
          return;
        }

        // setup
        globalIndex dofIndices[2]{};
        real64 localResidual[2]{};
        real64 localJacobian[2][2]{};

        real64 const poreVolume = volume[ei] * porosity[ei][0];
        real64 const dPoreVolume_dPres = volume[ei] * dPoro_dPres[ei][0];

        localIndex localRow = dofNumber[ei] - rankOffset;

        for( integer idof = 0; idof < 2; ++idof )
        {
          dofIndices[idof] = dofNumber[ei] + idof;
        }

        // compute accumulation
        int signPotDiff[2] = {1, -1};

        for( integer ip = 0; ip < numofPhases; ++ip )
        {
          real64 const phaseMass = poreVolume * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip];
          real64 const phaseMass_n = PhaseMass_n[ei][ip];

          localResidual[ip] += phaseMass - phaseMass_n;

          real64 const dPhaseMass_dP =  dPoreVolume_dPres * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip]
                                       + poreVolume * phaseVolFrac[ei][ip] * dPhaseDens[ei][0][ip][0];
          localJacobian[ip][0] += dPhaseMass_dP;

          real64 const dPhaseMass_dS = poreVolume * phaseDens[ei][0][ip];


          // if ( ip == 0)
          // {
          //   localJacobian[ip][1] += dPhaseMass_dS;
          // } else {
          //   localJacobian[ip][1] -= dPhaseMass_dS;
          // }

          localJacobian[ip][1] += signPotDiff[ip] * dPhaseMass_dS;

        }

        // complete

        integer const numRows = numofPhases;

        for( integer i = 0; i < numRows; ++i )
        {
          localRhs[localRow + i] += localResidual[i];
          localMatrix.addToRow< serialAtomic >( localRow + i,
                                                dofIndices,
                                                localJacobian[i],
                                                2 );
        }

      } );

    } );
  } );


}

void ImmiscibleMultiphaseFlow::assembleFluxTerms( real64 const dt,
                                                  DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
      immiscibleMultiphaseKernels::
        FaceBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                                                   dofManager.rankOffset(),
                                                                                   dofKey,
                                                                                   m_hasCapPressure,
                                                                                   m_useTotalMassEquation,
                                                                                   getName(),
                                                                                   mesh.getElemManager(),
                                                                                   stencilWrapper,
                                                                                   dt,
                                                                                   localMatrix.toViewConstSizes(),
                                                                                   localRhs.toView() );
    } );
  } );
}

// Ryan: Looks like this will need to be overwritten as well...
// I have left the CompositionalMultiphaseFVM implementation for reference
void ImmiscibleMultiphaseFlow::setupDofs( DomainPartition const & domain,
                                          DofManager & dofManager ) const
{
  GEOS_UNUSED_VAR( domain, dofManager );
  // add a field for the cell-centered degrees of freedom
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  //// this call with instruct GEOS to reorder the dof numbers
  //dofManager.setLocalReorderingType( viewKeyStruct::elemDofFieldString(),
  //                                   DofManager::LocalReorderingType::ReverseCutHillMcKee );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}

void ImmiscibleMultiphaseFlow::applyBoundaryConditions( real64 const time_n,
                                                        real64 const dt,
                                                        DomainPartition & domain,
                                                        DofManager const & dofManager,
                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                        arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // apply pressure boundary conditions.
  applyDirichletBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

// apply flux boundary conditions
  applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

  //   for( localIndex row = 0; row < localMatrix.toViewConstSizes().numRows(); ++row )
  // {
  //   std::cout << "row " << row << std::endl;
  //   std::cout << "\tcolumns: " << localMatrix.toViewConstSizes().getColumns( row ) << std::endl;
  //   std::cout << "\tvalues: " << localMatrix.toViewConstSizes().getEntries( row ) << std::endl;
  // }

}

namespace
{
char const bcLogMessage[] =
  "CompositionalMultiphaseBase {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";
}


void ImmiscibleMultiphaseFlow::applyDirichletBC( real64 const time_n,
                                                 real64 const dt,
                                                 DofManager const & dofManager,
                                                 DomainPartition & domain,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  // // Only validate BC at the beginning of Newton loop
  // if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  // {
  //   bool const bcConsistent = validateDirichletBC( domain, time_n + dt );
  //   GEOS_ERROR_IF( !bcConsistent, GEOS_FMT( "CompositionalMultiphaseBase {}: inconsistent boundary conditions", getDataContext() ) );
  // }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {

    // 1. Apply pressure Dirichlet BCs, store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::pressure::key(), fields::flow::bcPressure::key() );
    // 2. Apply saturation BC (phase volume fraction) and store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key(),
                                             fields::immiscibleMultiphaseFlow::bcPhaseVolumeFraction::key() );

    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    // 3. Call constitutive update, back-calculate target global component densities and apply to the system
    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {

      arrayView1d< real64 const > const bcPres =
        subRegion.getReference< array1d< real64 > >( fields::flow::bcPressure::key() );
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const bcPhaseVolFraction =
        subRegion.getReference< array2d< real64, immiscibleFlow::LAYOUT_PHASE > >(
          fields::immiscibleMultiphaseFlow::bcPhaseVolumeFraction::key() );

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< real64 const > const pres =
        subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFraction =
        subRegion.getReference< array2d< real64, immiscibleFlow::LAYOUT_PHASE > >(
          fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() );

      integer const numPhase = m_numPhases;


      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 3.1. Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    bcPres[ei],
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 3.2. For each phase, apply target saturation value
        for( integer ip = 0; ip < numPhase-1; ++ip )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ip + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcPhaseVolFraction[ei][ip],
                                                      phaseVolFraction[ei][ip] );
          localRhs[localRow + ip + 1] = rhsValue;
        }
      } );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::applySourceFluxBC( real64 const time,
                                                  real64 const dt,
                                                  DofManager const & dofManager,
                                                  DomainPartition & domain,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Step 1: count individual source flux boundary conditions

  std::map< string, localIndex > bcNameToBcId;
  localIndex bcCounter = 0;

  fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&] ( SourceFluxBoundaryCondition const & bc )
  {
    // collect all the bc names to idx
    bcNameToBcId[bc.getName()] = bcCounter;
    bcCounter++;
  } );

  if( bcCounter == 0 )
  {
    return;
  }

  // Step 2: count the set size for each source flux (each source flux may have multiple target sets)

  array1d< globalIndex > bcAllSetsSize( bcNameToBcId.size() );

  computeSourceFluxSizeScalingFactor( time,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {

    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&]( SourceFluxBoundaryCondition const & fs,
                                                         string const & setName,
                                                         SortedArrayView< localIndex const > const & targetSet,
                                                         ElementSubRegionBase & subRegion,
                                                         string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
                                   getName(), time+dt, fs.getCatalogName(), fs.getName(),
                                   setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      if( targetSet.size() == 0 )
      {
        return;
      }
      if( !subRegion.hasWrapper( dofKey ) )
      {
        if( fs.getLogLevel() >= 1 )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: trying to apply SourceFlux, but its targetSet named '{}' intersects with non-simulated region named '{}'.",
                                   getDataContext(), setName, subRegion.getName() ) );
        }
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // Step 3.1: get the values of the source boundary condition that need to be added to the rhs

      array1d< globalIndex > dofArray( targetSet.size() );
      array1d< real64 > rhsContributionArray( targetSet.size() );
      arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
      localIndex const rankOffset = dofManager.rankOffset();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd( 0.0 );

      // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
      fs.computeRhsContribution< FieldSpecificationAdd,
                                 parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                           time + dt,
                                                           dt,
                                                           subRegion,
                                                           dofNumber,
                                                           rankOffset,
                                                           localMatrix,
                                                           dofArray.toView(),
                                                           rhsContributionArrayView,
                                                           [] GEOS_HOST_DEVICE ( localIndex const )
      {
        return 0.0;
      } );

      // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];
      integer const fluidPhaseId = fs.getComponent();
      integer const numFluidPhases = m_numPhases;
      integer useTotalMassEquation = 0;
      forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                           targetSet,
                                                           rankOffset,
                                                           ghostRank,
                                                           fluidPhaseId,
                                                           numFluidPhases,
                                                           useTotalMassEquation,
                                                           dofNumber,
                                                           rhsContributionArrayView,
                                                           localRhs,
                                                           massProd] GEOS_HOST_DEVICE ( localIndex const a )
      {
        // we need to filter out ghosts here, because targetSet may contain them
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor here!
        massProd += rhsValue;
        if( useTotalMassEquation > 0 )
        {
          // for all "fluid components", we add the value to the total mass balance equation
          globalIndex const totalMassBalanceRow = dofNumber[ei] - rankOffset;
          localRhs[totalMassBalanceRow] += rhsValue;
          if( fluidPhaseId < numFluidPhases - 1 )
          {
            globalIndex const compMassBalanceRow = totalMassBalanceRow + fluidPhaseId + 1; // component mass bal equations are shifted
            localRhs[compMassBalanceRow] += rhsValue;
          }
        }
        else
        {
          globalIndex const compMassBalanceRow = dofNumber[ei] - rankOffset + fluidPhaseId;
          localRhs[compMassBalanceRow] += rhsValue;
        }
      } );

      // SourceFluxStatsAggregator::forAllFluxStatWrappers( subRegion, fs.getName(),
      //                                                    [&]( SourceFluxStatsAggregator::WrappedStats & wrapper )
      // {
      //   // set the new sub-region statistics for this timestep
      //   array1d< real64 > massProdArr{ m_numPhases };
      //   massProdArr[fluidPhaseId] = massProd.get();
      //   wrapper.gatherTimeStepStats( time, dt, massProdArr.toViewConst(), targetSet.size() );
      // } );
    } );
  } );
}


real64 ImmiscibleMultiphaseFlow::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                        real64 const & GEOS_UNUSED_PARAM( dt ),
                                                        DomainPartition const & domain,
                                                        DofManager const & dofManager,
                                                        arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  solverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[numNorm]{};
      real64 subRegionResidualNormalizer[numNorm]{};

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // step 1: compute the norm in the subRegion

      real64 subRegionFlowResidualNorm[1]{};
      real64 subRegionFlowResidualNormalizer[1]{};

      immiscibleMultiphaseKernels::
        ResidualNormKernelFactory::createAndLaunch< parallelDevicePolicy<> >( normType,
                                                                              2,
                                                                              rankOffset,
                                                                              dofKey,
                                                                              localRhs,
                                                                              subRegion,
                                                                              solid,
                                                                              m_nonlinearSolverParameters.m_minNormalizer,
                                                                              subRegionFlowResidualNorm,
                                                                              subRegionFlowResidualNormalizer );
      subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
      subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];

      // step 2: first reduction across meshBodies/regions/subRegions
      if( normType == solverBaseKernels::NormType::Linf )
      {
        solverBaseKernels::LinfResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, localResidualNorm );
      }
      else
      {
        solverBaseKernels::L2ResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, subRegionResidualNormalizer, localResidualNorm, localResidualNormalizer );
      }
    } );
  } );

  real64 residualNorm = 0.0;
  residualNorm = localResidualNorm[0];
  if( normType == solverBaseKernels::NormType::Linf )
  {
    solverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm[0], residualNorm );
  }
  else
  {
    solverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm[0], localResidualNormalizer[0], residualNorm );
  }

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residualNorm );
  }


  return residualNorm;
}

void ImmiscibleMultiphaseFlow::applySystemSolution( DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localSolution,
                                                    real64 const scalingFactor,
                                                    real64 const dt,
                                                    DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );

  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
  //DofManager::CompMask PhaseMask( m_numDofPerCell, 1, m_numPhases );

  // 1. apply the pressure update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               fields::flow::pressure::key(),
                               scalingFactor,
                               pressureMask );

  // 2. apply the phaseVolumeFraction update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key(),
                               scalingFactor,
                               ~pressureMask );

  // 3. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    std::vector< string > fields{ fields::flow::pressure::key(), fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() };

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );

}


void ImmiscibleMultiphaseFlow::updatePhaseVolumeFraction( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolumeFraction = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    phaseVolumeFraction[ei][1] = 1.0 - phaseVolumeFraction[ei][0];
  } );
}


void ImmiscibleMultiphaseFlow::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      arrayView1d< real64 > const & pres =
        subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const & pres_n =
        subRegion.template getField< fields::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      // after the update, save the new saturation
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      phaseVolFrac.setValues< parallelDevicePolicy<> >( phaseVolFrac_n );

      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const & phaseMass_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const & phaseMass =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass >();
      phaseMass.setValues< parallelDevicePolicy<> >( phaseMass_n );

      if( m_isThermal )
      {
        arrayView1d< real64 > const & temp =
          subRegion.template getField< fields::flow::temperature >();
        arrayView1d< real64 const > const & temp_n =
          subRegion.template getField< fields::flow::temperature_n >();
        temp.setValues< parallelDevicePolicy<> >( temp_n );
      }

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateFluidState( subRegion );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::implicitStepComplete( real64 const & time,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  // Step 1: save the converged aquifer state
  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // Step 3: save the converged solid state
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      if( m_keepFlowVariablesConstantDuringInitStep )
      {
        porousMaterial.ignoreConvergedState(); // newPorosity <- porosity_n
      }
      else
      {
        porousMaterial.saveConvergedState(); // porosity_n <- porosity
      }

      // Step 4: save converged state for the relperm model to handle hysteresis
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      string const & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPermMaterial =
        getConstitutiveModel< RelativePermeabilityBase >( subRegion, relPermName );
      relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac );

      // Step 5: if capillary pressure is supported, send the converged porosity and permeability to the capillary pressure model
      // note: this is needed when the capillary pressure depends on porosity and permeability (Leverett J-function for instance)
      if( m_hasCapPressure )
      {
        arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

        string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
        PermeabilityBase const & permeabilityMaterial =
          getConstitutiveModel< PermeabilityBase >( subRegion, permName );
        arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

        string const & capPressName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        CapillaryPressureBase const & capPressureMaterial =
          getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressName );
        capPressureMaterial.saveConvergedRockState( porosity, permeability );
      }
    } );
  } );
}

void ImmiscibleMultiphaseFlow::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::saveConvergedState( subRegion );

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
  phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseMass =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass >();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseMass_n =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();
  phaseMass_n.setValues< parallelDevicePolicy<> >( phaseMass );

}

void ImmiscibleMultiphaseFlow::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             auto & subRegion )
    {
      // update porosity, permeability, and solid internal energy
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updatePhaseVolumeFraction( subRegion );
      updateFluidState( subRegion );
    } );
  } );
}

real64 ImmiscibleMultiphaseFlow::setNextDt( const geos::real64 & currentDt, geos::DomainPartition & domain )
{
  return SolverBase::setNextDt( currentDt, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, ImmiscibleMultiphaseFlow, string const &, Group * const )

} // namespace geos
