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
 * @file ProppantTransport.cpp
 */

#include "ProppantTransport.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidSelector.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidFields.hpp"
#include "constitutive/fluid/singlefluid/ParticleFluidSelector.hpp"
#include "constitutive/fluid/singlefluid/ParticleFluidFields.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/permeability/ProppantPermeability.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportFields.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportKernels.hpp"
#include "mesh/MeshFields.hpp"


/**
 * @namespace the geos namespace that encapsulates the majority of the code
 */
namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace proppantTransportKernels;

ProppantTransport::ProppantTransport( const string & name,
                                      Group * const parent ):
  FlowSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::bridgingFactorString(), &m_bridgingFactor ).setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bridging factor used for bridging/screen-out calculation" );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).setApplyDefaultValue( 0.6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum proppant concentration" );

  registerWrapper( viewKeyStruct::proppantDiameterString(), &m_proppantDiameter ).setApplyDefaultValue( 0.4e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant diameter" );

  registerWrapper( viewKeyStruct::proppantDensityString(), &m_proppantDensity ).setApplyDefaultValue( 2500.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant density" );

  registerWrapper( viewKeyStruct::criticalShieldsNumberString(), &m_criticalShieldsNumber ).setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Critical Shields number" );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).setApplyDefaultValue( 0.03 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Friction coefficient" );

  registerWrapper( viewKeyStruct::updateProppantPackingString(), &m_updateProppantPacking ).setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag that enables/disables proppant-packing update" );

}

void ProppantTransport::postInputInitialization()
{
  FlowSolverBase::postInputInitialization();
}

void ProppantTransport::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::proppant;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             CellElementSubRegion & subRegion )
    {
      subRegion.registerField< proppantConcentration >( getName() );
      subRegion.registerField< proppantConcentration_n >( getName() );
      subRegion.registerField< componentConcentration >( getName() );
      subRegion.registerField< componentConcentration_n >( getName() );
      subRegion.registerField< bcComponentConcentration >( getName() );
      subRegion.registerField< cellBasedFlux >( getName() ).
        reference().resizeDimension< 1 >( 3 );
      subRegion.registerField< isProppantBoundary >( getName() );

      setConstitutiveNames( subRegion );
    } );

    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion & subRegion )
    {
      subRegion.registerField< proppantConcentration >( getName() );
      subRegion.registerField< proppantConcentration_n >( getName() );
      subRegion.registerField< componentConcentration >( getName() );
      subRegion.registerField< componentConcentration_n >( getName() );
      subRegion.registerField< bcComponentConcentration >( getName() );
      subRegion.registerField< componentDensity_n >( getName() );
      subRegion.registerField< cellBasedFlux >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< isProppantBoundary >( getName() );
      subRegion.registerField< isProppantMobile >( getName() );
      subRegion.registerField< proppantPackVolumeFraction >( getName() );
      subRegion.registerField< proppantExcessPackVolume >( getName() );
      subRegion.registerField< proppantLiftFlux >( getName() );

      setConstitutiveNames( subRegion );

    } );
  } );
}


void ProppantTransport::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< SlurryFluidBase >( subRegion );
  GEOS_THROW_IF( fluidName.empty(),
                 GEOS_FMT( "{}: Fluid model not found on subregion {}",
                           getDataContext(), subRegion.getName() ),
                 InputError );

  subRegion.registerWrapper< string >( viewKeyStruct::proppantNamesString() );
  string & proppantName = subRegion.getReference< string >( viewKeyStruct::proppantNamesString() );
  proppantName = getConstitutiveName< ParticleFluidBase >( subRegion );
  GEOS_THROW_IF( proppantName.empty(),
                 GEOS_FMT( "{}: Proppant model not found on subregion {}",
                           getDataContext(), subRegion.getName() ),
                 InputError );

}


void ProppantTransport::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager & cm = domain.getConstitutiveManager();

  // Validate proppant models in regions
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )

    {
      if( m_numDofPerCell < 1 )
      {
        SlurryFluidBase const & fluid0 = cm.getConstitutiveRelation< SlurryFluidBase >( subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        m_numComponents = fluid0.numFluidComponents();

        m_numDofPerCell = m_numComponents + 1;
      }
    } );

    if( m_numComponents > 0 )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            CellElementSubRegion & subRegion )

      {
        subRegion.getField< fields::proppant::componentConcentration >().resizeDimension< 1 >( m_numComponents );
        subRegion.getField< fields::proppant::componentConcentration_n >().resizeDimension< 1 >( m_numComponents );
      } );
    }
  } );
}

void ProppantTransport::resizeFractureFields( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  if( m_numComponents > 0 )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion & subRegion )
    {
      subRegion.getField< fields::proppant::componentConcentration >().resizeDimension< 1 >( m_numComponents );
      subRegion.getField< fields::proppant::componentConcentration_n >().resizeDimension< 1 >( m_numComponents );
      subRegion.getField< fields::proppant::componentDensity_n >().resizeDimension< 1 >( m_numComponents );
      subRegion.getField< fields::proppant::bcComponentConcentration >().resizeDimension< 1 >( m_numComponents );
    } );
  }
}

void ProppantTransport::updateFluidModel( ObjectManagerBase & dataGroup )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres  = dataGroup.getField< fields::flow::pressure >();

  arrayView2d< real64 const > const componentConc  = dataGroup.getField< fields::proppant::componentConcentration >();

  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    FluidUpdateKernel::launch( fluidWrapper,
                               pres,
                               componentConc );
  } );
}

void ProppantTransport::updateComponentDensity( ObjectManagerBase & dataGroup )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres  = dataGroup.getField< fields::flow::pressure >();

  arrayView2d< real64 const > const componentConc  = dataGroup.getField< fields::proppant::componentConcentration >();

  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    ComponentDensityUpdateKernel::launch( fluidWrapper,
                                          pres,
                                          componentConc );
  } );
}


void ProppantTransport::updateProppantModel( ObjectManagerBase & dataGroup )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const proppantConc  = dataGroup.getField< fields::proppant::proppantConcentration >();

  SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  arrayView2d< real64 const > const fluidDens            = fluid.fluidDensity();
  arrayView2d< real64 const > const dFluidDens_dPres     = fluid.dFluidDensity_dPressure();
  arrayView3d< real64 const > const dFluidDens_dCompConc = fluid.dFluidDensity_dComponentConcentration();
  arrayView2d< real64 const > const fluidVisc            = fluid.fluidViscosity();
  arrayView2d< real64 const > const dFluidVisc_dPres     = fluid.dFluidViscosity_dPressure();
  arrayView3d< real64 const > const dFluidVisc_dCompConc = fluid.dFluidViscosity_dComponentConcentration();

  ParticleFluidBase & proppant = getConstitutiveModel< ParticleFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::proppantNamesString() ) );

  constitutiveUpdatePassThru( proppant, [&]( auto & castedProppant )
  {
    typename TYPEOFREF( castedProppant ) ::KernelWrapper proppantWrapper = castedProppant.createKernelWrapper();
    ProppantUpdateKernel::launch( proppantWrapper,
                                  proppantConc,
                                  fluidDens,
                                  dFluidDens_dPres,
                                  dFluidDens_dCompConc,
                                  fluidVisc,
                                  dFluidVisc_dPres,
                                  dFluidVisc_dCompConc );
  } );
}

void ProppantTransport::updateProppantMobility( ObjectManagerBase & dataGroup )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const conc = dataGroup.getField< fields::proppant::proppantConcentration >();
  arrayView1d< real64 const > const aperture = dataGroup.getReference< array1d< real64 > >( fields::elementAperture::key() );
  arrayView1d< integer > const isProppantMobile = dataGroup.getField< fields::proppant::isProppantMobile >();

  real64 const minAperture = m_minAperture;
  real64 const maxProppantConcentration = m_maxProppantConcentration;

  forAll< parallelDevicePolicy<> >( dataGroup.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
  {
    isProppantMobile[a] = aperture[a] > minAperture && conc[a] < maxProppantConcentration;
  } );

}

void ProppantTransport::updateState( ObjectManagerBase & dataGroup )
{
  GEOS_MARK_FUNCTION;

  updateFluidModel( dataGroup );
  updateProppantModel( dataGroup );
}

void ProppantTransport::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );


  integer const numComponents = m_numComponents;

  // We have to redo the below loop after fractures are generated
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { fields::proppant::proppantConcentration::key(),
                                       fields::proppant::componentConcentration::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );

    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      // We have to redo the below loop after fractures are generated
      updateState( subRegion );

      SlurryFluidBase const & fluid =
        getConstitutiveModel< SlurryFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
      arrayView3d< real64 const > const componentDens = fluid.componentDensity();
      arrayView2d< real64 > const componentDens_n = subRegion.getField< fields::proppant::componentDensity_n >();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        for( localIndex c = 0; c < numComponents; ++c )
        {
          componentDens_n[ei][c] = componentDens[ei][0][c];
        }
      } );
    } );
  } );
  m_minAperture = m_bridgingFactor * m_proppantDiameter;
}

void ProppantTransport::preStepUpdate( real64 const & time,
                                       real64 const & GEOS_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {

    FlowSolverBase::precomputeData( mesh, regionNames );

    if( time <= 0 )
    {
      mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                    ElementSubRegionBase & subRegion )
      {
        updateProppantMobility( subRegion );
      } );
    }

    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

      arrayView3d< real64 const > const componentDens = fluid.componentDensity();
      arrayView2d< real64 > const componentDens_n = subRegion.getField< fields::proppant::componentDensity_n >();

      arrayView1d< real64 > const excessPackVolume = subRegion.getField< fields::proppant::proppantExcessPackVolume >();
      arrayView2d< real64 > const cellBasedFlux = subRegion.getField< fields::proppant::cellBasedFlux >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        for( localIndex c = 0; c < m_numComponents; ++c )
        {
          componentDens_n[ei][c] = componentDens[ei][0][c];
        }
        excessPackVolume[ei] = 0.0;
        LvArray::tensorOps::fill< 3 >( cellBasedFlux[ei], 0.0 );
      } );
    } );

  } );

  updateCellBasedFlux( time, domain );
}

void ProppantTransport::postStepUpdate( real64 const & time_n,
                                        real64 const & dt_return,
                                        DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      updateProppantMobility( subRegion );
    } );

    real64 const maxProppantConcentration = m_maxProppantConcentration;
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const & packVolFrac = subRegion.getField< fields::proppant::proppantPackVolumeFraction >();
      arrayView1d< real64 > const & proppantConc = subRegion.getField< fields::proppant::proppantConcentration >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( proppantConc[ei] >= maxProppantConcentration || packVolFrac[ei] >= 1.0 )
        {
          packVolFrac[ei] = 1.0;
          proppantConc[ei] = maxProppantConcentration;
        }
      } );
    } );
  } );
  if( m_updateProppantPacking == 1 )
  {
    updateProppantPackVolume( time_n, dt_return, domain );
  }
}

void ProppantTransport::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                           real64 const & GEOS_UNUSED_PARAM( dt ),
                                           DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const proppantConc = subRegion.getField< fields::proppant::proppantConcentration >();
      arrayView1d< real64 > const proppantConc_n = subRegion.getField< fields::proppant::proppantConcentration_n >();
      proppantConc_n.setValues< parallelDevicePolicy<> >( proppantConc );

      arrayView2d< real64 const > const componentConc = subRegion.getField< fields::proppant::componentConcentration >();
      arrayView2d< real64 > const componentConc_n = subRegion.getField< fields::proppant::componentConcentration_n >();
      componentConc_n.setValues< parallelDevicePolicy<> >( componentConc );
    } );
  } );
}

void ProppantTransport::implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                              real64 const & GEOS_UNUSED_PARAM( dt ),
                                              DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const proppantLiftFlux =
        subRegion.getField< fields::proppant::proppantLiftFlux >();
      proppantLiftFlux.zero();
    } );
  } );
}

void ProppantTransport::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                   DofManager & dofManager ) const
{
  for( auto const & meshTarget: getMeshTargets())
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: MeshBody = ({},{}) - target region = {}",
                               getName(), meshTarget.first.first.c_str(), meshTarget.first.second.c_str(), meshTarget.second ));
  }

  dofManager.addField( fields::proppant::proppantConcentration::key(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  dofManager.addCoupling( fields::proppant::proppantConcentration::key(),
                          fields::proppant::proppantConcentration::key(),
                          DofManager::Connector::Face );
}


void ProppantTransport::assembleSystem( real64 const GEOS_UNUSED_PARAM( time ),
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationTerms( dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void ProppantTransport::assembleAccumulationTerms( real64 const dt,
                                                   DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( fields::proppant::proppantConcentration::key() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      arrayView2d< real64 const > const componentDens_n = subRegion.getField< fields::proppant::componentDensity_n >();
      arrayView1d< real64 const > const proppantConc = subRegion.getField< fields::proppant::proppantConcentration >();
      arrayView1d< real64 const > const proppantConc_n = subRegion.getField< fields::proppant::proppantConcentration_n >();
      arrayView1d< real64 const > const proppantPackVolFrac = subRegion.getField< fields::proppant::proppantPackVolumeFraction >();
      arrayView1d< real64 const > const proppantLiftFlux = subRegion.getField< fields::proppant::proppantLiftFlux >();

      SlurryFluidBase const & fluid =
        getConstitutiveModel< SlurryFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );

      arrayView3d< real64 const > const componentDens = fluid.componentDensity();
      arrayView3d< real64 const > const dCompDens_dPres = fluid.dComponentDensity_dPressure();
      arrayView4d< real64 const > const dCompDens_dCompConc = fluid.dComponentDensity_dComponentConcentration();

      AccumulationKernel::launch( subRegion.size(),
                                  m_numComponents,
                                  m_numDofPerCell,
                                  dofManager.rankOffset(),
                                  dofNumber,
                                  elemGhostRank,
                                  proppantConc_n,
                                  proppantConc,
                                  componentDens_n,
                                  componentDens,
                                  dCompDens_dPres,
                                  dCompDens_dCompConc,
                                  volume,
                                  proppantPackVolFrac,
                                  proppantLiftFlux,
                                  dt,
                                  m_maxProppantConcentration,
                                  localMatrix,
                                  localRhs );
    } );
  } );
}


void ProppantTransport::assembleFluxTerms( real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );


    string const dofKey = dofManager.getKey( fields::proppant::proppantConcentration::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofKey );

    typename FluxKernel::FlowAccessors flowAccessors( elemManager, getName() );
    typename FluxKernel::ParticleFluidAccessors particleFluidAccessors( elemManager, getName() );
    typename FluxKernel::SlurryFluidAccessors slurryFluidAccessors( elemManager, getName() );
    typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, getName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
    {

      SurfaceElementStencilWrapper stencilWrapper = stencil.createKernelWrapper();

      FluxKernel::launch( stencilWrapper,
                          m_numDofPerCell,
                          dt,
                          dofManager.rankOffset(),
                          m_updateProppantPacking,
                          downVector,
                          dofNumberAccessor.toNestedViewConst(),
                          flowAccessors.get< fields::ghostRank >(),
                          flowAccessors.get< fields::flow::pressure >(),
                          flowAccessors.get< fields::proppant::proppantConcentration >(),
                          slurryFluidAccessors.get< fields::slurryfluid::componentDensity >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dComponentDensity_dPressure >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dComponentDensity_dComponentConcentration >(),
                          flowAccessors.get< fields::flow::gravityCoefficient >(),
                          slurryFluidAccessors.get< fields::singlefluid::density >(),
                          slurryFluidAccessors.get< fields::singlefluid::dDensity_dPressure >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dDensity_dProppantConcentration >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dDensity_dComponentConcentration >(),
                          slurryFluidAccessors.get< fields::singlefluid::viscosity >(),
                          slurryFluidAccessors.get< fields::singlefluid::dViscosity_dPressure >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dViscosity_dProppantConcentration >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dViscosity_dComponentConcentration >(),
                          slurryFluidAccessors.get< fields::slurryfluid::fluidDensity >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dFluidDensity_dPressure >(),
                          slurryFluidAccessors.get< fields::slurryfluid::dFluidDensity_dComponentConcentration >(),
                          particleFluidAccessors.get< fields::particlefluid::settlingFactor >(),
                          particleFluidAccessors.get< fields::particlefluid::dSettlingFactor_dPressure >(),
                          particleFluidAccessors.get< fields::particlefluid::dSettlingFactor_dProppantConcentration >(),
                          particleFluidAccessors.get< fields::particlefluid::dSettlingFactor_dComponentConcentration >(),
                          particleFluidAccessors.get< fields::particlefluid::collisionFactor >(),
                          particleFluidAccessors.get< fields::particlefluid::dCollisionFactor_dProppantConcentration >(),
                          flowAccessors.get< fields::proppant::isProppantMobile >(),
                          permAccessors.get< fields::permeability::permeability >(),
                          permAccessors.get< fields::permeability::permeabilityMultiplier >(),
                          flowAccessors.get< fields::elementAperture >(),
                          localMatrix,
                          localRhs );
    } );
  } );
}

void ProppantTransport::applyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  string const dofKey = dofManager.getKey( fields::proppant::proppantConcentration::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    //  Apply Dirichlet BC for proppant concentration

    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::proppant::proppantConcentration::key(),
                                             [&]( FieldSpecificationBase const & fs,
                                                  string const &,
                                                  SortedArrayView< localIndex const > const & lset,
                                                  ElementSubRegionBase & subRegion,
                                                  string const & )
    {
      arrayView1d< globalIndex const > const
      dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 const > const
      proppantConc = subRegion.getReference< array1d< real64 > >( fields::proppant::proppantConcentration::key() );

      fs.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                         parallelDevicePolicy<> >( lset,
                                                                   time_n + dt,
                                                                   subRegion,
                                                                   dofNumber,
                                                                   rankOffset,
                                                                   localMatrix,
                                                                   localRhs,
                                                                   proppantConc );
    } );

    //  Apply Dirichlet BC for component concentration
    if( m_numComponents > 0 )
    {
      map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

      fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                               mesh,
                                               fields::proppant::proppantConcentration::key(),
                                               [&]( FieldSpecificationBase const &,
                                                    string const & setName,
                                                    SortedArrayView< localIndex const > const &,
                                                    ElementSubRegionBase & subRegion,
                                                    string const & )
      {

        string const & subRegionName = subRegion.getName();
        GEOS_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0,
                       getDataContext() << ": Conflicting proppant boundary conditions on set " << setName );
        bcStatusMap[subRegionName][setName].resize( m_numComponents );
        bcStatusMap[subRegionName][setName].setValues< serialPolicy >( false );

      } );

      fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                               mesh,
                                               fields::proppant::componentConcentration::key(),
                                               [&] ( FieldSpecificationBase const & fs,
                                                     string const & setName,
                                                     SortedArrayView< localIndex const > const & targetSet,
                                                     ElementSubRegionBase & subRegion,
                                                     string const & )
      {

        string const & subRegionName = subRegion.getName();
        localIndex const comp = fs.getComponent();

        GEOS_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0,
                       getDataContext() << ": Proppant boundary condition not prescribed on set '" << setName << "'" );
        GEOS_ERROR_IF( bcStatusMap[subRegionName][setName][comp],
                       getDataContext() << ": Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
        bcStatusMap[subRegionName][setName][comp] = true;

        fs.applyFieldValue< FieldSpecificationEqual >( targetSet,
                                                       time_n + dt,
                                                       subRegion,
                                                       fields::proppant::bcComponentConcentration::key() );

      } );

      bool bcConsistent = true;
      for( auto const & bcStatusEntryOuter : bcStatusMap )
      {
        for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
        {
          for( localIndex ic = 0; ic < m_numComponents; ++ic )
          {
            bcConsistent &= bcStatusEntryInner.second[ic];
            GEOS_WARNING_IF( !bcConsistent,
                             getDataContext() << ": Composition boundary condition not applied to component " <<
                             ic << " on region '" << bcStatusEntryOuter.first << "'," <<
                             " set '" << bcStatusEntryInner.first << "'" );
          }
        }
      }

      GEOS_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

      fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                               mesh,
                                               fields::proppant::proppantConcentration::key(),
                                               [&] ( FieldSpecificationBase const &,
                                                     string const &,
                                                     SortedArrayView< localIndex const > const & targetSet,
                                                     ElementSubRegionBase & subRegion,
                                                     string const & )
      {
        arrayView1d< integer const > const ghostRank =
          subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
        arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

        arrayView2d< real64 const > const compConc =
          subRegion.getReference< array2d< real64 > >( fields::proppant::componentConcentration::key() );
        arrayView2d< real64 const > const bcCompConc =
          subRegion.getReference< array2d< real64 > >( fields::proppant::bcComponentConcentration::key() );

        forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
        {
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 )
            return;

          globalIndex const dofIndex = dofNumber[ei];
          localIndex const localRow = dofIndex - rankOffset;
          real64 rhsValue;

          for( localIndex ic = 0; ic < m_numComponents; ++ic )
          {
            FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                        rankOffset,
                                                        localMatrix,
                                                        rhsValue,
                                                        bcCompConc[ei][ic],
                                                        compConc[ei][ic] );
            localRhs[localRow + ic + 1] = rhsValue;
          }
        } );
      } );
    }
  } );
}

real64
ProppantTransport::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                          real64 const & GEOS_UNUSED_PARAM( dt ),
                                          DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 localResidualNorm = 0.0;
  real64 localResidualNormalizer = 0.0;

  solverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  localIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( fields::proppant::proppantConcentration::key() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[1]{};
      real64 subRegionResidualNormalizer[1]{};

      // step 1: compute the norm in the subRegion

      proppantTransportKernels::
        ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( normType,
                                                   m_numDofPerCell,
                                                   rankOffset,
                                                   dofKey,
                                                   localRhs,
                                                   subRegion,
                                                   m_nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm,
                                                   subRegionResidualNormalizer );

      // step 2: first reduction across meshBodies/regions/subRegions

      if( normType == solverBaseKernels::NormType::Linf )
      {
        if( subRegionResidualNorm[0] > localResidualNorm )
        {
          localResidualNorm = subRegionResidualNorm[0];
        }
      }
      else
      {
        localResidualNorm += subRegionResidualNorm[0];
        localResidualNormalizer += subRegionResidualNormalizer[0];
      }
    } );
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  if( normType == solverBaseKernels::NormType::Linf )
  {
    solverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm, residualNorm );
  }
  else
  {
    solverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm, localResidualNormalizer, residualNorm );
  }

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", ProppantTransport::coupledSolverAttributePrefix(), residualNorm );
  }

  return residualNorm;
}

void ProppantTransport::applySystemSolution( DofManager const & dofManager,
                                             arrayView1d< real64 const > const & localSolution,
                                             real64 const scalingFactor,
                                             real64 const dt,
                                             DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  dofManager.addVectorToField( localSolution,
                               fields::proppant::proppantConcentration::key(),
                               fields::proppant::proppantConcentration::key(),
                               scalingFactor,
                               { m_numDofPerCell, 0, 1 } );


  if( m_numDofPerCell > 1 )
  {
    dofManager.addVectorToField( localSolution,
                                 fields::proppant::proppantConcentration::key(),
                                 fields::proppant::componentConcentration::key(),
                                 scalingFactor,
                                 { m_numDofPerCell, 1, m_numDofPerCell } );
  }


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::proppant::proppantConcentration::key(),
                                       fields::proppant::componentConcentration::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );

    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      updateComponentDensity( subRegion );
    } );

  } );

}

void ProppantTransport::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const & proppantConc = subRegion.getField< fields::proppant::proppantConcentration >();
      arrayView1d< real64 const > const & proppantConc_n = subRegion.getField< fields::proppant::proppantConcentration_n >();
      proppantConc.setValues< parallelDevicePolicy<> >( proppantConc_n );

      arrayView2d< real64 > const & componentConc = subRegion.getField< fields::proppant::componentConcentration >();
      arrayView2d< real64 const > const & componentConc_n = subRegion.getField< fields::proppant::componentConcentration_n >();
      componentConc.setValues< parallelDevicePolicy<> >( componentConc_n );

      updateState( subRegion );
    } );
  } );
}



void ProppantTransport::updateCellBasedFlux( real64 const GEOS_UNUSED_PARAM( time_n ),
                                             DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    string const meshBodyName = mesh.getParent().getParent().getName();
    string const meshLevelName = mesh.getName();

    ElementRegionManager & elemManager = mesh.getElemManager();

    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & cellBasedFluxAccessor =
      elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( fields::proppant::cellBasedFlux::key() );

    typename FluxKernel::CellBasedFluxFlowAccessors flowAccessors( elemManager, getName() );
    typename FluxKernel::CellBasedFluxSlurryFluidAccessors slurryFluidAccessors( elemManager, getName() );
    typename FluxKernel::PermeabilityAccessors permAccessors( elemManager, getName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
    {
      SurfaceElementStencilWrapper stencilWrapper = stencil.createKernelWrapper();

      FluxKernel::launchCellBasedFluxCalculation( stencilWrapper,
                                                  downVector,
                                                  flowAccessors.get< fields::flow::pressure >(),
                                                  flowAccessors.get< fields::flow::gravityCoefficient >(),
                                                  slurryFluidAccessors.get< fields::singlefluid::density >(),
                                                  slurryFluidAccessors.get< fields::singlefluid::viscosity >(),
                                                  permAccessors.get< fields::permeability::permeability >(),
                                                  permAccessors.get< fields::permeability::permeabilityMultiplier >(),
                                                  flowAccessors.get< fields::elementAperture >(),
                                                  cellBasedFluxAccessor.toNestedView() );
    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::proppant::cellBasedFlux::key() }, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void ProppantTransport::updateProppantPackVolume( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    // For data modified through an accessor, we must create the view accessor
    // every time in order to ensure the data gets properly touched on device
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantConc =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( fields::proppant::proppantConcentration::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantPackVolFrac =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( fields::proppant::proppantPackVolumeFraction::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantExcessPackVolume =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( fields::proppant::proppantExcessPackVolume::key() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantLiftFlux =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( fields::proppant::proppantLiftFlux::key() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const
    aperture = elemManager.constructArrayViewAccessor< real64, 1 >( fields::elementAperture::key() );

    typename ProppantPackVolumeKernel::FlowAccessors flowAccessors( elemManager, getName() );
    typename ProppantPackVolumeKernel::SlurryFluidAccessors slurryFluidAccessors( elemManager, getName() );
    typename ProppantPackVolumeKernel::ParticleFluidAccessors particleFluidAccessors( elemManager, getName() );

    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
    {
      ProppantPackVolumeKernel::launchProppantPackVolumeCalculation( stencil,
                                                                     dt,
                                                                     m_proppantDensity,
                                                                     m_proppantDiameter,
                                                                     m_maxProppantConcentration,
                                                                     downVector,
                                                                     m_criticalShieldsNumber,
                                                                     m_frictionCoefficient,
                                                                     particleFluidAccessors.get< fields::particlefluid::settlingFactor >(),
                                                                     slurryFluidAccessors.get< fields::singlefluid::density >(),
                                                                     slurryFluidAccessors.get< fields::slurryfluid::fluidDensity >(),
                                                                     slurryFluidAccessors.get< fields::slurryfluid::fluidViscosity >(),
                                                                     flowAccessors.get< fields::proppant::isProppantMobile >(),
                                                                     flowAccessors.get< fields::proppant::isProppantBoundary >(),
                                                                     flowAccessors.get< fields::elementAperture >(),
                                                                     flowAccessors.get< fields::elementVolume >(),
                                                                     flowAccessors.get< fields::ghostRank >(),
                                                                     flowAccessors.get< fields::proppant::cellBasedFlux >(),
                                                                     proppantConc.toNestedView(),
                                                                     proppantPackVolFrac.toNestedView(),
                                                                     proppantExcessPackVolume.toNestedView(),
                                                                     proppantLiftFlux.toNestedView() );
    } );

    {
      FieldIdentifiers fieldsToBeSync;
      fieldsToBeSync.addElementFields( { fields::proppant::proppantConcentration::key(),
                                         fields::proppant::proppantPackVolumeFraction::key(),
                                         fields::proppant::proppantExcessPackVolume::key(),
                                         fields::proppant::proppantLiftFlux::key() },
                                       regionNames );

      CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
    }

    elemManager.forElementSubRegions( regionNames,
                                      [&]( localIndex const,
                                           ElementSubRegionBase & subRegion )
    {
      updateProppantMobility( subRegion );
    } );


    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
    {
      ProppantPackVolumeKernel::launchProppantPackVolumeUpdate( stencil,
                                                                downVector,
                                                                m_maxProppantConcentration,
                                                                flowAccessors.get< fields::proppant::isProppantMobile >(),
                                                                proppantExcessPackVolume.toNestedViewConst(),
                                                                proppantConc.toNestedView(),
                                                                proppantPackVolFrac.toNestedView() );
    } );

    {
      FieldIdentifiers fieldsToBeSync;

      fieldsToBeSync.addElementFields( { fields::proppant::proppantConcentration::key(),
                                         fields::proppant::proppantPackVolumeFraction::key() },
                                       regionNames );

      CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
    }

    elemManager.forElementSubRegions( regionNames,
                                      [&]( localIndex const,
                                           ElementSubRegionBase & subRegion )
    {
      updateProppantMobility( subRegion );
    } );

  } );

}


REGISTER_CATALOG_ENTRY( SolverBase, ProppantTransport, string const &, Group * const )
} /* namespace geos */
