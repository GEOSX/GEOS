/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProppantTransport.cpp
 */

#include "ProppantTransport.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SlurryFluidBase.hpp"
#include "constitutive/fluid/ParticleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/FaceElementRegion.hpp"

#include "physicsSolvers/fluidFlow/ProppantTransportKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace ProppantTransportKernels;

ProppantTransport::ProppantTransport( const std::string & name,
                                      Group * const parent ):
  FlowSolverBase( name, parent )
{
  this->registerWrapper( viewKeyStruct::proppantNamesString, &m_proppantModelNames )->setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of proppant constitutive object to use for this solver." );

  registerWrapper( viewKeyStruct::bridgingFactorString, &m_bridgingFactor )->setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Bridging factor used for bridging/screen-out calculation" );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString, &m_maxProppantConcentration )->setApplyDefaultValue( 0.6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum proppant concentration" );

  registerWrapper( viewKeyStruct::proppantDiameterString, &m_proppantDiameter )->setApplyDefaultValue( 0.4e-3 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Proppant diameter" );

  registerWrapper( viewKeyStruct::proppantDensityString, &m_proppantDensity )->setApplyDefaultValue( 2500.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Proppant density" );

  registerWrapper( viewKeyStruct::criticalShieldsNumberString, &m_criticalShieldsNumber )->setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Critical Shields number" );

  registerWrapper( viewKeyStruct::frictionCoefficientString, &m_frictionCoefficient )->setApplyDefaultValue( 0.03 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Friction coefficient" );

  registerWrapper( viewKeyStruct::updateProppantPackingString, &m_updateProppantPacking )->setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Flag that enables/disables proppant-packing update" );

}

void ProppantTransport::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();
  CheckModelNames( m_proppantModelNames, viewKeyStruct::proppantNamesString );
}

void ProppantTransport::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    forTargetSubRegions< CellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantConcentrationString )->setDefaultValue( 0.0 )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString )->setDefaultValue( 0.0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::componentConcentrationString )->setDefaultValue( 0.0 )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString )->setDefaultValue( 0.0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::updatedComponentConcentrationString )->setDefaultValue( 0.0 );
      subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::cellBasedFluxString )->setDefaultValue( { 0.0, 0.0, 0.0 } );
      subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::isProppantBoundaryString )->setDefaultValue( 0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString )->setDefaultValue( 0.0 );
    } );


    forTargetSubRegions< FaceElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 FaceElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantConcentrationString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::componentConcentrationString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::updatedComponentConcentrationString );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::oldProppantConcentrationString );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::oldComponentDensityString );
      subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::cellBasedFluxString );
      subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::isInterfaceElementString );
      subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::isProppantBoundaryString );
      subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::isProppantMobileString );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::poroMultiplierString )->setDefaultValue( 1.0 );
      subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::transTMultiplierString )->setDefaultValue( { 1.0, 1.0, 1.0 } );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString )->setDefaultValue( 0.0 );
    } );

  }
}

void ProppantTransport::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  ConstitutiveManager & cm = *domain->getConstitutiveManager();

  // Validate proppant models in regions
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping< SlurryFluidBase >( *meshLevel.getElemManager(), m_fluidModelNames );
    ValidateModelMapping< ParticleFluidBase >( *meshLevel.getElemManager(), m_proppantModelNames );
  }

  SlurryFluidBase const * fluid0 = cm.GetConstitutiveRelation< SlurryFluidBase >( m_fluidModelNames[0] );

  m_numComponents = fluid0->numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    forTargetSubRegions< CellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
    {
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString ).resizeDimension< 1 >( NC );
    } );
  }


}

void ProppantTransport::ResizeFractureFields( MeshLevel & mesh )
{
  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const,
                                                            FaceElementSubRegion & subRegion )
    {
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::updatedComponentConcentrationString ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString ).resizeDimension< 1 >( NC );
    } );
  }
}

void ProppantTransport::UpdateFluidModel( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  SlurryFluidBase & fluid = GetConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView1d< real64 const > const & pres  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView2d< real64 const > const & componentConc  = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString );
  arrayView2d< real64 const > const & dComponentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

  arrayView2d< real64 > const & updatedComponentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::updatedComponentConcentrationString );

  forAll< serialPolicy >( dataGroup.size(), [&] ( localIndex const a )
  {
    for( localIndex c = 0; c < m_numComponents; ++c )
    {
      updatedComponentConc[a][c] = componentConc[a][c] + dComponentConc[a][c];
    }

    fluid.PointUpdateFluidProperty( pres[a] + dPres[a], updatedComponentConc[a], 0.0, a, 0 );
  } );

}

void ProppantTransport::UpdateComponentDensity( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  SlurryFluidBase & fluid = GetConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView1d< real64 const > const & pres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView2d< real64 const > const & componentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString );
  arrayView2d< real64 const > const & dComponentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

  arrayView2d< real64 > const & updatedComponentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::updatedComponentConcentrationString );

  forAll< serialPolicy >( dataGroup.size(), [&] ( localIndex const a )
  {
    for( localIndex c = 0; c < m_numComponents; ++c )
    {
      updatedComponentConc[a][c] = componentConc[a][c] + dComponentConc[a][c];
    }

    fluid.PointUpdateComponentDensity( pres[a] + dPres[a], updatedComponentConc[a], a, 0 );
  } );

}


void ProppantTransport::UpdateProppantModel( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  localIndex const NC = m_numComponents;

  arrayView1d< real64 const > const & proppantConc  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );
  arrayView1d< real64 const > const & dProppantConc = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );

  SlurryFluidBase const & fluid = GetConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView2d< real64 const > const & fluidDens            = fluid.fluidDensity();
  arrayView2d< real64 const > const & dFluidDens_dPres     = fluid.dFluidDensity_dPressure();
  arrayView3d< real64 const > const & dFluidDens_dCompConc = fluid.dFluidDensity_dComponentConcentration();
  arrayView2d< real64 const > const & fluidVisc            = fluid.fluidViscosity();
  arrayView2d< real64 const > const & dFluidVisc_dPres     = fluid.dFluidViscosity_dPressure();
  arrayView3d< real64 const > const & dFluidVisc_dCompConc = fluid.dFluidViscosity_dComponentConcentration();

  ParticleFluidBase & particle = GetConstitutiveModel< ParticleFluidBase >( dataGroup, m_proppantModelNames[targetIndex] );

  forAll< serialPolicy >( dataGroup.size(), [&] ( localIndex const a )
  {
    particle.PointUpdate( NC,
                          proppantConc[a] + dProppantConc[a],
                          fluidDens[a][0],
                          dFluidDens_dPres[a][0],
                          dFluidDens_dCompConc[a][0],
                          fluidVisc[a][0],
                          dFluidVisc_dPres[a][0],
                          dFluidVisc_dCompConc[a][0], a );
  } );

}

void ProppantTransport::UpdateProppantMobility( Group & dataGroup )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & conc = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );
  arrayView1d< real64 const > const & aperture = dataGroup.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementApertureString );
  arrayView1d< integer > const & isProppantMobile = dataGroup.getReference< array1d< integer > >( viewKeyStruct::isProppantMobileString );

  real64 const minAperture = m_minAperture;
  real64 const maxProppantConcentration = m_maxProppantConcentration;

  forAll< parallelHostPolicy >( dataGroup.size(), [=]( localIndex const a )
  {
    isProppantMobile[a] = aperture[a] > minAperture && conc[a] < maxProppantConcentration;
  } );

}

void ProppantTransport::UpdateState( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup, targetIndex );
  UpdateProppantModel( dataGroup, targetIndex );
}

void ProppantTransport::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantConcentrationString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::componentConcentrationString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );

  ResetViews( mesh );

  // We have to redo the below loop after fractures are generated

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    UpdateState( subRegion, targetIndex );

    arrayView1d< real64 > const & proppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );
    arrayView1d< real64 > const & proppantConcOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldProppantConcentrationString );

    SlurryFluidBase const & fluid = GetConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );
    arrayView3d< real64 const > const & componentDens = fluid.componentDensity();
    arrayView2d< real64 > const & componentDensOld = subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString );

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      proppantConcOld[ei] = proppantConc[ei];
      for( localIndex c = 0; c < NC; ++c )
      {
        componentDensOld[ei][c] = componentDens[ei][0][c];
      }
    } );
  } );

  m_downVector = gravityVector();
  m_downVector.Normalize();

  m_minAperture = m_bridgingFactor * m_proppantDiameter;

  real64 const oneMinusMaxConc = (1.0 - m_maxProppantConcentration);
  m_proppantPackPermeability = m_proppantDiameter * m_proppantDiameter / 180.0 *
                               (oneMinusMaxConc * oneMinusMaxConc * oneMinusMaxConc) / (m_maxProppantConcentration * m_maxProppantConcentration);

}

real64 ProppantTransport::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      const int cycleNumber,
                                      DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FlowSolverBase::PrecomputeData( mesh );

  NodeManager const & nodeManager = *mesh.getNodeManager();
  FaceManager const & faceManager = *mesh.getFaceManager();

  real64 dt_return = dt;

  ImplicitStepSetup( time_n, dt, domain );

  if( cycleNumber == 0 )
  {
    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager.ApplyInitialConditions( &domain );

    localIndex const NC = m_numComponents;

    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
    {
      subRegion.CalculateElementGeometricQuantities( nodeManager, faceManager );

      UpdateState( subRegion, targetIndex );

      arrayView1d< real64 > const & dProppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );
      arrayView2d< real64 > const & dComponentConc = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

      arrayView1d< real64 > const & proppantConcOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldProppantConcentrationString );
      arrayView1d< real64 const > const & proppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );

      SlurryFluidBase const & fluid = GetConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );

      arrayView2d< real64 > const & componentDensOld = subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString );
      arrayView3d< real64 const > const & componentDens = fluid.componentDensity();

      arrayView1d< R1Tensor > const & cellBasedFlux = subRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::cellBasedFluxString );;
      arrayView1d< real64 > const & proppantLiftFlux = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString );
      arrayView1d< real64 > const & packVf = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString );

      arrayView1d< real64 > const & poroMultiplier = subRegion.getReference< array1d< real64 > >( viewKeyStruct::poroMultiplierString );
      arrayView1d< R1Tensor > const & transTMultiplier = subRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::transTMultiplierString );
      arrayView1d< real64 > const & excessPackV = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString );

      arrayView1d< integer > const & isInterfaceElement = subRegion.getReference< array1d< integer > >( viewKeyStruct::isInterfaceElementString );
      arrayView1d< integer > const & isProppantMobile = subRegion.getReference< array1d< integer > >( viewKeyStruct::isProppantMobileString );


      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
      {
        dProppantConc[ei] = 0.0;
        proppantConcOld[ei] = proppantConc[ei];

        for( localIndex c = 0; c < NC; ++c )
        {
          dComponentConc[ei][c] = 0.0;
          componentDensOld[ei][c] = componentDens[ei][0][c];
        }

        cellBasedFlux[ei] = 0.0;
        proppantLiftFlux[ei] = 0.0;

        packVf[ei] = 0.0;
        excessPackV[ei] = 0.0;

        poroMultiplier[ei] = 1.0;
        transTMultiplier[ei] = 1.0;

        isInterfaceElement[ei] = 0;
        isProppantMobile[ei] = 1;
      } );

    } );
  }

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateProppantMobility( subRegion );
  } );

  UpdateCellBasedFlux( time_n, domain );

  // currently the only method is implicit time integration
  dt_return = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateProppantMobility( subRegion );
  } );

  if( m_updateProppantPacking == 1 )
  {
    UpdateProppantPackVolume( time_n, dt_return, domain );
  }

  return dt_return;
}


void ProppantTransport::PreStepUpdate( real64 const & time,
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FlowSolverBase::PrecomputeData( mesh );

  NodeManager const & nodeManager = *mesh.getNodeManager();
  FaceManager const & faceManager = *mesh.getFaceManager();

  if( time <= 0 )
  {
    /* Below must be called after ImplicitStepSetup */

    ResizeFractureFields( mesh );

    localIndex const NC = m_numComponents;

    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
    {
      subRegion.CalculateElementGeometricQuantities( nodeManager, faceManager );

      UpdateState( subRegion, targetIndex );

      arrayView1d< real64 > const & dProppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );
      arrayView2d< real64 > const & dComponentConc = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

      arrayView1d< real64 > const & proppantConcOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldProppantConcentrationString );
      arrayView1d< real64 const > const & proppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );

      SlurryFluidBase const & fluid = GetConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );

      arrayView2d< real64 > const & componentDensOld = subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString );
      arrayView3d< real64 const > const & componentDens = fluid.componentDensity();

      arrayView1d< R1Tensor > const & cellBasedFlux = subRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::cellBasedFluxString );
      arrayView1d< real64 > const & proppantLiftFlux = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString );
      arrayView1d< real64 > const & packVf = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString );

      arrayView1d< real64 > const & poroMultiplier = subRegion.getReference< array1d< real64 > >( viewKeyStruct::poroMultiplierString );
      arrayView1d< R1Tensor > const & transTMultiplier = subRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::transTMultiplierString );
      arrayView1d< real64 > const & excessPackV = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString );

      arrayView1d< integer > const & isInterfaceElement = subRegion.getReference< array1d< integer > >( viewKeyStruct::isInterfaceElementString );
      arrayView1d< integer > const & isProppantMobile = subRegion.getReference< array1d< integer > >( viewKeyStruct::isProppantMobileString );

      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
      {
        dProppantConc[ei] = 0.0;
        proppantConcOld[ei] = proppantConc[ei];

        for( localIndex c = 0; c < NC; ++c )
        {
          dComponentConc[ei][c] = 0.0;
          componentDensOld[ei][c] = componentDens[ei][0][c];
        }

        cellBasedFlux[ei] = 0.0;
        proppantLiftFlux[ei] = 0.0;

        packVf[ei] = 0.0;
        excessPackV[ei] = 0.0;

        poroMultiplier[ei] = 1.0;
        transTMultiplier[ei] = 1.0;

        isInterfaceElement[ei] = 0;
        isProppantMobile[ei] = 1;
      } );
    } );
  }

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateProppantMobility( subRegion );
    UpdateState( subRegion, targetIndex );
  } );

  UpdateCellBasedFlux( time, domain );
}

void ProppantTransport::PostStepUpdate( real64 const & time_n,
                                        real64 const & dt_return,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateProppantMobility( subRegion );
  } );

  if( m_updateProppantPacking == 1 )
  {
    UpdateProppantPackVolume( time_n, dt_return, domain );
  }
}

void ProppantTransport::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const & GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition & domain )
{

  localIndex const NC = m_numComponents;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  ResetViews( mesh );

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dProppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );
    arrayView2d< real64 > const & dComponentConc = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

    arrayView1d< real64 > const & proppantConcOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldProppantConcentrationString );
    arrayView1d< real64 const > const & proppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );

    SlurryFluidBase const & fluid = GetConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );

    arrayView2d< real64 > const & componentDensOld = subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString );
    arrayView3d< real64 const > const & componentDens = fluid.componentDensity();

    arrayView1d< R1Tensor > const & cellBasedFlux = subRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::cellBasedFluxString );;
    arrayView1d< real64 > const & proppantLiftFlux = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString );
    arrayView1d< real64 > const & excessPackV = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString );

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      dProppantConc[ei] = 0.0;
      proppantConcOld[ei] = proppantConc[ei];

      for( localIndex c = 0; c < NC; ++c )
      {
        dComponentConc[ei][c] = 0.0;
        componentDensOld[ei][c] = componentDens[ei][0][c];
      }

      excessPackV[ei] = 0.0;
      proppantLiftFlux[ei] = 0.0;
      cellBasedFlux[ei] = 0.0;
    } );
  } );

}

void ProppantTransport::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & GEOSX_UNUSED_PARAM( dt ),
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & proppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );
    arrayView1d< real64 const > const & dProppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );

    arrayView2d< real64 > const & componentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString );
    arrayView2d< real64 const > const & dComponentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

    arrayView1d< real64 > const & proppantLiftFlux =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString );

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      proppantConc[ei] += dProppantConc[ei];
      proppantLiftFlux[ei] = 0.0;

      for( localIndex c = 0; c < NC; ++c )
      {
        componentConc[ei][c] += dComponentConc[ei][c];
      }
    } );
  } );

}

void ProppantTransport::SetupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                   DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::proppantConcentrationString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::proppantConcentrationString,
                          viewKeyStruct::proppantConcentrationString,
                          DofManager::Connector::Face );
}


void ProppantTransport::AssembleSystem( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  AssembleAccumulationTerms( domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  AssembleFluxTerms( time,
                     dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void ProppantTransport::AssembleAccumulationTerms( DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  localIndex const NC = m_numComponents;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
  {
    string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString );
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

    arrayView1d< real64 const > const & proppantConcOld =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldProppantConcentrationString );
    arrayView2d< real64 const > const & componentDensOld =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString );
    arrayView1d< real64 const > const & proppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );
    arrayView1d< real64 const > const & dProppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );
    arrayView1d< real64 const > const & proppantPackVf =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString );

    SlurryFluidBase const & fluid = GetConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );

    arrayView3d< real64 const > const & componentDens = fluid.componentDensity();
    arrayView3d< real64 const > const & dCompDens_dPres = fluid.dComponentDensity_dPressure();
    arrayView4d< real64 const > const & dCompDens_dCompConc = fluid.dComponentDensity_dComponentConcentration();

    globalIndex const rankOffset = dofManager.rankOffset();
    localIndex const NDOF = m_numDofPerCell;

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        stackArray1d< globalIndex, MAX_NUM_COMPONENTS > localAccumDOF( m_numDofPerCell );
        stackArray1d< real64, MAX_NUM_COMPONENTS > localAccum( m_numDofPerCell );
        stackArray2d< real64, MAX_NUM_COMPONENTS * MAX_NUM_COMPONENTS > localAccumJacobian( m_numDofPerCell, m_numDofPerCell );

        real64 effectiveVolume = volume[ei];
        real64 packPoreVolume = 0.0;

        if( proppantPackVf[ei] < 1.0 )
        {
          effectiveVolume = volume[ei] * ( 1.0 - proppantPackVf[ei] );
          packPoreVolume = volume[ei] * proppantPackVf[ei] * ( 1.0 - m_maxProppantConcentration );
        }

        AccumulationKernel::Compute( NC,
                                     proppantConcOld[ei],
                                     proppantConc[ei] + dProppantConc[ei],
                                     componentDensOld[ei],
                                     componentDens[ei][0],
                                     dCompDens_dPres[ei][0],
                                     dCompDens_dCompConc[ei][0],
                                     effectiveVolume,
                                     packPoreVolume,
                                     localAccum,
                                     localAccumJacobian );

        globalIndex const elemDOF = dofNumber[ei];

        for( localIndex idof = 0; idof < NDOF; ++idof )
        {
          localAccumDOF[idof] = elemDOF + idof;
        }
        // add contribution to global residual and dRdP

        localIndex const localRow = dofNumber[ei] - rankOffset;
        for( localIndex idof = 0; idof < NDOF; ++idof )
        {
          localRhs[localRow + idof] += localAccum[idof];
          localMatrix.addToRow< serialAtomic >( localRow + idof,
                                                localAccumDOF.data(),
                                                localAccumJacobian[idof].dataIfContiguous(),
                                                NDOF );
        }
      }
    } );
  } );

}


void ProppantTransport::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  FluxApproximationBase const & fluxApprox = *fvManager.getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = elemManager.ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofKey );

  FluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & pres  = m_pressure.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dPres = m_deltaPressure.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantConc   = m_proppantConcentration.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dProppantConc  = m_deltaProppantConcentration.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & gravCoef = m_gravCoef.toViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dens        = m_density.toViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres = m_dDensity_dPressure.toViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dDens_dProppantConc = m_dDensity_dProppantConcentration.toViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dDens_dComponentConc = m_dDensity_dComponentConcentration.toViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & visc        = m_viscosity.toViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dVisc_dPres = m_dViscosity_dPressure.toViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dVisc_dProppantConc = m_dViscosity_dProppantConcentration.toViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dVisc_dComponentConc = m_dViscosity_dComponentConcentration.toViewConst();

  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & componentDens = m_componentDensity.toViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dComponentDens_dPres = m_dComponentDensity_dPressure.toViewConst();
  FluxKernel::ElementViewConst< arrayView4d< real64 const > > const & dComponentDens_dComponentConc = m_dComponentDensity_dComponentConcentration.toViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & fluidDensity = m_fluidDensity.toViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dFluidDens_dPres = m_dFluidDensity_dPressure.toViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const & dFluidDens_dComponentConc = m_dFluidDensity_dComponentConcentration.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & settlingFactor = m_settlingFactor.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dPres = m_dSettlingFactor_dPressure.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dProppantConc = m_dSettlingFactor_dProppantConcentration.toViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dSettlingFactor_dComponentConc = m_dSettlingFactor_dComponentConcentration.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & collisionFactor = m_collisionFactor.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & dCollisionFactor_dProppantConc = m_dCollisionFactor_dProppantConcentration.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< integer const > > const & isProppantMobile  = m_isProppantMobile.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const & transTMultiplier  = m_transTMultiplier.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & aperture  = m_elementAperture.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf  = m_proppantPackVolumeFraction.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantLiftFlux  = m_proppantLiftFlux.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< integer const > > const & isInterfaceElement  = m_isInterfaceElement.toViewConst();

  FluxKernel::ElementViewConst< arrayView1d< integer const > > const & elemGhostRank = m_elemGhostRank.toViewConst();

  fluxApprox.forStencils< FaceElementStencil >( [&]( auto const & stencil )
  {

    FluxKernel::Launch( stencil,
                        m_numDofPerCell,
                        dt,
                        dofManager.rankOffset(),
                        transTMultiplier,
                        m_updateProppantPacking,
                        m_downVector,
                        dofNumber,
                        elemGhostRank,
                        pres,
                        dPres,
                        proppantConc,
                        dProppantConc,
                        componentDens,
                        dComponentDens_dPres,
                        dComponentDens_dComponentConc,
                        gravCoef,
                        dens,
                        dDens_dPres,
                        dDens_dProppantConc,
                        dDens_dComponentConc,
                        visc,
                        dVisc_dPres,
                        dVisc_dProppantConc,
                        dVisc_dComponentConc,
                        fluidDensity,
                        dFluidDens_dPres,
                        dFluidDens_dComponentConc,
                        settlingFactor,
                        dSettlingFactor_dPres,
                        dSettlingFactor_dProppantConc,
                        dSettlingFactor_dComponentConc,
                        collisionFactor,
                        dCollisionFactor_dProppantConc,
                        isProppantMobile,
                        proppantPackVf,
                        aperture,
                        proppantLiftFlux,
                        isInterfaceElement,
                        localMatrix,
                        localRhs );
  } );
}

void ProppantTransport::ApplyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString );
  globalIndex const rankOffset = dofManager.rankOffset();

  //  Apply Dirichlet BC for proppant concentration

  fsManager.Apply( time_n + dt,
                   &domain,
                   "ElementRegions",
                   viewKeyStruct::proppantConcentrationString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & )
  {
    arrayView1d< globalIndex const > const &
    dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const &
    proppantConc = subRegion->getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString );

    arrayView1d< real64 const > const &
    dProppantConc = subRegion->getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );

    fs->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                        parallelHostPolicy >( lset,
                                                              time_n + dt,
                                                              subRegion,
                                                              dofNumber,
                                                              rankOffset,
                                                              localMatrix,
                                                              localRhs,
                                                              [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      return proppantConc[a] + dProppantConc[a];
    } );
  } );

  //  Apply Dirichlet BC for component concentration

  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     viewKeyStruct::proppantConcentrationString,
                     [&]( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( fs ),
                          string const & setName,
                          SortedArrayView< localIndex const > const & GEOSX_UNUSED_PARAM( targetSet ),
                          Group * subRegion,
                          string const & )
    {

      string const & subRegionName = subRegion->getName();
      GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0, "Conflicting proppant boundary conditions on set " << setName );
      bcStatusMap[subRegionName][setName].resize( NC );
      bcStatusMap[subRegionName][setName] = false;

    } );

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     viewKeyStruct::componentConcentrationString,
                     [&] ( FieldSpecificationBase const * const fs,
                           string const & setName,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group * subRegion,
                           string const & )
    {

      string const & subRegionName = subRegion->getName();
      localIndex const comp = fs->GetComponent();

      GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0, "Proppant boundary condition not prescribed on set '" << setName << "'" );
      GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
      bcStatusMap[subRegionName][setName][comp] = true;

      fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                      time_n + dt,
                                                      subRegion,
                                                      viewKeyStruct::bcComponentConcentrationString );

    } );

    bool bcConsistent = true;
    for( auto const & bcStatusEntryOuter : bcStatusMap )
    {
      for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          bcConsistent &= bcStatusEntryInner.second[ic];
          GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic
                                                                                                      << " on region '" << bcStatusEntryOuter.first << "',"
                                                                                                      << " set '" << bcStatusEntryInner.first << "'" );
        }
      }
    }

    GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     viewKeyStruct::proppantConcentrationString,
                     [&] ( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( bc ),
                           string const & GEOSX_UNUSED_PARAM( setName ),
                           SortedArrayView< localIndex const > const & targetSet,
                           Group * subRegion,
                           string const & )
    {
      arrayView1d< integer const > const ghostRank =
        subRegion->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
      arrayView1d< globalIndex const > const & dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

      arrayView2d< real64 const > const & compConc =
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString );
      arrayView2d< real64 const > const & deltaCompConc =
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );
      arrayView2d< real64 const > const & bcCompConc =
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString );

      forAll< parallelHostPolicy >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
          return;

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcCompConc[a][ic],
                                                      compConc[a][ic] + deltaCompConc[a][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }
      } );
    } );
  }
}

real64
ProppantTransport::CalculateResidualNorm( DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localRhs )
{
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString );

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

    RAJA::ReduceSum< parallelHostReduce, real64 > localSum( 0.0 );

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;
        for( localIndex idof = 0; idof < m_numDofPerCell; ++idof )
        {
          real64 const val = localRhs[lid] / volume[ei];
          localSum += val * val;
        }
      }
    } );

    localResidualNorm += localSum.get();
  } );

  // compute global residual norm
  real64 const globalResidualNorm = MpiWrapper::Sum( localResidualNorm, MPI_COMM_GEOSX );
  return sqrt( globalResidualNorm );
}

void ProppantTransport::ApplySystemSolution( DofManager const & dofManager,
                                             arrayView1d< real64 const > const & localSolution,
                                             real64 const scalingFactor,
                                             DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::proppantConcentrationString,
                               viewKeyStruct::deltaProppantConcentrationString,
                               scalingFactor,
                               0, 1 );


  if( m_numDofPerCell > 1 )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::proppantConcentrationString,
                                 viewKeyStruct::deltaComponentConcentrationString,
                                 scalingFactor,
                                 1, m_numDofPerCell );
  }

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaProppantConcentrationString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaComponentConcentrationString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    //UpdateState( subRegion );
    UpdateComponentDensity( subRegion, targetIndex );
  } );

}

void ProppantTransport::SolveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void ProppantTransport::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dProppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );
    arrayView2d< real64 > const & dComponentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      dProppantConc[ei] = 0.0;
      for( localIndex c = 0; c < NC; ++c )
      {
        dComponentConc[ei][c] = 0.0;
      }
    } );

    UpdateState( subRegion, targetIndex );
  } );

}

void ProppantTransport::ResetViews( MeshLevel & mesh )
{
  FlowSolverBase::ResetViews( mesh );
  ElementRegionManager * const elemManager = mesh.getElemManager();

  m_pressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaPressureString );

  m_proppantConcentration =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantConcentrationString );
  m_deltaProppantConcentration =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString );

  m_cellBasedFlux =
    elemManager->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor > >( viewKeyStruct::cellBasedFluxString );

  m_proppantPackVolumeFraction =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString );

  m_proppantExcessPackVolume =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString );

  m_proppantLiftFlux =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantLiftFluxString );

  m_isProppantMobile =
    elemManager->ConstructViewAccessor< array1d< integer >, arrayView1d< integer > >( viewKeyStruct::isProppantMobileString );

  m_isInterfaceElement =
    elemManager->ConstructViewAccessor< array1d< integer >, arrayView1d< integer > >( viewKeyStruct::isInterfaceElementString );

  m_isProppantBoundaryElement =
    elemManager->ConstructViewAccessor< array1d< integer >, arrayView1d< integer > >( viewKeyStruct::isProppantBoundaryString );

  {
    using keys = SlurryFluidBase::viewKeyStruct;

    m_density =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::densityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dDensity_dPressure =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dDens_dPresString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dDensity_dProppantConcentration =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dDens_dProppantConcString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dDensity_dComponentConcentration =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dDens_dCompConcString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_componentDensity =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::componentDensityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dComponentDensity_dPressure =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dCompDens_dPresString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dComponentDensity_dComponentConcentration =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dCompDens_dCompConcString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_fluidDensity =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::fluidDensityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dFluidDensity_dPressure =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dFluidDens_dPresString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dFluidDensity_dComponentConcentration =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dFluidDens_dCompConcString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_fluidViscosity =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::fluidViscosityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_viscosity =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::viscosityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dViscosity_dPressure =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dVisc_dPresString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dViscosity_dProppantConcentration =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dVisc_dProppantConcString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dViscosity_dComponentConcentration =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dVisc_dCompConcString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
  }
  {
    using keys = ParticleFluidBase::viewKeyStruct;

    m_settlingFactor =
      elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::settlingFactorString,
                                                                                              targetRegionNames(),
                                                                                              proppantModelNames() );

    m_dSettlingFactor_dPressure =
      elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::dSettlingFactor_dPressureString,
                                                                                              targetRegionNames(),
                                                                                              proppantModelNames() );

    m_dSettlingFactor_dProppantConcentration =
      elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::dSettlingFactor_dProppantConcentrationString,
                                                                                              targetRegionNames(),
                                                                                              proppantModelNames() );

    m_dSettlingFactor_dComponentConcentration =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dSettlingFactor_dComponentConcentrationString,
                                                                                              targetRegionNames(),
                                                                                              proppantModelNames() );

    m_collisionFactor =
      elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::collisionFactorString,
                                                                                              targetRegionNames(),
                                                                                              proppantModelNames() );

    m_dCollisionFactor_dProppantConcentration =
      elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::dCollisionFactor_dProppantConcentrationString,
                                                                                              targetRegionNames(),
                                                                                              proppantModelNames() );
  }

  m_transTMultiplier =
    elemManager->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor > >( viewKeyStruct::transTMultiplierString );
}



void ProppantTransport::UpdateCellBasedFlux( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = *fvManager.getFluxApproximation( m_discretizationName );

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & pres               = m_pressure.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & gravCoef           = m_gravCoef.toViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dens               = m_density.toViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & visc               = m_viscosity.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & aperture           = m_elementAperture.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf     = m_proppantPackVolumeFraction.toViewConst();
  FluxKernel::ElementViewConst< arrayView1d< R1Tensor const > > const & transTMultiplier = m_transTMultiplier.toViewConst();

  FluxKernel::ElementView< arrayView1d< R1Tensor > > & cellBasedFlux  = m_cellBasedFlux.toView();

  fluxApprox.forAllStencils( [&]( auto const & stencil )
  {

    FluxKernel::LaunchCellBasedFluxCalculation( stencil,
                                                transTMultiplier,
                                                m_downVector,
                                                pres,
                                                gravCoef,
                                                dens,
                                                visc,
                                                aperture,
                                                proppantPackVf,
                                                cellBasedFlux );
  } );


  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::cellBasedFluxString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );

}

void ProppantTransport::UpdateProppantPackVolume( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain )
{

  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = *fvManager.getFluxApproximation( m_discretizationName );

  ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & conc = m_proppantConcentration.toView();
  ProppantPackVolumeKernel::ElementView< arrayView1d< real64 const > > const & settlingFactor = m_settlingFactor.toViewConst();
  ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & density = m_density.toViewConst();
  ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & fluidDensity = m_fluidDensity.toViewConst();
  ProppantPackVolumeKernel::ElementView< arrayView2d< real64 const > > const & fluidViscosity = m_fluidViscosity.toViewConst();

  ProppantPackVolumeKernel::ElementView< arrayView1d< integer > > const & isProppantMobile = m_isProppantMobile.toView();
  ProppantPackVolumeKernel::ElementViewConst< arrayView1d< integer > > const & isProppantBoundaryElement = m_isProppantBoundaryElement.toViewConst();
  ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantPackVf  = m_proppantPackVolumeFraction.toView();
  ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantExcessPackV  = m_proppantExcessPackVolume.toView();
  ProppantPackVolumeKernel::ElementView< arrayView1d< integer > > const & isInterfaceElement = m_isInterfaceElement.toView();
  ProppantPackVolumeKernel::ElementView< arrayView1d< R1Tensor const > > const & cellBasedFlux  = m_cellBasedFlux.toViewConst();
  ProppantPackVolumeKernel::ElementView< arrayView1d< real64 > > const & proppantLiftFlux  = m_proppantLiftFlux.toView();
  ProppantPackVolumeKernel::ElementViewConst< arrayView1d< real64 const > > const & volume  = m_volume.toViewConst();
  ProppantPackVolumeKernel::ElementViewConst< arrayView1d< real64 const > > const & aperture  = m_elementAperture.toViewConst();
  ProppantPackVolumeKernel::ElementViewConst< arrayView1d< integer const > > const & elemGhostRank = m_elemGhostRank.toViewConst();

  fluxApprox.forAllStencils( [&]( auto const & stencil )
  {
    ProppantPackVolumeKernel::LaunchProppantPackVolumeCalculation( stencil,
                                                                   dt,
                                                                   m_proppantDensity,
                                                                   m_proppantDiameter,
                                                                   m_maxProppantConcentration,
                                                                   m_downVector,
                                                                   m_criticalShieldsNumber,
                                                                   m_frictionCoefficient,
                                                                   conc,
                                                                   settlingFactor,
                                                                   density,
                                                                   fluidDensity,
                                                                   fluidViscosity,
                                                                   isProppantMobile,
                                                                   isProppantBoundaryElement,
                                                                   proppantPackVf,
                                                                   proppantExcessPackV,
                                                                   aperture,
                                                                   volume,
                                                                   elemGhostRank,
                                                                   cellBasedFlux,
                                                                   proppantLiftFlux );
  } );

  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantConcentrationString ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantPackVolumeFractionString ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantExcessPackVolumeString ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantLiftFluxString ) );

    CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );
  }

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateProppantMobility( subRegion );
  } );


  fluxApprox.forAllStencils( [&]( auto const & stencil )
  {
    ProppantPackVolumeKernel::LaunchProppantPackVolumeUpdate( stencil,
                                                              m_downVector,
                                                              m_maxProppantConcentration,
                                                              conc,
                                                              isProppantMobile,
                                                              proppantPackVf,
                                                              proppantExcessPackV );
  } );


  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantConcentrationString ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantPackVolumeFractionString ) );

    CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );
  }

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateProppantMobility( subRegion );
  } );



  fluxApprox.forAllStencils( [&]( auto const & stencil )
  {
    ProppantPackVolumeKernel::LaunchInterfaceElementUpdate( stencil,
                                                            m_downVector,
                                                            isProppantMobile,
                                                            isInterfaceElement );
  } );

  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( string( viewKeyStruct::isInterfaceElementString ) );

    CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );
  }

  // update poroMultiplier and transTMultiplier

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {

    arrayView1d< real64 const > const & proppantPackVfNew =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString );
    arrayView1d< real64 const > const & aperture0 =
      subRegion.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementApertureString );

    arrayView1d< real64 > const & poroMultiplier =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::poroMultiplierString );
    arrayView1d< R1Tensor > const & transTMultiplier =
      subRegion.getReference< array1d< R1Tensor > >( viewKeyStruct::transTMultiplierString );

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
    {

      poroMultiplier[ei] = 1.0 - m_maxProppantConcentration * proppantPackVfNew[ei];

      //K0 horizontal

      transTMultiplier[ei][0] = ( 1.0 - proppantPackVfNew[ei] )
                                + proppantPackVfNew[ei] * m_proppantPackPermeability * 12.0 / ( aperture0[ei] * aperture0[ei] );

      //K1 vertical

      transTMultiplier[ei][1] = 1.0 /
                                ( proppantPackVfNew[ei] * aperture0[ei] * aperture0[ei] / 12.0 / m_proppantPackPermeability
                                  + ( 1.0 - proppantPackVfNew[ei] ) );

    } );
  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, ProppantTransport, std::string const &, Group * const )
} /* namespace geosx */
