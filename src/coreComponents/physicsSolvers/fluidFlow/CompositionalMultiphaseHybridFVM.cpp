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
 * @file CompositionalMultiphaseHybridFVM.cpp
 */

#include "CompositionalMultiphaseHybridFVM.hpp"

#include "mesh/DomainPartition.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "finiteVolume/MimeticInnerProductDispatch.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CompositionalMultiphaseHybridFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/AccumulationKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ResidualNormKernel.hpp"

/**
 * @namespace the geos namespace that encapsulates the majority of the code
 */
namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace isothermalCompositionalMultiphaseBaseKernels;
using namespace compositionalMultiphaseHybridFVMKernels;
using namespace mimeticInnerProduct;

CompositionalMultiphaseHybridFVM::CompositionalMultiphaseHybridFVM( const std::string & name,
                                                                    Group * const parent ):
  CompositionalMultiphaseBase( name, parent ),
  m_lengthTolerance( 0 )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM;
}

void CompositionalMultiphaseHybridFVM::registerDataOnMesh( Group & meshBodies )
{
  GEOS_MARK_FUNCTION;

  // 1) Register the elem-centered data
  CompositionalMultiphaseBase::registerDataOnMesh( meshBodies );

  // 2) Register the face data
  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getBaseDiscretization();

    FaceManager & faceManager = meshLevel.getFaceManager();

    // primary variables: face pressure changes

    faceManager.registerField< fields::flow::facePressure_n >( getName() );

    // auxiliary data for the buoyancy coefficient
    faceManager.registerField< fields::flow::mimGravityCoefficient >( getName() );
  } );
}

void CompositionalMultiphaseHybridFVM::initializePreSubGroups()
{
  CompositionalMultiphaseBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  GEOS_THROW_IF( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ),
                 getCatalogName() << " " << getDataContext() <<
                 ": the HybridMimeticDiscretization must be selected with CompositionalMultiphaseHybridFVM",
                 InputError );

  GEOS_THROW_IF( m_hasCapPressure,
                 getCatalogName() << " " << getDataContext() <<
                 ": capillary pressure is not yet supported by CompositionalMultiphaseHybridFVM",
                 InputError );
}

void CompositionalMultiphaseHybridFVM::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );
  if( dynamicCast< QuasiRTInnerProduct const * >( &mimeticInnerProductBase )  ||
      dynamicCast< QuasiTPFAInnerProduct const * >( &mimeticInnerProductBase )  ||
      dynamicCast< SimpleInnerProduct const * >( &mimeticInnerProductBase ) )
  {
    GEOS_ERROR( getCatalogName() << " " << getDataContext() <<
                "The QuasiRT, QuasiTPFA, and Simple inner products are only available in SinglePhaseHybridFVM" );
  }

  m_lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * 1e-8;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    CompositionalMultiphaseBase::initializePostInitialConditionsPreSubGroups();

    // in the flux kernel, we need to make sure that we act only on the target regions
    // for that, we need the following region filter
    for( string const & regionName : regionNames )
    {
      m_regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
    }

    // check that multipliers are stricly larger than 0, which would work with SinglePhaseFVM, but not with SinglePhaseHybridFVM.
    // To deal with a 0 multiplier, we would just have to skip the corresponding face in the FluxKernel
    arrayView1d< real64 const > const & transMultiplier = faceManager.getField< fields::flow::transMultiplier >();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );
    forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      minVal.min( transMultiplier[iface] );
    } );

    GEOS_THROW_IF( minVal.get() <= 0.0,
                   getCatalogName() << " " << getDataContext() <<
                   ": the transmissibility multipliers used in SinglePhaseHybridFVM must strictly larger than 0.0",
                   std::runtime_error );

    FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
    fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
    {
      GEOS_LOG_RANK_0( getCatalogName() << " " << getDataContext() << ": An aquifer boundary condition named " <<
                       bc.getName() << " was requested in the XML file. \n" <<
                       "This type of boundary condition is not yet supported by CompositionalMultiphaseHybridFVM and will be ignored" );
    } );
  } );

}

void CompositionalMultiphaseHybridFVM::precomputeData( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  FlowSolverBase::precomputeData( mesh, regionNames );

  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  array1d< RAJA::ReduceSum< serialReduce, real64 > > mimFaceGravCoefNumerator;
  array1d< RAJA::ReduceSum< serialReduce, real64 > > mimFaceGravCoefDenominator;
  mimFaceGravCoefNumerator.resize( faceManager.size() );
  mimFaceGravCoefDenominator.resize( faceManager.size() );

  // node data

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // face data

  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getField< fields::flow::transMultiplier >();

  arrayView1d< real64 > const mimFaceGravCoef =
    faceManager.getField< fields::flow::mimGravityCoefficient >();

  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

  real64 const lengthTolerance = m_lengthTolerance;

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & subRegion )
  {
    arrayView2d< real64 const > const & elemCenter =
      subRegion.template getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
    string & permModelName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
    arrayView3d< real64 const > const & elemPerm =
      getConstitutiveModel< PermeabilityBase >( subRegion, permModelName ).permeability();
    arrayView1d< real64 const > const elemGravCoef =
      subRegion.template getReference< array1d< real64 > >( fields::flow::gravityCoefficient::key() );
    arrayView1d< real64 const > const & elemVolume = subRegion.getElementVolume();
    arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

    // here we precompute some quantities (mimFaceFracCoef) used in the FluxKernel to assemble the one-sided gravity term in the transport
    // scheme
    // This one-sided gravity term is currently always treated with TPFA, as in MRST.
    // In the future, I will change that (here and in the FluxKernel) to have a consistent inner product for the gravity term as well
    compositionalMultiphaseHybridFVMKernels::
      simpleKernelLaunchSelector< PrecomputeKernel,
                                  mimeticInnerProduct::TPFAInnerProduct >( subRegion.numFacesPerElement(),
                                                                           subRegion.size(),
                                                                           faceManager.size(),
                                                                           nodePosition,
                                                                           faceToNodes,
                                                                           elemCenter,
                                                                           elemVolume,
                                                                           elemPerm,
                                                                           elemGravCoef,
                                                                           elemToFaces,
                                                                           transMultiplier,
                                                                           lengthTolerance,
                                                                           mimFaceGravCoefNumerator.toView(),
                                                                           mimFaceGravCoefDenominator.toView(),
                                                                           mimFaceGravCoef );

  } );

}

void CompositionalMultiphaseHybridFVM::implicitStepSetup( real64 const & time_n,
                                                          real64 const & dt,
                                                          DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // setup the elem-centered fields
  CompositionalMultiphaseBase::implicitStepSetup( time_n, dt, domain );

  // setup the face fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 > const & facePres_n =
      faceManager.getField< fields::flow::facePressure_n >();
    arrayView1d< real64 const > const & facePres =
      faceManager.getField< fields::flow::facePressure >();
    facePres_n.setValues< parallelDevicePolicy<> >( facePres );
  } );

}


void CompositionalMultiphaseHybridFVM::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                  DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Face );

  // this call with instruct GEOS to reorder the dof numbers
  //dofManager.setLocalReorderingType( viewKeyStruct::elemDofFieldString(),
  //                                   DofManager::LocalReorderingType::ReverseCutHillMcKee );

  // for the volume balance equation, disable global coupling
  // this equation is purely local (not coupled to neighbors or other physics)
  dofManager.disableGlobalCouplingForEquation( viewKeyStruct::elemDofFieldString(),
                                               m_numComponents );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::faceDofFieldString(),
                       FieldLocation::Face,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::faceDofFieldString(),
                          viewKeyStruct::faceDofFieldString(),
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::faceDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );

}


void CompositionalMultiphaseHybridFVM::assembleFluxTerms( real64 const dt,
                                                          DomainPartition const & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // node data (for transmissibility computation)

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

    // face data

    // get the face-based DOF numbers for the assembly
    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

    // get the element dof numbers for the assembly
    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    // get the face-centered pressures
    arrayView1d< real64 const > const & facePres =
      faceManager.getField< fields::flow::facePressure >();

    // get the face-centered depth
    arrayView1d< real64 const > const & faceGravCoef =
      faceManager.getField< fields::flow::gravityCoefficient >();
    arrayView1d< real64 const > const & mimFaceGravCoef =
      faceManager.getField< fields::flow::mimGravityCoefficient >();

    // get the face-centered transMultiplier
    arrayView1d< real64 const > const & transMultiplier =
      faceManager.getField< fields::flow::transMultiplier >();

    // get the face-to-nodes connectivity for the transmissibility calculation
    ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

    arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList().toViewConst();
    arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList().toViewConst();
    arrayView2d< localIndex const > const & elemList          = faceManager.elementList().toViewConst();


    // tolerance for transmissibility calculation
    real64 const lengthTolerance = m_lengthTolerance;

    FluxKernel::CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getName() );
    FluxKernel::MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getName() );

    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames,
                                                                                [&]( localIndex const,
                                                                                     localIndex const er,
                                                                                     localIndex const esr,
                                                                                     ElementRegionBase const &,
                                                                                     CellElementSubRegion const & subRegion )
    {
      PermeabilityBase const & permeabilityModel =
        getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() ) );

      mimeticInnerProductReducedDispatch( mimeticInnerProductBase,
                                          [&] ( auto const mimeticInnerProduct )
      {
        using IP_TYPE = TYPEOFREF( mimeticInnerProduct );
        kernelLaunchSelector< FluxKernel,
                              IP_TYPE >( subRegion.numFacesPerElement(),
                                         m_numComponents, m_numPhases,
                                         er, esr, subRegion,
                                         permeabilityModel,
                                         m_regionFilter.toViewConst(),
                                         nodePosition,
                                         elemRegionList,
                                         elemSubRegionList,
                                         elemList,
                                         faceToNodes,
                                         faceDofNumber,
                                         faceGhostRank,
                                         facePres,
                                         faceGravCoef,
                                         mimFaceGravCoef,
                                         transMultiplier,
                                         compFlowAccessors.get( fields::flow::phaseMobility{} ),
                                         compFlowAccessors.get( fields::flow::dPhaseMobility{} ),
                                         compFlowAccessors.get( fields::flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::phaseDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::dPhaseDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::phaseMassDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::dPhaseMassDensity{} ),
                                         multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ),
                                         multiFluidAccessors.get( fields::multifluid::dPhaseCompFraction{} ),
                                         elemDofNumber.toNestedViewConst(),
                                         dofManager.rankOffset(),
                                         lengthTolerance,
                                         dt,
                                         m_useTotalMassEquation,
                                         localMatrix,
                                         localRhs );

      } );
    } );

  } );
}

void CompositionalMultiphaseHybridFVM::assembleStabilizedFluxTerms( real64 const dt,
                                                                    DomainPartition const & domain,
                                                                    DofManager const & dofManager,
                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                    arrayView1d< real64 > const & localRhs ) const
{
  // stab not implemented
  GEOS_UNUSED_VAR( dt, domain, dofManager, localMatrix, localRhs );
  GEOS_ERROR( "Stabilized flux not available for this flow solver" );
}

real64 CompositionalMultiphaseHybridFVM::scalingForSystemSolution( DomainPartition & domain,
                                                                   DofManager const & dofManager,
                                                                   arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  real64 scalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::flow::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();
      auto const subRegionData =
        isothermalCompositionalMultiphaseBaseKernels::
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
                                                     dofKey,
                                                     subRegion,
                                                     localSolution );

      scalingFactor = std::min( scalingFactor, subRegionData.localMinVal );
    } );

    FaceManager const & faceManager = mesh.getFaceManager();

    string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );
    real64 const maxRelativePresChange = m_maxRelativePresChange;

    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    arrayView1d< real64 const > const & facePressure =
      faceManager.getField< fields::flow::facePressure >();
    globalIndex const rankOffset = dofManager.rankOffset();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minFaceVal( 1.0 );
    forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      if( faceGhostRank[iface] < 0 && faceDofNumber[iface] >= 0 )
      {
        real64 const facePres = facePressure[iface];
        real64 const absPresChange = LvArray::math::abs( localSolution[faceDofNumber[iface] - rankOffset] );
        if( facePres > isothermalCompositionalMultiphaseBaseKernels::minDensForDivision )
        {
          real64 const relativePresChange = absPresChange / facePres;
          if( relativePresChange > maxRelativePresChange )
          {
            minFaceVal.min( maxRelativePresChange / relativePresChange );
          }
        }
      }
    } );

    if( minFaceVal.get() < scalingFactor )
    {
      scalingFactor = minFaceVal.get();
    }
  } );
  return LvArray::math::max( MpiWrapper::min( scalingFactor ), m_minScalingFactor );
}


bool CompositionalMultiphaseHybridFVM::checkSystemSolution( DomainPartition & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer localCheck = 1;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    // check cell-centered variables for each region
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure =
        subRegion.getField< fields::flow::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens =
        subRegion.getField< fields::flow::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::flow::temperatureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::flow::globalCompDensityScalingFactor >();
      // check that pressure and component densities are non-negative
      auto const subRegionData =
        isothermalCompositionalMultiphaseBaseKernels::
          SolutionCheckKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                     m_allowNegativePressure,
                                                     compositionalMultiphaseUtilities::ScalingType::Global,
                                                     scalingFactor,
                                                     pressure,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     elemDofKey,
                                                     subRegion,
                                                     localSolution );

      localCheck = std::min( localCheck, subRegionData.localMinVal );
    } );
  } );

  return MpiWrapper::min( localCheck );
}

void CompositionalMultiphaseHybridFVM::applyBoundaryConditions( real64 const time_n,
                                                                real64 const dt,
                                                                DomainPartition & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  CompositionalMultiphaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // TODO: implement face boundary conditions here
}

void CompositionalMultiphaseHybridFVM::applyAquiferBC( real64 const time,
                                                       real64 const dt,
                                                       DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time, dt, dofManager, domain, localMatrix, localRhs );

}

void CompositionalMultiphaseHybridFVM::saveAquiferConvergedState( real64 const & time,
                                                                  real64 const & dt,
                                                                  DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time, dt, domain );
}


real64 CompositionalMultiphaseHybridFVM::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                                real64 const & dt,
                                                                DomainPartition const & domain,
                                                                DofManager const & dofManager,
                                                                arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 localResidualNorm = 0.0;
  real64 localResidualNormalizer = 0.0;

  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  // local residual
  globalIndex const rankOffset = dofManager.rankOffset();
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::faceDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // here we compute the cell-centered residual norm in the derived class
    // to avoid duplicating a synchronization point

    // get a view into local residual vector

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[2]{};
      real64 subRegionResidualNormalizer[2]{};

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // step 1.1: compute the norm in the subRegion

      isothermalCompositionalMultiphaseBaseKernels::
        ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( normType,
                                                   numFluidComponents(),
                                                   rankOffset,
                                                   elemDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   fluid,
                                                   solid,
                                                   m_nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm,
                                                   subRegionResidualNormalizer );

      // step 1.2: reduction across meshBodies/regions/subRegions

      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        // take max between mass and volume residual
        subRegionResidualNorm[0] = LvArray::math::max( subRegionResidualNorm[0], subRegionResidualNorm[1] );
        if( subRegionResidualNorm[0] > localResidualNorm )
        {
          localResidualNorm = subRegionResidualNorm[0];
        }
      }
      else
      {
        // sum up mass and volume residual
        subRegionResidualNorm[0] = subRegionResidualNorm[0] + subRegionResidualNorm[1];
        localResidualNorm += subRegionResidualNorm[0];
        localResidualNormalizer += subRegionResidualNormalizer[0];
      }
    } );

    // step 2: compute the residual for the face-based constraints

    real64 faceResidualNorm[1]{};
    real64 faceResidualNormalizer[1]{};

    // step 2.1: compute the norm for the local faces

    compositionalMultiphaseHybridFVMKernels::
      ResidualNormKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( normType,
                                                 rankOffset,
                                                 faceDofKey,
                                                 localRhs,
                                                 m_regionFilter.toViewConst(),
                                                 getName(),
                                                 elemManager,
                                                 faceManager,
                                                 dt,
                                                 m_nonlinearSolverParameters.m_minNormalizer,
                                                 faceResidualNorm,
                                                 faceResidualNormalizer );

    // step 2.2: reduction across meshBodies/regions/subRegions

    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      if( faceResidualNorm[0] > localResidualNorm )
      {
        localResidualNorm = faceResidualNorm[0];
      }
    }
    else
    {
      localResidualNorm += faceResidualNorm[0];
      localResidualNormalizer += faceResidualNormalizer[0];
    }
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  if( normType == physicsSolverBaseKernels::NormType::Linf )
  {
    physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm, residualNorm );
  }
  else
  {
    physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm, localResidualNormalizer, residualNorm );
  }

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residualNorm );
  }

  return residualNorm;
}

void CompositionalMultiphaseHybridFVM::applySystemSolution( DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localSolution,
                                                            real64 const scalingFactor,
                                                            real64 const GEOS_UNUSED_PARAM( dt ),
                                                            DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // 1. apply the elem-based update
  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               fields::flow::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               fields::flow::globalCompDensity::key(),
                               scalingFactor,
                               ~pressureMask );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::faceDofFieldString(),
                               fields::flow::facePressure::key(),
                               scalingFactor );

  // 3. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    {
      fieldsToBeSync.addElementFields( { fields::flow::pressure::key(),
                                         fields::flow::globalCompDensity::key() },
                                       regionNames );

      fieldsToBeSync.addFields( FieldLocation::Face, { fields::flow::facePressure::key() } );
    };

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}


void CompositionalMultiphaseHybridFVM::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // 1. Reset the cell-centered fields
  CompositionalMultiphaseBase::resetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const & facePres_n =
      faceManager.getField< fields::flow::facePressure_n >();
    arrayView1d< real64 > const & facePres =
      faceManager.getField< fields::flow::facePressure >();
    facePres.setValues< parallelDevicePolicy<> >( facePres_n );
  } );
}

void CompositionalMultiphaseHybridFVM::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  MultiFluidBase const & fluid =
    getConstitutiveModel< MultiFluidBase >( dataGroup,
                                            dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  RelativePermeabilityBase const & relperm =
    getConstitutiveModel< RelativePermeabilityBase >( dataGroup,
                                                      dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() ) );

  compositionalMultiphaseHybridFVMKernels::
    PhaseMobilityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               dataGroup,
                                               fluid,
                                               relperm );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, CompositionalMultiphaseHybridFVM, std::string const &, Group * const )
} /* namespace geos */
