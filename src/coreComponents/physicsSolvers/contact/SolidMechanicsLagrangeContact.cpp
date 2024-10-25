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
 * @file SolidMechanicsLagrangeContact.cpp
 *
 */
#define GEOS_DISPATCH_VEM
#include "SolidMechanicsLagrangeContact.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/FrictionSelector.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp" // needed to register pressure(_n)
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "physicsSolvers/contact/LogLevelsInfo.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/solvers/PreconditionerJacobi.hpp"
#include "linearAlgebra/solvers/PreconditionerBlockJacobi.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "finiteElement/elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp"
#include "finiteElement/elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;
using namespace finiteElement;
const localIndex geos::SolidMechanicsLagrangeContact::m_maxFaceNodes = 11;

SolidMechanicsLagrangeContact::SolidMechanicsLagrangeContact( const string & name,
                                                              Group * const parent ):
  ContactSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::stabilizationNameString(), &m_stabilizationName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the stabilization to use in the lagrangian contact solver" );

  registerWrapper( viewKeyStruct::stabilizationScalingCoefficientString(), &m_stabilitzationScalingCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "It be used to increase the scale of the stabilization entries. A value < 1.0 results in larger entries in the stabilization matrix." );

  addLogLevel< logInfo::Configuration >();
}

void SolidMechanicsLagrangeContact::postInputInitialization()
{
  ContactSolverBase::postInputInitialization();

  setMGRStrategy();
}

void SolidMechanicsLagrangeContact::setMGRStrategy()
{
  LinearSolverParameters & linearSolverParameters = m_linearSolverParameters.get();

  if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    return;

  linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::lagrangianContactMechanics;
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "{}: MGR strategy set to {}", getName(),
                                      EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy )));
  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
}

void SolidMechanicsLagrangeContact::registerDataOnMesh( Group & meshBodies )
{
  ContactSolverBase::registerDataOnMesh( meshBodies );

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::rotationMatrixString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the rotation matrices on the fracture." ).
        reference().resizeDimension< 1, 2 >( 3, 3 );

      subRegion.registerField< fields::contact::deltaTraction >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the normal traction tolerance." );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the normal displacement tolerance." );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::slidingToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName()).
        setDescription( "An array that holds the sliding tolerance." );

      // Needed just because SurfaceGenerator initialize the field "pressure" (NEEDED!!!)
      // It is used in "TwoPointFluxApproximation.cpp", called by "SurfaceGenerator.cpp"
      subRegion.registerField< flow::pressure >( getName() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName());
      subRegion.registerField< flow::pressure_n >( getName() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( this->getName());

    } );

    forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                      MeshLevel & mesh,
                                                      arrayView1d< string const > const & )
    {
      FaceManager & faceManager = mesh.getFaceManager();

      faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::transMultiplierString() ).
        setApplyDefaultValue( 1.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the permeability transmissibility multipliers" );
    } );

  } );
}

void SolidMechanicsLagrangeContact::initializePreSubGroups()
{
  ContactSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_stabilizationName ) )
  {

    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_stabilizationName );
    fluxApprox.addFieldName( contact::traction::key() );
    fluxApprox.setCoeffName( "penaltyStiffness" );


    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            SurfaceElementRegion const & region )
      {
        array1d< string > & stencilTargetRegions = fluxApprox.targetRegions( meshBodyName );
        stencilTargetRegions.emplace_back( region.getName() );
      } );
    } );
  }

}

void SolidMechanicsLagrangeContact::setupSystem( DomainPartition & domain,
                                                 DofManager & dofManager,
                                                 CRSMatrix< real64, globalIndex > & localMatrix,
                                                 ParallelVector & rhs,
                                                 ParallelVector & solution,
                                                 bool const GEOS_UNUSED_PARAM( setSparsity ) )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, true ); // "true" is to force setSparsity

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner( domain );
  }
}

void SolidMechanicsLagrangeContact::implicitStepSetup( real64 const & time_n,
                                                       real64 const & dt,
                                                       DomainPartition & domain )
{
  computeRotationMatrices( domain );
  computeTolerances( domain );
  computeFaceDisplacementJump( domain );

  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );
}

void SolidMechanicsLagrangeContact::implicitStepComplete( real64 const & time_n,
                                                          real64 const & dt,
                                                          DomainPartition & domain )
{
  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 > const & deltaTraction = subRegion.getField< contact::deltaTraction >();
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 > const & oldDispJump = subRegion.getField< contact::oldDispJump >();
      arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< integer > const & oldFractureState = subRegion.getField< contact::oldFractureState >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          deltaTraction[kfe][i] = 0.0;
          oldDispJump[kfe][i] = dispJump[kfe][i];
        }
        oldFractureState[kfe] = fractureState[kfe];
      } );
    } );

    // Need a synchronization of deltaTraction as will be used in AssembleStabilization
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { contact::deltaTraction::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );

  } );
}

SolidMechanicsLagrangeContact::~SolidMechanicsLagrangeContact()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsLagrangeContact::computeTolerances( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  real64 minNormalTractionTolerance( 1e10 );
  real64 maxNormalTractionTolerance( -1e10 );
  real64 minNormalDisplacementTolerance( 1e10 );
  real64 maxNormalDisplacementTolerance( -1e10 );
  real64 minSlidingTolerance( 1e10 );
  real64 maxSlidingTolerance( -1e10 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // Get the "face to element" map (valid for the entire mesh)
    FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
    arrayView2d< localIndex const > const & faceToElemRegion = faceToElem.m_toElementRegion;
    arrayView2d< localIndex const > const & faceToElemSubRegion = faceToElem.m_toElementSubRegion;
    arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex;

    // Get the volume for all elements
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

    // Bulk modulus accessor
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
      elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::bulkModulusString() );
    // Shear modulus accessor
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
      elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::shearModulusString() );

    using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
    ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
      elemManager.constructViewAccessor< CellElementSubRegion::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );
    ElementRegionManager::ElementViewConst< NodeMapViewType > const elemToNodeView = elemToNode.toNestedViewConst();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
        arrayView1d< real64 const > const & faceArea = subRegion.getElementArea().toViewConst();
        arrayView3d< real64 const > const & faceRotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
        ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

        arrayView1d< real64 > const & normalTractionTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
        arrayView1d< real64 > const & normalDisplacementTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );
        arrayView1d< real64 > const & slidingTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

        RAJA::ReduceMin< ReducePolicy< parallelHostPolicy >, real64 > minSubRegionNormalTractionTolerance( 1e10 );
        RAJA::ReduceMax< ReducePolicy< parallelHostPolicy >, real64 > maxSubRegionNormalTractionTolerance( -1e10 );
        RAJA::ReduceMin< ReducePolicy< parallelHostPolicy >, real64 > minSubRegionNormalDisplacementTolerance( 1e10 );
        RAJA::ReduceMax< ReducePolicy< parallelHostPolicy >, real64 > maxSubRegionNormalDisplacementTolerance( -1e10 );
        RAJA::ReduceMin< ReducePolicy< parallelHostPolicy >, real64 > minSubRegionSlidingTolerance( 1e10 );
        RAJA::ReduceMax< ReducePolicy< parallelHostPolicy >, real64 > maxSubRegionSlidingTolerance( -1e10 );

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            real64 const area = faceArea[kfe];
            // approximation of the stiffness along coordinate directions
            // ( first, second ) index -> ( element index, direction )
            // 1. T -> top (index 0), B -> bottom (index 1)
            // 2. the coordinate direction (x, y, z)
            real64 stiffDiagApprox[ 2 ][ 3 ];
            real64 averageYoungModulus = 0.0;
            real64 averageConstrainedModulus = 0.0;
            real64 averageBoxSize0 = 0.0;

            for( localIndex i = 0; i < elemsToFaces.sizeOfArray( kfe ); ++i )
            {
              localIndex const faceIndex = elemsToFaces[kfe][i];
              localIndex const er = faceToElemRegion[faceIndex][0];
              localIndex const esr = faceToElemSubRegion[faceIndex][0];
              localIndex const ei = faceToElemIndex[faceIndex][0];

              real64 const volume = elemVolume[er][esr][ei];

              // Get the "element to node" map for the specific region/subregion
              NodeMapViewType const & cellElemsToNodes = elemToNodeView[er][esr];
              localIndex const numNodesPerElem = cellElemsToNodes.size( 1 );

              // Compute the box size
              real64 maxSize[3];
              real64 minSize[3];
              for( localIndex j = 0; j < 3; ++j )
              {
                maxSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
                minSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
              }
              for( localIndex a = 1; a < numNodesPerElem; ++a )
              {
                for( localIndex j = 0; j < 3; ++j )
                {
                  maxSize[j] = fmax( maxSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
                  minSize[j] = fmin( minSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
                }
              }

              real64 boxSize[3];
              for( localIndex j = 0; j < 3; ++j )
              {
                boxSize[j] = maxSize[j] - minSize[j];
              }

              // Get linear elastic isotropic constitutive parameters for the element
              real64 const K = bulkModulus[er][esr][ei];
              real64 const G = shearModulus[er][esr][ei];
              real64 const E = 9.0 * K * G / ( 3.0 * K + G );
              real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );
              real64 const M = K + 4.0 / 3.0 * G;

              // Combine E and nu to obtain a stiffness approximation (like it was an hexahedron)
              for( localIndex j = 0; j < 3; ++j )
              {
                stiffDiagApprox[ i ][ j ] = E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 4.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] );
              }

              averageYoungModulus += 0.5*E;
              averageConstrainedModulus += 0.5*M;
              averageBoxSize0 += 0.5*boxSize[0];
            }

            // Average the stiffness and compute the inverse
            real64 invStiffApprox[ 3 ][ 3 ] = { { 0 } };
            for( localIndex j = 0; j < 3; ++j )
            {
              invStiffApprox[ j ][ j ] = ( stiffDiagApprox[ 0 ][ j ] + stiffDiagApprox[ 1 ][ j ] ) / ( stiffDiagApprox[ 0 ][ j ] * stiffDiagApprox[ 1 ][ j ] );
            }

            // Rotate in the local reference system, computing R^T * (invK) * R
            real64 temp[ 3 ][ 3 ];
            LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( temp, faceRotationMatrix[ kfe ], invStiffApprox );
            real64 rotatedInvStiffApprox[ 3 ][ 3 ];
            LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, faceRotationMatrix[ kfe ] );
            LvArray::tensorOps::scale< 3, 3 >( rotatedInvStiffApprox, area );

            // Finally, compute tolerances for the given fracture element

            normalDisplacementTolerance[kfe] = rotatedInvStiffApprox[ 0 ][ 0 ] * averageYoungModulus / 2.e+7;
            minSubRegionNormalDisplacementTolerance.min( normalDisplacementTolerance[kfe] );
            maxSubRegionNormalDisplacementTolerance.max( normalDisplacementTolerance[kfe] );

            slidingTolerance[kfe] = sqrt( rotatedInvStiffApprox[ 1 ][ 1 ] * rotatedInvStiffApprox[ 1 ][ 1 ] +
                                          rotatedInvStiffApprox[ 2 ][ 2 ] * rotatedInvStiffApprox[ 2 ][ 2 ] ) * averageYoungModulus / 2.e+7;
            minSubRegionSlidingTolerance.min( slidingTolerance[kfe] );
            maxSubRegionSlidingTolerance.max( slidingTolerance[kfe] );

            normalTractionTolerance[kfe] = 1.0 / 2.0 * averageConstrainedModulus / averageBoxSize0 * normalDisplacementTolerance[kfe];
            minSubRegionNormalTractionTolerance.min( normalTractionTolerance[kfe] );
            maxSubRegionNormalTractionTolerance.max( normalTractionTolerance[kfe] );
          }
        } );

        minNormalDisplacementTolerance = std::min( minNormalDisplacementTolerance, minSubRegionNormalDisplacementTolerance.get() );
        maxNormalDisplacementTolerance = std::max( maxNormalDisplacementTolerance, maxSubRegionNormalDisplacementTolerance.get() );
        minSlidingTolerance = std::min( minSlidingTolerance, minSubRegionSlidingTolerance.get() );
        maxSlidingTolerance = std::max( maxSlidingTolerance, maxSubRegionSlidingTolerance.get() );
        minNormalTractionTolerance = std::min( minNormalTractionTolerance, minSubRegionNormalTractionTolerance.get() );
        maxNormalTractionTolerance = std::max( maxNormalTractionTolerance, maxSubRegionNormalTractionTolerance.get() );
      }
    } );
  } );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Configuration,
                              GEOS_FMT( "{}: normal displacement tolerance = [{}, {}], sliding tolerance = [{}, {}], normal traction tolerance = [{}, {}]",
                                        this->getName(), minNormalDisplacementTolerance, maxNormalDisplacementTolerance,
                                        minSlidingTolerance, maxSlidingTolerance,
                                        minNormalTractionTolerance, maxNormalTractionTolerance ) );
}

void SolidMechanicsLagrangeContact::resetStateToBeginningOfStep( DomainPartition & domain )
{
  SolidMechanicsLagrangianFEM::resetStateToBeginningOfStep( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 > const & traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 > const & deltaTraction = subRegion.getField< contact::deltaTraction >();
      arrayView2d< real64 > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 const > const & oldDispJump = subRegion.getField< contact::oldDispJump >();

      arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< integer const > const & oldFractureState = subRegion.getField< contact::oldFractureState >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          traction[kfe][i] -= deltaTraction[kfe][i];
          deltaTraction[kfe][i] = 0.0;

          dispJump[kfe][i] = oldDispJump[kfe][i];
        }
        fractureState[kfe] = oldFractureState[kfe];
      } );
    } );
  } );
}

void SolidMechanicsLagrangeContact::computeFaceDisplacementJump( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    EdgeManager const & edgeManager = mesh.getEdgeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
    ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
    arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
    arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
    arrayView2d< real64 const > faceNormals = faceManager.faceNormal();
    arrayView1d< real64 const > faceAreas = faceManager.faceArea();

    // Get the coordinates for all nodes
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u =
      nodeManager.getField< solidMechanics::totalDisplacement >();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView3d< real64 > const &
        rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
        ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();
        arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

        arrayView2d< real64 > const & dispJump = subRegion.getField< contact::dispJump >();
        arrayView1d< real64 > const & slip = subRegion.getField< fields::contact::slip >();
        arrayView1d< real64 > const & aperture = subRegion.getField< fields::elementAperture >();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( elemsToFaces.sizeOfArray( kfe ) != 2 )
          {
            return;
          }

          // Contact constraints
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );

          stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea0;
          stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea1;
          computeFaceNodalArea( elemsToFaces[kfe][0],
                                nodePosition,
                                faceToNodeMap,
                                faceToEdgeMap,
                                edgeToNodeMap,
                                faceCenters,
                                faceNormals,
                                faceAreas,
                                nodalArea0 );

          computeFaceNodalArea( elemsToFaces[kfe][1],
                                nodePosition,
                                faceToNodeMap,
                                faceToEdgeMap,
                                edgeToNodeMap,
                                faceCenters,
                                faceNormals,
                                faceAreas,
                                nodalArea1 );

          real64 globalJumpTemp[ 3 ] = { 0 };
          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            for( localIndex i = 0; i < 3; ++i )
            {
              globalJumpTemp[ i ] +=
                ( -u[faceToNodeMap( elemsToFaces[kfe][0], a )][i] * nodalArea0[a]
                  +u[faceToNodeMap( elemsToFaces[kfe][1], a )][i] * nodalArea1[a] ) / area[kfe];
            }
          }

          real64 dispJumpTemp[ 3 ];
          LvArray::tensorOps::Ri_eq_AjiBj< 3, 3 >( dispJumpTemp, rotationMatrix[ kfe ], globalJumpTemp );
          LvArray::tensorOps::copy< 3 >( dispJump[ kfe ], dispJumpTemp );

          slip[ kfe ] = LvArray::math::sqrt( LvArray::math::square( dispJump( kfe, 1 ) ) +
                                             LvArray::math::square( dispJump( kfe, 2 ) ) );
          aperture[ kfe ] = dispJump[ kfe ][ 0 ];
        } );
      }
    } );
  } );
}

void SolidMechanicsLagrangeContact::setupDofs( DomainPartition const & domain,
                                               DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );

  // restrict coupling to fracture regions only
  map< std::pair< string, string >, array1d< string > > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    array1d< string > regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< SurfaceElementRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addField( contact::traction::key(),
                       FieldLocation::Elem,
                       3,
                       meshTargets );

  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  dofManager.addCoupling( contact::traction::key(), stabilizationMethod );
}

void SolidMechanicsLagrangeContact::assembleSystem( real64 const time,
                                                    real64 const dt,
                                                    DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  synchronizeFractureState( domain );

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  assembleContact( domain, dofManager, localMatrix, localRhs );

  // for sequential: add (fixed) pressure force contribution into residual (no derivatives)
  if( m_isFixedStressPoromechanicsUpdate )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel const & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      assembleForceResidualPressureContribution( mesh, regionNames, dofManager, localMatrix, localRhs );
    } );
  }
}

void SolidMechanicsLagrangeContact::assembleContact( DomainPartition & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    /// assemble Kut
    assembleForceResidualDerivativeWrtTraction( mesh, regionNames, dofManager, localMatrix, localRhs );
    /// assemble Ktu, Ktt blocks.
    assembleTractionResidualDerivativeWrtDisplacementAndTraction( mesh, regionNames, dofManager, localMatrix, localRhs );
    /// assemble stabilization
    assembleStabilization( mesh, domain.getNumericalMethodManager(), dofManager, localMatrix, localRhs );
  } );
}

void SolidMechanicsLagrangeContact::
  assembleForceResidualPressureContribution( MeshLevel const & mesh,
                                             arrayView1d< string const > const & regionNames,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > faceNormal = faceManager.faceNormal();
  arrayView1d< real64 const > faceAreas = faceManager.faceArea();

  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

  arrayView1d< globalIndex const > const &
  dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    arrayView1d< real64 const > const & pressure = subRegion.getReference< array1d< real64 > >( flow::pressure::key() );
    ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

      real64 Nbar[3];
      Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
      Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
      Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
      LvArray::tensorOps::normalize< 3 >( Nbar );

      globalIndex rowDOF[3 * m_maxFaceNodes];
      real64 nodeRHS[3 * m_maxFaceNodes];
      stackArray1d< real64, 3 * m_maxFaceNodes > dRdP( 3*m_maxFaceNodes );

      for( localIndex kf=0; kf<2; ++kf )
      {
        localIndex const faceIndex = elemsToFaces[kfe][kf];
        // Compute local area contribution for each node
        stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
        computeFaceNodalArea( elemsToFaces[kfe][kf],
                              nodePosition,
                              faceToNodeMap,
                              faceToEdgeMap,
                              edgeToNodeMap,
                              faceCenters,
                              faceNormal,
                              faceAreas,
                              nodalArea );

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          real64 const nodalForceMag = -( pressure[kfe] ) * nodalArea[a];
          real64 globalNodalForce[ 3 ];
          LvArray::tensorOps::scaledCopy< 3 >( globalNodalForce, Nbar, nodalForceMag );

          for( localIndex i=0; i<3; ++i )
          {
            rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + LvArray::integerConversion< globalIndex >( i );
            // Opposite sign w.r.t. theory because of minus sign in stiffness matrix definition (K < 0)
            nodeRHS[3*a+i] = +globalNodalForce[i] * pow( -1, kf );
          }
        }

        for( localIndex idof = 0; idof < numNodesPerFace * 3; ++idof )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], nodeRHS[idof] );
          }
        }
      }
    } );
  } );
}

real64 SolidMechanicsLagrangeContact::calculateResidualNorm( real64 const & time,
                                                             real64 const & dt,
                                                             DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 const solidResidual = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  real64 const contactResidual = calculateContactResidualNorm( domain, dofManager, localRhs );

  return sqrt( solidResidual * solidResidual + contactResidual * contactResidual );
}

real64 SolidMechanicsLagrangeContact::calculateContactResidualNorm( DomainPartition const & domain,
                                                                    DofManager const & dofManager,
                                                                    arrayView1d< real64 const > const & localRhs )
{
  string const & dofKey = dofManager.getKey( contact::traction::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  real64 stickResidual = 0.0;
  real64 slipResidual = 0.0;
  real64 slipNormalizer = 0.0;
  real64 openResidual = 0.0;
  real64 openNormalizer = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                        [&]( localIndex const, FaceElementSubRegion const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();

      RAJA::ReduceSum< parallelHostReduce, real64 > stickSum( 0.0 );
      RAJA::ReduceSum< parallelHostReduce, real64 > slipSum( 0.0 );
      RAJA::ReduceMax< parallelHostReduce, real64 > slipMax( 0.0 );
      RAJA::ReduceSum< parallelHostReduce, real64 > openSum( 0.0 );
      RAJA::ReduceMax< parallelHostReduce, real64 > openMax( 0.0 );
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const k )
      {
        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
          switch( fractureState[k] )
          {
            case contact::FractureState::Stick:
              {
                for( localIndex dim = 0; dim < 3; ++dim )
                {
                  real64 const norm = localRhs[localRow + dim] / area[k];
                  stickSum += norm * norm;
                }
                break;
              }
            case contact::FractureState::Slip:
            case contact::FractureState::NewSlip:
              {
                for( localIndex dim = 0; dim < 3; ++dim )
                {
                  slipSum += localRhs[localRow + dim] * localRhs[localRow + dim];
                  slipMax.max( LvArray::math::abs( localRhs[localRow + dim] ) );
                }
                break;
              }
            case contact::FractureState::Open:
              {
                for( localIndex dim = 0; dim < 3; ++dim )
                {
                  openSum += localRhs[localRow + dim] * localRhs[localRow + dim];
                  openMax.max( LvArray::math::abs( localRhs[localRow + dim] ) );
                }
                break;
              }
          }
        }
      } );

      stickResidual += stickSum.get();
      slipResidual += slipSum.get();
      slipNormalizer = LvArray::math::max( slipNormalizer, slipMax.get());
      openResidual += openSum.get();
      openNormalizer = LvArray::math::max( openNormalizer, openMax.get());
    } );
  } );

  stickResidual = MpiWrapper::sum( stickResidual );
  stickResidual = sqrt( stickResidual );

  slipResidual = MpiWrapper::sum( slipResidual );
  slipNormalizer = MpiWrapper::max( slipNormalizer );
  slipResidual = sqrt( slipResidual ) / ( slipNormalizer + 1.0 );

  openResidual = MpiWrapper::sum( openResidual );
  openNormalizer = MpiWrapper::max( openNormalizer );
  openResidual = sqrt( openResidual ) / ( openNormalizer + 1.0 );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    std::cout << GEOS_FMT( "        ( Rstick Rslip Ropen ) = ( {:15.6e} {:15.6e} {:15.6e} )", stickResidual, slipResidual, openResidual );
  }

  return sqrt( stickResidual * stickResidual + slipResidual * slipResidual + openResidual * openResidual );
}

void SolidMechanicsLagrangeContact::createPreconditioner( DomainPartition const & domain )
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    // TODO: move among inputs (xml)
    string const leadingBlockApproximation = "blockJacobi";

    LinearSolverParameters mechParams = getLinearSolverParameters();
    // Because of boundary conditions
    mechParams.isSymmetric = false;

    std::unique_ptr< BlockPreconditioner< LAInterface > > precond;
    std::unique_ptr< PreconditionerBase< LAInterface > > tracPrecond;

    if( leadingBlockApproximation == "jacobi" )
    {
      precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::LowerUpperTriangular,
                                                                        SchurComplementOption::FirstBlockDiagonal,
                                                                        BlockScalingOption::UserProvided );
      // Using GEOSX implementation of Jacobi preconditioner
      // tracPrecond = std::make_unique< PreconditionerJacobi< LAInterface > >();

      // Using LAI implementation of Jacobi preconditioner
      LinearSolverParameters tracParams;
      tracParams.preconditionerType = LinearSolverParameters::PreconditionerType::jacobi;
      tracPrecond = LAInterface::createPreconditioner( tracParams );
    }
    else if( leadingBlockApproximation == "blockJacobi" )
    {
      precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::LowerUpperTriangular,
                                                                        SchurComplementOption::FirstBlockUserDefined,
                                                                        BlockScalingOption::UserProvided );
      tracPrecond = std::make_unique< PreconditionerBlockJacobi< LAInterface > >( mechParams.dofsPerNode );
    }
    else
    {
      GEOS_ERROR( "SolidMechanicsLagrangeContact::CreatePreconditioner leadingBlockApproximation option " << leadingBlockApproximation << " not supported" );
    }

    // Preconditioner for the leading block: tracPrecond
    precond->setupBlock( 0,
                         { { contact::traction::key(), { 3, true } } },
                         std::move( tracPrecond ) );

    if( mechParams.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
    {
      if( getRigidBodyModes().empty() )
      {
        MeshLevel const & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();
        LAIHelperFunctions::computeRigidBodyModes( mesh,
                                                   m_dofManager,
                                                   { solidMechanics::totalDisplacement::key() },
                                                   getRigidBodyModes() );
      }
    }

    // Preconditioner for the Schur complement: mechPrecond
    std::unique_ptr< PreconditionerBase< LAInterface > > mechPrecond = LAInterface::createPreconditioner( mechParams, getRigidBodyModes() );
    precond->setupBlock( 1,
                         { { solidMechanics::totalDisplacement::key(), { 3, true } } },
                         std::move( mechPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void SolidMechanicsLagrangeContact::computeRotationMatrices( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();

    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             FaceElementSubRegion & subRegion )
    {
      ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

      arrayView3d< real64 > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
      arrayView2d< real64 > const unitNormal   = subRegion.getNormalVector();
      arrayView2d< real64 > const unitTangent1 = subRegion.getTangentVector1();
      arrayView2d< real64 > const unitTangent2 = subRegion.getTangentVector2();

      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        if( elemsToFaces.sizeOfArray( kfe ) != 2 )
        {
          return;
        }

        localIndex const f0 = elemsToFaces[kfe][0];
        localIndex const f1 = elemsToFaces[kfe][1];

        stackArray1d< real64, 3 > Nbar( 3 );
        Nbar[ 0 ] = faceNormal[f0][0] - faceNormal[f1][0];
        Nbar[ 1 ] = faceNormal[f0][1] - faceNormal[f1][1];
        Nbar[ 2 ] = faceNormal[f0][2] - faceNormal[f1][2];
        LvArray::tensorOps::normalize< 3 >( Nbar );

        computationalGeometry::RotationMatrix_3D( Nbar.toSliceConst(), rotationMatrix[kfe] );
        real64 const columnVector1[3] = { rotationMatrix[kfe][ 0 ][ 1 ],
                                          rotationMatrix[kfe][ 1 ][ 1 ],
                                          rotationMatrix[kfe][ 2 ][ 1 ] };

        real64 const columnVector2[3] = { rotationMatrix[kfe][ 0 ][ 2 ],
                                          rotationMatrix[kfe][ 1 ][ 2 ],
                                          rotationMatrix[kfe][ 2 ][ 2 ] };

        LvArray::tensorOps::copy< 3 >( unitNormal[kfe], Nbar );
        LvArray::tensorOps::copy< 3 >( unitTangent1[kfe], columnVector1 );
        LvArray::tensorOps::copy< 3 >( unitTangent2[kfe], columnVector2 );
      } );
    } );
  } );
}

void SolidMechanicsLagrangeContact::computeFaceIntegrals( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                                                          localIndex const (&faceToNodes)[11],
                                                          localIndex const (&faceToEdges)[11],
                                                          localIndex const & numFaceVertices,
                                                          real64 const & faceArea,
                                                          real64 const (&faceCenter)[3],
                                                          real64 const (&faceNormal)[3],
                                                          arrayView2d< localIndex const > const & edgeToNodes,
                                                          real64 const & invCellDiameter,
                                                          real64 const (&cellCenter)[3],
                                                          stackArray1d< real64, FaceManager::maxFaceNodes() > & basisIntegrals,
                                                          real64 (& threeDMonomialIntegrals)[3] ) const
{
  GEOS_MARK_FUNCTION;
  localIndex const MFN = m_maxFaceNodes; // Max number of face vertices.
  basisIntegrals.resize( numFaceVertices );
  // Rotate the face.
  //  - compute rotation matrix.
  real64 faceRotationMatrix[ 3 ][ 3 ];
  computationalGeometry::RotationMatrix_3D( faceNormal, faceRotationMatrix );
  //  - below we compute the diameter, the rotated vertices and the rotated center.
  real64 faceRotatedVertices[ MFN ][ 2 ];
  real64 faceDiameter = 0;

  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    // apply the transpose (that is the inverse) of the rotation matrix to face vertices.
    // NOTE:
    // the second and third rows of the transpose of the rotation matrix rotate on the 2D face.
    faceRotatedVertices[numVertex][0] =
      faceRotationMatrix[ 0 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 0 ) +
      faceRotationMatrix[ 1 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 1 ) +
      faceRotationMatrix[ 2 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 2 );
    faceRotatedVertices[numVertex][1] =
      faceRotationMatrix[ 0 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 0 ) +
      faceRotationMatrix[ 1 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 1 ) +
      faceRotationMatrix[ 2 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 2 );
  }

  faceDiameter = computationalGeometry::computeDiameter< 2 >( faceRotatedVertices,
                                                              numFaceVertices );
  real64 const invFaceDiameter = 1.0/faceDiameter;
  // - rotate the face centroid as done for the vertices.
  real64 faceRotatedCentroid[2];
  faceRotatedCentroid[0] =
    faceRotationMatrix[ 0 ][ 1 ]*faceCenter[0] +
    faceRotationMatrix[ 1 ][ 1 ]*faceCenter[1] +
    faceRotationMatrix[ 2 ][ 1 ]*faceCenter[2];
  faceRotatedCentroid[1] =
    faceRotationMatrix[ 0 ][ 2 ]*faceCenter[0] +
    faceRotationMatrix[ 1 ][ 2 ]*faceCenter[1] +
    faceRotationMatrix[ 2 ][ 2 ]*faceCenter[2];
  // - compute edges' lengths, outward pointing normals and local edge-to-nodes map.
  real64 edgeOutwardNormals[ MFN ][ 2 ];
  real64 edgeLengths[ MFN ];
  localIndex localEdgeToNodes[ MFN ][ 2 ];

  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    if( edgeToNodes( faceToEdges[numEdge], 0 ) == faceToNodes[ numEdge ] )
    {
      localEdgeToNodes[ numEdge ][ 0 ] = numEdge;
      localEdgeToNodes[ numEdge ][ 1 ] = (numEdge+1)%numFaceVertices;
    }
    else
    {
      localEdgeToNodes[ numEdge ][ 0 ] = (numEdge+1)%numFaceVertices;
      localEdgeToNodes[ numEdge ][ 1 ] = numEdge;
    }
    real64 edgeTangent[2];
    edgeTangent[0] = faceRotatedVertices[(numEdge+1)%numFaceVertices][0] -
                     faceRotatedVertices[numEdge][0];
    edgeTangent[1] = faceRotatedVertices[(numEdge+1)%numFaceVertices][1] -
                     faceRotatedVertices[numEdge][1];
    edgeOutwardNormals[numEdge][0] = edgeTangent[1];
    edgeOutwardNormals[numEdge][1] = -edgeTangent[0];
    real64 signTestVector[2];
    signTestVector[0] = faceRotatedVertices[numEdge][0] - faceRotatedCentroid[0];
    signTestVector[1] = faceRotatedVertices[numEdge][1] - faceRotatedCentroid[1];
    if( signTestVector[0]*edgeOutwardNormals[numEdge][0] +
        signTestVector[1]*edgeOutwardNormals[numEdge][1] < 0 )
    {
      edgeOutwardNormals[numEdge][0] = -edgeOutwardNormals[numEdge][0];
      edgeOutwardNormals[numEdge][1] = -edgeOutwardNormals[numEdge][1];
    }
    edgeLengths[numEdge] = LvArray::math::sqrt< real64 >( edgeTangent[0]*edgeTangent[0] +
                                                          edgeTangent[1]*edgeTangent[1] );
    edgeOutwardNormals[numEdge][0] /= edgeLengths[numEdge];
    edgeOutwardNormals[numEdge][1] /= edgeLengths[numEdge];
  }

  // Compute boundary quadrature weights (also equal to the integrals of basis functions on the
  // boundary).
  real64 boundaryQuadratureWeights[ MFN ];
  for( localIndex numWeight = 0; numWeight < numFaceVertices; ++numWeight )
    boundaryQuadratureWeights[numWeight] = 0.0;
  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    boundaryQuadratureWeights[ localEdgeToNodes[ numEdge ][ 0 ] ] += 0.5*edgeLengths[numEdge];
    boundaryQuadratureWeights[ localEdgeToNodes[ numEdge ][ 1 ] ] += 0.5*edgeLengths[numEdge];
  }

  // Compute scaled monomials' integrals on edges.
  real64 monomBoundaryIntegrals[3] = { 0.0 };
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    monomBoundaryIntegrals[0] += boundaryQuadratureWeights[ numVertex ];
    monomBoundaryIntegrals[1] += (faceRotatedVertices[ numVertex ][ 0 ] - faceRotatedCentroid[0]) *
                                 invFaceDiameter*boundaryQuadratureWeights[ numVertex ];
    monomBoundaryIntegrals[2] += (faceRotatedVertices[ numVertex ][ 1 ] - faceRotatedCentroid[1]) *
                                 invFaceDiameter*boundaryQuadratureWeights[ numVertex ];
  }

  // Compute non constant 2D and 3D scaled monomials' integrals on the face.
  real64 monomInternalIntegrals[2] = { 0.0 };
  for( localIndex numSubTriangle = 0; numSubTriangle < numFaceVertices; ++numSubTriangle )
  {
    localIndex const nextVertex = (numSubTriangle+1)%numFaceVertices;
    // - compute value of 2D monomials at the quadrature point on the sub-triangle (the
    //   barycenter).
    //   The result is ((v(0)+v(1)+faceCenter)/3 - faceCenter) / faceDiameter =
    //   = (v(0) + v(1) - 2*faceCenter)/(3*faceDiameter).
    real64 monomialValues[2];
    for( localIndex i = 0; i < 2; ++i )
    {
      monomialValues[i] = (faceRotatedVertices[numSubTriangle][i] +
                           faceRotatedVertices[nextVertex][i] -
                           2.0*faceRotatedCentroid[i]) / (3.0*faceDiameter);
    }
    // compute value of 3D monomials at the quadrature point on the sub-triangle (the
    // barycenter).  The result is
    // ((v(0) + v(1) + faceCenter)/3 - cellCenter)/cellDiameter.
    real64 threeDMonomialValues[3];
    for( localIndex i = 0; i < 3; ++i )
    {
      threeDMonomialValues[i] = ( (faceCenter[i] +
                                   nodesCoords[faceToNodes[ numSubTriangle ]][i] +
                                   nodesCoords[faceToNodes[ nextVertex ]][i]) / 3.0 -
                                  cellCenter[i] ) * invCellDiameter;
    }
    // compute quadrature weight associated to the quadrature point (the area of the
    // sub-triangle).
    real64 edgesTangents[2][2];               // used to compute the area of the sub-triangle
    for( localIndex i = 0; i < 2; ++i )
    {
      edgesTangents[0][i] = faceRotatedVertices[numSubTriangle][i] - faceRotatedCentroid[i];
    }
    for( localIndex i = 0; i < 2; ++i )
    {
      edgesTangents[1][i] = faceRotatedVertices[nextVertex][i] - faceRotatedCentroid[i];
    }
    real64 subTriangleArea = 0.5*LvArray::math::abs
                               ( edgesTangents[0][0]*edgesTangents[1][1] -
                               edgesTangents[0][1]*edgesTangents[1][0] );
    // compute the integrals on the sub-triangle and add it to the global integrals
    for( localIndex i = 0; i < 2; ++i )
    {
      monomInternalIntegrals[ i ] += monomialValues[ i ]*subTriangleArea;
    }
    for( localIndex i = 0; i < 3; ++i )
    {
      // threeDMonomialIntegrals is assumed to be initialized to 0 by the caller
      threeDMonomialIntegrals[ i ] += threeDMonomialValues[ i ]*subTriangleArea;
    }
  }

  // Compute integral of basis functions times normal derivative of monomials on the boundary.
  real64 basisTimesMonomNormalDerBoundaryInt[ MFN ][ 2 ];
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      basisTimesMonomNormalDerBoundaryInt[ numVertex ][ i ] = 0.0;
    }
  }
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      real64 thisEdgeIntTimesNormal_i = edgeOutwardNormals[numVertex][i]*edgeLengths[numVertex];
      basisTimesMonomNormalDerBoundaryInt[ localEdgeToNodes[ numVertex ][ 0 ] ][i] += thisEdgeIntTimesNormal_i;
      basisTimesMonomNormalDerBoundaryInt[ localEdgeToNodes[ numVertex ][ 1 ] ][i] += thisEdgeIntTimesNormal_i;
    }
  }
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      basisTimesMonomNormalDerBoundaryInt[ numVertex ][ i ] *= 0.5*invFaceDiameter;
    }
  }

  // Compute integral mean of basis functions on this face.
  real64 const invFaceArea = 1.0/faceArea;
  real64 const monomialDerivativeInverse = (faceDiameter*faceDiameter)*invFaceArea;
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    real64 piNablaDofs[ 3 ];
    piNablaDofs[ 1 ] = monomialDerivativeInverse *
                       basisTimesMonomNormalDerBoundaryInt[ numVertex ][ 0 ];
    piNablaDofs[ 2 ] = monomialDerivativeInverse *
                       basisTimesMonomNormalDerBoundaryInt[ numVertex ][ 1 ];
    piNablaDofs[ 0 ] = (boundaryQuadratureWeights[ numVertex ] -
                        piNablaDofs[ 1 ]*monomBoundaryIntegrals[ 1 ] -
                        piNablaDofs[ 2 ]*monomBoundaryIntegrals[ 2 ])/monomBoundaryIntegrals[ 0 ];
    basisIntegrals[ numVertex ] = piNablaDofs[ 0 ]*faceArea +
                                  (piNablaDofs[ 1 ]*monomInternalIntegrals[ 0 ] +
                                   piNablaDofs[ 2 ]*monomInternalIntegrals[ 1 ]);
  }
}

void SolidMechanicsLagrangeContact::computeFaceNodalArea( localIndex const kf0,
                                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                                                          ArrayOfArraysView< localIndex const > const & faceToEdgeMap,
                                                          arrayView2d< localIndex const > const & edgeToNodeMap,
                                                          arrayView2d< real64 const > const faceCenters,
                                                          arrayView2d< real64 const > const faceNormals,
                                                          arrayView1d< real64 const > const faceAreas,
                                                          stackArray1d< real64, FaceManager::maxFaceNodes() > & basisIntegrals ) const
{
  GEOS_MARK_FUNCTION;
  localIndex const TriangularPermutation[3] = { 0, 1, 2 };
  localIndex const QuadrilateralPermutation[4] = { 0, 1, 3, 2 };
  localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

  basisIntegrals.resize( numNodesPerFace );
  for( localIndex a = 0; a < numNodesPerFace; ++a )
  {
    basisIntegrals[a] = 0.0;
  }
  localIndex const * const permutation = ( numNodesPerFace == 3 ) ? TriangularPermutation : QuadrilateralPermutation;
  if( numNodesPerFace == 3 )
  {
    real64 xLocal[3][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[3];
    for( localIndex q=0; q<H1_TriangleFace_Lagrange1_Gauss1::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, xLocal );
      H1_TriangleFace_Lagrange1_Gauss1::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        basisIntegrals[a] += detJ * N[permutation[a]];
      }
    }
  }
  else if( numNodesPerFace == 4 )
  {
    real64 xLocal[4][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[4];
    for( localIndex q=0; q<H1_QuadrilateralFace_Lagrange1_GaussLegendre2::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight( q, xLocal );
      H1_QuadrilateralFace_Lagrange1_GaussLegendre2::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        basisIntegrals[a] += detJ * N[permutation[a]];
      }
    }
  }
  else if( numNodesPerFace > 4 && numNodesPerFace <= m_maxFaceNodes )
  {
    // we need to L2 projector based on VEM to approximate the quadrature weights
    // we need to use extra geometry information to computing L2 projector

    localIndex const MFN = m_maxFaceNodes; // Max number of face vertices.
    localIndex const faceIndex = kf0;
    localIndex const numFaceNodes = faceToNodeMap[ faceIndex ].size();

    // get the face center and normal.
    real64 const faceArea = faceAreas[ faceIndex ];
    localIndex faceToNodes[ MFN ];
    localIndex faceToEdges[ MFN ];
    for( localIndex i = 0; i < numFaceNodes; ++i )
    {
      faceToNodes[i] = faceToNodeMap[ faceIndex ][ i ];
      faceToEdges[i] = faceToEdgeMap[ faceIndex ][ i ];
    }
    // - get outward face normal and center
    real64 faceNormal[3] = { faceNormals[faceIndex][0],
                             faceNormals[faceIndex][1],
                             faceNormals[faceIndex][2] };
    real64 const faceCenter[3] { faceCenters[faceIndex][0],
                                 faceCenters[faceIndex][1],
                                 faceCenters[faceIndex][2] };
    // - compute integrals calling auxiliary method
    real64 threeDMonomialIntegrals[3] = { 0.0 };
    real64 const invCellDiameter = 0.0;
    real64 const cellCenter[3] { 0.0, 0.0, 0.0 };
    computeFaceIntegrals( nodePosition,
                          faceToNodes,
                          faceToEdges,
                          numFaceNodes,
                          faceArea,
                          faceCenter,
                          faceNormal,
                          edgeToNodeMap,
                          invCellDiameter,
                          cellCenter,
                          basisIntegrals,
                          threeDMonomialIntegrals );
  }
  else
  {
    GEOS_ERROR( "SolidMechanicsLagrangeContact " << getDataContext() << ": face with " << numNodesPerFace <<
                " nodes. Only triangles and quadrilaterals and PEBI prisms up to 11 sides are supported." );
  }
}

void SolidMechanicsLagrangeContact::
  assembleForceResidualDerivativeWrtTraction( MeshLevel const & mesh,
                                              arrayView1d< string const > const & regionNames,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > faceNormals = faceManager.faceNormal();
  arrayView1d< real64 const > faceAreas = faceManager.faceArea();

  string const & tracDofKey = dofManager.getKey( contact::traction::key() );
  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    arrayView1d< globalIndex const > const & tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
    arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( contact::traction::key() );
    arrayView3d< real64 const > const & rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
    ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
    {
      if( elemsToFaces.sizeOfArray( kfe ) != 2 )
      {
        return;
      }
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );

      globalIndex rowDOF[3 * m_maxFaceNodes]; // this needs to be changed when dealing with arbitrary element types
      real64 nodeRHS[3 * m_maxFaceNodes];
      stackArray2d< real64, 3 * m_maxFaceNodes * 3 > dRdT( 3 * m_maxFaceNodes, 3 );
      globalIndex colDOF[3];
      for( localIndex i = 0; i < 3; ++i )
      {
        colDOF[i] = tracDofNumber[kfe] + i;
      }

      for( localIndex kf = 0; kf < 2; ++kf )
      {
        constexpr int normalSign[2] = { 1, -1 };
        // Testing the face integral on polygonal faces
        stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
        localIndex const faceIndex = elemsToFaces[kfe][kf];
        computeFaceNodalArea( faceIndex,
                              nodePosition,
                              faceToNodeMap,
                              faceToEdgeMap,
                              edgeToNodeMap,
                              faceCenters,
                              faceNormals,
                              faceAreas,
                              nodalArea );

        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          real64 const localNodalForce[ 3 ] = { traction( kfe, 0 ) * nodalArea[a],
                                                traction( kfe, 1 ) * nodalArea[a],
                                                traction( kfe, 2 ) * nodalArea[a] };
          real64 globalNodalForce[ 3 ];
          LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( globalNodalForce, rotationMatrix[ kfe ], localNodalForce );

          for( localIndex i = 0; i < 3; ++i )
          {
            rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
            // Opposite sign w.r.t. to formulation presented in
            // Algebraically Stabilized Lagrange Multiplier Method for Frictional Contact Mechanics with
            // Hydraulically Active Fractures
            // Franceschini, A., Castelletto, N., White, J. A., Tchelepi, H. A.
            // Computer Methods in Applied Mechanics and Engineering (2020) 368, 113161
            // doi: 10.1016/j.cma.2020.113161
            nodeRHS[3*a+i] = +globalNodalForce[i] * normalSign[ kf ];

            // Opposite sign w.r.t. to the same formulation as above
            dRdT( 3*a+i, 0 ) = rotationMatrix( kfe, i, 0 ) * normalSign[ kf ] * nodalArea[a];
            dRdT( 3*a+i, 1 ) = rotationMatrix( kfe, i, 1 ) * normalSign[ kf ] * nodalArea[a];
            dRdT( 3*a+i, 2 ) = rotationMatrix( kfe, i, 2 ) * normalSign[ kf ] * nodalArea[a];
          }
        }

        for( localIndex idof = 0; idof < numNodesPerFace * 3; ++idof )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                        colDOF,
                                                        dRdT[idof].dataIfContiguous(),
                                                        3 );
            RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], nodeRHS[idof] );
          }
        }
      }
    } );
  } );
}


void SolidMechanicsLagrangeContact::
  assembleTractionResidualDerivativeWrtDisplacementAndTraction( MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > faceNormals = faceManager.faceNormal();
  arrayView1d< real64 const > faceAreas = faceManager.faceArea();

  string const & tracDofKey = dofManager.getKey( contact::traction::key() );
  string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 FaceElementSubRegion const & subRegion )
  {
    string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
    FrictionBase const & frictionLaw = getConstitutiveModel< FrictionBase >( subRegion, frictionLawName );

    arrayView1d< globalIndex const > const & tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & area = subRegion.getElementArea();
    arrayView3d< real64 const > const &
    rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
    ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();
    arrayView2d< real64 const > const & traction = subRegion.getField< contact::traction >();
    arrayView1d< integer const > const & fractureState = subRegion.getField< contact::fractureState >();
    arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
    arrayView2d< real64 const > const & previousDispJump = subRegion.getField< contact::oldDispJump >();
    arrayView1d< real64 const > const & slidingTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

    constitutiveUpdatePassThru( frictionLaw, [&] ( auto & castedFrictionLaw )
    {
      using FrictionType = TYPEOFREF( castedFrictionLaw );
      typename FrictionType::KernelWrapper frictionWrapper = castedFrictionLaw.createKernelUpdates();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        if( elemsToFaces.sizeOfArray( kfe ) != 2 )
        {
          return;
        }

        if( ghostRank[kfe] < 0 )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          globalIndex nodeDOF[2 * 3 * m_maxFaceNodes];
          globalIndex elemDOF[3];
          for( localIndex i = 0; i < 3; ++i )
          {
            elemDOF[i] = tracDofNumber[kfe] + i;
          }

          real64 elemRHS[3] = {0.0, 0.0, 0.0};
          real64 const Ja = area[kfe];

          stackArray2d< real64, 2 * 3 * m_maxFaceNodes * 3 > dRdU( 3, 2 * 3 * m_maxFaceNodes );
          stackArray2d< real64, 3 * 3 > dRdT( 3, 3 );

          switch( fractureState[kfe] )
          {
            case contact::FractureState::Stick:
              {
                for( localIndex i = 0; i < 3; ++i )
                {
                  if( i == 0 )
                  {
                    elemRHS[i] = Ja * dispJump[kfe][i];
                  }
                  else
                  {
                    elemRHS[i] = Ja * ( dispJump[kfe][i] - previousDispJump[kfe][i] );
                  }
                }

                for( localIndex kf = 0; kf < 2; ++kf )
                {
                  // Compute local area contribution for each node
                  stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
                  computeFaceNodalArea( elemsToFaces[kfe][kf],
                                        nodePosition,
                                        faceToNodeMap,
                                        faceToEdgeMap,
                                        edgeToNodeMap,
                                        faceCenters,
                                        faceNormals,
                                        faceAreas,
                                        nodalArea );

                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + i;

                      for( localIndex j = 0; j < 3; ++j )
                      {
                        dRdU( j, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, j ) * pow( -1, kf );
                      }
                    }
                  }
                }
                break;
              }
            case contact::FractureState::Slip:
            case contact::FractureState::NewSlip:
              {
                elemRHS[0] = Ja * dispJump[kfe][0];

                for( localIndex kf = 0; kf < 2; ++kf )
                {
                  // Compute local area contribution for each node
                  stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
                  computeFaceNodalArea( elemsToFaces[kfe][kf],
                                        nodePosition,
                                        faceToNodeMap,
                                        faceToEdgeMap,
                                        edgeToNodeMap,
                                        faceCenters,
                                        faceNormals,
                                        faceAreas,
                                        nodalArea );

                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] +
                                                                      LvArray::integerConversion< globalIndex >( i );
                      dRdU( 0, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, 0 ) * pow( -1, kf );
                    }
                  }
                }

                real64 dLimitTau_dNormalTraction = 0;
                real64 const limitTau = frictionWrapper.computeLimitTangentialTractionNorm( traction[kfe][0],
                                                                                            dLimitTau_dNormalTraction );

                real64 sliding[ 2 ] = { dispJump[kfe][1] - previousDispJump[kfe][1], dispJump[kfe][2] - previousDispJump[kfe][2] };
                real64 slidingNorm = sqrt( sliding[ 0 ]*sliding[ 0 ] + sliding[ 1 ]*sliding[ 1 ] );

                if( !( ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 ) && ( fractureState[kfe] == contact::FractureState::NewSlip ) )
                    && slidingNorm > slidingTolerance[kfe] )
                {
                  for( localIndex i = 1; i < 3; ++i )
                  {
                    elemRHS[i] = Ja * ( traction[kfe][i] - limitTau * sliding[ i-1 ] / slidingNorm );

                    dRdT( i, 0 ) = Ja * dLimitTau_dNormalTraction * sliding[ i-1 ] / slidingNorm;
                    dRdT( i, i ) = Ja;
                  }

                  // A symmetric 2x2 matrix.
                  real64 dUdgT[ 3 ];
                  dUdgT[ 0 ] = (slidingNorm * slidingNorm - sliding[ 0 ] * sliding[ 0 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 1 ] = (slidingNorm * slidingNorm - sliding[ 1 ] * sliding[ 1 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 2 ] = -sliding[ 0 ] * sliding[ 1 ] * limitTau / std::pow( slidingNorm, 3 );

                  for( localIndex kf = 0; kf < 2; ++kf )
                  {
                    // Compute local area contribution for each node
                    stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea;
                    computeFaceNodalArea( elemsToFaces[kfe][kf],
                                          nodePosition,
                                          faceToNodeMap,
                                          faceToEdgeMap,
                                          edgeToNodeMap,
                                          faceCenters,
                                          faceNormals,
                                          faceAreas,
                                          nodalArea );

                    for( localIndex a = 0; a < numNodesPerFace; ++a )
                    {
                      for( localIndex i = 0; i < 3; ++i )
                      {
                        real64 const localRowB[ 2 ] = { rotationMatrix( kfe, i, 1 ), rotationMatrix( kfe, i, 2 ) };
                        real64 localRowE[ 2 ];
                        LvArray::tensorOps::Ri_eq_symAijBj< 2 >( localRowE, dUdgT, localRowB );

                        dRdU( 1, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea[a] * localRowE[ 0 ] * pow( -1, kf );
                        dRdU( 2, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea[a] * localRowE[ 1 ] * pow( -1, kf );
                      }
                    }
                  }
                }
                else
                {
                  real64 vaux[ 2 ] = { traction[kfe][1], traction[kfe][2] };
                  real64 vauxNorm = sqrt( vaux[ 0 ]*vaux[ 0 ] + vaux[ 1 ]*vaux[ 1 ] );
                  if( vauxNorm > 0.0 )
                  {
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      elemRHS[i] = Ja * traction[kfe][i] * ( 1.0 - limitTau / vauxNorm );

                      dRdT( i, 0 ) = Ja * traction[kfe][i] * dLimitTau_dNormalTraction / vauxNorm;
                      dRdT( i, i ) = Ja;
                    }
                  }
                  else
                  {
                    for( int i = 1; i < 3; ++i )
                    {
                      elemRHS[i] = 0.0;

                      dRdT( i, i ) = Ja;
                    }
                  }
                }
                break;
              }
            case contact::FractureState::Open:
              {
                for( int i = 0; i < 3; ++i )
                {
                  elemRHS[i] = Ja * traction[kfe][i];

                  dRdT( i, i ) = Ja;
                }
                break;
              }
          }

          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

          for( localIndex idof = 0; idof < 3; ++idof )
          {
            localRhs[localRow + idof] += elemRHS[idof];

            if( fractureState[kfe] != contact::FractureState::Open )
            {

              localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow + idof,
                                                                        nodeDOF,
                                                                        dRdU[idof].dataIfContiguous(),
                                                                        2 * 3 * numNodesPerFace );
            }

            if( fractureState[kfe] != contact::FractureState::Stick )
            {
              localMatrix.addToRow< serialAtomic >( localRow + idof,
                                                    elemDOF,
                                                    dRdT[idof].dataIfContiguous(),
                                                    3 );
            }
          }
        }
      } );
    } );
  } );
}

void SolidMechanicsLagrangeContact::assembleStabilization( MeshLevel const & mesh,
                                                           NumericalMethodsManager const & numericalMethodManager,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & tracDofKey = dofManager.getKey( contact::traction::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the finite volume method used to compute the stabilization
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemRegion = faceToElem.m_toElementRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemSubRegion = faceToElem.m_toElementSubRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex.toViewConst();

  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
  FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

  GEOS_ERROR_IF( !fractureSubRegion.hasField< contact::traction >(), "The fracture subregion must contain traction field." );
  ArrayOfArraysView< localIndex const > const elem2dToFaces = fractureSubRegion.faceList().toViewConst();

  // Get the state of fracture elements
  arrayView1d< integer const > const & fractureState =
    fractureSubRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );

  // Get the tractions and stabilization contribution to the local jump
  arrayView2d< real64 const > const & traction = fractureSubRegion.getField< contact::traction >();
  arrayView2d< real64 const > const & deltaTraction = fractureSubRegion.getField< contact::deltaTraction >();

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager.constructViewAccessor< real64_array, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // Get area and rotation matrix for all faces
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList().toViewConst();
  arrayView2d< real64 const > faceCenters = faceManager.faceCenter();
  arrayView2d< real64 const > faceNormals = faceManager.faceNormal();

  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView3d< real64 const > const &
  faceRotationMatrix = fractureSubRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );

  // Bulk modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::bulkModulusString() );
  // Shear modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::shearModulusString() );

  using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
  ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
    elemManager.constructViewAccessor< CellElementSubRegion::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );
  ElementRegionManager::ElementViewConst< NodeMapViewType > const elemToNodeView = elemToNode.toNestedViewConst();

  arrayView1d< globalIndex const > const & tracDofNumber = fractureSubRegion.getReference< globalIndex_array >( tracDofKey );

  stabilizationMethod.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

    forAll< serialPolicy >( stencil.size(), [=] ( localIndex const iconn )
    {
      localIndex const numFluxElems = sei.sizeOfArray( iconn );

      // A fracture connector has to be an edge shared by two faces
      if( numFluxElems == 2 )
      {
        // Find shared edge (pair of nodes)
        real64 Nbar0[3];
        real64 Nbar1[3];
        Nbar0[ 0 ] = faceRotationMatrix[ sei[iconn][0] ][0][0];
        Nbar0[ 1 ] = faceRotationMatrix[ sei[iconn][0] ][1][0];
        Nbar0[ 2 ] = faceRotationMatrix[ sei[iconn][0] ][2][0];
        Nbar1[ 0 ] = faceRotationMatrix[ sei[iconn][1] ][0][0];
        Nbar1[ 1 ] = faceRotationMatrix[ sei[iconn][1] ][1][0];
        Nbar1[ 2 ] = faceRotationMatrix[ sei[iconn][1] ][2][0];

        real64 normalProduct = LvArray::tensorOps::AiBi< 3 >( Nbar0, Nbar1 );

        localIndex const id1 = ( normalProduct > 0.0 ) ? 0 : 1;

        localIndex const numNodesPerFace0 = faceToNodeMap.sizeOfArray( elem2dToFaces[sei[iconn][0]][0] );
        array1d< localIndex > nodes0( numNodesPerFace0 );
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          nodes0[i] = faceToNodeMap( elem2dToFaces[sei[iconn][0]][0], i );
        }
        localIndex const numNodesPerFace1 = faceToNodeMap.sizeOfArray( elem2dToFaces[sei[iconn][1]][0] );
        array1d< localIndex > nodes1( numNodesPerFace1 );
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          nodes1[i] = faceToNodeMap( elem2dToFaces[sei[iconn][1]][id1], i );
        }
        std::sort( nodes0.begin(), nodes0.end() );
        std::sort( nodes1.begin(), nodes1.end() );
        array1d< localIndex > edge( std::max( numNodesPerFace0, numNodesPerFace1 ) );
        edge.setValues< serialPolicy >( -1 );
        std::set_intersection( nodes0.begin(), nodes0.end(), nodes1.begin(), nodes1.end(), edge.begin() );
        localIndex realNodes = 0;
        for( localIndex i = 0; i < edge.size(); ++i )
        {
          if( edge[i] > -1 )
          {
            realNodes++;
          }
        }
        GEOS_ERROR_IF( realNodes != 2,
                       getDataContext() << ": An edge shared by two fracture elements must have 2 nodes." );
        edge.resize( realNodes );

        // Compute nodal area factor
        localIndex node0index0 = -1;
        localIndex node1index0 = -1;
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          if( edge[0] == faceToNodeMap( elem2dToFaces[sei[iconn][0]][0], i ) )
          {
            node0index0 = i;
          }
          if( edge[1] == faceToNodeMap( elem2dToFaces[sei[iconn][0]][0], i ) )
          {
            node1index0 = i;
          }
        }
        localIndex node0index1 = -1;
        localIndex node1index1 = -1;
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          if( edge[0] == faceToNodeMap( elem2dToFaces[sei[iconn][1]][id1], i ) )
          {
            node0index1 = i;
          }
          if( edge[1] == faceToNodeMap( elem2dToFaces[sei[iconn][1]][id1], i ) )
          {
            node1index1 = i;
          }
        }
        stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea0;
        stackArray1d< real64, FaceManager::maxFaceNodes() > nodalArea1;
        localIndex const faceIndex0 = elem2dToFaces[sei[iconn][0]][0];
        localIndex const faceIndex1 = elem2dToFaces[sei[iconn][1]][id1];

        computeFaceNodalArea( faceIndex0,
                              nodePosition,
                              faceToNodeMap,
                              faceToEdgeMap,
                              edgeToNodeMap,
                              faceCenters,
                              faceNormals,
                              faceArea,
                              nodalArea0 );

        computeFaceNodalArea( faceIndex1,
                              nodePosition,
                              faceToNodeMap,
                              faceToEdgeMap,
                              edgeToNodeMap,
                              faceCenters,
                              faceNormals,
                              faceArea,
                              nodalArea1 );

        real64 const areafac = nodalArea0[node0index0] * nodalArea1[node0index1] + nodalArea0[node1index0] * nodalArea1[node1index1];

        // first index: face, second index: element (T/B), third index: dof (x, y, z)
        real64 stiffDiagApprox[ 2 ][ 2 ][ 3 ];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Get fracture, face and region/subregion/element indices (for elements on both sides)
          localIndex const fractureIndex = sei[iconn][kf];

          for( localIndex i = 0; i < 2; ++i )
          {
            localIndex const faceIndex = ( kf == 0 || id1 == 0 ) ? elem2dToFaces[fractureIndex][i] : elem2dToFaces[fractureIndex][1 - i];
            localIndex const ke = faceToElemIndex[faceIndex][0] >= 0 ? 0 : 1;

            localIndex const er  = faceToElemRegion[faceIndex][ke];
            localIndex const esr = faceToElemSubRegion[faceIndex][ke];
            localIndex const ei  = faceToElemIndex[faceIndex][ke];

            real64 const volume = elemVolume[er][esr][ei];

            // Get the "element to node" map for the specific region/subregion
            NodeMapViewType const & cellElemsToNodes = elemToNodeView[er][esr];
            localIndex const numNodesPerElem = cellElemsToNodes.size( 1 );

            // Compute the box size
            real64 maxSize[3];
            real64 minSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              maxSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
              minSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
            }
            for( localIndex a = 1; a < numNodesPerElem; ++a )
            {
              for( localIndex j = 0; j < 3; ++j )
              {
                maxSize[j] = fmax( maxSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
                minSize[j] = fmin( minSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
              }
            }
            real64 boxSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              boxSize[j] = maxSize[j] - minSize[j];
            }

            // Get linear elastic isotropic constitutive parameters for the element
            real64 const K = bulkModulus[er][esr][ei];
            real64 const G = shearModulus[er][esr][ei];
            real64 const E = 9.0 * K * G / ( 3.0 * K + G );
            real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );

            // Combine E and nu to obtain a stiffness approximation (like it was an hexahedron)
            for( localIndex j = 0; j < 3; ++j )
            {
              stiffDiagApprox[ kf ][ i ][ j ] = m_stabilitzationScalingCoefficient * E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 2.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] );
            }
          }
        }
        real64 invTotStiffApprox[ 3 ][ 3 ] = { { 0 } };
        for( localIndex i = 0; i < 3; ++i )
        {
          // K(i,i)^-1 = Ka(i,i)^-1 + Kb(i,i)^-1
          // T -> top (index 0), B -> bottom (index 1)
          // Ka(i,i) = KT(i,i) + KB(i,i)
          // Kb(i,i) = KT(i,i) + KB(i,i)
          invTotStiffApprox[ i ][ i ] = 1.0 / ( stiffDiagApprox[ 0 ][ 0 ][ i ] + stiffDiagApprox[ 1 ][ 0 ][ i ] )
                                        + 1.0 / ( stiffDiagApprox[ 0 ][ 1 ][ i ] + stiffDiagApprox[ 1 ][ 1 ][ i ] );
        }

        array2d< real64 > avgRotationMatrix( 3, 3 );

        // To be able to compute an average rotation matrix, normal has to point in the same direction.
        if( normalProduct < 0.0 )
        {
          LvArray::tensorOps::scale< 3 >( Nbar1, -1.0 );
          normalProduct *= -1.0;
        }
        // If the surfaces are co-planar, then use the first rotation matrix
        if( std::abs( normalProduct - 1.0 ) < 1.e+2*machinePrecision )
        {
          LvArray::tensorOps::copy< 3, 3 >( avgRotationMatrix, faceRotationMatrix[ sei[iconn][0] ] );
        }
        // otherwise, compute the average rotation matrix
        else
        {
          real64 avgNbar[3];
          avgNbar[ 0 ] = faceArea[elem2dToFaces[ sei[iconn][0] ][0]] * Nbar0[0] + faceArea[elem2dToFaces[ sei[iconn][1] ][0]] * Nbar1[0];
          avgNbar[ 1 ] = faceArea[elem2dToFaces[ sei[iconn][0] ][0]] * Nbar0[1] + faceArea[elem2dToFaces[ sei[iconn][1] ][0]] * Nbar1[1];
          avgNbar[ 2 ] = faceArea[elem2dToFaces[ sei[iconn][0] ][0]] * Nbar0[2] + faceArea[elem2dToFaces[ sei[iconn][1] ][0]] * Nbar1[2];
          LvArray::tensorOps::normalize< 3 >( avgNbar );

          computationalGeometry::RotationMatrix_3D( avgNbar, avgRotationMatrix );
        }

        // Compute R^T * (invK) * R
        real64 temp[ 3 ][ 3 ];
        real64 rotatedInvStiffApprox[ 3 ][ 3 ];
        LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( temp, avgRotationMatrix, invTotStiffApprox );
        LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, avgRotationMatrix );

        // Add nodal area contribution
        stackArray2d< real64, 3*3 > totalInvStiffApprox( 3, 3 );
        for( localIndex i = 0; i < 3; ++i )
        {
          for( localIndex j = 0; j < 3; ++j )
          {
            totalInvStiffApprox( i, j ) = -rotatedInvStiffApprox[ i ][ j ] * areafac;
          }
        }

        // Get DOF numbering
        localIndex fractureIndex[2];
        localIndex nDof[2];
        globalIndex elemDOF[2][3];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          fractureIndex[kf] = sei[iconn][kf];
          for( localIndex i = 0; i < 3; ++i )
          {
            elemDOF[kf][i] = tracDofNumber[fractureIndex[kf]] + i;
          }
          nDof[kf] = 0;
          switch( fractureState[fractureIndex[kf]] )
          {
            case ( contact::FractureState::Stick ):
              {
                nDof[kf] = 3;
                break;
              }
            case ( contact::FractureState::NewSlip ):
            case ( contact::FractureState::Slip ):
              {
                nDof[kf] = 1;
                break;
              }
            case ( contact::FractureState::Open ):
              {
                nDof[kf] = 0;
                break;
              }
          }
        }

        // Define local "transmissibility" matrices
        stackArray2d< real64, 3*3 > totalInvStiffApprox00( nDof[0], nDof[0] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox01( nDof[0], nDof[1] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox10( nDof[1], nDof[0] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox11( nDof[1], nDof[1] );
        for( localIndex i = 0; i < nDof[0]; ++i )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            totalInvStiffApprox00( i, j ) = totalInvStiffApprox( i, j );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            totalInvStiffApprox01( i, j ) = -totalInvStiffApprox( i, j );
          }
        }

        for( localIndex i = 0; i < nDof[1]; ++i )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            totalInvStiffApprox10( i, j ) = -totalInvStiffApprox( i, j );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            totalInvStiffApprox11( i, j ) = totalInvStiffApprox( i, j );
          }
        }

        // Compute rhs
        stackArray1d< real64, 3 > rhs0( 3 );
        if( nDof[0] > 0 )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            rhs0( 0 ) += totalInvStiffApprox00( 0, j ) * ( traction[fractureIndex[0]][j] );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            rhs0( 0 ) += totalInvStiffApprox01( 0, j ) * ( traction[fractureIndex[1]][j] );
          }
          for( localIndex i = 1; i < nDof[0]; ++i )
          {
            for( localIndex j = 0; j < nDof[0]; ++j )
            {
              rhs0( i ) += totalInvStiffApprox00( i, j ) * ( deltaTraction[fractureIndex[0]][j] );
            }
            for( localIndex j = 0; j < nDof[1]; ++j )
            {
              rhs0( i ) += totalInvStiffApprox01( i, j ) * ( deltaTraction[fractureIndex[1]][j] );
            }
          }
        }

        stackArray1d< real64, 3 > rhs1( 3 );
        if( nDof[1] > 0 )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            rhs1( 0 ) += totalInvStiffApprox10( 0, j ) * ( traction[fractureIndex[0]][j] );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            rhs1( 0 ) += totalInvStiffApprox11( 0, j ) * ( traction[fractureIndex[1]][j] );
          }
          for( localIndex i = 1; i < nDof[1]; ++i )
          {
            for( localIndex j = 0; j < nDof[0]; ++j )
            {
              rhs1( i ) += totalInvStiffApprox10( i, j ) * ( deltaTraction[fractureIndex[0]][j] );
            }
            for( localIndex j = 0; j < nDof[1]; ++j )
            {
              rhs1( i ) += totalInvStiffApprox11( i, j ) * ( deltaTraction[fractureIndex[1]][j] );
            }
          }
        }

        // Global matrix and rhs assembly
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[kf][0] - rankOffset );

          stackArray2d< real64, 3*3 > const & totalInvStiffApproxDiag = ( kf == 0 ) ? totalInvStiffApprox00 : totalInvStiffApprox11;
          stackArray2d< real64, 3*3 > const & totalInvStiffApproxOffDiag = ( kf == 0 ) ? totalInvStiffApprox01 : totalInvStiffApprox10;
          stackArray1d< real64, 3 > const & rhs = ( kf == 0 ) ? rhs0 : rhs1;

          // Only assemble contribution if "row" fracture element is local
          // TODO: use parallel atomics
          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            for( localIndex idof = 0; idof < nDof[kf]; ++idof )
            {
              // (i,i)-block
              localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow + idof,
                                                                              elemDOF[kf],
                                                                              totalInvStiffApproxDiag[idof].dataIfContiguous(),
                                                                              nDof[kf] );
              // (i,j)-block
              if( nDof[1-kf] > 0 )
              {
                localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow + idof,
                                                                                elemDOF[1 - kf],
                                                                                totalInvStiffApproxOffDiag[idof].dataIfContiguous(),
                                                                                nDof[1 - kf] );
              }

              // residual
              RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow + idof], rhs[idof] );
            }
          }
        }
      }
    } );
  } );
}

void SolidMechanicsLagrangeContact::applySystemSolution( DofManager const & dofManager,
                                                         arrayView1d< real64 const > const & localSolution,
                                                         real64 const scalingFactor,
                                                         real64 const dt,
                                                         DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::deltaTraction::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::traction::key(),
                               scalingFactor );

  // fractureStateString is synchronized in UpdateFractureState
  // oldFractureStateString and oldDispJumpString used locally only

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::traction::key(),
                                       contact::deltaTraction::key(),
                                       contact::dispJump::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void SolidMechanicsLagrangeContact::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  computeFaceDisplacementJump( domain );
}

bool SolidMechanicsLagrangeContact::resetConfigurationToDefault( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  using namespace fields::contact;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();
        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( fractureState[kfe] != FractureState::Open )
          {
            fractureState[kfe] = FractureState::Stick;
          }
        } );
      }
    } );
  } );
  return false;
}

bool SolidMechanicsLagrangeContact::updateConfiguration( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  using namespace fields::contact;

  real64 changedArea = 0;
  real64 totalArea = 0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {
      string const & fricitonLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
      FrictionBase const & frictionLaw = getConstitutiveModel< FrictionBase >( subRegion, fricitonLawName );

      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView1d< integer > const & fractureState = subRegion.getField< contact::fractureState >();
      arrayView1d< real64 const > const & faceArea = subRegion.getElementArea().toViewConst();

      arrayView1d< real64 const > const & normalTractionTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
      arrayView1d< real64 const > const & normalDisplacementTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );

      RAJA::ReduceSum< parallelHostReduce, real64 > changed( 0 );
      RAJA::ReduceSum< parallelHostReduce, real64 > total( 0 );

      constitutiveUpdatePassThru( frictionLaw, [&] ( auto & castedFrictionLaw )
      {
        using FrictionType = TYPEOFREF( castedFrictionLaw );
        typename FrictionType::KernelWrapper frictionWrapper = castedFrictionLaw.createKernelUpdates();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            integer const originalFractureState = fractureState[kfe];
            if( originalFractureState == FractureState::Open )
            {
              if( dispJump[kfe][0] <= -normalDisplacementTolerance[kfe] )
              {
                fractureState[kfe] = FractureState::Stick;
                if( getLogLevel() >= 10 )
                  GEOS_LOG( GEOS_FMT( "{}: {} -> {}: dispJump = {}, normalDisplacementTolerance = {}",
                                      kfe, originalFractureState, fractureState[kfe],
                                      dispJump[kfe][0], normalDisplacementTolerance[kfe] ) );
              }
            }
            else if( traction[kfe][0] > normalTractionTolerance[kfe] )
            {
              fractureState[kfe] = FractureState::Open;
              if( getLogLevel() >= 10 )
                GEOS_LOG( GEOS_FMT( "{}: {} -> {}: traction = {}, normalTractionTolerance = {}",
                                    kfe, originalFractureState, fractureState[kfe],
                                    traction[kfe][0], normalTractionTolerance[kfe] ) );
            }
            else
            {
              real64 currentTau = sqrt( traction[kfe][1]*traction[kfe][1] + traction[kfe][2]*traction[kfe][2] );

              real64 dLimitTangentialTractionNorm_dTraction = 0.0;
              real64 const limitTau =
                frictionWrapper.computeLimitTangentialTractionNorm( traction[kfe][0],
                                                                    dLimitTangentialTractionNorm_dTraction );

              if( originalFractureState == FractureState::Stick && currentTau >= limitTau )
              {
                currentTau *= (1.0 - m_slidingCheckTolerance);
              }
              else if( originalFractureState != FractureState::Stick && currentTau <= limitTau )
              {
                currentTau *= (1.0 + m_slidingCheckTolerance);
              }
              if( currentTau > limitTau )
              {
                if( originalFractureState == FractureState::Stick )
                {
                  fractureState[kfe] = FractureState::NewSlip;
                }
                else
                {
                  fractureState[kfe] = FractureState::Slip;
                }
              }
              else
              {
                fractureState[kfe] = FractureState::Stick;
              }
              if( getLogLevel() >= 10 && fractureState[kfe] != originalFractureState )
                GEOS_LOG( GEOS_FMT( "{}: {} -> {}: currentTau = {}, limitTau = {}",
                                    kfe, originalFractureState, fractureState[kfe],
                                    currentTau, limitTau ) );
            }

            changed += faceArea[kfe] * !compareFractureStates( originalFractureState, fractureState[kfe] );
            total += faceArea[kfe];
          }
        } );
      } );

      changedArea += changed.get();
      totalArea += total.get();
    } );
  } );

  // Need to synchronize the fracture state due to the use will be made of in AssemblyStabilization
  synchronizeFractureState( domain );

  // Compute global area of changed elements
  changedArea = MpiWrapper::sum( changedArea );
  // and total area of fracture elements
  totalArea = MpiWrapper::sum( totalArea );

  GEOS_LOG_LEVEL_INFO_RANK_0( logInfo::Configuration, GEOS_FMT( "  {}: changed area {} out of {}", getName(), changedArea, totalArea ) );

  // Assume converged if changed area is below certain fraction of total area
  return changedArea <= m_nonlinearSolverParameters.m_configurationTolerance * totalArea;
}

bool SolidMechanicsLagrangeContact::isFractureAllInStickCondition( DomainPartition const & domain ) const
{
  globalIndex numStick, numNewSlip, numSlip, numOpen;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    computeFractureStateStatistics( mesh, numStick, numNewSlip, numSlip, numOpen );
  } );

  return ( ( numNewSlip + numSlip + numOpen ) == 0 );
}

real64 SolidMechanicsLagrangeContact::setNextDt( real64 const & currentDt,
                                                 DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
  return currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangeContact, string const &, Group * const )

} /* namespace geos */
