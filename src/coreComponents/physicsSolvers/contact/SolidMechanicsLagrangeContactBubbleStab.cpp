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
 * @file SolidMechanicsLagrangeContactBubbleStab.cpp
 *
 */
#include "SolidMechanicsLagrangeContactBubbleStab.hpp"

#include "physicsSolvers/contact/kernels/RotationMatrixKernel.hpp"
#include "physicsSolvers/contact/kernels/SolidMechanicsLagrangeContactKernels.hpp"
#include "physicsSolvers/contact/kernels/SolidMechanicsALMSimultaneousKernels.hpp"
#include "physicsSolvers/contact/kernels/SolidMechanicsALMJumpUpdateKernels.hpp"
#include "physicsSolvers/contact/kernels/SolidMechanicsALMBubbleKernels.hpp"
#include "physicsSolvers/contact/LogLevelsInfo.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/FrictionSelector.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsLagrangeContactBubbleStab::SolidMechanicsLagrangeContactBubbleStab( const string & name,
                                                                                  Group * const parent ):
  ContactSolverBase( name, parent )
{
  m_faceTypeToFiniteElements["Quadrilateral"] =  std::make_unique< finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre2 >();
  m_faceTypeToFiniteElements["Triangle"] =  std::make_unique< finiteElement::H1_TriangleFace_Lagrange1_Gauss1 >();

}

SolidMechanicsLagrangeContactBubbleStab::~SolidMechanicsLagrangeContactBubbleStab()
{
  // TODO Auto-generated destructor stub
}

void SolidMechanicsLagrangeContactBubbleStab::registerDataOnMesh( Group & meshBodies )
{
  ContactSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & )
  {
    FaceManager & faceManager = meshLevel.getFaceManager();

    // Register the total bubble displacement
    faceManager.registerField< solidMechanics::totalBubbleDisplacement >( this->getName() ).
      reference().resizeDimension< 1 >( 3 );

    // Register the incremental bubble displacement
    faceManager.registerField< solidMechanics::incrementalBubbleDisplacement >( this->getName() ).
      reference().resizeDimension< 1 >( 3 );
  } );

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      // Register the rotation matrix
      subRegion.registerField< contact::rotationMatrix >( this->getName() ).
        reference().resizeDimension< 1, 2 >( 3, 3 );
    } );
  } );
}

void SolidMechanicsLagrangeContactBubbleStab::setupDofs( DomainPartition const & domain,
                                                         DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );

  map< std::pair< string, string >, array1d< string > > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                arrayView1d< string const > const & )
  {
    array1d< string > regions;
    regions.emplace_back( getUniqueFractureRegionName() );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addField( solidMechanics::totalBubbleDisplacement::key(),
                       FieldLocation::Face,
                       3,
                       meshTargets );

  // Add coupling between bubble
  // Useful to create connection between bubble dofs for Augmented Lagrangian formulation
  dofManager.addCoupling( solidMechanics::totalBubbleDisplacement::key(),
                          solidMechanics::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem );

  dofManager.addField( contact::traction::key(),
                       FieldLocation::Elem,
                       3,
                       meshTargets );

  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );
}

void SolidMechanicsLagrangeContactBubbleStab::setupSystem( DomainPartition & domain,
                                                           DofManager & dofManager,
                                                           CRSMatrix< real64, globalIndex > & localMatrix,
                                                           ParallelVector & rhs,
                                                           ParallelVector & solution,
                                                           bool const GEOS_UNUSED_PARAM( setSparsity ) )
{


  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, true ); // "true" is to force setSparsity

  // Create the lists of interface elements that have same type.
  createFaceTypeList( domain );

  // Create the lists of interface elements that have same type and same fracture state.
  updateStickSlipList( domain );

  // Create the list of cell elements that they are enriched with bubble functions.
  createBubbleCellList( domain );

  dofManager.setDomain( domain );
  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  // Set the sparsity pattern without the Abu and Aub blocks.
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  this->addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  this->addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/localMatrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

}

void SolidMechanicsLagrangeContactBubbleStab::implicitStepSetup( real64 const & time_n,
                                                                 real64 const & dt,
                                                                 DomainPartition & domain )
{

  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    arrayView2d< real64 > const incrBubbleDisp =
      faceManager.getField< fields::solidMechanics::incrementalBubbleDisplacement >();

    arrayView3d< real64 > const rotationMatrix =
      subRegion.getField< fields::contact::rotationMatrix >().toView();

    // Compute rotation matrices
    rotationMatrixKernel::ComputeRotationMatricesKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                                           faceNormal,
                                                                                           elemsToFaces,
                                                                                           rotationMatrix );

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [ = ]
                                      GEOS_HOST_DEVICE ( localIndex const k )
    {
      localIndex const kf0 = elemsToFaces[k][0];
      localIndex const kf1 = elemsToFaces[k][1];
      LvArray::tensorOps::fill< 3 >( incrBubbleDisp[kf0], 0.0 );
      LvArray::tensorOps::fill< 3 >( incrBubbleDisp[kf1], 0.0 );
    } );
  } );
}

void SolidMechanicsLagrangeContactBubbleStab::assembleSystem( real64 const time,
                                                              real64 const dt,
                                                              DomainPartition & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  assembleContact( domain, dofManager, localMatrix, localRhs );
}

void SolidMechanicsLagrangeContactBubbleStab::assembleContact( DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = this->getUniqueFractureRegionName();

    forFiniteElementOnStickFractureSubRegions( meshName, [&] ( string const &,
                                                               finiteElement::FiniteElementBase const & subRegionFE,
                                                               arrayView1d< localIndex const > const & faceElementList,
                                                               bool const )
    {
      solidMechanicsALMKernels::ALMSimultaneousFactory kernelFactory( dispDofNumber,
                                                                      bubbleDofNumber,
                                                                      dofManager.rankOffset(),
                                                                      localMatrix,
                                                                      localRhs,
                                                                      dt,
                                                                      faceElementList );

      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             constitutive::CoulombFriction >( mesh,
                                                              fractureRegionName,
                                                              faceElementList,
                                                              subRegionFE,
                                                              viewKeyStruct::frictionLawNameString(),
                                                              kernelFactory );

      GEOS_UNUSED_VAR( maxTraction );
    } );
  } );

  // Loop for assembling contributes of bubble elements (Abb, Abu, Aub)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );


    solidMechanicsALMKernels::ALMBubbleFactory kernelFactory( dispDofNumber,
                                                              bubbleDofNumber,
                                                              dofManager.rankOffset(),
                                                              localMatrix,
                                                              localRhs,
                                                              dt,
                                                              gravityVectorData );

    real64 maxTraction = finiteElement::
                           regionBasedKernelApplication
                         < parallelDevicePolicy< >,
                           constitutive::ElasticIsotropic,
                           CellElementSubRegion >( mesh,
                                                   regionNames,
                                                   getDiscretizationName(),
                                                   SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                   kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );

  } );

}

void SolidMechanicsLagrangeContactBubbleStab::implicitStepComplete( real64 const & time_n,
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
      arrayView2d< real64 > const deltaDispJump  = subRegion.getField< contact::deltaDispJump >();
      arrayView2d< real64 const > const & dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 > const & oldDispJump = subRegion.getField< contact::oldDispJump >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        LvArray::tensorOps::fill< 3 >( deltaDispJump[kfe], 0.0 );
        LvArray::tensorOps::fill< 3 >( deltaTraction[kfe], 0.0 );
        LvArray::tensorOps::copy< 3 >( oldDispJump[kfe], dispJump[kfe] );
      } );
    } );
  } );
}

real64 SolidMechanicsLagrangeContactBubbleStab::calculateResidualNorm( real64 const & time,
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

real64 SolidMechanicsLagrangeContactBubbleStab::calculateContactResidualNorm( DomainPartition const & domain,
                                                                              DofManager const & dofManager,
                                                                              arrayView1d< real64 const > const & localRhs )
{
  string const & dofKey = dofManager.getKey( contact::traction::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  real64 stickResidual = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                        [&]( localIndex const, FaceElementSubRegion const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();

      RAJA::ReduceSum< parallelHostReduce, real64 > stickSum( 0.0 );
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const k )
      {
        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
          for( localIndex dim = 0; dim < 3; ++dim )
          {
            real64 const norm = localRhs[localRow + dim] / area[k];
            stickSum += norm * norm;
          }
        }
      } );

      stickResidual += stickSum.get();
    } );
  } );

  stickResidual = MpiWrapper::sum( stickResidual );
  stickResidual = sqrt( stickResidual );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    std::cout << GEOS_FMT( "        ( R  ) = ( {:15.6e}  )", stickResidual );
  }

  return sqrt( stickResidual * stickResidual );
}


void SolidMechanicsLagrangeContactBubbleStab::applySystemSolution( DofManager const & dofManager,
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

void SolidMechanicsLagrangeContactBubbleStab::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  computeFaceDisplacementJump( domain );
}

real64 SolidMechanicsLagrangeContactBubbleStab::setNextDt( real64 const & currentDt,
                                                           DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
  return currentDt;
}

void SolidMechanicsLagrangeContactBubbleStab::addCouplingNumNonzeros( DomainPartition & domain,
                                                                      DofManager & dofManager,
                                                                      arrayView1d< localIndex > const & rowLengths ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );
    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    arrayView1d< globalIndex const > const dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
    {

      arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
      arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      for( localIndex bi=0; bi<bubbleElemsList.size(); ++bi )
      {
        localIndex const cellIndex = bubbleElemsList[bi];
        localIndex const k = faceElemsList[bi][0];

        localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[k] - rankOffset );

        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          for( localIndex i=0; i<3; ++i )
          {
            rowLengths[localRow + i] += numDispDof;
          }
        }

        for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
        {
          const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
          localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );

          if( localDispRow >= 0 && localDispRow < rowLengths.size() )
          {
            for( int d=0; d<3; ++d )
            {
              rowLengths[localDispRow + d] += 3;
            }
          }
        }
      }

    } );

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
    {
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
      localIndex const numDispDof = 3*numNodesPerFace;

      for( int k=0; k<2; ++k )
      {
        localIndex const kf = elemsToFaces[kfe][k];

        localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[kf] - rankOffset );

        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          for( localIndex i=0; i<3; ++i )
          {
            rowLengths[localRow + i] += numDispDof;
          }
        }

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          const localIndex & node = faceToNodeMap( kf, a );
          localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );

          if( localDispRow >= 0 && localDispRow < rowLengths.size() )
          {
            for( int d=0; d<3; ++d )
            {
              rowLengths[localDispRow + d] += 3;
            }
          }
        }
      }

    }

  } );
}

void SolidMechanicsLagrangeContactBubbleStab::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                          DofManager const & dofManager,
                                                                          SparsityPatternView< globalIndex > const & pattern ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );
    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    arrayView1d< globalIndex const > const dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    static constexpr int maxNumDispDof = 3 * 8;

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
    {

      arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
      arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      for( localIndex bi=0; bi<bubbleElemsList.size(); ++bi )
      {
        localIndex const cellIndex = bubbleElemsList[bi];
        localIndex const k = faceElemsList[bi][0];

        // working arrays
        stackArray1d< globalIndex, maxNumDispDof > eqnRowIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > eqnRowIndicesBubble( 3 );
        stackArray1d< globalIndex, maxNumDispDof > dofColIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > dofColIndicesBubble( 3 );

        for( localIndex idof = 0; idof < 3; ++idof )
        {
          eqnRowIndicesBubble[idof] = bubbleDofNumber[k] + idof - rankOffset;
          dofColIndicesBubble[idof] = bubbleDofNumber[k] + idof;
        }

        for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
        {
          const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
          for( localIndex idof = 0; idof < 3; ++idof )
          {
            eqnRowIndicesDisp[3*a + idof] = dispDofNumber[node] + idof - rankOffset;
            dofColIndicesDisp[3*a + idof] = dispDofNumber[node] + idof;
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
        {
          if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesBubble.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesBubble[j] );
            }
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesBubble.size(); ++i )
        {
          if( eqnRowIndicesBubble[i] >= 0 && eqnRowIndicesBubble[i] < pattern.numRows() )
          {
            for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesBubble[i], dofColIndicesDisp[j] );
            }
          }
        }

      }

    } );

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    static constexpr int maxNumDispFaceDof = 3 * 4;

    for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
    {

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
      localIndex const numDispDof = 3*numNodesPerFace;

      for( int k=0; k<2; ++k )
      {
        localIndex const kf = elemsToFaces[kfe][k];
        localIndex const kf_other = elemsToFaces[kfe][(1+k)%2];

        // working arrays
        stackArray1d< globalIndex, maxNumDispFaceDof > eqnRowIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > eqnRowIndicesBubble( 3 );
        stackArray1d< globalIndex, maxNumDispFaceDof > dofColIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > dofColIndicesBubble( 3 );

        for( localIndex idof = 0; idof < 3; ++idof )
        {
          eqnRowIndicesBubble[idof] = bubbleDofNumber[kf] + idof - rankOffset;
          dofColIndicesBubble[idof] = bubbleDofNumber[kf] + idof;
        }

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          const localIndex & node = faceToNodeMap( kf_other, a );
          for( localIndex idof = 0; idof < 3; ++idof )
          {
            eqnRowIndicesDisp[3*a + idof] = dispDofNumber[node] + idof - rankOffset;
            dofColIndicesDisp[3*a + idof] = dispDofNumber[node] + idof;
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
        {
          if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesBubble.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesBubble[j] );
            }
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesBubble.size(); ++i )
        {
          if( eqnRowIndicesBubble[i] >= 0 && eqnRowIndicesBubble[i] < pattern.numRows() )
          {
            for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesBubble[i], dofColIndicesDisp[j] );
            }
          }
        }

      }
    }
  } );

}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangeContactBubbleStab, string const &, Group * const )

} /* namespace geos */
