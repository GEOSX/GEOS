/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * SolidMechanicsAugmentedLagrangianContact.cpp
 */

#include "mesh/DomainPartition.hpp"
#include "SolidMechanicsAugmentedLagrangianContact.hpp"

#include "physicsSolvers/contact/SolidMechanicsALMKernels.hpp"
#include "physicsSolvers/contact/SolidMechanicsALMJumpUpdateKernels.hpp"
#include "physicsSolvers/contact/SolidMechanicsALMBubbleKernels.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsAugmentedLagrangianContact::SolidMechanicsAugmentedLagrangianContact( const string & name,
                                                                                    Group * const parent ):
  ContactSolverBase( name, parent )
{

  m_faceTypeToFiniteElements["Quadrilateral"] =  std::make_unique< finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre2 >();
  m_faceTypeToFiniteElements["Triangle"] =  std::make_unique< finiteElement::H1_TriangleFace_Lagrange1_Gauss1 >();

  // TODO Implement the MGR strategy

  // Set the default linear solver parameters
  //LinearSolverParameters & linParams = m_linearSolverParameters.get();
  //linParams.dofsPerNode = 3;
  //linParams.isSymmetric = true;
  //linParams.amg.separateComponents = true;
}

SolidMechanicsAugmentedLagrangianContact::~SolidMechanicsAugmentedLagrangianContact()
{}

void SolidMechanicsAugmentedLagrangianContact::registerDataOnMesh( dataRepository::Group & meshBodies )
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

      // Register the traction field
      subRegion.registerField< contact::traction >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the displacement jump
      subRegion.registerField< contact::dispJump >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the delta displacement jump
      subRegion.registerField< contact::deltaDispJump >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the displacement jump old
      subRegion.registerField< contact::oldDispJump >( this->getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the displacement jump
      subRegion.registerField< contact::penalty >( this->getName() ).
        reference().resizeDimension< 1 >( 2 );

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

    } );
  } );

}

void SolidMechanicsAugmentedLagrangianContact::setupDofs( DomainPartition const & domain,
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

  // The dofManager can not create the connection due to the coupling
  // between totalDisplacement and totalBubbleDisplacement
  // These connection are created using the two functions
  // addCouplingNumNonzeros and addCouplingSparsityPattern
  // dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
  //                         solidMechanics::totalBubbleDisplacement::key(),
  //                         DofManager::Connector::Elem);

}

void SolidMechanicsAugmentedLagrangianContact::setupSystem( DomainPartition & domain,
                                                            DofManager & dofManager,
                                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                                            ParallelVector & rhs,
                                                            ParallelVector & solution,
                                                            bool const setSparsity )
{

  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( setSparsity );

  // Create the list of interface elements that have same type.
  createFaceTypeList( domain );

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
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

}

void SolidMechanicsAugmentedLagrangianContact::implicitStepSetup( real64 const & time_n,
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

    arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    arrayView3d< real64 > const
    rotationMatrix = subRegion.getField< fields::contact::rotationMatrix >().toView();

    // Compute rotation matrices
    solidMechanicsALMKernels::ComputeRotationMatricesKernel::
      launch< parallelDevicePolicy<> >( subRegion.size(),
                                        faceNormal,
                                        elemsToFaces,
                                        rotationMatrix );

    arrayView2d< real64 > const
    penalty = subRegion.getField< fields::contact::penalty >().toView();

    // Set the penalty coefficients
    // TODO: This is only temporary. The setting of penalty will be adaptive in the final version.
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      penalty[k] [0] = 1.e+7;
      penalty[k] [1] = 1.e+7;
    } );

  } );

  // Set the tollerances
  computeTolerances( domain );
}

void SolidMechanicsAugmentedLagrangianContact::assembleSystem( real64 const time,
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

  //ParallelMatrix parallel_matrix;
  //parallel_matrix.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOSX );
  //parallel_matrix.write("mech.mtx");
  //abort();

  // Loop for assembling contributes from interface elements (Aut*eps^-1*Atu and Aub*eps^-1*Abu)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = this->getUniqueFractureRegionName();

    forFiniteElementOnFractureSubRegions( meshName, [&] ( string const & finiteElementName,
                                                          arrayView1d< localIndex const > const & faceElementList )
    {

      finiteElement::FiniteElementBase & subRegionFE = *(m_faceTypeToFiniteElements[finiteElementName]);

      solidMechanicsALMKernels::ALMFactory kernelFactory( dispDofNumber,
                                                          bubbleDofNumber,
                                                          dofManager.rankOffset(),
                                                          localMatrix,
                                                          localRhs,
                                                          dt,
                                                          faceElementList );

      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             constitutive::NullModel >( mesh,
                                                        fractureRegionName,
                                                        faceElementList,
                                                        subRegionFE,
                                                        "",
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

void SolidMechanicsAugmentedLagrangianContact::implicitStepComplete( real64 const & time_n,
                                                                     real64 const & dt,
                                                                     DomainPartition & domain )
{


  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();
    arrayView2d< real64 > const oldDispJump = subRegion.getField< contact::oldDispJump >();
    arrayView2d< real64 > const deltaDispJump  = subRegion.getField< contact::deltaDispJump >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      LvArray::tensorOps::fill< 3 >( deltaDispJump[kfe], 0.0 );
      LvArray::tensorOps::copy< 3 >( oldDispJump[kfe], dispJump[kfe] );
    } );

  } );

}

real64 SolidMechanicsAugmentedLagrangianContact::calculateResidualNorm( real64 const & time,
                                                                        real64 const & dt,
                                                                        DomainPartition const & domain,
                                                                        DofManager const & dofManager,
                                                                        arrayView1d< real64 const > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  real64 const solidResidualNorm = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  string const bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

  globalIndex const rankOffset = dofManager.rankOffset();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  // Bubble residual
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [ elemsToFaces, localRhs, localSum, bubbleDofNumber, rankOffset, ghostRank] GEOS_HOST_DEVICE ( localIndex const kfe )
    {

      if( ghostRank[kfe] < 0 )
      {
        for( int kk=0; kk<2; ++kk )
        {
          localIndex const k = elemsToFaces[kfe][kk];
          localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[k] - rankOffset );
          for( localIndex i = 0; i < 3; ++i )
          {
            localSum += localRhs[localRow + i] * localRhs[localRow + i];
          }
        }
      }

    } );
    real64 const localResidualNorm[2] = { localSum.get(), SolidMechanicsLagrangianFEM::getMaxForce() };

    int const rank     = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int const numRanks = MpiWrapper::commSize( MPI_COMM_GEOSX );
    array1d< real64 > globalValues( numRanks * 2 );

    // Everything is done on rank 0
    MpiWrapper::gather( localResidualNorm,
                        2,
                        globalValues.data(),
                        2,
                        0,
                        MPI_COMM_GEOSX );

    if( rank==0 )
    {
      for( int r=0; r<numRanks; ++r )
      {
        // sum/max across all ranks
        globalResidualNorm[0] += globalValues[r*2];
        globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
      }
    }

    MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );
  } );

  real64 const bubbleResidualNorm = sqrt( globalResidualNorm[0] )/(globalResidualNorm[1]+1);  // the + 1 is for the first
  // time-step when maxForce = 0;

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    std::cout << GEOS_FMT( "        ( RBubbleDisp ) = ( {:4.2e} )", bubbleResidualNorm );
  }

  return sqrt( solidResidualNorm * solidResidualNorm + bubbleResidualNorm * bubbleResidualNorm );

}

void SolidMechanicsAugmentedLagrangianContact::applySystemSolution( DofManager const & dofManager,
                                                                    arrayView1d< real64 const > const & localSolution,
                                                                    real64 const scalingFactor,
                                                                    real64 const dt,
                                                                    DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager,
                                                    localSolution,
                                                    scalingFactor,
                                                    dt,
                                                    domain );

  dofManager.addVectorToField( localSolution,
                               solidMechanics::totalBubbleDisplacement::key(),
                               solidMechanics::totalBubbleDisplacement::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               solidMechanics::totalBubbleDisplacement::key(),
                               solidMechanics::incrementalBubbleDisplacement::key(),
                               scalingFactor );


  // Loop for updating the displacement jump
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( solidMechanics::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = this->getUniqueFractureRegionName();

    CRSMatrix< real64, globalIndex > const voidMatrix;
    array1d< real64 > const voidRhs;

    forFiniteElementOnFractureSubRegions( meshName, [&] ( string const & finiteElementName,
                                                          arrayView1d< localIndex const > const & faceElementList )
    {

      finiteElement::FiniteElementBase & subRegionFE = *(m_faceTypeToFiniteElements[finiteElementName]);

      solidMechanicsALMKernels::ALMJumpUpdateFactory kernelFactory( dispDofNumber,
                                                                    bubbleDofNumber,
                                                                    dofManager.rankOffset(),
                                                                    voidMatrix.toViewConstSizes(),
                                                                    voidRhs.toView(),
                                                                    dt,
                                                                    faceElementList );

      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             constitutive::NullModel >( mesh,
                                                        fractureRegionName,
                                                        faceElementList,
                                                        subRegionFE,
                                                        "",
                                                        kernelFactory );

      GEOS_UNUSED_VAR( maxTraction );

    } );
  } );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addFields( FieldLocation::Face,
                              { solidMechanics::incrementalBubbleDisplacement::key(),
                                solidMechanics::totalBubbleDisplacement::key() } );

    fieldsToBeSync.addElementFields( { contact::dispJump::key(),
                                       contact::deltaDispJump::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );

}

void SolidMechanicsAugmentedLagrangianContact::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
}

bool SolidMechanicsAugmentedLagrangianContact::updateConfiguration( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  int hasConfigurationConverged = true;
  int condConv = true;
  // TODO: This function's design is temporary and intended solely for testing the stick mode.
  // In the final version the parallelHostPolicy will be substitute with the parallelDevicePolicy<>.
  array2d< real64 > traction_new;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {

      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();

      arrayView2d< real64 const > const deltaDispJump = subRegion.getField< contact::deltaDispJump >();
      arrayView2d< real64 const > const penalty = subRegion.getField< contact::penalty >();

      arrayView1d< real64 const > const normalDisplacementTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );
      arrayView1d< real64 > const & slidingTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );
      arrayView1d< real64 const > const & normalTractionTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );

      std::ptrdiff_t const sizes[ 2 ] = {subRegion.size(), 3};
      traction_new.resize( 2, sizes );
      arrayView2d< real64 > const traction_new_v = traction_new.toView();

      // TODO: Create a Kernel to update Traction
      forAll< parallelHostPolicy >( subRegion.size(),
                                    [ traction_new_v, traction, penalty, dispJump, deltaDispJump ] ( localIndex const kfe )
      {
        real64 eps_N = penalty[kfe][0];
        real64 eps_T = penalty[kfe][1];
        traction_new_v[kfe][0] = traction[kfe][0] + eps_N * dispJump[kfe][0];
        traction_new_v[kfe][1] = traction[kfe][1] + eps_T * deltaDispJump[kfe][1];
        traction_new_v[kfe][2] = traction[kfe][2] + eps_T * deltaDispJump[kfe][2];
      } );

      forAll< parallelHostPolicy >( subRegion.size(),
                                    [ &condConv, ghostRank, traction_new_v, normalTractionTolerance, normalDisplacementTolerance, slidingTolerance, deltaDispJump, dispJump ] ( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          if( traction_new_v[kfe][0] >= normalTractionTolerance[kfe] )
          {
            GEOS_ERROR( "Unsuported open mode! Only stick mode is supported with ALM" );
          }
          else
          {
            real64 deltaDisp = sqrt( pow( deltaDispJump[kfe][1], 2 ) + pow( deltaDispJump[kfe][2], 2 ));
            if( std::abs( dispJump[kfe][0] ) > normalDisplacementTolerance[kfe] )
            {
              //GEOS_LOG_LEVEL(2,
              //GEOS_FMT( " Element: {}, Stick mode and g_n > tol1 => compenetration, g_n: {:4.2e} tol1: {:4.2e}",
              //kfe,std::abs(dispJump[kfe][0]), normalDisplacementTolerance[kfe] ))
              condConv = false;
            }
            if( deltaDisp > slidingTolerance[kfe] )
            {
              //GEOS_LOG_LEVEL(2,
              //GEOS_FMT( " Element: {}, Stick and dg_t > tol2, dg_t: {:4.2e} tol2: {:4.2e}",
              //kfe, deltaDisp, slidingTolerance[kfe] ))
              condConv = false;
            }
          }
        }
      } );

    } );
  } );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {

      arrayView2d< real64 > const traction_new_v = traction_new.toView();
      arrayView2d< real64 > const traction = subRegion.getField< contact::traction >();
      forAll< parallelHostPolicy >( subRegion.size(), [ traction, traction_new_v ] ( localIndex const kfe )
      {
        traction[kfe][0] = traction_new_v( kfe, 0 );
        traction[kfe][1] = traction_new_v( kfe, 1 );
        traction[kfe][2] = traction_new_v( kfe, 2 );
      } );
    } );
  } );

  if( !condConv )
  {
    hasConfigurationConverged = false;
  }

  // Need to synchronize the fracture state due to the use will be made of in AssemblyStabilization
  synchronizeFractureState( domain );

  // Compute if globally the fracture state has changed
  int hasConfigurationConvergedGlobally;
  MpiWrapper::allReduce( &hasConfigurationConverged,
                         &hasConfigurationConvergedGlobally,
                         1,
                         MPI_LAND,
                         MPI_COMM_GEOSX );

  return hasConfigurationConvergedGlobally;

}

void SolidMechanicsAugmentedLagrangianContact::createFaceTypeList( DomainPartition const & domain )
{

  // Generate lists containing elements of various face types
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    array1d< localIndex > keys( subRegion.size());
    array1d< localIndex > vals( subRegion.size());
    array1d< localIndex > quadList;
    array1d< localIndex > triList;
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nTri_r( 0 );
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nQuad_r( 0 );

    arrayView1d< localIndex > const keys_v = keys.toView();
    arrayView1d< localIndex > const vals_v = vals.toView();
    // Determine the size of the lists and generate the vector keys and vals for parallel indexing into lists.
    // (With RAJA, parallelizing this operation seems the most viable approach.)
    forAll< parallelDevicePolicy<> >( subRegion.size(), [ nTri_r, nQuad_r, faceToNodeMap, keys_v, vals_v ] GEOS_HOST_DEVICE ( localIndex const kfe )
    {

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
      if( numNodesPerFace == 3 )
      {
        keys_v[kfe]=0;
        vals_v[kfe]=kfe;
        nTri_r += 1;
      }
      else if( numNodesPerFace == 4 )
      {
        keys_v[kfe]=1;
        vals_v[kfe]=kfe;
        nQuad_r += 1;
      }
      else
      {
        GEOS_ERROR( "SolidMechanicsAugmentedLagrangianContact:: invalid face type" );
      }
    } );

    localIndex nQuad = static_cast< localIndex >(nQuad_r.get());
    localIndex nTri = static_cast< localIndex >(nTri_r.get());

    // Sort vals according to keys to ensure that
    // elements of the same type are adjacent in the vals list.
    // This arrangement allows for efficient copying into the container
    // by leveraging parallelism.
    RAJA::sort_pairs< parallelDevicePolicy<> >( keys_v, vals_v );

    quadList.resize( nQuad );
    triList.resize( nTri );
    arrayView1d< localIndex > const quadList_v = quadList.toView();
    arrayView1d< localIndex > const triList_v = triList.toView();

    forAll< parallelDevicePolicy<> >( nTri, [triList_v, vals_v] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      triList_v[kfe] = vals_v[kfe];
    } );

    forAll< parallelDevicePolicy<> >( nQuad, [quadList_v, vals_v, nTri] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      quadList_v[kfe] = vals_v[nTri+kfe];
    } );

    this->m_faceTypesToFaceElements[meshName]["Quadrilateral"] =  quadList;
    this->m_faceTypesToFaceElements[meshName]["Triangle"] =  triList;
  } );

}

void SolidMechanicsAugmentedLagrangianContact::createBubbleCellList( DomainPartition & domain ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    array1d< localIndex > tmpSpace( 2*subRegion.size());
    SortedArray< localIndex > faceIdList;

    arrayView1d< localIndex > const tmpSpace_v = tmpSpace.toView();
    // Store indexes of faces in the temporany array.
    {
      ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [ tmpSpace_v, elemsToFaces ] GEOS_HOST_DEVICE ( localIndex const kfe )
      {

        localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
        tmpSpace_v[2*kfe] = kf0, tmpSpace_v[2*kfe+1] = kf1;

      } );
    }

    // Sort indexes to enable efficient searching using binary search.
    RAJA::stable_sort< parallelDevicePolicy<> >( tmpSpace_v );
    faceIdList.insert( tmpSpace_v.begin(), tmpSpace_v.end());

    // Search for bubble element on each CellElementSubRegion and
    // store element indexes, global and local face indexes.
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion & cellElementSubRegion )
    {

      arrayView2d< localIndex const > const elemsToFaces = cellElementSubRegion.faceList().toViewConst();

      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nBubElems_r( 0 );

      localIndex const n_max = cellElementSubRegion.size() * elemsToFaces.size( 1 );
      array1d< localIndex > keys( n_max );
      array1d< localIndex > perms( n_max );
      array1d< localIndex > vals( n_max );
      array1d< localIndex > localFaceIds( n_max );

      arrayView1d< localIndex > const keys_v = keys.toView();
      arrayView1d< localIndex > const perms_v = perms.toView();
      arrayView1d< localIndex > const vals_v = vals.toView();
      arrayView1d< localIndex > const localFaceIds_v = localFaceIds.toView();
      SortedArrayView< localIndex const > const faceIdList_v = faceIdList.toViewConst();

      forAll< parallelDevicePolicy<> >( cellElementSubRegion.size(),
                                        [ keys_v, vals_v, perms_v, localFaceIds_v, nBubElems_r, elemsToFaces, faceIdList_v ]
                                        GEOS_HOST_DEVICE ( localIndex const kfe )
      {
        for( int i=0; i < elemsToFaces.size( 1 ); ++i )
        {
          perms_v[kfe*elemsToFaces.size( 1 )+i] = kfe*elemsToFaces.size( 1 )+i;
          if( faceIdList_v.contains( elemsToFaces[kfe][i] ))
          {
            keys_v[kfe*elemsToFaces.size( 1 )+i] = 0;
            vals_v[kfe*elemsToFaces.size( 1 )+i] = kfe;
            localFaceIds_v[kfe*elemsToFaces.size( 1 )+i] = i;
            nBubElems_r += 1;
          }
          else
          {
            keys_v[kfe*elemsToFaces.size( 1 )+i] = 1;
            vals_v[kfe*elemsToFaces.size( 1 )+i] = -1;
            localFaceIds_v[kfe*elemsToFaces.size( 1 )+i] = -1;
          }
        }
      } );

      // Sort perms according to keys to ensure that bubble elements are adjacent
      // and occupy the first positions of the list.
      // This arrangement allows for efficient copying into the container
      // by leveraging parallelism.
      localIndex nBubElems = static_cast< localIndex >(nBubElems_r.get());
      RAJA::sort_pairs< parallelDevicePolicy<> >( keys_v, perms_v );

      array1d< localIndex > bubbleElemsList;
      bubbleElemsList.resize( nBubElems );

      arrayView1d< localIndex > const bubbleElemsList_v = bubbleElemsList.toView();

      forAll< parallelDevicePolicy<> >( n_max, [ keys_v, vals_v, perms_v ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys_v[k] = vals_v[perms_v[k]];
      } );

      forAll< parallelDevicePolicy<> >( nBubElems, [ bubbleElemsList_v, keys_v ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        bubbleElemsList_v[k] = keys_v[k];
      } );
      cellElementSubRegion.setBubbleElementsList( bubbleElemsList.toViewConst());

      forAll< parallelDevicePolicy<> >( n_max, [ keys_v, localFaceIds_v, perms_v ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys_v[k] = localFaceIds_v[perms_v[k]];
      } );

      array2d< localIndex > faceElemsList;
      faceElemsList.resize( nBubElems, 2 );

      arrayView2d< localIndex > const faceElemsList_v = faceElemsList.toView();

      forAll< parallelDevicePolicy<> >( nBubElems,
                                        [ faceElemsList_v, keys_v, elemsToFaces, bubbleElemsList_v ]
                                        GEOS_HOST_DEVICE ( localIndex const k )
      {
        localIndex const kfe =  bubbleElemsList_v[k];
        faceElemsList_v[k][0] = elemsToFaces[kfe][keys_v[k]];
        faceElemsList_v[k][1] = keys_v[k];
      } );
      cellElementSubRegion.setFaceElementsList( faceElemsList.toViewConst());

    } );

  } );

}

void SolidMechanicsAugmentedLagrangianContact::addCouplingNumNonzeros( DomainPartition & domain,
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

void SolidMechanicsAugmentedLagrangianContact::addCouplingSparsityPattern( DomainPartition const & domain,
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

// TODO: Is it possible to define this method once? Similar to SolidMechanicsLagrangeContact
void SolidMechanicsAugmentedLagrangianContact::computeTolerances( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

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
        arrayView3d< real64 const > const & faceRotationMatrix = subRegion.getField< fields::contact::rotationMatrix >().toView();
        ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

        arrayView1d< real64 > const & normalTractionTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
        arrayView1d< real64 > const & normalDisplacementTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );
        arrayView1d< real64 > const & slidingTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

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
            normalDisplacementTolerance[kfe] = rotatedInvStiffApprox[ 0 ][ 0 ] * averageYoungModulus / 2.e+8;
            slidingTolerance[kfe] = sqrt( rotatedInvStiffApprox[ 1 ][ 1 ] * rotatedInvStiffApprox[ 1 ][ 1 ] +
                                          rotatedInvStiffApprox[ 2 ][ 2 ] * rotatedInvStiffApprox[ 2 ][ 2 ] ) * averageYoungModulus / 2.e+8;
            normalTractionTolerance[kfe] = 1.0 / 2.0 * averageConstrainedModulus / averageBoxSize0 * normalDisplacementTolerance[kfe];
          }
        } );
      }
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsAugmentedLagrangianContact, string const &, Group * const )
} /* namespace geos */
