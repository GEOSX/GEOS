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


/**
 * @file AcousticPOD.cpp
 */

#include "AcousticPOD.hpp"
#include "AcousticPODKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticTimeSchemeSEMKernel.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp"
#include "events/EventManager.hpp"
#include "physicsSolvers/wavePropagation/shared/PrecomputeSourcesAndReceiversKernel.hpp"

namespace geos
{

using namespace dataRepository;

AcousticPOD::AcousticPOD( const std::string & name,
                          Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::computePODmatrixString(), &m_computePODmatrix ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Bool for POD matrices computation" );

  registerWrapper( viewKeyStruct::computeSourceValueString(), &m_computeSourceValue ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 1 ).
    setDescription( "Bool for POD source and receivers computation" );

  registerWrapper( viewKeyStruct::massPODString(), &m_massPOD ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Mass POD forward matrix" );

  registerWrapper( viewKeyStruct::massGradientPODString(), &m_massGradientPOD ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Mass gradient POD forward matrix" );

  registerWrapper( viewKeyStruct::dampingPODString(), &m_dampingPOD ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Damping POD forward matrix" );

  registerWrapper( viewKeyStruct::invAPODString(), &m_invAPOD ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Inverse scalar product POD matrix" );

  registerWrapper( viewKeyStruct::sourceConstantsPODString(), &m_sourceConstantsPOD ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "sourceConstants POD vector" );

  registerWrapper( viewKeyStruct::a_np1String(), &m_a_np1 ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "POD basis coefficients at time n+1" );

  registerWrapper( viewKeyStruct::a_nString(), &m_a_n ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "POD basis coefficients at time n" );

  registerWrapper( viewKeyStruct::a_nm1String(), &m_a_nm1 ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "POD basis coefficients at time n-1" );

  registerWrapper( viewKeyStruct::perturbationString(), &m_perturbation ).
    setInputFlag( InputFlags::FALSE ).
    setDefaultValue( 0 ).
    setDescription( "POD model perturbation coefficient" );

  registerWrapper( viewKeyStruct::sizePODString(), &m_sizePOD ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "POD basis size" );

  registerWrapper( viewKeyStruct::orderInitString(), &m_orderInit ).
    setInputFlag( InputFlags::FALSE ).
    setDefaultValue( 0 ).
    setDescription( "Taylor order for initial condition computation" );
}

AcousticPOD::~AcousticPOD()
{
  // TODO Auto-generated destructor stub
}

void AcousticPOD::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void AcousticPOD::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< fields::ForcingRHS,
                               fields::AcousticMassVector,
			       fields::AcousticMassGradientVector,
                               fields::AcousticDampingVector,
                               fields::AcousticFreeSurfaceNodeIndicator >( getName() );


    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::AcousticFreeSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::AcousticVelocity >( getName() );
      subRegion.registerField< fields::AcousticDensity >( getName() );
      subRegion.registerField< fields::PartialGradient >( getName() );
    } );

  } );
}

void AcousticPOD::postInputInitialization()
{
  WaveSolverBase::postInputInitialization();

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );

  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
				  [&] ( string const &,
					MeshLevel & mesh,
					arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< fields::AcousticFreeSurfaceFaceIndicator >();
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const
      nodeCoords32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< fields::AcousticMassVector >();
    {
      GEOS_MARK_SCOPE( mass_zero );
      mass.zero();
    }
    arrayView1d< real32 > const massGradient = nodeManager.getField< fields::AcousticMassGradientVector >();
    {
      GEOS_MARK_SCOPE( massGradient_zero );
      massGradient.zero();
    }
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping = nodeManager.getField< fields::AcousticDampingVector >();
    {
      GEOS_MARK_SCOPE( damping_zero );
      damping.zero();
    }


    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
											  CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      arrayView1d< real32 const > const velocity = elementSubRegion.getField< fields::AcousticVelocity >();
      arrayView1d< real32 const > const gradient = elementSubRegion.getField< fields::PartialGradient >();

      finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
	using FE_TYPE = TYPEOFREF( finiteElement );

	acousticPODKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );
	kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
							       nodeCoords32,
							       elemsToNodes,
							       velocity,
							       gradient,
							       mass,
							       massGradient);

	acousticPODKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );
	kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
							       nodeCoords32,
							       elemsToFaces,
							       facesToNodes,
							       facesDomainBoundaryIndicator,
							       freeSurfaceFaceIndicator,
							       velocity,
							       damping);
      } );
    } );
  } );


  if( m_computePODmatrix )
  {
    DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
    /// Compute Stiffness matrix POD
    forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                    [&] ( string const &,
                                          MeshLevel & mesh,
                                          arrayView1d< string const > const & regionNames )
    {
      NodeManager const & nodeManager = mesh.getNodeManager();
      FaceManager const & faceManager = mesh.getFaceManager();

      arrayView1d< real32 const > const mass = nodeManager.getField< fields::AcousticMassVector >();
      arrayView1d< real32 const > const massGradient = nodeManager.getField< fields::AcousticMassGradientVector >();
      arrayView1d< real32 const > const damping = nodeManager.getField< fields::AcousticDampingVector >();

      int sizePOD = m_sizePOD;
      arrayView1d< localIndex const > const nodesGhostRank = nodeManager.ghostRank();

      if( m_forward )
      {
        m_massPOD.resize( sizePOD, sizePOD );
	m_massGradientPOD.resize( sizePOD, sizePOD );
        m_dampingPOD.resize( sizePOD, sizePOD );
        m_massPOD.zero();
        m_dampingPOD.zero();
      }

      m_a_np1.resize( sizePOD );
      m_a_n.resize( sizePOD );
      m_a_nm1.resize( sizePOD );
      m_rhsPOD.resize( sizePOD );

      m_a_n.zero();
      m_a_nm1.zero();
      m_a_np1.zero();
      m_rhsPOD.zero();

      m_invAPOD.resize( sizePOD, sizePOD );
      m_invAPOD.zero();

      arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const
	nodeCoords32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            CellElementSubRegion & elementSubRegion )
      {
        GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                        "Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8), using the SEM formulation",
                        InputError );

        if( m_forward )
        {
          arrayView2d< real32 > const dampingPOD = m_dampingPOD.toView();
          arrayView2d< real32 > const massPOD = m_massPOD.toView();
	  arrayView2d< real32 > const massGradientPOD = m_massGradientPOD.toView();
	  int const countPhi = m_sizePOD;
	  int const shotIndex = m_shotIndex;

	  computeMassAndDampingPOD( massPOD,
				    massGradientPOD,
				    dampingPOD,
				    mass,
				    massGradient,
				    damping,
				    countPhi,
				    shotIndex,
				    nodesGhostRank );
	}
      } );
    } );
  }
}

void AcousticPOD::precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodeCoords = baseMesh.getNodeManager().referencePosition();

  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  arrayView2d< real64 const > const faceCenter  = faceManager.faceCenter();

  int sizePOD = m_sizePOD;

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();

  m_sourceConstantsPOD.resize( sourceCoordinates.size( 0 ), sizePOD );

  arrayView2d< real64 > const sourceConstantsPOD = m_sourceConstantsPOD.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  localIndex const computeSourceValue = m_computeSourceValue;
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstantsPOD.zero();
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();

  m_receiverConstantsPOD.resize( receiverCoordinates.size( 0 ), sizePOD );

  arrayView2d< real64 > const receiverConstantsPOD = m_receiverConstantsPOD.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< EXEC_POLICY >( -1 );
  receiverConstantsPOD.zero();
  receiverIsLocal.zero();

  arrayView2d< real32 > const sourceValue = m_sourceValue.toView();
  if( computeSourceValue )
  {
    sourceValue.zero();
  }

  real64 dt = 0;
  EventManager const & event = this->getGroupByPath< EventManager >( "/Problem/Events" );
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + this->getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & elementSubRegion )
  {
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                    "Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8), using the SEM formulation",
                    InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();
    int const countPhi = m_sizePOD;
    int const shotIndex = m_shotIndex;

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      acousticPODKernels::
        PrecomputeSourceAndReceiverKernel::
        launch< EXEC_POLICY, FE_TYPE >
        ( elementSubRegion.size(),
	  numNodesPerElem,
	  numFacesPerElem,
	  nodeCoords,
	  elemGhostRank,
	  elemsToNodes,
	  elemsToFaces,
	  elemCenter,
	  faceNormal,
	  faceCenter,
	  sourceCoordinates,
	  sourceIsAccessible,
	  sourceNodeIds,
	  sourceConstantsPOD,
	  computeSourceValue,
	  receiverCoordinates,
	  receiverIsLocal,
	  receiverNodeIds,
	  receiverConstantsPOD,
	  sourceValue,
	  countPhi,
	  shotIndex,
	  dt,
	  m_timeSourceFrequency,
	  m_timeSourceDelay,
	  m_rickerOrder );
    } );
  } );

  sourceConstantsPOD.move( LvArray::MemorySpace::host, true );
  MpiWrapper::allReduce( sourceConstantsPOD.data(),
                         sourceConstantsPOD.data(),
                         sourceConstantsPOD.size( 0 )*sourceConstantsPOD.size( 1 ),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOS );

  if( computeSourceValue )
  {
    sourceValue.move( LvArray::MemorySpace::host, true );
    MpiWrapper::allReduce( sourceValue.data(),
			   sourceValue.data(),
			   sourceValue.size( 0 )*sourceValue.size( 1 ),
			   MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
			   MPI_COMM_GEOS );
  }
}

void AcousticPOD::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstantsPOD.toViewConst();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real32 const > const sourceValue   = m_sourceValue.toViewConst();

  GEOS_THROW_IF( cycleNumber > sourceValue.size( 0 ), "Too many steps compared to array size", std::runtime_error );

  forAll< EXEC_POLICY >( sourceConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const isrc )
    {
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        real32 const localIncrement = sourceConstants[isrc][inode] * sourceValue[cycleNumber][isrc];
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhs[inode], localIncrement );
      }
    } );
}

void AcousticPOD::initializePostInitialConditionsPreSubGroups()
{

  GEOS_MARK_FUNCTION;
  {
    GEOS_MARK_SCOPE( WaveSolverBase::initializePostInitialConditionsPreSubGroups );
    WaveSolverBase::initializePostInitialConditionsPreSubGroups();
  }
  if( m_usePML )
  {
    AcousticPOD::initializePML();
  }

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    if(m_computePODmatrix)
    {
      MeshLevel & baseMesh = domain.getMeshBodies().getGroup< MeshBody >( meshBodyName ).getBaseDiscretization();
      precomputeSourceAndReceiverTerm( baseMesh, mesh, regionNames );
    }
  } );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );
}


void AcousticPOD::applyFreeSurfaceBC( real64 time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< fields::AcousticFreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< fields::AcousticFreeSurfaceNodeIndicator >();

  //freeSurfaceFaceIndicator.zero();
  //freeSurfaceNodeIndicator.zero();

  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "FreeSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        freeSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          freeSurfaceNodeIndicator[dof] = 1;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}


/*** Checks if a directory exists.
* @param dirName Directory name to check existence of.
* @return true is dirName exists and is a directory.
*/
bool dirExistsPOD( const std::string & dirName )
{
  struct stat buffer;
  return stat( dirName.c_str(), &buffer ) == 0;
}


real64 AcousticPOD::explicitStepForward( real64 const & time_n,
                                         real64 const & dt,
                                         integer cycleNumber,
                                         DomainPartition & domain,
                                         bool computeGradient )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & GEOS_UNUSED_PARAM ( regionNames ) )
  {
    arrayView1d< real32 > const a_nm1 = m_a_nm1.toView();
    arrayView1d< real32 > const a_n = m_a_n.toView();
    arrayView1d< real32 > const a_np1 = m_a_np1.toView();

    forAll< EXEC_POLICY >( a_n.size(), [=] GEOS_HOST_DEVICE ( localIndex const m )
      {
        a_nm1[m] = a_n[m];
        a_n[m]   = a_np1[m];
      } );
  } );

  return dtOut;
}

real64 AcousticPOD::explicitStepInternal( real64 const & time_n,
                                          real64 const & dt,
                                          integer cycleNumber,
                                          DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    computeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );
    synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );
  } );
  return dt;
}

void AcousticPOD::computeUnknowns( real64 const & time_n,
				   real64 const & dt,
				   integer cycleNumber,
				   DomainPartition & domain,
				   MeshLevel & mesh,
				   arrayView1d< string const > const & regionNames )
{
  arrayView1d< real32 > const rhs = m_rhsPOD.toView();

  arrayView1d< real32 > const a_nm1 = m_a_nm1.toView();
  arrayView1d< real32 > const a_n = m_a_n.toView();
  arrayView1d< real32 > const a_np1 = m_a_np1.toView();

  EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & minTime = event.getReference< real64 >( EventManager::viewKeyStruct::minTimeString() );
  integer const cycleForSource = int(round( -minTime / dt + cycleNumber ));
  addSourceToRightHandSide( cycleForSource, rhs );

  /// calculate your time integrators
  real64 const dt2 = dt*dt;

  GEOS_MARK_SCOPE ( updateP );

  //if( cycleNumber == 0 )
  //{
  //  computeInitialConditions();
  //}
  //else
  //{
    arrayView2d< real32 const > const damping = m_dampingPOD.toViewConst();
    arrayView2d< real32 const > const mass = m_massPOD.toViewConst();
    arrayView2d< real32 const > const invA = m_invAPOD.toViewConst();
    array1d< real32 > b( a_n.size());
    arrayView1d< real32 > bV = b.toView();

    forAll< EXEC_POLICY >( a_n.size(), [=] GEOS_HOST_DEVICE ( localIndex const m )
    {
      bV[m] = dt2 * ( rhs[m] - a_n[m]);
      for( localIndex n = 0; n < a_n.size(); ++n )
      {
	bV[m] += 2 * mass[m][n] * a_n[n];
	bV[m] -= (mass[m][n] - dt * 0.5 * damping[m][n]) * a_nm1[n];
      }
      rhs[m] = 0.0;
    } );

    forAll< EXEC_POLICY >( a_n.size(), [=] GEOS_HOST_DEVICE ( localIndex const m )
    {
      a_np1[m] = 0.0;
      for( localIndex n = 0; n < a_n.size(); ++n )
      {
	a_np1[m] += invA[m][n] * bV[n];
      }
    } );
    //}
}

void AcousticPOD::computeInitialConditions()
{
  arrayView1d< real32 > const a_n = m_a_n;
  arrayView1d< real32 > const a_np1 = m_a_np1;
  int const size = a_np1.size();
  array1d< real32 > b(size);
  arrayView1d< real32 > bV = b.toView();

  localIndex const shotIndex = m_shotIndex;
  real32 const alpha = m_perturbation;
  localIndex const ord = m_orderInit;

  real32 fact = 1;
  for(int i=0; i<ord+1; ++i)
  {
    if( i>0 )
    {
      fact *= i;
    }

    std::string fileName1 = GEOS_FMT( "phi/shot_{:05}/initialConditions/order_{:02}/vector_00000.dat", shotIndex, i);
    std::ifstream wf1( fileName1, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf1,
		   getDataContext() << ": Could not open file "<< fileName1 << " for reading",
		   InputError );
    bV.move( LvArray::MemorySpace::host, true );
    wf1.read( (char *)&bV[0], size*sizeof( real32 ) );
    wf1.close( );

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const j )
    {
      a_n[j] += pow(alpha,i) * bV[j] / fact;
    } );

    std::string fileName2 = GEOS_FMT( "phi/shot_{:05}/initialConditions/order_{:02}/vector_00001.dat", shotIndex, i);
    std::ifstream wf2( fileName2, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf2,
                   getDataContext() << ": Could not open file "<< fileName2 << " for reading",
                   InputError );
    bV.move( LvArray::MemorySpace::host, true );
    wf2.read( (char *)&bV[0], size*sizeof( real32 ) );
    wf2.close( );

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const j )
    {
      a_np1[j] += pow(alpha,i) * bV[j] / fact;
    } );
  }
}

void AcousticPOD::synchronizeUnknowns( real64 const & time_n,
				       real64 const & dt,
				       integer const,
				       DomainPartition & domain,
				       MeshLevel & mesh,
				       arrayView1d< string const > const & )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const a_n = m_a_n.toView();
  arrayView1d< real32 > const a_np1 = m_a_np1.toView();

  /// synchronize pressure fields
  FieldIdentifiers fieldsToBeSync;

  if( m_usePML )
  {
    fieldsToBeSync.addFields( FieldLocation::Node, {
        fields::AuxiliaryVar1PML::key(),
        fields::AuxiliaryVar4PML::key() } );
  }

  /*CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldsToBeSync,
                                mesh,
                                domain.getNeighbors(),
                                true );
  */
  /// compute the seismic traces since last step.
  arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();

  computeAllSeismoTraces( time_n, dt, a_np1, a_n, pReceivers );
  incrementIndexSeismoTrace( time_n );

  if( m_usePML )
  {
    arrayView2d< real32 > const grad_n = nodeManager.getField< fields::AuxiliaryVar2PML >();
    arrayView1d< real32 > const divV_n = nodeManager.getField< fields::AuxiliaryVar3PML >();
    grad_n.zero();
    divV_n.zero();
  }
}


void AcousticPOD::initializePML()
{
  GEOS_MARK_FUNCTION;

  registerWrapper< parametersPML >( viewKeyStruct::parametersPMLString() ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Parameters needed to compute damping in the PML region" );

  parametersPML & param = getReference< parametersPML >( viewKeyStruct::parametersPMLString() );

  /// Get the default thicknesses and wave speeds in the PML regions from the PerfectlyMatchedLayer
  /// field specification parameters (from the xml)
  real32 minThicknessPML=0;
  real32 smallestXMinPML=0;
  real32 largestXMaxPML=0;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  fsManager.forSubGroups< PerfectlyMatchedLayer >( [&] ( PerfectlyMatchedLayer const & fs )
  {
    param.xMinPML=fs.getMin();
    param.xMaxPML=fs.getMax();
    param.thicknessMinXYZPML=fs.getThicknessMinXYZ();
    param.thicknessMaxXYZPML=fs.getThicknessMaxXYZ();
    param.reflectivityPML = fs.getReflectivity();
    param.waveSpeedMinXYZPML=fs.getWaveSpeedMinXYZ();
    param.waveSpeedMaxXYZPML=fs.getWaveSpeedMaxXYZ();
    minThicknessPML=fs.minThickness;
    smallestXMinPML=fs.smallestXMin;
    largestXMaxPML=fs.largestXMax;
  } );

  /// Now compute the PML parameters above internally
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    NodeManager & nodeManager = mesh.getNodeManager();
    /// WARNING: the array below is one of the PML auxiliary variables
    arrayView1d< real32 > const indicatorPML = nodeManager.getField< fields::AuxiliaryVar4PML >();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const
      X = nodeManager.getField< fields::referencePosition32 >().toViewConst();
    indicatorPML.zero();

    real32 xInteriorMin[3]{};
    real32 xInteriorMax[3]{};
    real32 xGlobalMin[3]{};
    real32 xGlobalMax[3]{};
    real32 cMin[3]{};
    real32 cMax[3]{};
    integer counterMin[3]{};
    integer counterMax[3]{};

    /// Set a node-based flag in the PML regions
    /// WARNING: the array used as a flag is one of the PML
    /// auxiliary variables to save memory
    fsManager.apply< ElementSubRegionBase,
                     PerfectlyMatchedLayer >( 0.0,
                                              mesh,
                                              PerfectlyMatchedLayer::catalogName(),
                                              [&]( PerfectlyMatchedLayer const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )

    {
      CellElementSubRegion::NodeMapType const & elemToNodes =
        subRegion.getReference< CellElementSubRegion::NodeMapType >( CellElementSubRegion::viewKeyStruct::nodeListString() );
      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst = elemToNodes.toViewConst();

      forAll< EXEC_POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const l )
        {
          localIndex const k = targetSet[ l ];
          localIndex const numNodesPerElem = elemToNodesViewConst[k].size();

          for( localIndex i=0; i<numNodesPerElem; ++i )
          {
            indicatorPML[elemToNodesViewConst[k][i]]=1.0;
          }
        } );
    } );


    /// find the interior and global coordinates limits
    RAJA::ReduceMin< parallelDeviceReduce, real32 > xMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > yMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > zMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > xMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > yMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > zMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > xMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > yMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > zMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > xMaxInterior( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > yMaxInterior( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > zMaxInterior( -LvArray::NumericLimits< real32 >::max );

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        xMinGlobal.min( X[a][0] );
        yMinGlobal.min( X[a][1] );
        zMinGlobal.min( X[a][2] );
        xMaxGlobal.max( X[a][0] );
        yMaxGlobal.max( X[a][1] );
        zMaxGlobal.max( X[a][2] );
        if( !isZero( indicatorPML[a] - 1.0 ))
        {
          xMinInterior.min( X[a][0] );
          yMinInterior.min( X[a][1] );
          zMinInterior.min( X[a][2] );
          xMaxInterior.max( X[a][0] );
          yMaxInterior.max( X[a][1] );
          zMaxInterior.max( X[a][2] );
        }
      } );

    xGlobalMin[0] = xMinGlobal.get();
    xGlobalMin[1] = yMinGlobal.get();
    xGlobalMin[2] = zMinGlobal.get();
    xGlobalMax[0] = xMaxGlobal.get();
    xGlobalMax[1] = yMaxGlobal.get();
    xGlobalMax[2] = zMaxGlobal.get();
    xInteriorMin[0] = xMinInterior.get();
    xInteriorMin[1] = yMinInterior.get();
    xInteriorMin[2] = zMinInterior.get();
    xInteriorMax[0] = xMaxInterior.get();
    xInteriorMax[1] = yMaxInterior.get();
    xInteriorMax[2] = zMaxInterior.get();

    for( integer i=0; i<3; ++i )
    {
      xGlobalMin[i] = MpiWrapper::min( xGlobalMin[i] );
      xGlobalMax[i] = MpiWrapper::max( xGlobalMax[i] );
      xInteriorMin[i] = MpiWrapper::min( xInteriorMin[i] );
      xInteriorMax[i] = MpiWrapper::max( xInteriorMax[i] );
    }


    /// if the coordinates limits and PML thicknesses are not provided
    /// from the xml, replace them with the above
    for( integer i=0; i<3; ++i )
    {
      if( param.xMinPML[i]<smallestXMinPML )
        param.xMinPML[i] = xInteriorMin[i];
      if( param.xMaxPML[i]>largestXMaxPML )
        param.xMaxPML[i] = xInteriorMax[i];
      if( param.thicknessMinXYZPML[i]<0 )
        param.thicknessMinXYZPML[i] = xInteriorMin[i]-xGlobalMin[i];
      if( param.thicknessMaxXYZPML[i]<0 )
        param.thicknessMaxXYZPML[i] = xGlobalMax[i]-xInteriorMax[i];
    }

    /// Compute the average wave speeds in the PML regions internally
    /// using the actual velocity field
    fsManager.apply< ElementSubRegionBase,
                     PerfectlyMatchedLayer >( 0.0,
                                              mesh,
                                              PerfectlyMatchedLayer::catalogName(),
                                              [&]( PerfectlyMatchedLayer const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )

    {
      CellElementSubRegion::NodeMapType const & elemToNodes =
        subRegion.getReference< CellElementSubRegion::NodeMapType >( CellElementSubRegion::viewKeyStruct::nodeListString() );
      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst = elemToNodes.toViewConst();
      arrayView1d< real32 const > const vel = subRegion.getReference< array1d< real32 > >( fields::AcousticVelocity::key());
      finiteElement::FiniteElementBase const &
      fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      real32 const xMin[3]{param.xMinPML[0], param.xMinPML[1], param.xMinPML[2]};
      real32 const xMax[3]{param.xMaxPML[0], param.xMaxPML[1], param.xMaxPML[2]};

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticPODKernels::
          waveSpeedPMLKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( targetSet,
          X,
          elemToNodesViewConst,
          vel,
          xMin,
          xMax,
          cMin,
          cMax,
          counterMin,
          counterMax );
      } );
    } );

    for( integer i=0; i<3; ++i )
    {
      cMin[i] = MpiWrapper::sum( cMin[i] );
      cMax[i] = MpiWrapper::sum( cMax[i] );
      counterMin[i] = MpiWrapper::sum( counterMin[i] );
      counterMax[i] = MpiWrapper::sum( counterMax[i] );
    }
    for( integer i=0; i<3; ++i )
    {
      cMin[i] /= std::max( 1, counterMin[i] );
      cMax[i] /= std::max( 1, counterMax[i] );
    }

    /// if the PML wave speeds are not provided from the xml
    /// replace them with the above
    for( integer i=0; i<3; ++i )
    {
      if( param.waveSpeedMinXYZPML[i]<0 )
        param.waveSpeedMinXYZPML[i] = cMin[i];
      if( param.waveSpeedMaxXYZPML[i]<0 )
        param.waveSpeedMaxXYZPML[i] = cMax[i];
    }

    /// add safeguards when PML thickness is negative or too small
    for( integer i=0; i<3; ++i )
    {
      if( param.thicknessMinXYZPML[i]<=minThicknessPML )
      {
        param.thicknessMinXYZPML[i]=LvArray::NumericLimits< real32 >::max;
        param.waveSpeedMinXYZPML[i]=0;
      }
      if( param.thicknessMaxXYZPML[i]<=minThicknessPML )
      {
        param.thicknessMaxXYZPML[i]=LvArray::NumericLimits< real32 >::max;
        param.waveSpeedMaxXYZPML[i]=0;
      }
    }

    /// WARNING: don't forget to reset the indicator to zero
    /// so it can be used by the PML application
    indicatorPML.zero();

    GEOS_LOG_LEVEL_RANK_0( 1,
                            "PML parameters are: \n"
                            << "\t inner boundaries xMin = "<<param.xMinPML<<"\n"
                            << "\t inner boundaries xMax = "<<param.xMaxPML<<"\n"
                            << "\t left, front, top max PML thicknesses  = "<<param.thicknessMinXYZPML<<"\n"
                            << "\t right, back, bottom max PML thicknesses  = "<<param.thicknessMaxXYZPML<<"\n"
                            << "\t left, front, top average wave speed  = "<<param.waveSpeedMinXYZPML<<"\n"
                            << "\t right, back, bottom average wave speed  = "<<param.waveSpeedMaxXYZPML<<"\n"
                            << "\t theoretical reflectivity = "<< param.reflectivityPML );

  } );
}



void AcousticPOD::applyPML( real64 const time, DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  parametersPML const & param = getReference< parametersPML >( viewKeyStruct::parametersPMLString() );

  /// Loop over the different mesh bodies; for wave propagation, there is only one mesh body
  /// which is the whole mesh
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    NodeManager & nodeManager = mesh.getNodeManager();

    /// Array views of the pressure p, PML auxiliary variables, and node coordinates
    arrayView1d< real32 const > const p_n = nodeManager.getField< fields::Pressure_n >();
    arrayView2d< real32 const > const v_n = nodeManager.getField< fields::AuxiliaryVar1PML >();
    arrayView2d< real32 > const grad_n = nodeManager.getField< fields::AuxiliaryVar2PML >();
    arrayView1d< real32 > const divV_n = nodeManager.getField< fields::AuxiliaryVar3PML >();
    arrayView1d< real32 const > const u_n = nodeManager.getField< fields::AuxiliaryVar4PML >();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const
      X = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// Select the subregions concerned by the PML (specified in the xml by the Field Specification)
    /// 'targetSet' contains the indices of the elements in a given subregion
    fsManager.apply< ElementSubRegionBase,
                     PerfectlyMatchedLayer >( time,
                                              mesh,
                                              PerfectlyMatchedLayer::catalogName(),
                                              [&]( PerfectlyMatchedLayer const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )

    {

      /// Get the element to nodes mapping in the subregion
      CellElementSubRegion::NodeMapType const & elemToNodes =
        subRegion.getReference< CellElementSubRegion::NodeMapType >( CellElementSubRegion::viewKeyStruct::nodeListString() );

      /// Get a const ArrayView of the mapping above
      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst = elemToNodes.toViewConst();

      /// Array view of the wave speed
      arrayView1d< real32 const > const vel = subRegion.getReference< array1d< real32 > >( fields::AcousticVelocity::key());

      /// Get the object needed to determine the type of the element in the subregion
      finiteElement::FiniteElementBase const &
      fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      real32 xMin[3];
      real32 xMax[3];
      real32 dMin[3];
      real32 dMax[3];
      real32 cMin[3];
      real32 cMax[3];
      for( integer i=0; i<3; ++i )
      {
        xMin[i] = param.xMinPML[i];
        xMax[i] = param.xMaxPML[i];
        dMin[i] = param.thicknessMinXYZPML[i];
        dMax[i] = param.thicknessMaxXYZPML[i];
        cMin[i] = param.waveSpeedMinXYZPML[i];
        cMax[i] = param.waveSpeedMaxXYZPML[i];
      }
      real32 const r = param.reflectivityPML;

      /// Get the type of the elements in the subregion
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        /// apply the PML kernel
        acousticPODKernels::
          PMLKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( targetSet,
          X,
          elemToNodesViewConst,
          vel,
          p_n,
          v_n,
          u_n,
          xMin,
          xMax,
          dMin,
          dMax,
          cMin,
          cMax,
          r,
          grad_n,
          divV_n );
      } );
    } );
  } );

}

void AcousticPOD::cleanup( real64 const time_n,
                           integer const cycleNumber,
                           integer const eventCounter,
                           real64 const eventProgress,
                           DomainPartition & domain )
{
  // call the base class cleanup (for reporting purposes)
  PhysicsSolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    arrayView1d< real32 const > const a_n = m_a_n.toView();
    arrayView1d< real32 const > const a_np1 =  m_a_np1.toView();
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0, a_np1, a_n, pReceivers );
  } );
}

void AcousticPOD::computeAllSeismoTraces( real64 const time_n,
                                          real64 const dt,
                                          arrayView1d< real32 const > const var_np1,
                                          arrayView1d< real32 const > const var_n,
                                          arrayView2d< real32 > varAtReceivers,
					  arrayView1d< real32 > coeffs,
					  bool add)
{

  /*
   * In forward case we compute seismo if time_n + dt  is the first time
   * step after the timeSeismo to write.
   *
   *  time_n        timeSeismo    time_n + dt
   *   ---|--------------|-------------|
   *
   * In backward (time_n goes decreasing) case we compute seismo if
   * time_n is the last time step before the timeSeismo to write.
   *
   *  time_n - dt    timeSeismo    time_n
   *   ---|--------------|-------------|
   */
  if( m_nsamplesSeismoTrace == 0 )
    return;
  integer const dir = m_forward ? +1 : -1;
  for( localIndex iSeismo = m_indexSeismoTrace; iSeismo < m_nsamplesSeismoTrace; iSeismo++ )
  {
    real64 const timeSeismo = m_dtSeismoTrace * (m_forward ? iSeismo : (m_nsamplesSeismoTrace - 1) - iSeismo);
    if( dir * timeSeismo > dir * (time_n + epsilonLoc))
      break;
    computeSeismoTracePOD( time_n,
                           (m_forward)?dt:-dt,
                           timeSeismo,
                           (m_forward)?m_indexSeismoTrace:(m_nsamplesSeismoTrace-m_indexSeismoTrace-1),
                           m_receiverConstantsPOD,
                           m_receiverIsLocal,
                           m_nsamplesSeismoTrace,
                           m_outputSeismoTrace,
                           var_np1,
                           var_n,
                           varAtReceivers );
  }
}

void AcousticPOD::computeMassAndDampingPOD( arrayView2d< real32 > const massPOD,
					    arrayView2d< real32 > const massGradientPOD,
					    arrayView2d< real32 > const dampingPOD,
					    arrayView1d< real32 const > const mass,
					    arrayView1d< real32 const > const massGradient,
					    arrayView1d< real32 const > const damping,
					    int const countPhi,
					    int const shotIndex,
					    arrayView1d< localIndex const > const nodesGhostRank )

{
  GEOS_LOG_RANK_0("Computing mass and damping POD...");
  int size = mass.size();
  array1d< real32 > phim( size );
  array1d< real32 > phin( size );
  arrayView1d< real32 > phimV = phim.toView();
  arrayView1d< real32 > phinV	= phin.toView();

  GEOS_MARK_SCOPE ( DirectRead );
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  for( localIndex m=0; m<countPhi; ++m )
  {
    std::string fileName1 = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", shotIndex, rank, m+1);
    std::ifstream wf1( fileName1, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf1,
		   ": Could not open file "<< fileName1 << " for reading",
		   InputError );
    phimV.move( LvArray::MemorySpace::host, true );
    wf1.read( (char *)&phimV[0], size*sizeof( real32 ) );
    wf1.close( );

    for( localIndex n=0; n<m+1; ++n )
    {
      std::string fileName2 = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", shotIndex, rank, n+1);
      std::ifstream wf2( fileName2, std::ios::in | std::ios::binary );
      GEOS_THROW_IF( !wf2,
		     ": Could not open file "<< fileName2 << " for reading",
		     InputError );
      phinV.move( LvArray::MemorySpace::host, true );
      wf2.read( (char *)&phinV[0], size*sizeof( real32 ) );
      wf2.close( );

      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	if( nodesGhostRank[a] < 0 )
        {
	  if( abs(phimV[a]) > 1e-12 && abs(phinV[a]) > 1e-12 )
	  {
	    real32 localIncrement = phimV[a] * mass[a] * phinV[a];
	    RAJA::atomicAdd< ATOMIC_POLICY >( &massPOD[m][n], localIncrement );
	    if( m!=n )
	    {
	      RAJA::atomicAdd< ATOMIC_POLICY >( &massPOD[n][m], localIncrement );
	    }

	    if( massGradient[a] != 0 )
	    {
	      localIncrement = phimV[a] * massGradient[a] * phinV[a];
              RAJA::atomicAdd< ATOMIC_POLICY >( &massGradientPOD[m][n], localIncrement );
              if( m!=n )
              {
                RAJA::atomicAdd< ATOMIC_POLICY >( &massGradientPOD[n][m], localIncrement );
              }
	    }

	    if( damping[a] != 0 )
	    {
	      localIncrement = phimV[a] * damping[a] * phinV[a];
	      RAJA::atomicAdd< ATOMIC_POLICY >( &dampingPOD[m][n], localIncrement );
	      if( m!=n )
	      {
		RAJA::atomicAdd< ATOMIC_POLICY >( &dampingPOD[n][m], localIncrement );
	      }
	    }
	  }
	}
      } ); // end loop over element
    }
  }
  massPOD.move( LvArray::MemorySpace::host, true );
  MpiWrapper::allReduce( massPOD.data(),
			 massPOD.data(),
			 massPOD.size( 0 )*massPOD.size( 1 ),
			 MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
			 MPI_COMM_GEOS );
  massGradientPOD.move( LvArray::MemorySpace::host, true );
  MpiWrapper::allReduce( massGradientPOD.data(),
                         massGradientPOD.data(),
                         massGradientPOD.size( 0 )*massGradientPOD.size( 1 ),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOS );
  dampingPOD.move( LvArray::MemorySpace::host, true );
  MpiWrapper::allReduce( dampingPOD.data(),
                         dampingPOD.data(),
                         dampingPOD.size( 0 )*dampingPOD.size( 1 ),
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                         MPI_COMM_GEOS );
}

void AcousticPOD::computeSeismoTracePOD( real64 const time_n,
					 real64 const dt,
					 real64 const timeSeismo,
					 localIndex iSeismo,
					 arrayView2d< real64 const > const receiverConstants,
					 arrayView1d< localIndex const > const receiverIsLocal,
					 localIndex const nsamplesSeismoTrace,
					 localIndex const outputSeismoTrace,
					 arrayView1d< real32 const > const var_np1,
					 arrayView1d< real32 const > const var_n,
					 arrayView2d< real32 > varAtReceivers )
{
  real64 const time_np1 = time_n + dt;

  real32 const a1 = (LvArray::math::abs( dt ) < WaveSolverBase::epsilonLoc ) ? 1.0 : (time_np1 - timeSeismo)/dt;
  real32 const a2 = 1.0 - a1;
  if( nsamplesSeismoTrace > 0 )
  {
    forAll< WaveSolverBase::EXEC_POLICY >( receiverConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
	varAtReceivers[iSeismo][ircv] = 0.0;
	real32 vtmp_np1 = 0.0;
	real32 vtmp_n = 0.0;

	for( localIndex m = 0; m<var_np1.size( 0 ); ++m )
	{
	  vtmp_np1 += var_np1[m] * receiverConstants[ircv][m];
	  vtmp_n += var_n[m] * receiverConstants[ircv][m];
	}
	// linear interpolation between the pressure value at time_n and time_(n+1)
	varAtReceivers[iSeismo][ircv] = a1*vtmp_n + a2*vtmp_np1;
      }
    } );
  }
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticPOD, string const &, dataRepository::Group * const )

} /* namespace geos */
