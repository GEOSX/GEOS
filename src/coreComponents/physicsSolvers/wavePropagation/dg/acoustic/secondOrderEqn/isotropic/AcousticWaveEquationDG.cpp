/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file AcousticWaveEquationDG.cpp
 */

#include "AcousticWaveEquationDG.hpp"
//#include "AcousticWaveEquationDGKernel.hpp"


#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "events/EventManager.hpp"
#include "physicsSolvers/wavePropagation/shared/PrecomputeSourcesAndReceiversKernel.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

AcousticWaveEquationDG::AcousticWaveEquationDG( const std::string & name,
                                                Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::sourceElemString(), &m_sourceElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the sources" );

  registerWrapper( viewKeyStruct::receiverElemString(), &m_rcvElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the receivers" );

  registerWrapper( viewKeyStruct::receiverRegionString(), &m_receiverRegion ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Region containing the receivers" );

}

AcousticWaveEquationDG::~AcousticWaveEquationDG()
{
  // TODO Auto-generated destructor stub
}

void AcousticWaveEquationDG::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void AcousticWaveEquationDG::registerDataOnMesh( Group & meshBodies )
{

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< acousticfieldsdg::AcousticFreeSurfaceNodeIndicator >( this->getName() );

    nodeManager.registerField<acousticfieldsdg::AcousticMassMatrix>( getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< acousticfieldsdg::AcousticFreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< acousticfieldsdg::AcousticVelocity >( this->getName() );

      subRegion.registerField< acousticfieldsdg::Pressure_nm1 >( this->getName() );
      subRegion.registerField< acousticfieldsdg::Pressure_n >( this->getName() );
      subRegion.registerField< acousticfieldsdg::Pressure_np1 >( this->getName() );
      subRegion.registerField< acousticfieldsdg::StiffnessVector >( this->getName() );
      subRegion.registerField< acousticfieldsdg::AcousticMassMatrix >( this->getName() );

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< DG_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {

        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

        subRegion.getField< acousticfieldsdg::Pressure_nm1 >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< acousticfieldsdg::Pressure_n >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< acousticfieldsdg::Pressure_np1 >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< acousticfieldsdg::StiffnessVector >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< acousticfieldsdg::AcousticMassMatrix >().resizeDimension< 1 >( numNodesPerElem );

      } );

    } );
  } );
}


void AcousticWaveEquationDG::precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh,
                                                              arrayView1d< string const > const & regionNames )
{
  arrayView1d< globalIndex const > const nodeLocalToGlobalMap = baseMesh.getNodeManager().localToGlobalMap().toViewConst();
  ArrayOfArraysView< localIndex const > const nodesToElements = baseMesh.getNodeManager().elementList().toViewConst();
  ArrayOfArraysView< localIndex const > const facesToNodes = baseMesh.getFaceManager().nodeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodeCoords = baseMesh.getNodeManager().referencePosition();


  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex > const sourceElem = m_sourceElem.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  arrayView1d< localIndex > const rcvElem = m_rcvElem.toView();
  receiverNodeIds.setValues< EXEC_POLICY >( -1 );
  receiverConstants.setValues< EXEC_POLICY >( -1 );
  receiverIsLocal.zero();

  mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                 localIndex const er,
                                                                                                 localIndex const esr,
                                                                                                 ElementRegionBase &,
                                                                                                 CellElementSubRegion & elementSubRegion )
  {
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Tetrahedron,
                   "Invalid type of element, the acoustic DG solver is designed for tetrahedral meshes only  ",
                   InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes = baseMesh.getElemManager().getRegion( er ).getSubRegion< CellElementSubRegion >( esr ).nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();
    arrayView1d< globalIndex const > const elemLocalToGlobalMap = elementSubRegion.localToGlobalMap().toViewConst();

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< DG_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      // AcousticWaveEquationDGKernels::
      //   PrecomputeSourceAndReceiverKernel::
      //   launch< EXEC_POLICY, FE_TYPE >
      //   ( elementSubRegion.size(),
      //   numNodesPerElem,
      //   numFacesPerElem,
      //   X,
      //   elemGhostRank,
      //   elemsToNodes,
      //   elemsToFaces,
      //   elemCenter,
      //   faceNormal,
      //   faceCenter,
      //   sourceCoordinates,
      //   sourceIsAccessible,
      //   sourceElem,
      //   sourceNodeIds,
      //   sourceConstants,
      //   receiverCoordinates,
      //   receiverIsLocal,
      //   rcvElem,
      //   receiverNodeIds,
      //   receiverConstants,
      //   sourceValue,
      //   dt,
      //   timeSourceFrequency,
      //   rickerOrder );
    } );
  } );

}

void AcousticWaveEquationDG::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;
  {
    GEOS_MARK_SCOPE( WaveSolverBase::initializePostInitialConditionsPreSubGroups );
    WaveSolverBase::initializePostInitialConditionsPreSubGroups();
  }

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  applyFreeSurfaceBC( 0.0, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    MeshLevel & baseMesh = domain.getMeshBodies().getGroup< MeshBody >( meshBodyName ).getBaseDiscretization();
    precomputeSourceAndReceiverTerm( baseMesh, mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function

    /// damping matrix to be computed for each dof in the boundary of the mesh

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< acousticfieldsdg::AcousticFreeSurfaceFaceIndicator >();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      computeTargetNodeSet( elemsToNodes, elementSubRegion.size(), fe.getNumQuadraturePoints() );

      //  arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfieldsdgdgdg::AcousticVelocity >();
      //  arrayView1d< real32 const > const density = elementSubRegion.getField< acousticfieldsdgdgdg::AcousticDensity >();

      /// Partial gradient if gradient as to be computed

      finiteElement::FiniteElementDispatchHandler< DG_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        // AcousticMatricesSEM::MassMatrix< FE_TYPE > kernelM( finiteElement );
        // kernelM.template computeMassMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
        //                                                                   nodeCoords,
        //                                                                   elemsToNodes,
        //                                                                   velocity,
        //                                                                   density,
        //                                                                   mass );

        // AcousticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );
        // kernelD.template computeDampingMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
        //                                                                      nodeCoords,
        //                                                                      elemsToFaces,
        //                                                                      facesToNodes,
        //                                                                      facesDomainBoundaryIndicator,
        //                                                                      freeSurfaceFaceIndicator,
        //                                                                      velocity,
        //                                                                      density,
        //                                                                      damping );


      } );
    } );
  } );

  // if( m_timestepStabilityLimit==1 )
  // {
  //   real64 dtOut = 0.0;
  //   computeTimeStep( dtOut );
  //   m_timestepStabilityLimit = 0;
  //   m_timeStep=dtOut;
  // }


  WaveSolverUtils::initTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ), m_receiverIsLocal );
}


real64 AcousticWaveEquationDG::computeTimeStep( real64 & dtOut )
{
  GEOS_ERROR( getDataContext() << ":  Time-Step computation for  acoustic dg wave propagator not yet implemented" );
  return dtOut;
}

//TODO: Modify to use on discontinuous variable
void AcousticWaveEquationDG::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView2d< real32 > const p_np1 = nodeManager.getField< acousticfieldsdg::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< acousticfieldsdg::AcousticFreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< acousticfieldsdg::AcousticFreeSurfaceNodeIndicator >();


  freeSurfaceFaceIndicator.zero();
  freeSurfaceNodeIndicator.zero();

  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "FreeSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    // string const & functionName = bc.getFunctionName();

    // if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    // {
    //   real64 const value = bc.getScale();

    //   for( localIndex i = 0; i < targetSet.size(); ++i )
    //   {
    //     localIndex const kf = targetSet[ i ];
    //     freeSurfaceFaceIndicator[kf] = 1;

    //     localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
    //     for( localIndex a=0; a < numNodes; ++a )
    //     {
    //       localIndex const dof = faceToNodeMap( kf, a );
    //       freeSurfaceNodeIndicator[dof] = 1;

    //       p_np1[dof] = value;
    //     }
    //   }
    // }
    // else
    // {
    //   GEOS_ERROR( "This option is not supported yet" );
    // }
  } );
}

// Here for retrocompatibily
real64 AcousticWaveEquationDG::explicitStepForward( real64 const & time_n,
                                                    real64 const & dt,
                                                    integer cycleNumber,
                                                    DomainPartition & domain,
                                                    bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}



real64 AcousticWaveEquationDG::explicitStepBackward( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer cycleNumber,
                                                     DomainPartition & domain,
                                                     bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  GEOS_ERROR( "Backward propagation for the first-order wave propagator not yet implemented" );
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}


void AcousticWaveEquationDG::prepareNextTimestep( MeshLevel & mesh )           
{                                                                               
  NodeManager & nodeManager = mesh.getNodeManager();                            
                                                                                
  arrayView2d< real32 > const p_nm1 = nodeManager.getField< acousticfieldsdg::Pressure_nm1 >();
  arrayView2d< real32 > const p_n   = nodeManager.getField< acousticfieldsdg::Pressure_n >();
  arrayView2d< real32 > const p_np1 = nodeManager.getField< acousticfieldsdg::Pressure_np1 >();
                                                                                
  arrayView2d< real32 > const stiffnessVector = nodeManager.getField< acousticfieldsdg::StiffnessVector >();
  //arrayView2d< real32 > const rhs = nodeManager.getField< acousticfieldsdgdgdg::ForcingRHS >();
                                                                                
  SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();
                                                                                
  // forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
  // {                                                                             
  //   localIndex const a = solverTargetNodesSet[n];                               
                                                                                
  //   p_nm1[a] = p_n[a];                                                          
  //   p_n[a]   = p_np1[a];                                                        
                                                                                
  //   stiffnessVector[a] = rhs[a] = 0.0;                                          
  // } );  

}  

void AcousticWaveEquationDG::computeUnknowns( real64 const & time_n,
                                               real64 const & dt,
                                               DomainPartition & domain,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  // arrayView1d< real32 const > const mass = nodeManager.getField< acousticfieldsdgdg::AcousticAcousticMassMatrix >();
  // arrayView1d< real32 const > const damping = nodeManager.getField< acousticfieldsdgdg::DampingVector >();

  arrayView2d< real32 > const p_nm1 = nodeManager.getField< acousticfieldsdg::Pressure_nm1 >();
  arrayView2d< real32 > const p_n = nodeManager.getField< acousticfieldsdg::Pressure_n >();
  arrayView2d< real32 > const p_np1 = nodeManager.getField< acousticfieldsdg::Pressure_np1 >();

  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< acousticfieldsdg::AcousticFreeSurfaceNodeIndicator >();
  arrayView2d< real32 > const stiffnessVector = nodeManager.getField< acousticfieldsdg::StiffnessVector >();
  // arrayView1d< real32 > const rhs = nodeManager.getField< acousticfieldsdgdg::ForcingRHS >();

  // auto kernelFactory = acousticWaveEquationSEMKernels::ExplicitAcousticSEMFactory( dt );

  // finiteElement::
  //   regionBasedKernelApplication< EXEC_POLICY,
  //                                 constitutive::NullModel,
  //                                 CellElementSubRegion >( mesh,
  //                                                         regionNames,
  //                                                         getDiscretizationName(),
  //                                                         "",
  //                                                         kernelFactory );
  //Modification of cycleNember useful when minTime < 0
  // addSourceToRightHandSide( time_n, rhs );

  // /// calculate your time integrators
  // real64 const dt2 = pow( dt, 2 );

  // SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();
  //   GEOS_MARK_SCOPE ( updateP );
  //   AcousticTimeSchemeDG::LeapFrogWithoutPML( dt, p_np1, p_n, p_nm1, mass, stiffnessVector, damping,
  //                                              rhs, freeSurfaceNodeIndicator, solverTargetNodesSet );
}

void AcousticWaveEquationDG::synchronizeUnknowns( real64 const & time_n,
                                                   real64 const & dt,
                                                   DomainPartition & domain,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView2d< real32 > const p_n = nodeManager.getField< acousticfieldsdg::Pressure_n >();
  arrayView2d< real32 > const p_np1 = nodeManager.getField< acousticfieldsdg::Pressure_np1 >();

  arrayView2d< real32 > const stiffnessVector = nodeManager.getField< acousticfieldsdg::StiffnessVector >();
  //arrayView1d< real32 > const rhs = nodeManager.getField< acousticfieldsdg::ForcingRHS >();

  /// synchronize pressure fields
  FieldIdentifiers fieldsToBeSync;
  fieldsToBeSync.addFields( FieldLocation::Node, { acousticfieldsdg::Pressure_np1::key() } );

  CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldsToBeSync,
                                mesh,
                                domain.getNeighbors(),
                                true );
  /// compute the seismic traces since last step.
  arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();

  //computeAllSeismoTraces( time_n, dt, p_np1, p_n, pReceivers );
  //incrementIndexSeismoTrace( time_n );

}



real64 AcousticWaveEquationDG::explicitStepInternal( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer const cycleNumber,
                                                     DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  real64 dtCompute;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    computeUnknowns( time_n, dt, domain, mesh, regionNames );
    synchronizeUnknowns( time_n, dt, domain, mesh, regionNames );
  } );

  return dt;


}

void AcousticWaveEquationDG::cleanup( real64 const time_n, integer const, integer const, real64 const, DomainPartition & domain )
{
  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< real32 const > const p_np1 = elementSubRegion.getField< acousticfieldsdg::Pressure_np1 >();
      arrayView2d< real32 const > const p_n = elementSubRegion.getField< acousticfieldsdg::Pressure_n >();

      arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, p_np1, p_n, pReceivers );

    } );
  } );

  // increment m_indexSeismoTrace
  while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
  {
    m_indexSeismoTrace++;
  }
}

void AcousticWaveEquationDG::initializePML()
{
  GEOS_ERROR( "PML for the first order acoustic wave propagator not yet implemented" );
}

void AcousticWaveEquationDG::applyPML( real64 const, DomainPartition & )
{
  GEOS_ERROR( "PML for the first order acoustic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticWaveEquationDG, string const &, dataRepository::Group * const )

} /* namespace geos */
