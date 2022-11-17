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
 * @file SolidMechanicsLagrangianFEM.cpp
 */

#define GEOSX_DISPATCH_VEM

#include "SolidMechanicsLagrangianFEM.hpp"
#include "kernels/ImplicitSmallStrainNewmark.hpp"
#include "kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "kernels/ExplicitSmallStrain.hpp"
#include "kernels/ExplicitFiniteStrain.hpp"

#include "codingUtilities/Utilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/TractionBoundaryCondition.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "LvArray/src/output.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;

SolidMechanicsLagrangianFEM::SolidMechanicsLagrangianFEM( const string & name,
                                                          Group * const parent ):
  SolverBase( name, parent ),
  m_newmarkGamma( 0.5 ),
  m_newmarkBeta( 0.25 ),
  m_massDamping( 0.0 ),
  m_stiffnessDamping( 0.0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_maxForce( 0.0 ),
  m_maxNumResolves( 10 ),
  m_strainTheory( 0 ),
  m_iComm( CommunicationTools::getInstance().getCommID() )
{

  registerWrapper( viewKeyStruct::newmarkGammaString(), &m_newmarkGamma ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of :math:`\\gamma` in the Newmark Method for Implicit Dynamic time integration option" );

  registerWrapper( viewKeyStruct::newmarkBetaString(), &m_newmarkBeta ).
    setApplyDefaultValue( 0.25 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of :math:`\\beta` in the Newmark Method for Implicit Dynamic time integration option. "
                    "This should be pow(newmarkGamma+0.5,2.0)/4.0 unless you know what you are doing." );

  registerWrapper( viewKeyStruct::massDampingString(), &m_massDamping ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of mass based damping coefficient. " );

  registerWrapper( viewKeyStruct::stiffnessDampingString(), &m_stiffnessDamping ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of stiffness based damping coefficient. " );

  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed after some other event is executed. "
                    "For example, if a SurfaceGenerator is specified, it will be executed after the mechanics solve. "
                    "However if a new surface is generated, then the mechanics solve must be executed again due to the "
                    "change in topology." );

  registerWrapper( viewKeyStruct::strainTheoryString(), &m_strainTheory ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Indicates whether or not to use "
                    "`Infinitesimal Strain Theory <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_, or "
                    "`Finite Strain Theory <https://en.wikipedia.org/wiki/Finite_strain_theory>`_. Valid Inputs are:\n"
                    " 0 - Infinitesimal Strain \n"
                    " 1 - Finite Strain" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setApplyDefaultValue( viewKeyStruct::noContactRelationNameString() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxForceString(), &m_maxForce ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "The maximum force contribution in the problem domain." );

}

void SolidMechanicsLagrangianFEM::postProcessInput()
{
  SolverBase::postProcessInput();

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.isSymmetric = true;
  linParams.dofsPerNode = 3;
  linParams.amg.separateComponents = true;

}

SolidMechanicsLagrangianFEM::~SolidMechanicsLagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanicsLagrangianFEM::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodes = meshLevel.getNodeManager();

    nodes.registerField< solidMechanics::totalDisplacement >( getName() ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerField< solidMechanics::incrementalDisplacement >( getName() ).
      reference().resizeDimension< 1 >( 3 );

    if( m_timeIntegrationOption != TimeIntegrationOption::QuasiStatic )
    {
      nodes.registerField< solidMechanics::velocity >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      nodes.registerField< solidMechanics::acceleration >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      nodes.registerField< solidMechanics::velocityTilde >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      nodes.registerField< solidMechanics::uhatTilde >( getName() ).
        reference().resizeDimension< 1 >( 3 );
    }

    nodes.registerField< solidMechanics::mass >( getName() );

    nodes.registerField< solidMechanics::externalForce >( getName() ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerField< solidMechanics::contactForce >( getName() ).
      reference().resizeDimension< 1 >( 3 );

    Group & nodeSets = nodes.sets();
    nodeSets.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::sendOrReceiveNodesString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE );

    nodeSets.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::nonSendOrReceiveNodesString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE );

    nodeSets.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::targetNodesString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE );

    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.excludeWrappersFromPacking( { viewKeyStruct::elemsAttachedToSendOrReceiveNodesString(),
                                              viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() } );
    } );
  } );
}

void SolidMechanicsLagrangianFEM::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  SolverBase::setConstitutiveNamesCallSuper( subRegion );

  subRegion.registerWrapper< string >( viewKeyStruct::solidMaterialNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );

  string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
  solidMaterialName = SolverBase::getConstitutiveName< SolidBase >( subRegion );
  GEOSX_ERROR_IF( solidMaterialName.empty(), GEOSX_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );
}

void SolidMechanicsLagrangianFEM::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  GEOSX_UNUSED_VAR( subRegion );
}


void SolidMechanicsLagrangianFEM::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
      solidMaterialName = SolverBase::getConstitutiveName< SolidBase >( subRegion );
    } );
  } );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const &
  feDiscretization = feDiscretizationManager.getGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_UNUSED_VAR( feDiscretization );
}



template< typename ... PARAMS >
real64 SolidMechanicsLagrangianFEM::explicitKernelDispatch( MeshLevel & mesh,
                                                            arrayView1d< string const > const & targetRegions,
                                                            string const & finiteElementName,
                                                            real64 const dt,
                                                            std::string const & elementListName )
{
  GEOSX_MARK_FUNCTION;
  real64 rval = 0;
  if( m_strainTheory==0 )
  {
    auto kernelFactory = solidMechanicsLagrangianFEMKernels::ExplicitSmallStrainFactory( dt, elementListName );
    rval = finiteElement::
             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                           constitutive::SolidBase,
                                           CellElementSubRegion >( mesh,
                                                                   targetRegions,
                                                                   finiteElementName,
                                                                   viewKeyStruct::solidMaterialNamesString(),
                                                                   kernelFactory );
  }
  else if( m_strainTheory==1 )
  {
    auto kernelFactory = solidMechanicsLagrangianFEMKernels::ExplicitFiniteStrainFactory( dt, elementListName );
    rval = finiteElement::
             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                           constitutive::SolidBase,
                                           CellElementSubRegion >( mesh,
                                                                   targetRegions,
                                                                   finiteElementName,
                                                                   viewKeyStruct::solidMaterialNamesString(),
                                                                   kernelFactory );
  }
  else
  {
    GEOSX_ERROR( "Invalid option for strain theory (0 = infinitesimal strain, 1 = finite strain" );
  }

  return rval;
}


void SolidMechanicsLagrangianFEM::initializePostInitialConditionsPreSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&]( string const &,
                                       MeshLevel & mesh,
                                       auto const & regionNames )
  {
    NodeManager & nodes = mesh.getNodeManager();
    Group & nodeSets = nodes.sets();

    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    EdgeManager const & edgeManager = mesh.getEdgeManager();
    arrayView1d< real64 > & mass = nodes.getField< solidMechanics::mass >();

    arrayView1d< integer const > const & nodeGhostRank = nodes.ghostRank();

    // to fill m_sendOrReceiveNodes and m_nonSendOrReceiveNodes, we first insert
    // the nodes one-by-one in the following std::sets. Then, when all the nodes
    // have been collected, we do a batch insertion into m_sendOrReceiveNodes and
    // m_nonSendOrReceiveNodes
    std::set< localIndex > tmpSendOrReceiveNodes;
    std::set< localIndex > tmpNonSendOrReceiveNodes;

    SortedArray< localIndex > & m_sendOrReceiveNodes = nodeSets.getReference< SortedArray< localIndex > >( viewKeyStruct::sendOrReceiveNodesString() );
    SortedArray< localIndex > & m_nonSendOrReceiveNodes = nodeSets.getReference< SortedArray< localIndex > >( viewKeyStruct::nonSendOrReceiveNodesString() );
    SortedArray< localIndex > & m_targetNodes = nodeSets.getReference< SortedArray< localIndex > >( viewKeyStruct::targetNodesString() );

    elementRegionManager.forElementRegionsComplete( regionNames,
                                                    [&]( localIndex const,
                                                         localIndex const er,
                                                         ElementRegionBase & elemRegion )
    {
      elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr, CellElementSubRegion & elementSubRegion )
      {
        string const & solidMaterialName = elementSubRegion.getReference< string >( viewKeyStruct::solidMaterialNamesString() );

        arrayView2d< real64 const > const
        rho = elementSubRegion.getConstitutiveModel( solidMaterialName ).getReference< array2d< real64 > >( SolidBase::viewKeyStruct::densityString() );

        SortedArray< localIndex > & elemsAttachedToSendOrReceiveNodes = getElemsAttachedToSendOrReceiveNodes( elementSubRegion );
        SortedArray< localIndex > & elemsNotAttachedToSendOrReceiveNodes = getElemsNotAttachedToSendOrReceiveNodes( elementSubRegion );

        std::set< localIndex > tmpElemsAttachedToSendOrReceiveNodes;
        std::set< localIndex > tmpElemsNotAttachedToSendOrReceiveNodes;

        elemsAttachedToSendOrReceiveNodes.setName(
          "SolidMechanicsLagrangianFEM::m_elemsAttachedToSendOrReceiveNodes["
          + std::to_string( er ) + "][" + std::to_string( esr ) + "]" );

        elemsNotAttachedToSendOrReceiveNodes.setName(
          "SolidMechanicsLagrangianFEM::m_elemsNotAttachedToSendOrReceiveNodes["
          + std::to_string( er ) + "][" + std::to_string( esr ) + "]" );

        arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
        finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe,
                                                                                 [&] ( auto const finiteElement )
        {
          using FE_TYPE = TYPEOFREF( finiteElement );
          using SUBREGION_TYPE = TYPEOFREF( elementSubRegion );

          typename FE_TYPE::template MeshData< SUBREGION_TYPE > meshData;
          finiteElement::FiniteElementBase::initialize< FE_TYPE, SUBREGION_TYPE >( nodes,
                                                                                   edgeManager,
                                                                                   faceManager,
                                                                                   elementSubRegion,
                                                                                   meshData );

          constexpr localIndex maxSupportPoints = FE_TYPE::maxSupportPoints;
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

          real64 N[maxSupportPoints];
          for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
          {
            typename FE_TYPE::StackVariables feStack;
            finiteElement.template setup< FE_TYPE >( k, meshData, feStack );
            localIndex const numSupportPoints =
              finiteElement.template numSupportPoints< FE_TYPE >( feStack );

            for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
            {
              FE_TYPE::calcN( q, feStack, N );

              for( localIndex a=0; a< numSupportPoints; ++a )
              {
                mass[elemsToNodes[k][a]] += rho[k][q] * detJ[k][q] * N[a];
              }
            }

            bool isAttachedToGhostNode = false;
            for( localIndex a=0; a<elementSubRegion.numNodesPerElement(); ++a )
            {
              if( nodeGhostRank[elemsToNodes[k][a]] >= -1 )
              {
                isAttachedToGhostNode = true;
                tmpSendOrReceiveNodes.insert( elemsToNodes[k][a] );
              }
              else
              {
                tmpNonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
              }
            }

            if( isAttachedToGhostNode )
            {
              tmpElemsAttachedToSendOrReceiveNodes.insert( k );
            }
            else
            {
              tmpElemsNotAttachedToSendOrReceiveNodes.insert( k );
            }
          }
        } );
        elemsAttachedToSendOrReceiveNodes.insert( tmpElemsAttachedToSendOrReceiveNodes.begin(),
                                                  tmpElemsAttachedToSendOrReceiveNodes.end() );
        elemsNotAttachedToSendOrReceiveNodes.insert( tmpElemsNotAttachedToSendOrReceiveNodes.begin(),
                                                     tmpElemsNotAttachedToSendOrReceiveNodes.end() );

        m_sendOrReceiveNodes.insert( tmpSendOrReceiveNodes.begin(),
                                     tmpSendOrReceiveNodes.end() );
        m_nonSendOrReceiveNodes.insert( tmpNonSendOrReceiveNodes.begin(),
                                        tmpNonSendOrReceiveNodes.end() );
        m_targetNodes = m_sendOrReceiveNodes;
        m_targetNodes.insert( m_nonSendOrReceiveNodes.begin(),
                              m_nonSendOrReceiveNodes.end() );

      } );
    } );

  } );
}



real64 SolidMechanicsLagrangianFEM::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                const int cycleNumber,
                                                DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator = this->getParent().getGroupPointer< SolverBase >( "SurfaceGen" );

  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );

    if( surfaceGenerator!=nullptr )
    {
      surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain );
    }
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic ||
           m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
  {
    int const maxNumResolves = m_maxNumResolves;
    int locallyFractured = 0;
    int globallyFractured = 0;
    implicitStepSetup( time_n, dt, domain );
    for( int solveIter=0; solveIter<maxNumResolves; ++solveIter )
    {
      setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );

      dtReturn = nonlinearImplicitStep( time_n,
                                        dt,
                                        cycleNumber,
                                        domain );

      if( surfaceGenerator!=nullptr )
      {
        if( surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain ) > 0 )
        {
          locallyFractured = 1;
        }
        MpiWrapper::allReduce( &locallyFractured,
                               &globallyFractured,
                               1,
                               MPI_MAX,
                               MPI_COMM_GEOSX );
      }
      if( globallyFractured == 0 )
      {
        break;
      }
      else
      {
        GEOSX_LOG_RANK_0( "Fracture Occurred. Resolve" );
      }
    }
    implicitStepComplete( time_n, dt, domain );
  }

  return dtReturn;
}

real64 SolidMechanicsLagrangianFEM::explicitStep( real64 const & time_n,
                                                  real64 const & dt,
                                                  const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodes = mesh.getNodeManager();
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    Group const & nodeSets = nodes.sets();

    SortedArrayView< localIndex const > const &
    m_sendOrReceiveNodes = nodeSets.getReference< SortedArray< localIndex > >( viewKeyStruct::sendOrReceiveNodesString() ).toViewConst();

    SortedArrayView< localIndex const > const &
    m_nonSendOrReceiveNodes = nodeSets.getReference< SortedArray< localIndex > >( viewKeyStruct::nonSendOrReceiveNodesString() ).toViewConst();

    // save previous constitutive state data in preparation for next timestep
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
      SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
      constitutiveRelation.saveConvergedState();
    } );

    FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

    arrayView1d< real64 const > const & mass = nodes.getField< solidMechanics::mass >();
    solidMechanics::arrayView2dLayoutVelocity const & vel = nodes.getField< solidMechanics::velocity >();
    solidMechanics::arrayView2dLayoutTotalDisplacement const & u = nodes.getField< solidMechanics::totalDisplacement >();
    solidMechanics::arrayView2dLayoutIncrDisplacement const & uhat = nodes.getField< solidMechanics::incrementalDisplacement >();
    solidMechanics::arrayView2dLayoutAcceleration const & acc = nodes.getField< solidMechanics::acceleration >();

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node,
                              { solidMechanics::velocity::key(),
                                solidMechanics::acceleration::key() } );
    m_iComm.resize( domain.getNeighbors().size() );
    CommunicationTools::getInstance().synchronizePackSendRecvSizes( fieldsToBeSync, mesh, domain.getNeighbors(), m_iComm, true );

    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, mesh, solidMechanics::acceleration::key() );

    //3: v^{n+1/2} = v^{n} + a^{n} dt/2
    solidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt / 2 );

    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, mesh, solidMechanics::velocity::key() );

    //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
    solidMechanicsLagrangianFEMKernels::displacementUpdate( vel, uhat, u, dt );

    fsManager.applyFieldValue( time_n + dt,
                               mesh,
                               solidMechanics::totalDisplacement::key(),
                               [&]( FieldSpecificationBase const & bc,
                                    SortedArrayView< localIndex const > const & targetSet )
    {
      integer const component = bc.getComponent();
      GEOSX_ERROR_IF_LT_MSG( component, 0, "Component index required for displacement BC " << bc.getName() );

      forAll< parallelDevicePolicy< 1024 > >( targetSet.size(),
                                              [=] GEOSX_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        vel( a, component ) = u( a, component );
      } );
    },
                               [&]( FieldSpecificationBase const & bc,
                                    SortedArrayView< localIndex const > const & targetSet )
    {
      integer const component = bc.getComponent();
      GEOSX_ERROR_IF_LT_MSG( component, 0, "Component index required for displacement BC " << bc.getName() );

      forAll< parallelDevicePolicy< 1024 > >( targetSet.size(),
                                              [=] GEOSX_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        uhat( a, component ) = u( a, component ) - vel( a, component );
        vel( a, component )  = uhat( a, component ) / dt;
      } );
    } );

    //Step 5. Calculate deformation input to constitutive model and update state to
    // Q^{n+1}
    explicitKernelDispatch( mesh,
                            regionNames,
                            this->getDiscretizationName(),
                            dt,
                            string( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() ) );

    // apply this over a set
    solidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_sendOrReceiveNodes.toViewConst() );

    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, mesh, solidMechanics::velocity::key() );

    parallelDeviceEvents packEvents;
    CommunicationTools::getInstance().asyncPack( fieldsToBeSync, mesh, domain.getNeighbors(), m_iComm, true, packEvents );

    waitAllDeviceEvents( packEvents );

    CommunicationTools::getInstance().asyncSendRecv( domain.getNeighbors(), m_iComm, true, packEvents );

    explicitKernelDispatch( mesh,
                            regionNames,
                            this->getDiscretizationName(),
                            dt,
                            string( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() ) );

    // apply this over a set
    solidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_nonSendOrReceiveNodes.toViewConst() );
    fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, mesh, solidMechanics::velocity::key() );

    // this includes  a device sync after launching all the unpacking kernels
    parallelDeviceEvents unpackEvents;
    CommunicationTools::getInstance().finalizeUnpack( mesh, domain.getNeighbors(), m_iComm, true, unpackEvents );

  } );

  return dt;
}



void SolidMechanicsLagrangianFEM::applyDisplacementBCImplicit( real64 const time,
                                                               DofManager const & dofManager,
                                                               DomainPartition & domain,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  string const dofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    fsManager.apply< NodeManager >( time,
                                    mesh,
                                    solidMechanics::totalDisplacement::key(),
                                    [&]( FieldSpecificationBase const & bc,
                                         string const &,
                                         SortedArrayView< localIndex const > const & targetSet,
                                         NodeManager & targetGroup,
                                         string const fieldName )
    {
      bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                         parallelDevicePolicy< 32 > >( targetSet,
                                                                       time,
                                                                       targetGroup,
                                                                       fieldName,
                                                                       dofKey,
                                                                       dofManager.rankOffset(),
                                                                       localMatrix,
                                                                       localRhs );
    } );
  } );
}

void SolidMechanicsLagrangianFEM::applyTractionBC( real64 const time,
                                                   DofManager const & dofManager,
                                                   DomainPartition & domain,
                                                   arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const dofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const blockLocalDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
    globalIndex const dofRankOffset = dofManager.rankOffset();

    fsManager.template apply< FaceManager,
                              TractionBoundaryCondition >( time,
                                                           mesh,
                                                           TractionBoundaryCondition::catalogName(),
                                                           [&]( TractionBoundaryCondition const & bc,
                                                                string const &,
                                                                SortedArrayView< localIndex const > const & targetSet,
                                                                Group &,
                                                                string const & )
    {
      bc.launch( time,
                 blockLocalDofNumber,
                 dofRankOffset,
                 faceManager,
                 targetSet,
                 localRhs );
    } );
  } );
}

void SolidMechanicsLagrangianFEM::applyChomboPressure( DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    FaceManager & faceManager = mesh.getFaceManager();
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real64 const > const faceArea  = faceManager.faceArea();
    arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    string const dofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
    arrayView1d< real64 const > const facePressure = faceManager.getReference< array1d< real64 > >( "ChomboPressure" );

    forAll< serialPolicy >( faceManager.size(), [=] ( localIndex const kf )
    {
      int const numNodes = LvArray::integerConversion< int >( faceToNodeMap.sizeOfArray( kf ));
      for( int a=0; a<numNodes; ++a )
      {
        localIndex const dof = dofNumber[ faceToNodeMap( kf, a ) ];
        if( dof < 0 || dof >= localRhs.size() )
          continue;

        for( int component=0; component<3; ++component )
        {
          real64 const value = -facePressure[ kf ] * faceNormal( kf, component ) * faceArea[kf] / numNodes;
          localRhs[ dof + component ] += value;
        }
      }
    } );
  } );
}



void
SolidMechanicsLagrangianFEM::
  implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                     real64 const & dt,
                     DomainPartition & domain )
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    solidMechanics::arrayView2dLayoutIncrDisplacement const uhat =
      nodeManager.getField< solidMechanics::incrementalDisplacement >();
    solidMechanics::arrayView2dLayoutTotalDisplacement const disp =
      nodeManager.getField< solidMechanics::totalDisplacement >();


    localIndex const numNodes = nodeManager.size();

    if( this->m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
    {
      solidMechanics::arrayViewConst2dLayoutAcceleration const a_n = nodeManager.getField< solidMechanics::acceleration >();
      solidMechanics::arrayView2dLayoutVelocity const v_n = nodeManager.getField< solidMechanics::velocity >();
      arrayView2d< real64 > const vtilde = nodeManager.getField< solidMechanics::velocityTilde >();
      arrayView2d< real64 > const uhatTilde = nodeManager.getField< solidMechanics::uhatTilde >();

      real64 const newmarkGamma = this->getReference< real64 >( solidMechanicsViewKeys.newmarkGamma );
      real64 const newmarkBeta = this->getReference< real64 >( solidMechanicsViewKeys.newmarkBeta );

      forAll< parallelDevicePolicy< 32 > >( numNodes, [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        for( int i=0; i<3; ++i )
        {
          vtilde[a][i] = v_n( a, i ) + (1.0-newmarkGamma) * a_n( a, i ) * dt;
          uhatTilde[a][i] = ( v_n( a, i ) + 0.5 * ( 1.0 - 2.0*newmarkBeta ) * a_n( a, i ) * dt ) *dt;
          uhat( a, i ) = uhatTilde[a][i];
          disp( a, i ) += uhatTilde[a][i];
        }
      } );
    }
    else if( this->m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
    {
      forAll< parallelDevicePolicy< 32 > >( numNodes, [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        for( int i=0; i<3; ++i )
        {
          uhat( a, i ) = 0.0;
        }
      } );
    }

    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();
    ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase >
    constitutiveRelations = elementRegionManager.constructFullConstitutiveAccessor< ConstitutiveBase >( constitutiveManager );

    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
      SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
      constitutiveRelation.saveConvergedState();
    } );
  } );

}

void SolidMechanicsLagrangianFEM::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    localIndex const numNodes = nodeManager.size();
    ElementRegionManager & elementRegionManager = mesh.getElemManager();

    solidMechanics::arrayView2dLayoutIncrDisplacement const uhat =
      nodeManager.getField< solidMechanics::incrementalDisplacement >();

    if( this->m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
    {
      solidMechanics::arrayView2dLayoutAcceleration const a_n = nodeManager.getField< solidMechanics::acceleration >();
      solidMechanics::arrayView2dLayoutVelocity const v_n = nodeManager.getField< solidMechanics::velocity >();
      arrayView2d< real64 const > const vtilde    = nodeManager.getField< solidMechanics::velocityTilde >();
      arrayView2d< real64 const > const uhatTilde = nodeManager.getField< solidMechanics::uhatTilde >();

      real64 const newmarkGamma = this->getReference< real64 >( solidMechanicsViewKeys.newmarkGamma );
      real64 const newmarkBeta = this->getReference< real64 >( solidMechanicsViewKeys.newmarkBeta );

      RAJA::forall< parallelDevicePolicy<> >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                              [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        for( int i=0; i<3; ++i )
        {
          a_n( a, i ) = 1.0 / ( newmarkBeta * dt*dt) * ( uhat( a, i ) - uhatTilde[a][i] );
          v_n[a][i] = vtilde[a][i] + newmarkGamma * a_n( a, i ) * dt;
        }
      } );
    }

    // save (converged) constitutive state data
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            CellElementSubRegion & subRegion )
    {
      string const & solidMaterialName = subRegion.template getReference< string >( viewKeyStruct::solidMaterialNamesString() );
      SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, solidMaterialName );
      constitutiveRelation.saveConvergedState();
    } );
  } );

}

void SolidMechanicsLagrangianFEM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                             DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  dofManager.addField( solidMechanics::totalDisplacement::key(),
                       FieldLocation::Node,
                       3,
                       getMeshTargets() );

  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          solidMechanics::totalDisplacement::key(),
                          DofManager::Connector::Elem );
}


void SolidMechanicsLagrangianFEM::setupSystem( DomainPartition & domain,
                                               DofManager & dofManager,
                                               CRSMatrix< real64, globalIndex > & localMatrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution,
                                               bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                  dofManager.numGlobalDofs(),
                                                  8*8*3*1.2 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    arrayView1d< globalIndex const > const
    dofNumber = nodeManager.getReference< globalIndex_array >( dofManager.getKey( solidMechanics::totalDisplacement::key() ) );


    if( m_contactRelationName != viewKeyStruct::noContactRelationNameString() )
    {
      ElementRegionManager const & elemManager = mesh.getElemManager();
      array1d< string > allFaceElementRegions;
      elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elemRegion )
      {
        allFaceElementRegions.emplace_back( elemRegion.getName() );
      } );

      finiteElement::
        fillSparsity< FaceElementSubRegion,
                      solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic >( mesh,
                                                                                            allFaceElementRegions,
                                                                                            this->getDiscretizationName(),
                                                                                            dofNumber,
                                                                                            dofManager.rankOffset(),
                                                                                            sparsityPattern );

    }
    finiteElement::
      fillSparsity< CellElementSubRegion,
                    solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic >( mesh,
                                                                                          regionNames,
                                                                                          this->getDiscretizationName(),
                                                                                          dofNumber,
                                                                                          dofManager.rankOffset(),
                                                                                          sparsityPattern );


  } );

  sparsityPattern.compress();
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );


}

void SolidMechanicsLagrangianFEM::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  localMatrix.zero();
  localRhs.zero();

  if( m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
  {
    GEOSX_UNUSED_VAR( dt );
    assemblyLaunch< constitutive::SolidBase,
                    solidMechanicsLagrangianFEMKernels::QuasiStaticFactory >( domain,
                                                                              dofManager,
                                                                              localMatrix,
                                                                              localRhs );
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
  {
    assemblyLaunch< constitutive::SolidBase,
                    solidMechanicsLagrangianFEMKernels::ImplicitNewmarkFactory >( domain,
                                                                                  dofManager,
                                                                                  localMatrix,
                                                                                  localRhs,
                                                                                  m_newmarkGamma,
                                                                                  m_newmarkBeta,
                                                                                  m_massDamping,
                                                                                  m_stiffnessDamping,
                                                                                  dt );
  }
}

void
SolidMechanicsLagrangianFEM::
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    string const dofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    fsManager.apply< NodeManager >( time_n + dt,
                                    mesh,
                                    viewKeyStruct::forceString(),
                                    [&]( FieldSpecificationBase const & bc,
                                         string const &,
                                         SortedArrayView< localIndex const > const & targetSet,
                                         NodeManager & targetGroup,
                                         string const & GEOSX_UNUSED_PARAM( fieldName ) )
    {
      // TODO: fix use of dummy name
      bc.applyBoundaryConditionToSystem< FieldSpecificationAdd,
                                         parallelDevicePolicy< 32 > >( targetSet,
                                                                       time_n + dt,
                                                                       targetGroup,
                                                                       solidMechanics::totalDisplacement::key(),
                                                                       dofKey,
                                                                       dofManager.rankOffset(),
                                                                       localMatrix,
                                                                       localRhs );
    } );

  } );

  applyTractionBC( time_n + dt, dofManager, domain, localRhs );

  FaceManager const & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();

  if( faceManager.hasWrapper( "ChomboPressure" ) )
  {
    fsManager.applyFieldValue( time_n, domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ), "ChomboPressure" );
    applyChomboPressure( dofManager, domain, localRhs );
  }

  applyDisplacementBCImplicit( time_n + dt, dofManager, domain, localMatrix, localRhs );
}

real64
SolidMechanicsLagrangianFEM::
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  real64 totalResidualNorm = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > const
    dofNumber = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( solidMechanics::totalDisplacement::key() ) );
    globalIndex const rankOffset = dofManager.rankOffset();

    arrayView1d< integer const > const ghostRank = nodeManager.ghostRank();

    RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

    SortedArrayView< localIndex const > const &
    targetNodes = nodeManager.sets().getReference< SortedArray< localIndex > >( viewKeyStruct::targetNodesString() ).toViewConst();

    forAll< parallelDevicePolicy<> >( targetNodes.size(),
                                      [localRhs, localSum, dofNumber, rankOffset, ghostRank, targetNodes] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      localIndex const nodeIndex = targetNodes[k];
      if( ghostRank[nodeIndex] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[nodeIndex] - rankOffset );

        for( localIndex dim = 0; dim < 3; ++dim )
        {
          localSum += localRhs[localRow + dim] * localRhs[localRow + dim];
        }
      }
    } );
    real64 const localResidualNorm[2] = { localSum.get(), this->m_maxForce };

    // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
    // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
    real64 globalResidualNorm[2] = {0, 0};

    int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
    array1d< real64 > globalValues( size * 2 );

    // Everything is done on rank 0
    MpiWrapper::gather( localResidualNorm,
                        2,
                        globalValues.data(),
                        2,
                        0,
                        MPI_COMM_GEOSX );

    if( rank==0 )
    {
      for( int r=0; r<size; ++r )
      {
        // sum/max across all ranks
        globalResidualNorm[0] += globalValues[r*2];
        globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
      }
    }

    MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );


    real64 const residual = sqrt( globalResidualNorm[0] )/(globalResidualNorm[1]+1); // the + 1 is for the first
                                                                                     // time-step when maxForce = 0;
    totalResidualNorm = std::max( residual, totalResidualNorm );
  } );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    std::cout << GEOSX_FMT( "( R{} ) = ( {:4.2e} ) ; ", coupledSolverAttributePrefix(), totalResidualNorm );
  }

  return totalResidualNorm;
}



void
SolidMechanicsLagrangianFEM::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  dofManager.addVectorToField( localSolution,
                               solidMechanics::totalDisplacement::key(),
                               solidMechanics::incrementalDisplacement::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               solidMechanics::totalDisplacement::key(),
                               solidMechanics::totalDisplacement::key(),
                               scalingFactor );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )

  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addFields( FieldLocation::Node,
                              { solidMechanics::incrementalDisplacement::key(),
                                solidMechanics::totalDisplacement::key() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void SolidMechanicsLagrangianFEM::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    solidMechanics::arrayView2dLayoutIncrDisplacement const & incdisp =
      nodeManager.getField< solidMechanics::incrementalDisplacement >();
    solidMechanics::arrayView2dLayoutTotalDisplacement const & disp =
      nodeManager.getField< solidMechanics::totalDisplacement >();

    // TODO need to finish this rewind
    forAll< parallelDevicePolicy< 32 > >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      for( localIndex i = 0; i < 3; ++i )
      {
        disp( a, i ) -= incdisp( a, i );
        incdisp( a, i ) = 0.0;
      }
    } );
  } );
}


void SolidMechanicsLagrangianFEM::applyContactConstraint( DofManager const & dofManager,
                                                          DomainPartition & domain,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  if( m_contactRelationName != viewKeyStruct::noContactRelationNameString() )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & )
    {
      FaceManager const & faceManager = mesh.getFaceManager();
      NodeManager & nodeManager = mesh.getNodeManager();
      ElementRegionManager & elemManager = mesh.getElemManager();

      solidMechanics::arrayViewConst2dLayoutTotalDisplacement const u =
        nodeManager.getField< solidMechanics::totalDisplacement >();
      arrayView2d< real64 > const fc = nodeManager.getField< solidMechanics::contactForce >();
      fc.zero();

      arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
      ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

      string const dofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
      arrayView1d< globalIndex > const nodeDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
      globalIndex const rankOffset = dofManager.rankOffset();

      // TODO: this bound may need to change
      constexpr localIndex maxNodexPerFace = 4;
      constexpr localIndex maxDofPerElem = maxNodexPerFace * 3 * 2;

      elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
      {
        ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

        real64 const contactStiffness = contact.stiffness();

        arrayView1d< real64 > const area = subRegion.getElementArea();
        arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();

        // TODO: use parallel policy?
        forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {
          real64 Nbar[ 3 ] = { faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0],
                               faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1],
                               faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2] };

          LvArray::tensorOps::normalize< 3 >( Nbar );


          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const kf1 = elemsToFaces[kfe][1];
          localIndex const numNodesPerFace=facesToNodes.sizeOfArray( kf0 );
          real64 const Ja = area[kfe] / numNodesPerFace;

          stackArray1d< globalIndex, maxDofPerElem > rowDOF( numNodesPerFace*3*2 );
          stackArray1d< real64, maxDofPerElem > nodeRHS( numNodesPerFace*3*2 );
          stackArray2d< real64, maxDofPerElem *maxDofPerElem > dRdP( numNodesPerFace*3*2, numNodesPerFace*3*2 );

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            real64 penaltyForce[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( Nbar );
            localIndex const node0 = facesToNodes[kf0][a];
            localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
            real64 gap[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( u[node1] );
            LvArray::tensorOps::subtract< 3 >( gap, u[node0] );
            real64 const gapNormal = LvArray::tensorOps::AiBi< 3 >( gap, Nbar );

            for( int i=0; i<3; ++i )
            {
              rowDOF[3*a+i]                     = nodeDofNumber[node0]+i;
              rowDOF[3*(numNodesPerFace + a)+i] = nodeDofNumber[node1]+i;
            }

            if( gapNormal < 0 )
            {
              LvArray::tensorOps::scale< 3 >( penaltyForce, -contactStiffness * gapNormal * Ja );
              for( int i=0; i<3; ++i )
              {
                LvArray::tensorOps::subtract< 3 >( fc[node0], penaltyForce );
                LvArray::tensorOps::add< 3 >( fc[node1], penaltyForce );
                nodeRHS[3*a+i]                     -= penaltyForce[i];
                nodeRHS[3*(numNodesPerFace + a)+i] += penaltyForce[i];

                dRdP( 3*a+i, 3*a+i )                                         -= contactStiffness * Ja * Nbar[i] * Nbar[i];
                dRdP( 3*a+i, 3*(numNodesPerFace + a)+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
                dRdP( 3*(numNodesPerFace + a)+i, 3*a+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
                dRdP( 3*(numNodesPerFace + a)+i, 3*(numNodesPerFace + a)+i ) -= contactStiffness * Ja * Nbar[i] * Nbar[i];
              }
            }
          }

          for( localIndex idof = 0; idof < numNodesPerFace*3*2; ++idof )
          {
            localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

            if( localRow >= 0 && localRow < localMatrix.numRows() )
            {
              localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                        rowDOF.data(),
                                                                        dRdP[idof].dataIfContiguous(),
                                                                        numNodesPerFace*3*2 );
              RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow], nodeRHS[idof] );
            }
          }
        } );
      } );
    } );
  }
}

real64
SolidMechanicsLagrangianFEM::scalingForSystemSolution( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( domain, dofManager, localSolution );

  return 1.0;
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianFEM, string const &, dataRepository::Group * const )
}
