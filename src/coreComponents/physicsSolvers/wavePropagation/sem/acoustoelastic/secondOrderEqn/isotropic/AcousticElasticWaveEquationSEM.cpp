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
 * @file AcousticElasticWaveEquationSEM.cpp
 */

#include "AcousticElasticWaveEquationSEM.hpp"
#include "AcousticElasticWaveEquationSEMKernel.hpp"
#include "AcoustoElasticTimeSchemeSEMKernel.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include <typeinfo>
#include <limits>

namespace geos
{
using namespace dataRepository;
using namespace fields;

void AcousticElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    nodeManager.registerField< acoustoelasticfields::CouplingVectorx >( getName() );
    nodeManager.registerField< acoustoelasticfields::CouplingVectory >( getName() );
    nodeManager.registerField< acoustoelasticfields::CouplingVectorz >( getName() );
  } );
}

void AcousticElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  auto acousSolver = acousticSolver();
  auto elasSolver = elasticSolver();

  auto acousNodesSet = acousSolver->getSolverNodesSet();
  auto elasNodesSet = elasSolver->getSolverNodesSet();

  for( auto val : acousNodesSet )
  {
    if( elasNodesSet.contains( val ) )
      m_interfaceNodesSet.insert( val );
  }
  localIndex const numInterfaceNodes = MpiWrapper::sum( m_interfaceNodesSet.size() );
  GEOS_THROW_IF( numInterfaceNodes == 0, "Failed to compute interface: check xml input (solver order)", std::runtime_error );

  m_acousRegions = acousSolver->getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  m_elasRegions = elasSolver->getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    arrayView2d< real64 const > const faceNormals          = faceManager.faceNormal().toViewConst();
    arrayView2d< real64 const > const faceCenters          = faceManager.faceCenter().toViewConst();
    ArrayOfArraysView< localIndex const > const faceToNode = faceManager.nodeList().toViewConst();
    arrayView2d< localIndex const > const faceToSubRegion  = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const faceToRegion     = faceManager.elementRegionList();
    arrayView2d< localIndex const > const faceToElement    = faceManager.elementList();

    arrayView1d< real32 > const couplingVectorx = nodeManager.getField< acoustoelasticfields::CouplingVectorx >();
    couplingVectorx.zero();

    arrayView1d< real32 > const couplingVectory = nodeManager.getField< acoustoelasticfields::CouplingVectory >();
    couplingVectory.zero();

    arrayView1d< real32 > const couplingVectorz = nodeManager.getField< acoustoelasticfields::CouplingVectorz >();
    couplingVectorz.zero();

    elemManager.forElementRegions( m_acousRegions, [&] ( localIndex const regionIndex, ElementRegionBase const & elemRegion )
    {
      elemRegion.forElementSubRegionsIndex( [&]( localIndex const subRegionIndex, ElementSubRegionBase const & elementSubRegion )
      {
        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
        arrayView2d< real64 const > const elemCenters = elementSubRegion.getElementCenter().toViewConst();

        finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
        {
          using FE_TYPE = TYPEOFREF( finiteElement );

          acousticElasticWaveEquationSEMKernels::CouplingKernel< FE_TYPE > kernelC;
          kernelC.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                                 nodeCoords,
                                                                 regionIndex,
                                                                 subRegionIndex,
                                                                 faceToSubRegion,
                                                                 faceToRegion,
                                                                 faceToElement,
                                                                 faceToNode,
                                                                 faceNormals,
                                                                 faceCenters,
                                                                 elemCenters,
                                                                 couplingVectorx,
                                                                 couplingVectory,
                                                                 couplingVectorz );
        } );
      } );
    } );
  } );
}

real64 AcousticElasticWaveEquationSEM::solverStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   int const,
                                                   DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  auto acousSolver = acousticSolver();
  auto elasSolver = elasticSolver();

  SortedArrayView< localIndex const > const interfaceNodesSet = m_interfaceNodesSet.toViewConst();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const acousticMass = nodeManager.getField< acousticfields::AcousticMassVector >();
    arrayView1d< real32 const > const elasticMass = nodeManager.getField< elasticfields::ElasticMassVector >();
    arrayView1d< localIndex const > const acousticFSNodeIndicator = nodeManager.getField< acousticfields::AcousticFreeSurfaceNodeIndicator >();
    arrayView1d< localIndex const > const elasticFSNodeIndicator = nodeManager.getField< elasticfields::ElasticFreeSurfaceNodeIndicator >();

    arrayView1d< real32 const > const p_n    = nodeManager.getField< acousticfields::Pressure_n >();
    arrayView1d< real32 const > const ux_nm1 = nodeManager.getField< elasticfields::Displacementx_nm1 >();
    arrayView1d< real32 const > const uy_nm1 = nodeManager.getField< elasticfields::Displacementy_nm1 >();
    arrayView1d< real32 const > const uz_nm1 = nodeManager.getField< elasticfields::Displacementz_nm1 >();
    arrayView1d< real32 const > const ux_n   = nodeManager.getField< elasticfields::Displacementx_n >();
    arrayView1d< real32 const > const uy_n   = nodeManager.getField< elasticfields::Displacementy_n >();
    arrayView1d< real32 const > const uz_n   = nodeManager.getField< elasticfields::Displacementz_n >();
    // acoutic -> elastic coupling vectors
    arrayView1d< real32 const > const atoex  = nodeManager.getField< acoustoelasticfields::CouplingVectorx >();
    arrayView1d< real32 const > const atoey  = nodeManager.getField< acoustoelasticfields::CouplingVectory >();
    arrayView1d< real32 const > const atoez  = nodeManager.getField< acoustoelasticfields::CouplingVectorz >();

    arrayView1d< real32 > const p_np1  = nodeManager.getField< acousticfields::Pressure_np1 >();
    arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

    elasSolver->computeUnknowns( time_n, dt, domain, mesh, m_elasRegions );

    AcoustoElasticTimeSchemeSEM::LeapFrog( dt, ux_np1, uy_np1, uz_np1, p_n, elasticMass, atoex, atoey, atoez,
                                           elasticFSNodeIndicator, interfaceNodesSet );

    elasSolver->synchronizeUnknowns( time_n, dt, domain, mesh, m_elasRegions );

    acousSolver->computeUnknowns( time_n, dt, domain, mesh, m_acousRegions );

    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const in )
    {
      localIndex const n = interfaceNodesSet[in];
      if( acousticFSNodeIndicator[n] == 1 )
        return;

      real32 const localIncrement = (
        atoex[n] * ( ux_np1[n] - 2.0 * ux_n[n] + ux_nm1[n] ) +
        atoey[n] * ( uy_np1[n] - 2.0 * uy_n[n] + uy_nm1[n] ) +
        atoez[n] * ( uz_np1[n] - 2.0 * uz_n[n] + uz_nm1[n] )
        ) / acousticMass[n];

      RAJA::atomicAdd< ATOMIC_POLICY >( &p_np1[n], localIncrement );
    } );

    acousSolver->synchronizeUnknowns( time_n, dt, domain, mesh, m_acousRegions );

    acousSolver->prepareNextTimestep( mesh );
    elasSolver->prepareNextTimestep( mesh );
  } );

  return dt;
}

void AcousticElasticWaveEquationSEM::cleanup( real64 const time_n,
                                              integer const cycleNumber,
                                              integer const eventCounter,
                                              real64 const eventProgress,
                                              DomainPartition & domain )
{
  elasticSolver()->cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
  acousticSolver()->cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticElasticWaveEquationSEM, string const &, Group * const )

} /* namespace geos */
