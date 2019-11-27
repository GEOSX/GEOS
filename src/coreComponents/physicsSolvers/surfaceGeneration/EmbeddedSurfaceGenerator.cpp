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
 * @file EmbeddedSurfaceGenerator.cpp
 */

#include "EmbeddedSurfaceGenerator.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/EmbeddedSurfaceRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp"
#include "meshUtilities/SimpleGeometricObjects/GeometricObjectManager.hpp"
#include "meshUtilities/SimpleGeometricObjects/BoundedThickPlane.hpp"

#ifdef USE_GEOSX_PTP
#include "physicsSolvers/GEOSX_PTP/ParallelTopologyChange.hpp"
#endif
#include <set>

namespace geosx
{
  using namespace dataRepository;
  using namespace constitutive;

EmbeddedSurfaceGenerator::EmbeddedSurfaceGenerator( const std::string& name,
                                    Group * const parent ):
  SolverBase( name, parent ),
  m_solidMaterialName("")
{
  registerWrapper(viewKeyStruct::solidMaterialNameString, &m_solidMaterialName, 0)->
      setInputFlag(InputFlags::REQUIRED)->
      setDescription("Name of the solid material used in solid mechanic solver");

  registerWrapper( viewKeyStruct::fractureRegionNameString, &m_fractureRegionName, 0 )->
      setInputFlag(dataRepository::InputFlags::OPTIONAL)->
      setApplyDefaultValue("FractureRegion");
}

EmbeddedSurfaceGenerator::~EmbeddedSurfaceGenerator()
{
}

void EmbeddedSurfaceGenerator::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0);

    NodeManager * const nodeManager = meshLevel->getNodeManager();
    EdgeManager * const edgeManager = meshLevel->getEdgeManager();

    nodeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::parentIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Parent index of node.");

    nodeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::childIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Child index of node.");

    edgeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::parentIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Parent index of the edge.");

    edgeManager->registerWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::childIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Child index of the edge.");
  }
}

void EmbeddedSurfaceGenerator::InitializePostSubGroups( Group * const problemManager )
{
  /*
   * Here we generate embedded elements for fractures (or faults) that already exist in the domain and
   * were specified in the input file.
   */

  // Get domain
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>( dataRepository::keys::domain );
  // Get geometric object manager
  GeometricObjectManager * geometricObjManager = problemManager->GetGroup<GeometricObjectManager>( "Geometry");

  // Get meshLevel
  Group     * const meshBodies = domain->getMeshBodies();
  MeshBody  * const meshBody   = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel  = meshBody->GetGroup<MeshLevel>(0);

  // Get managers
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  EdgeManager * const edgeManager = meshLevel->getEdgeManager();
  array1d<R1Tensor> const & nodesCoord = nodeManager->referencePosition();

  // Get EmbeddedSurfaceSubRegions
  EmbeddedSurfaceRegion    * const    embeddedSurfaceRegion =
      elemManager->GetRegion<EmbeddedSurfaceRegion>(this->m_fractureRegionName);
  EmbeddedSurfaceSubRegion * const embeddedSurfaceSubRegion =
      embeddedSurfaceRegion->GetSubRegion<EmbeddedSurfaceSubRegion>(0);

  // Loop over all the fracture planes
  geometricObjManager->forSubGroups<BoundedThickPlane>( [&]( BoundedThickPlane * const fracture ) -> void
  {
    /* 1. Find out if an element is cut but the fracture or not.
     * Loop over all the elements and for each one of them loop over the nodes and compute the
     * dot product between the distance between the plane center and the node and the normal
     * vector defining the plane. If two scalar products have different signs the plane cuts the
     * cell. If a nodes gives a 0 dot product it has to be neglected or the method won't work.
     */
    R1Tensor planeCenter  = fracture->getCenter();
    R1Tensor normalVector = fracture->getNormal();
    // Initialize variables
    globalIndex nodeIndex;
    integer isPositive, isNegative;
    R1Tensor distVec;

    elemManager->forElementRegions( [&](ElementRegionBase * const region )->void
    {
      Group * subRegions = region->GetGroup(ElementRegionBase::viewKeyStruct::elementSubRegions);
      subRegions->forSubGroups<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion ) -> void
      {
        FixedOneToManyRelation const & cellToNodes = subRegion->nodeList();
        FixedOneToManyRelation const & cellToEdges = subRegion->edgeList();
        for(localIndex cellIndex =0; cellIndex<subRegion->size(); cellIndex++)
        {
          isPositive = 0;
          isNegative = 0;
          for(localIndex kn =0; kn<subRegion->numNodesPerElement(); kn++)
          {
            nodeIndex = cellToNodes[cellIndex][kn];
            distVec  = nodesCoord[nodeIndex];
            distVec -= planeCenter;
            // check if the dot product is zero
            if ( Dot(distVec, normalVector) > 0 )
            {
              isPositive = 1;
            } else if ( Dot(distVec, normalVector) < 0 )
            {
              isNegative = 1;
            }
          } // end loop over nodes
          if (isPositive * isNegative == 1)
            // TODO in reality this condition is not sufficient because the fracture is a bounded plane. I should also check
            // that the cell is inside the actual fracture plane. For now let us assume the plane is not bounded.
          {
            embeddedSurfaceSubRegion->AddNewEmbeddedSurface( cellIndex,
                                                             normalVector,
                                                             *nodeManager,
                                                             *edgeManager,
                                                             cellToEdges,
                                                             fracture );
          }
        } // end loop over cells
      });// end loop over subregions
    });// end loop over elementRegions
  });// end loop over thick planes

  std::cout << "Number of embedded surface elements: " << embeddedSurfaceSubRegion->size() << std::endl;
}

void EmbeddedSurfaceGenerator::InitializePostInitialConditions_PreSubGroups( Group * const  GEOSX_UNUSED_ARG ( problemManager ) )
{
  // I don't think there is  much to do here.
}


void EmbeddedSurfaceGenerator::postRestartInitialization( Group * const GEOSX_UNUSED_ARG( domain0 ) )
{
  // Not sure about this for now.
  std::cout << "postRestartInitialization \n";
}


real64 EmbeddedSurfaceGenerator::SolverStep( real64 const & GEOSX_UNUSED_ARG( time_n),
                                             real64 const & GEOSX_UNUSED_ARG( dt ),
                                             const int GEOSX_UNUSED_ARG( cycleNumber ),
                                             DomainPartition * const  GEOSX_UNUSED_ARG( domain ) )
{
  real64 rval = 0;
  /*
   * This should be method that generates new fracture elements based on the propagation criterion of choice.
   */
  return rval;
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        EmbeddedSurfaceGenerator,
                        std::string const &, dataRepository::Group * const )

} /* namespace geosx */
