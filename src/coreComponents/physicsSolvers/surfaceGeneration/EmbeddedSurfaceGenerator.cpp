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
#include "mesh/FaceElementRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp"


#ifdef USE_GEOSX_PTP
#include "GEOSX_PTP/ParallelTopologyChange.hpp"
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
  // TODO Auto-generated destructor stub
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

void EmbeddedSurfaceGenerator::InitializePostSubGroups( Group * const GEOSX_UNUSED_ARG ( problemManager ) )
{

}

void EmbeddedSurfaceGenerator::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_ARG( problemManager ) )
{
  std::cout << "2. InitializePostInitialConditions_PreSubGroups \n";
}


void EmbeddedSurfaceGenerator::postRestartInitialization( Group * const GEOSX_UNUSED_ARG( domain0 ) )
{
  std::cout << "postRestartInitialization \n";
}


real64 EmbeddedSurfaceGenerator::SolverStep( real64 const & GEOSX_UNUSED_ARG( time_n),
                                             real64 const & GEOSX_UNUSED_ARG( dt ),
                                             const int GEOSX_UNUSED_ARG( cycleNumber ),
                                             DomainPartition * const  domain )
{
  real64 rval = 0;

  // hardcoding the plane
  R1Tensor planeCenter  = {0, 0, 0};
  R1Tensor normalVector = {0, 1, 0};

  // Get meshLevel
  Group * const meshBodies = domain->getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);

  // Get managers
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  // NodeManager * const nodeManager = meshLevel->getNodeManager();
  // EdgeManager * const edgeManager = meshLevel->getEdgeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();
  array1d<R1Tensor> & faceCenter  = faceManager->faceCenter();

  // Initialize variables
  globalIndex faceIndex;
  real64 sumScalarProduct;
  integer numEmbeddedSurfaceElem = 0;
  R1Tensor distVec;

  // 1. Count the number of embedded surface elemetns
  elemManager->forElementRegions( [&](ElementRegionBase * const region )->void
      {
        Group * subRegions = region->GetGroup(ElementRegionBase::viewKeyStruct::elementSubRegions);
        subRegions->forSubGroups<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion ) -> void
        {
          FixedOneToManyRelation const & cellToFaces = subRegion->faceList();
          for(localIndex cellIndex =0; cellIndex<subRegion->size(); cellIndex++)
          {
            sumScalarProduct = 0;
            for(localIndex kf =0; kf<subRegion->numFacesPerElement(); kf++)
            {
              faceIndex = cellToFaces[cellIndex][kf];
              distVec  = planeCenter;
              distVec -= faceCenter[faceIndex];
              sumScalarProduct *= Dot(distVec, normalVector);
            }
            if (sumScalarProduct < 0)
            {
              numEmbeddedSurfaceElem += 1;
            }
          }
        });

      });
  std::cout << "number of embedded surface elements: " << numEmbeddedSurfaceElem;
  return rval;
}


REGISTER_CATALOG_ENTRY( SolverBase,
                        EmbeddedSurfaceGenerator,
                        std::string const &, dataRepository::Group * const )

} /* namespace geosx */
