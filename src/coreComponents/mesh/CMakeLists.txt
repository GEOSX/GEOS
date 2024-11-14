# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: mesh

Contains:
 - components allowing to import, store and access the mesh.
 - partitioning and communication tools (ParMETIS, Scotch and VTK interfaces).
 - GEOS components defining simple geometric objects.
 - basic geometric and mesh utilities.
#]]

#
# Specify all headers
#
set( mesh_headers
     BufferOps.hpp
     CellElementRegion.hpp
     CellElementRegionSelector.hpp
     CellElementSubRegion.hpp
     DomainPartition.hpp
     EdgeManager.hpp
     ElementRegionBase.hpp
     ElementRegionManager.hpp
     ElementSubRegionBase.hpp
     ElementType.hpp
     ExternalDataSourceBase.hpp
     ExternalDataSourceManager.hpp
     EmbeddedSurfaceNodeManager.hpp
     EmbeddedSurfaceSubRegion.hpp
     MeshFields.hpp
     FaceElementSubRegion.hpp
     FaceManager.hpp
     FieldIdentifiers.hpp
     InterObjectRelation.hpp
     MeshBody.hpp
     MeshForLoopInterface.hpp
     MeshLevel.hpp
     MeshManager.hpp
     MeshObjectPath.hpp
     NodeManager.hpp
     ObjectManagerBase.hpp
     generators/ParticleBlock.hpp
     generators/ParticleBlockManager.hpp
     ParticleManager.hpp
     ParticleRegion.hpp
     ParticleRegionBase.hpp
     ParticleSubRegion.hpp
     ParticleSubRegionBase.hpp
     generators/ParticleBlockABC.hpp
     ParticleType.hpp
     Perforation.hpp
     PerforationData.hpp
     PerforationFields.hpp
     SurfaceElementRegion.hpp
     SurfaceElementSubRegion.hpp
     ToElementRelation.hpp
     ToParticleRelation.hpp
     WellElementRegion.hpp
     WellElementSubRegion.hpp
     generators/CellBlock.hpp
     generators/CellBlockABC.hpp
     generators/CellBlockManager.hpp
     generators/CellBlockManagerABC.hpp
     generators/CellBlockUtilities.hpp
     generators/LineBlock.hpp
     generators/LineBlockABC.hpp
     generators/ExternalMeshGeneratorBase.hpp
     generators/FaceBlock.hpp
     generators/FaceBlockABC.hpp
     generators/InternalMeshGenerator.hpp
     generators/InternalWellGenerator.hpp
     generators/InternalWellboreGenerator.hpp
     generators/MeshComponentBase.hpp
     generators/MeshGeneratorBase.hpp
     generators/ParMETISInterface.hpp
     generators/ParticleMeshGenerator.hpp
     generators/PartitionDescriptor.hpp
     generators/PrismUtilities.hpp
     generators/Region.hpp
     generators/WellGeneratorBase.hpp
     mpiCommunications/CommID.hpp
     mpiCommunications/CommunicationTools.hpp
     mpiCommunications/MPI_iCommData.hpp
     mpiCommunications/NeighborCommunicator.hpp
     mpiCommunications/NeighborData.hpp
     mpiCommunications/PartitionBase.hpp
     mpiCommunications/SpatialPartition.hpp
     simpleGeometricObjects/Rectangle.hpp
     simpleGeometricObjects/Disc.hpp
     simpleGeometricObjects/CustomPolarObject.hpp
     simpleGeometricObjects/Box.hpp
     simpleGeometricObjects/Cylinder.hpp
     simpleGeometricObjects/GeometricObjectManager.hpp
     simpleGeometricObjects/SimpleGeometricObjectBase.hpp
     simpleGeometricObjects/PlanarGeometricObject.hpp
     simpleGeometricObjects/ThickPlane.hpp
     utilities/AverageOverQuadraturePointsKernel.hpp
     utilities/CIcomputationKernel.hpp
     utilities/ComputationalGeometry.hpp
     utilities/MeshMapUtilities.hpp
     utilities/StructuredGridUtilities.hpp )

#
# Specify all sources
#
set( mesh_sources
     BufferOps.cpp
     CellElementRegion.cpp
     CellElementRegionSelector.cpp
     CellElementSubRegion.cpp
     DomainPartition.cpp
     EdgeManager.cpp
     ElementRegionBase.cpp
     ElementRegionManager.cpp
     ElementSubRegionBase.cpp
     EmbeddedSurfaceNodeManager.cpp
     EmbeddedSurfaceSubRegion.cpp
     ExternalDataSourceBase.cpp
     ExternalDataSourceManager.cpp
     FaceElementSubRegion.cpp
     FaceManager.cpp
     MeshBody.cpp
     MeshLevel.cpp
     MeshManager.cpp
     MeshObjectPath.cpp
     NodeManager.cpp
     ObjectManagerBase.cpp
     generators/ParticleBlock.cpp
     generators/ParticleBlockManager.cpp
     ParticleManager.cpp
     ParticleRegion.cpp
     ParticleRegionBase.cpp
     ParticleSubRegion.cpp
     ParticleSubRegionBase.cpp
     Perforation.cpp
     PerforationData.cpp
     SurfaceElementRegion.cpp
     SurfaceElementSubRegion.cpp
     ToElementRelation.cpp
     ToParticleRelation.cpp
     WellElementRegion.cpp
     WellElementSubRegion.cpp
     generators/CellBlock.cpp
     generators/CellBlockManager.cpp
     generators/CellBlockUtilities.cpp
     generators/LineBlock.cpp
     generators/ExternalMeshGeneratorBase.cpp
     generators/FaceBlock.cpp
     generators/InternalMeshGenerator.cpp
     generators/InternalWellGenerator.cpp
     generators/InternalWellboreGenerator.cpp
     generators/MeshComponentBase.cpp
     generators/MeshGeneratorBase.cpp
     generators/ParMETISInterface.cpp
     generators/ParticleMeshGenerator.cpp
     generators/Region.cpp
     generators/WellGeneratorBase.cpp
     mpiCommunications/CommID.cpp
     mpiCommunications/CommunicationTools.cpp
     mpiCommunications/MPI_iCommData.cpp
     mpiCommunications/NeighborCommunicator.cpp
     mpiCommunications/PartitionBase.cpp
     mpiCommunications/SpatialPartition.cpp
     simpleGeometricObjects/Rectangle.cpp
     simpleGeometricObjects/Disc.cpp
     simpleGeometricObjects/CustomPolarObject.cpp
     simpleGeometricObjects/Box.cpp
     simpleGeometricObjects/Cylinder.cpp
     simpleGeometricObjects/GeometricObjectManager.cpp
     simpleGeometricObjects/SimpleGeometricObjectBase.cpp
     simpleGeometricObjects/PlanarGeometricObject.cpp
     simpleGeometricObjects/ThickPlane.cpp
     utilities/ComputationalGeometry.cpp )

set( dependencyList ${parallelDeps} schema constitutive finiteElement parmetis metis )

if( ENABLE_VTK )
    message(STATUS "Adding VTK readers")
    set( mesh_headers ${mesh_headers}
         generators/CollocatedNodes.hpp
         generators/VTKFaceBlockUtilities.hpp
         generators/VTKHierarchicalDataSource.hpp
         generators/VTKMeshGenerator.hpp
         generators/VTKMeshGeneratorTools.hpp
         generators/VTKWellGenerator.hpp
         generators/VTKUtilities.hpp         
         )
    set( mesh_sources ${mesh_sources}
         generators/CollocatedNodes.cpp
         generators/VTKFaceBlockUtilities.cpp
         generators/VTKHierarchicalDataSource.cpp
         generators/VTKMeshGenerator.cpp
         generators/VTKMeshGeneratorTools.cpp
         generators/VTKWellGenerator.cpp
         generators/VTKUtilities.cpp
         )
    list( APPEND dependencyList VTK::IOLegacy VTK::FiltersParallelDIY2 )
    if( ENABLE_MPI )
      list( APPEND dependencyList VTK::IOParallelXML VTK::ParallelMPI )
    endif()
endif()

if( ENABLE_SCOTCH )
    set( mesh_headers ${mesh_headers} generators/PTScotchInterface.hpp )
    set( mesh_sources ${mesh_sources} generators/PTScotchInterface.cpp )
    list( APPEND dependencyList ptscotch )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       mesh
                 SOURCES    ${mesh_sources}
                 HEADERS    ${mesh_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS mesh LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( ENABLE_VTK )
    # Avoid compiling with nvcc which sometimes segfaults in debug
    set_source_files_properties( generators/VTKMeshGeneratorTools.cpp PROPERTIES LANGUAGE CXX )
endif()

if( GEOS_ENABLE_TESTS )
    add_subdirectory( unitTests )
endif( )

