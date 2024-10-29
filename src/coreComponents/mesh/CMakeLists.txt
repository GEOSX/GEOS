# Specify all headers
set( mesh_headers
     BufferOps.hpp
     CellElementRegion.hpp
     CellElementSubRegion.hpp
     DomainPartition.hpp
     EdgeManager.hpp
     ElementRegionBase.hpp
     ElementRegionManager.hpp
     ElementSubRegionBase.hpp
     ElementType.hpp
     EmbeddedSurfaceNodeManager.hpp
     EmbeddedSurfaceSubRegion.hpp
     FaceElementSubRegion.hpp
     FaceManager.hpp
     FieldIdentifiers.hpp
     InterObjectRelation.hpp
     MeshBody.hpp
     MeshFields.hpp
     MeshForLoopInterface.hpp
     MeshLevel.hpp
     MeshManager.hpp
     MeshObjectPath.hpp
     NodeManager.hpp
     ObjectManagerBase.hpp
     ParticleManager.hpp
     ParticleRegion.hpp
     ParticleRegionBase.hpp
     ParticleSubRegion.hpp
     ParticleSubRegionBase.hpp
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
     generators/ExternalMeshGeneratorBase.hpp
     generators/FaceBlock.hpp
     generators/FaceBlockABC.hpp
     generators/InternalMeshGenerator.hpp
     generators/InternalWellGenerator.hpp
     generators/InternalWellboreGenerator.hpp
     generators/LineBlock.hpp
     generators/LineBlockABC.hpp
     generators/MeshGeneratorBase.hpp
     generators/ParMETISInterface.hpp
     generators/PartitionDescriptor.hpp
     generators/PartitionDescriptorABC.hpp
     generators/PrismUtilities.hpp
     generators/WellGeneratorABC.hpp
     generators/WellGeneratorBase.hpp
     mpiCommunications/CommID.hpp
     mpiCommunications/CommunicationTools.hpp
     mpiCommunications/MPI_iCommData.hpp
     mpiCommunications/NeighborCommunicator.hpp
     mpiCommunications/NeighborData.hpp
     mpiCommunications/PartitionBase.hpp
     mpiCommunications/SpatialPartition.hpp
     particleGenerators/ParticleBlock.hpp
     particleGenerators/ParticleBlockABC.hpp
     particleGenerators/ParticleBlockManager.hpp
     particleGenerators/ParticleMeshGenerator.hpp
     particleGenerators/ParticleMeshGeneratorBase.hpp
     simpleGeometricObjects/Box.hpp
     simpleGeometricObjects/CustomPolarObject.hpp
     simpleGeometricObjects/Cylinder.hpp
     simpleGeometricObjects/Disc.hpp
     simpleGeometricObjects/GeometricObjectManager.hpp
     simpleGeometricObjects/PlanarGeometricObject.hpp
     simpleGeometricObjects/Rectangle.hpp
     simpleGeometricObjects/SimpleGeometricObjectBase.hpp
     simpleGeometricObjects/ThickPlane.hpp
     utilities/AverageOverQuadraturePointsKernel.hpp     
     utilities/CIcomputationKernel.hpp
     utilities/ComputationalGeometry.hpp
     utilities/MeshMapUtilities.hpp
     utilities/StructuredGridUtilities.hpp )

# Specify all sources
set( mesh_sources
     BufferOps.cpp
     CellElementRegion.cpp
     CellElementSubRegion.cpp
     DomainPartition.cpp
     EdgeManager.cpp
     ElementRegionBase.cpp
     ElementRegionManager.cpp
     ElementSubRegionBase.cpp
     EmbeddedSurfaceNodeManager.cpp
     EmbeddedSurfaceSubRegion.cpp
     FaceElementSubRegion.cpp
     FaceManager.cpp
     MeshBody.cpp
     MeshLevel.cpp
     MeshManager.cpp
     MeshObjectPath.cpp
     NodeManager.cpp
     ObjectManagerBase.cpp
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
     generators/ExternalMeshGeneratorBase.cpp
     generators/FaceBlock.cpp
     generators/InternalMeshGenerator.cpp
     generators/InternalWellGenerator.cpp
     generators/InternalWellboreGenerator.cpp
     generators/LineBlock.cpp
     generators/MeshGeneratorBase.cpp
     generators/ParMETISInterface.cpp
     generators/PartitionDescriptor.cpp
     generators/WellGeneratorBase.cpp
     mpiCommunications/CommID.cpp
     mpiCommunications/CommunicationTools.cpp
     mpiCommunications/MPI_iCommData.cpp
     mpiCommunications/NeighborCommunicator.cpp
     mpiCommunications/PartitionBase.cpp
     mpiCommunications/SpatialPartition.cpp
     particleGenerators/ParticleBlock.cpp
     particleGenerators/ParticleBlockManager.cpp
     particleGenerators/ParticleMeshGenerator.cpp
     particleGenerators/ParticleMeshGeneratorBase.cpp
     simpleGeometricObjects/Box.cpp
     simpleGeometricObjects/CustomPolarObject.cpp
     simpleGeometricObjects/Cylinder.cpp
     simpleGeometricObjects/Disc.cpp
     simpleGeometricObjects/GeometricObjectManager.cpp
     simpleGeometricObjects/PlanarGeometricObject.cpp
     simpleGeometricObjects/Rectangle.cpp
     simpleGeometricObjects/SimpleGeometricObjectBase.cpp
     simpleGeometricObjects/ThickPlane.cpp
     utilities/ComputationalGeometry.cpp )

set( dependencyList ${parallelDeps} schema dataRepository constitutive finiteElement parmetis metis )

if( ENABLE_VTK )
    message(STATUS "Adding VTK readers")
    set( mesh_headers ${mesh_headers}
         generators/CollocatedNodes.hpp
         generators/VTKFaceBlockUtilities.hpp
         generators/VTKMeshGenerator.hpp
         generators/VTKMeshGeneratorTools.hpp
         generators/VTKWellGenerator.hpp
         generators/VTKUtilities.hpp
         )
    set( mesh_sources ${mesh_sources}
         generators/CollocatedNodes.cpp
         generators/VTKFaceBlockUtilities.cpp
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

blt_add_library( NAME       mesh
                 SOURCES    ${mesh_sources}
                 HEADERS    ${mesh_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS} )

target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

if( ENABLE_VTK )
    # Avoid compiling with nvcc which sometimes segfaults in debug
    set_source_files_properties( generators/VTKMeshGeneratorTools.cpp PROPERTIES LANGUAGE CXX )
endif()

if( GEOS_ENABLE_TESTS )
    add_subdirectory( unitTests )
endif( )

