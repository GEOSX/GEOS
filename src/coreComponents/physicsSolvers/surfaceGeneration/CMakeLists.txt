# Specify solver headers
set( physicsSolvers_headers
     ${physicsSolvers_headers}
     surfaceGeneration/EmbeddedSurfaceGenerator.hpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.hpp
     surfaceGeneration/ParallelTopologyChange.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     surfaceGeneration/SurfaceGeneratorFields.hpp
     surfaceGeneration/LogLevelsInfo.hpp
     surfaceGeneration/kernels/surfaceGenerationKernels.hpp
     surfaceGeneration/kernels/surfaceGenerationKernelsHelpers.hpp
     PARENT_SCOPE )

# Specify solver sources
set( physicsSolvers_sources
     ${physicsSolvers_sources}
     surfaceGeneration/EmbeddedSurfaceGenerator.cpp
     surfaceGeneration/EmbeddedSurfacesParallelSynchronization.cpp
     surfaceGeneration/ParallelTopologyChange.cpp
     surfaceGeneration/ParallelTopologyChangeNoFixup.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     PARENT_SCOPE )