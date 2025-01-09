include(${CMAKE_CURRENT_LIST_DIR}/lassen-clang-13-cuda-11.cmake)

set( GEOS_ENABLE_CONTACT OFF CACHE BOOL "" FORCE)
set( GEOS_ENABLE_FLUIDFLOW ON CACHE BOOL "" FORCE)
set( GEOS_ENABLE_INDUCEDSEISMICITY OFF CACHE BOOL "" FORCE)
set( GEOS_ENABLE_MULTIPHYSICS OFF CACHE BOOL "" FORCE)
set( GEOS_ENABLE_SIMPLEPDE OFF CACHE BOOL "" FORCE)
set( GEOS_ENABLE_SOLIDMECHANICS ON CACHE BOOL "" FORCE)
set( GEOS_ENABLE_SURFACEGENERATION OFF CACHE BOOL "" FORCE)
set( GEOS_ENABLE_WAVEPROPAGATION OFF CACHE BOOL "" FORCE)
# set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
set(ENABLE_CUDA_NVTOOLSEXT ON CACHE BOOL "")

# $<$<COMPILE_LANGUAGE:CUDA>:-ptx>