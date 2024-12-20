set(CCE_VERSION 15.0.0)

include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/frontier-cce@${CCE_VERSION}.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/frontier-base.cmake)

set( CONDUIT_DIR "${GEOS_TPL_DIR}/conduit-0.8.7" CACHE PATH "" )
set( HDF5_DIR "${GEOS_TPL_DIR}/hdf5-1.12.2" CACHE PATH "" )

set( ENABLE_SILO FALSE CACHE BOOL "" )
set( ENABLE_VTK FALSE CACHE BOOL "" )

set( BLAS_DIR "/opt/rocm-${ROCM_VERSION}/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOS_TPL_DIR}/pugixml-1.11.4" CACHE PATH "" )
set( FMT_DIR "${GEOS_TPL_DIR}/fmt-8.0.1" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOS_TPL_DIR}/suite-sparse-5.10.1" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )
set( ENABLE_HYPRE_MIXINT TRUE CACHE STRING "" )

set( HYPRE_DIR "/lustre/orion/geo127/world-shared/hypre/hypre_v2.27.0-218-ge2806c33d_cce-15.0.0_rocm-5.4.3_mixint_umpire-2022.03.0_caliper-2.8.0_rel/" CACHE PATH "" )

set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
set( CALIPER_DIR "${GEOS_TPL_DIR}/caliper-2.8.0" CACHE PATH "" )
