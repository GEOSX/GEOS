include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/tioga-cce-15.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/tioga-base.cmake)

set( CONDUIT_DIR "${GEOS_TPL_DIR}/conduit-0.8.7" CACHE PATH "" )
set( HDF5_DIR "${GEOS_TPL_DIR}/hdf5-1.14.1-2" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-5.4.3/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOS_TPL_DIR}/pugixml-1.13" CACHE PATH "" )
set( FMT_DIR "${GEOS_TPL_DIR}/fmt-10.0.0" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOS_TPL_DIR}/suite-sparse-5.10.1" CACHE PATH "" )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )
set( ENABLE_HYPRE_MIXINT FALSE CACHE STRING "" )
set( HYPRE_DIR "${GEOS_TPL_DIR}/hypre-develop" CACHE PATH "" )

set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
set( CALIPER_DIR "${GEOS_TPL_DIR}/caliper-2.8.0" CACHE PATH "" )

# haven't build I/O TPLs on tioga yet
set( ENABLE_SILO OFF CACHE BOOL "" FORCE )