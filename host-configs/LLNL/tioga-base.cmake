set( GEOS_ENABLE_VEM OFF CACHE BOOL "" FORCE )

set( ENABLE_MATHPRESSO OFF CACHE BOOL "" )
set( ENABLE_PAMELA OFF CACHE BOOL "" )
set( ENABLE_PVTPackage ON CACHE BOOL "" )
set( ENABLE_PETSC OFF CACHE BOOL "" FORCE )
set( ENABLE_CALIPER ON CACHE BOOL "" )
set( ENABLE_PAPI OFF CACHE BOOL "" )
set( ENABLE_ESSL OFF CACHE BOOL "" )
set( ENABLE_TRILINOS OFF CACHE BOOL "" )
set( ENABLE_VTK ON CACHE BOOL "" )
set( ENABLE_OPENMP OFF CACHE BOOL "" FORCE )

set( CAMP_STANDALONE TRUE CACHE BOOL "" )

# ROCM options
set( ENABLE_ROCM ON CACHE BOOL "" FORCE )
set( ROCM_ROOT "${HIP_ROOT}" CACHE PATH "" )

set( GEOS_BUILD_OBJ_LIBS OFF CACHE BOOL "" FORCE )
set( ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "" )
set( gtest_disable_pthreads ON CACHE BOOL "" )

set( ENABLE_TESTS ON CACHE BOOL "" FORCE )
set( ENABLE_EXAMPLES ON CACHE BOOL "" FORCE )
set( ENABLE_BENCHMARKS ON CACHE BOOL "" FORCE )
set( ENABLE_DOCS OFF CACHE BOOL "" FORCE )

set( ENABLE_SCOTCH OFF CACHE BOOL "" FORCE )
set( ENABLE_SUPERLU_DIST OFF CACHE BOOL "" FORCE )

# HYPRE options
set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )
set( ENABLE_HYPRE_MIXINT ON CACHE STRING "" )

# TPLs
include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
