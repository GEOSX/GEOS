#######################################
#
# Pangea2 - Base hostconfig
#
#
########################################

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )
set( ENABLE_UNCRUSTIFY ON CACHE BOOL "" FORCE )
set( ENABLE_VTK ON CACHE BOOL "" FORCE )
set( ENABLE_PVTPackage ON CACHE BOOL "" FORCE )
set( ENABLE_PETSC OFF CACHE BOOL "" FORCE )
set( ENABLE_MATHPRESSO ON CACHE BOOL "" FORCE )
set( ENABLE_HYPRE ON CACHE BOOL "" FORCE )
set( ENABLE_XML_UPDATES OFF CACHE BOOL "" FORCE )
set( ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE )
set( ENABLE_BENCHMARKS OFF CACHE BOOL "" FORCE )

#######################################
# RAJA/CHAI SETUP
#######################################
option( RAJA_ENABLE_TBB "" OFF )
option( ENABLE_CALIPER "Enables CALIPER" ON )

set( ENABLE_CUDA      "OFF"       CACHE PATH "" FORCE )
set( CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE )
set( CHAI_ARGS         ""          CACHE PATH "" FORCE )

set( ENABLE_OPENMP     "OFF"        CACHE PATH "" FORCE )
set( RAJA_ENABLE_OPENMP "OFF"        CACHE PATH "" FORCE )

include( ${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake )
