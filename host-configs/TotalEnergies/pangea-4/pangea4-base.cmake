#######################################
#
# Pangea4 - base config for CPU cluster
#
#   - RAJA   CPU
#   - CHAI   CPU
#   - CUDA   OFF
#   - OPENMP OFF
#   - HYPRE ON
#
#######################################

#######################################
# SCIENTIFIC LIBRARIES
#######################################

set( ENABLE_HYPRE       ON  CACHE BOOL "" FORCE )
set( ENABLE_MATHPRESSO  ON  CACHE BOOL "" FORCE )
set( ENABLE_PAMELA      ON  CACHE BOOL "" FORCE )
set( ENABLE_PETSC       OFF CACHE BOOL "" FORCE )
set( ENABLE_PVTPackage  ON  CACHE BOOL "" FORCE )
set( ENABLE_SCOTCH      ON  CACHE BOOL "" FORCE )
set( ENABLE_SUITESPARSE ON  CACHE BOOL "" FORCE )
set( ENABLE_TRILINOS    OFF CACHE BOOL "" FORCE )
set( ENABLE_VTK         ON  CACHE BOOL "" FORCE )

#######################################
# DEVELOPMENT TOOLS
#######################################

set( ENABLE_DOXYGEN           ON CACHE BOOL "" FORCE )
set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )
set( ENABLE_SPHINX            ON CACHE BOOL "" FORCE )
set( ENABLE_UNCRUSTIFY        ON CACHE BOOL "" FORCE )
set( ENABLE_XML_UPDATES       ON CACHE BOOL "" FORCE )

#######################################
# PERFORMANCE TOOLS
#######################################

set( ENABLE_BENCHMARKS ON CACHE BOOL "" FORCE )
set( ENABLE_CALIPER    ON CACHE BOOL "" FORCE )

#######################################
# RAJA/CHAI SETUP
#######################################

set( ENABLE_OPENMP      OFF         CACHE BOOL "" FORCE )
set( ENABLE_CUDA        OFF         CACHE BOOL "" FORCE )

set( ENABLE_CHAI        ON          CACHE BOOL   "" FORCE )
set( CHAI_BUILD_TYPE    "cpu-no-rm" CACHE STRING "" FORCE )
set( CHAI_ARGS          ""          CACHE STRING "" FORCE )

set( ENABLE_RAJA        ON          CACHE BOOL "" FORCE )
set( RAJA_ENABLE_HIP    OFF         CACHE BOOL "" FORCE )
set( RAJA_ENABLE_OPENMP OFF         CACHE BOOL "" FORCE )
set( RAJA_ENABLE_TBB    OFF         CACHE BOOL "" FORCE )

#######################################
# PYTHON SETUP
#######################################

set( ENABLE_PYGEOSX         ON CACHE BOOL "" )
set( ENABLE_VTK_WRAP_PYTHON ON CACHE BOOL "" )

#######################################
# TPLs
#######################################

set(GEOS_TPL_DIR "$ENV{GEOS_TPL_DIR}" CACHE PATH "" FORCE)
if (NOT DEFINED GEOS_TPL_DIR)
  message(FATAL_ERROR "You must set GEOS_TPL_DIR with -DGEOS_TPL_DIR=")
endif ()
