# CVX .cmake for Intel 1-API compilers and MPI 
#
# file: CPU-GCC_10.2.0-ompi_hpcx.cmake
# 
# Michael E. Thomadakis michael.thomadakis@chevron.com

# detect host and name the configuration file
site_name(HOST_NAME)

#set(CONFIG_NAME "CPU-OPTO3-Hypre-GCC_10.2.0-ompi_hpcx-OMP" CACHE PATH "")
# message("## CVX CONFIG_NAME = ${CONFIG_NAME}")

message("## CVX CMAKE_CURRENT_LIST_DIR = ${CMAKE_CURRENT_LIST_DIR}") 

# Default is ON
set(ENABLE_WARNINGS_AS_ERRORS OFF CACHE BOOL "")
   
# Local mtml settings 
set(GEOSX_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
set(GEOS_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE) 
set(GEOSX_DIR "$ENV{GEOSX_DIR}" CACHE PATH "" FORCE)
message("## GEOSX_TPL_DIR = ${GEOSX_TPL_DIR}")
message("## GEOS_TPL_DIR = ${GEOS_TPL_DIR}")
message("## GEOSX_DIR = ${GEOSX_DIR}")

# set paths to C, C++, and Fortran compilers. Note that while GEOSX does not contain any Fortran code,
# some of the third-party libraries do contain Fortran code. Thus a Fortran compiler must be specified.
set(CMAKE_C_COMPILER "$ENV{GEOSX_CC}" CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE " -Wno-error -pthread -O3 -fno-fast-math -mdaz-ftz -DNDEBUG " CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

set(CMAKE_CXX_COMPILER "$ENV{GEOSX_CXX}" CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

set(CMAKE_Fortran_COMPILER "$ENV{GEOSX_FORT}" CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# enable MPI and set paths to compilers and executable.
# Note that the MPI compilers are wrappers around standard serial compilers.
# Therefore, the MPI compilers must wrap the appropriate serial compilers specified
# in CMAKE_C_COMPILER, CMAKE_CXX_COMPILER, and CMAKE_Fortran_COMPILER.
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "$ENV{GEOSX_MPICC}" CACHE PATH "")
set(MPI_CXX_COMPILER "$ENV{GEOSX_MPICXX}" CACHE PATH "")
set(MPI_Fortran_COMPILER "$ENV{GEOSX_MPIFORT}" CACHE PATH "")
set(MPIEXEC "$ENV{GEOSX_MPIRUN}" CACHE PATH "")

set(BLAS_LIBRARIES "/data/saet/software/aarch64/lib/libopenblas.so" CACHE STRING "")
# set(BLAS_LIBRARIES "$ENV{NVHPC_ROOT}/math_libs/lib64/libnvblas.so" CACHE STRING "")
message("## BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE )

# MikeT : Check if there is value to enable these 
# set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)

# disable CUDA and OpenMP
set(CUDA_ENABLED OFF CACHE BOOL "" FORCE)
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)

set(ENABLE_SPHINX_EXECUTABLE ON CACHE BOOL "")
set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")

# enable PAMELA and PVTPackage
set(ENABLE_PAMELA OFF CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
set(ENABLE_VALGRIND OFF CACHE BOOL "")
set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_CALIPER_HYPRE "$ENV{ENABLE_CALIPER_HYPRE}" CACHE BOOL "")

# Hypre ON
if(NOT DEFINED ENABLE_HYPRE)
  set(ENABLE_HYPRE ON CACHE BOOL "" FORCE)
  # set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE) 
  set(ENABLE_TRILINOS "$ENV{ENABLE_TRILINOS}" CACHE BOOL "" FORCE)
endif()
if(ENABLE_HYPRE)
  set(GEOSX_LA_INTERFACE "Hypre" CACHE STRING "" FORCE)
else()
  set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
endif()

if( (ENABLE_HYPRE AND ENABLE_TRILINOS) OR (NOT ENABLE_TRILINOS AND NOT ENABLE_HYPRE))
  MESSAGE(SEND_ERROR "Exactly one of ENABLE_HYPRE and ENABLE_TRILINOS must be defined.")
  MESSAGE(SEND_ERROR "ENABLE_HYPRE = ${ENABLE_HYPRE}.")
  MESSAGE(SEND_ERROR "ENABLE_TRILINOS = ${ENABLE_TRILINOS}.")
endif()

MESSAGE(STATUS "## GEOSX_LA_INTERFACE = ${GEOSX_LA_INTERFACE}")

# disable Doxygen
set(ENABLE_DOXYGEN OFF CACHE PATH "")

# enable tests
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

# let GEOSX define some third party libraries information for you
message("## tpls.cmake = ${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake")
include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
message("## CMAKE_CURRENT_LIST_DIR = ${CMAKE_CURRENT_LIST_DIR}")
