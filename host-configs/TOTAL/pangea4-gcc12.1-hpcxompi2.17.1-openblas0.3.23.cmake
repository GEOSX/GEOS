#######################################
#
# Pangea4 - gcc - hpcxompi - openblas
#
# Uses :
#   - cray wrappers for gcc (cc, CC, ftn)
#   - OpenBLAS      for BLAS and LAPACK
#   - HPC-X OpenMPI for MPI
#
#######################################
#
# Requires modules :
#   - PrgEnv-gnu           = 8.4.0, which loads :
#     . gcc                = 12.1
#     . craype             = 2.7.23
#     . cray-dsmml         = 0.2.2
#     . craype-network-ofi = 1.0
#     . libfabric          = 1.13.1
#   - cmake                = 3.27.2
#   - cray-python          = 3.10.10
#   - craype-x86-milan     = 1.0
#     PrgEnv-gnu loads gcc 12 that does not support craype-x86-genoa
#   - hpcx                 = 2.17.1
#   - openblas             = 0.3.23
#
# Load modules this way :
#   - module purge
#   - module load PrgEnv-gnu/8.4.0 craype-x86-milan cmake/3.27.2 cray-python/3.10.10
#   - module unload cray-libsci/23.09.1.1 cray-mpich/8.1.27
#   - module load hpcx openblas/0.3.23
#
########################################

set( CONFIG_NAME "pangea4-gcc12.1-hpcxompi2.17.1-openblas0.3.23" CACHE PATH "" )

include(${CMAKE_CURRENT_LIST_DIR}/pangea4-base.cmake)

#######################################
# COMPILER SETUP
#######################################

# use :
#  - cray wrappers for gnu compilers so that we link properly with infiniband network
#  - explicit optimization flags even when using cray wrappers

if( NOT DEFINED ENV{GCC_PATH} )
    message( FATAL_ERROR "GCC is not loaded. Please load the PrgEnv-gnu/8.4.0 module." )
endif()

set( CMAKE_C_COMPILER       "cc"  CACHE PATH "" )
set( CMAKE_CXX_COMPILER     "CC"  CACHE PATH "" )
set( CMAKE_Fortran_COMPILER "ftn" CACHE PATH "" )

set( COMMON_FLAGS  "-march=native -mtune=native" )
set( RELEASE_FLAGS "-O3 -DNDEBUG"                )
set( DEBUG_FLAGS   "-O0 -g"                      )

set( CMAKE_C_FLAGS               ${COMMON_FLAGS}  CACHE STRING "" )
set( CMAKE_CXX_FLAGS             ${COMMON_FLAGS}  CACHE STRING "" )
set( CMAKE_Fortran_FLAGS         ${COMMON_FLAGS}  CACHE STRING "" )
set( CMAKE_CXX_FLAGS_RELEASE     ${RELEASE_FLAGS} CACHE STRING "" )
set( CMAKE_C_FLAGS_RELEASE       ${RELEASE_FLAGS} CACHE STRING "" )
set( CMAKE_Fortran_FLAGS_RELEASE ${RELEASE_FLAGS} CACHE STRING "" )
set( CMAKE_CXX_FLAGS_DEBUG       ${DEBUG_FLAGS}   CACHE STRING "" )
set( CMAKE_C_FLAGS_DEBUG         ${DEBUG_FLAGS}   CACHE STRING "" )
set( CMAKE_Fortran_FLAGS_DEBUG   ${DEBUG_FLAGS}   CACHE STRING "" )

#######################################
# MPI SETUP
#######################################

# use :
# - HPC-X OpenMPI library

set( ENABLE_MPI ON CACHE BOOL "" )

if( NOT DEFINED ENV{HPCX_MPI_DIR} )
    message( FATAL_ERROR "HPC-X OpenMPI is not loaded. Please load the hpcx module." )
endif()

#######################################                                                                                                                                           
# BLAS/LAPACK SETUP                                                                                                                                                               
#######################################

# use :
# - OpenBLAS library

find_library(OPENBLAS_LIB openblas)

if(NOT OPENBLAS_LIB)
    message(FATAL_ERROR "OpenBLAS is not loaded. Please load the openblas/0.3.23 module.")
endif()

include( ${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake )
