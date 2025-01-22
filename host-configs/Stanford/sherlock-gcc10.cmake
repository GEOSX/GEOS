# User defined versions
set(GCC_VERSION "10.1.0" CACHE STRING "Version of GCC")
set(OMPI_VERSION "4.1.2" CACHE STRING "Version of Open MPI")
set(OPENBLAS_VERSION "0.3.10" CACHE STRING "Version of OpenBLAS")

set(CONFIG_NAME "sherlock-gcc${GCC_VERSION}-ompi${OMPI_VERSION}-openblas${OPENBLAS_VERSION}" CACHE PATH "")

set(SOFTWARE_ROOT /share/software/user/open CACHE PATH "")
set(GCC_ROOT "${SOFTWARE_ROOT}/gcc/${GCC_VERSION}" CACHE PATH "")
set(MPI_ROOT "${SOFTWARE_ROOT}/openmpi/${OMPI_VERSION}" CACHE PATH "")

site_name(HOST_NAME)

# Compilers
set(CMAKE_C_COMPILER       "${GCC_ROOT}/bin/gcc"      CACHE PATH "")
set(CMAKE_CXX_COMPILER     "${GCC_ROOT}/bin/g++"      CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${GCC_ROOT}/bin/gfortran" CACHE PATH "")

# OpenMP options
#set(ENABLE_OPENMP ON CACHE BOOL "")

# MPI options
set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER       "${MPI_ROOT}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_ROOT}/bin/mpic++"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_ROOT}/bin/mpifort" CACHE PATH "")
set(MPIEXEC              "${MPI_ROOT}/bin/mpirun"  CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# CUDA options
set(ENABLE_CUDA OFF CACHE PATH "" FORCE)

# Blas/Lapack options
set(OPENBLAS_ROOT "${SOFTWARE_ROOT}/openblas/${OPENBLAS_VERSION}" CACHE STRING "")
set(BLAS_LIBRARIES "${OPENBLAS_ROOT}/lib/libopenblas.so" CACHE STRING "")
set(LAPACK_LIBRARIES "${OPENBLAS_ROOT}/lib/libopenblas.so" CACHE STRING "")

set(ENABLE_VALGRIND OFF CACHE BOOL "")
set(ENABLE_CALIPER ON CACHE BOOL "")

if (DEFINED ENV{GEOSX_TPL_DIR})
    set(GEOS_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
endif()
include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
