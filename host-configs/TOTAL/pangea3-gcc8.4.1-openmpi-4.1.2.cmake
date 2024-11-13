# hostconfig for pangea3
#
#
set(CONFIG_NAME "pangea3-gcc8.4.1-ompi-4.1.2" CACHE PATH "")

# Set up the tpls
set(GEOS_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
if (NOT DEFINED GEOS_TPL_DIR)
  message(FATAL_ERROR "You must set GEOS_TPL_DIR with -D GEOS_TPL_DIR=")
endif ()


# C options
set(CMAKE_C_COMPILER gcc CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

# C++ options
set(CMAKE_CXX_COMPILER g++ CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g" CACHE STRING "")

# Fortran options
set(CMAKE_Fortran_COMPILER gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -mcpu=power9 -mtune=power9" CACHE STRING "")
#set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")


# OpenMP options
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
#set(OpenMP_Fortran_FLAGS "-qsmp=omp" CACHE STRING "")
#set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")

# MPI options
if (DEFINED ENV{MPI_ROOT})
  set(ENABLE_MPI ON CACHE BOOL "")
  set(MPI_C_COMPILER         $ENV{MPI_ROOT}/bin/mpicc  CACHE PATH "")
  set(MPI_CXX_COMPILER       $ENV{MPI_ROOT}/bin/mpicxx CACHE PATH "")
  set(MPI_Fortran_COMPILER   $ENV{MPI_ROOT}/bin/mpifort CACHE PATH "")
  set(MPIEXEC                $ENV{MPI_ROOT}/bin/mpirun  CACHE STRING "")
  set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")
else()
  message(FATAL_ERROR "You must have MPI_ROOT variable set, we advise loading module ompi/4.1.2")
endif()

# Cuda options
if (DEFINED ENV{CUDA_ROOT})
  set(ENABLE_CUDA ON CACHE BOOL "")
  set(ENABLE_CUDA_NVTOOLSEXT ON CACHE BOOL "")
  set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_ROOT} CACHE PATH "")
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "")
  set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
  set(CUDA_ARCH sm_70 CACHE STRING "")
  set(CMAKE_CUDA_ARCHITECTURES 70 CACHE STRING "")
  set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-relaxed-constexpr --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3 -Xcompiler -mcpu=powerpc64le -Xcompiler -mtune=powerpc64le" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")

  # Uncomment this line to make nvcc output register usage for each kernel.
  # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --resource-usage" CACHE STRING "" FORCE)
else()
  message(FATAL_ERROR "You must have CUDA_ROOT environment variable set, we advise loading module cuda/11.5.0")
endif()

# GTEST options
set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)



# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

# Silo configure script doesn't recognize systype
set(SILO_BUILD_TYPE powerpc64-unknown-linux-gnu CACHE STRING "")

set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")

set(ENABLE_ESSL OFF CACHE BOOL "")

set(ENABLE_PYGEOSX ON CACHE BOOL "")

if (DEFINED ENV{OPENBLAS_ROOT})
  set(BLAS_LIBRARIES $ENV{OPENBLAS_ROOT}/lib/libopenblas.a)
  set(LAPACK_LIBRARIES $ENV{OPENBLAS_ROOT}/lib/libopenblas.a)
else()
  message(FATAL_ERROR "You must have OPENBLAS_ROOT environment variable set, we advise loading module openblas/0.3.18")
endif()

set(ENABLE_DOXYGEN OFF CACHE PATH "")

set(PETSC_OMP_DIR ${GEOS_TPL_ROOT_DIR}/omp-links-for-petsc CACHE STRING "")

# PETSc doesn't seem to work correctly with clang.
set(ENABLE_PETSC OFF CACHE BOOL "")

set(ENABLE_HYPRE ON CACHE BOOL "")
set(ENABLE_HYPRE_DEVICE "CUDA" CACHE BOOL "")

# disable benchmarks, they are incompatible with P3's nvcc version (cuda 11.5.0)
set(ENABLE_BENCHMARKS OFF CACHE BOOL "")

include( ${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake )
