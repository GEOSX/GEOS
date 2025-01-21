#------------------------------------------------------------------------------
# CUDA support
#------------------------------------------------------------------------------
# _blt_tutorial_surface_cuda_config_start
#
# CVX .cmake for Nvidia nvcc and OpenMPI MPI; OpenMP enabled
# For Nvidia LauncPad environment
#
# Michael E. Thomadakis michael.thomadakis@chevron.com
# Chevron Technical Center, Innovation and HPC R&D
# 

# detect host and name the configuration file
site_name(HOST_NAME)
message("## CVX CMAKE_CURRENT_LIST_DIR = ${CMAKE_CURRENT_LIST_DIR}") 

# Default is ON
set(ENABLE_WARNINGS_AS_ERRORS OFF CACHE BOOL "")

# Info
set(GEOSX_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
set(GEOS_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
set(GEOSX_DIR "$ENV{GEOSX_DIR}" CACHE PATH "" FORCE)
message("## GEOSX_TPL_DIR = ${GEOSX_TPL_DIR}")
message("## GEOS_TPL_DIR = ${GEOS_TPL_DIR}")
message("## GEOSX_DIR = ${GEOSX_DIR}")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "$ENV{GEOSX_MPICC}" CACHE PATH "")
set(MPI_CXX_COMPILER "$ENV{GEOSX_MPICXX}" CACHE PATH "")
set(MPI_Fortran_COMPILER "$ENV{GEOSX_MPIFORT}" CACHE PATH "")
set(MPIEXEC "$ENV{GEOSX_MPIRUN}" CACHE PATH "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

set(ENABLE_SPHINX_EXECUTABLE OFF CACHE BOOL "")
set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_CALIPER_HYPRE "$ENV{ENABLE_CALIPER_HYPRE}" CACHE BOOL "")
set(ENABLE_CUDA_NVTOOLSEXT "$ENV{ENABLE_CUDA_NVTOOLSEXT}" CACHE BOOL "")
set(ENABLE_ADIAK ON CACHE BOOL "")

set(ENABLE_CUDA ON CACHE BOOL "")
# set(ENABLE_CUDA ON CACHE PATH "" FORCE)

set(CUDA_TOOLKIT_ROOT_DIR "$ENV{CUDA_HOME}"  CACHE PATH "")
message("## CUDA_TOOLKIT_ROOT_DIR = ${CUDA_TOOLKIT_ROOT_DIR}")
if(ENABLE_CUDA)
  # GEOS_USE_DEVICE and CALC_FEM_SHAPE_IN_KERNEL must be defined and visible

  set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")
  message("## ENABLE_CUDA = ${ENABLE_CUDA}; ENABLE_HYPRE_DEVICE = ${ENABLE_HYPRE_DEVICE}")

  set(ENABLE_CUDA_NVTOOLSEXT "$ENV{ENABLE_CUDA_NVTOOLSEXT}" CACHE BOOL "")
  message("## ENABLE_CUDA_NVTOOLSEXT = ${ENABLE_CUDA_NVTOOLSEXT}")

  # H100, H200 sm_90 
  set(CUDA_ARCH sm_90 CACHE STRING "")
  set(CMAKE_CUDA_ARCHITECTURES 90 CACHE STRING "")

  set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
  set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "") 
  # set(CMAKE_CUDA_STANDARD 14 CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS " $ENV{GEOS_CUDA_FLAGS_CLI} -restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations -Xcompiler -DGEOS_USE_DEVICE  -Xcompiler -DCALC_FEM_SHAPE_IN_KERNEL " CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS_RELEASE " ${CMAKE_CUDA_FLAGS} -O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -Ofast  " CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "") 
  set(CMAKE_CUDA_FLAGS_DEBUG " ${CMAKE_CUDA_FLAGS} -g -G -O0 -Xcompiler -O0  " CACHE STRING "") 

  set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")

  # CHECK A100 
  # -D CUPTI_PREFIX=path-to-your-cupti-installation
  # set(ENABLE_CALIPER_WITH_CUPTI ON CACHE BOOL "")
  
endif()


# Compilers

# --library-path=${CUDA_TOOLKIT_LIB_DIR}

set(CMAKE_C_COMPILER "$ENV{GEOSX_CC}" CACHE PATH "")
set(CMAKE_C_FLAGS_RELEASE " -O3  -Wno-error -pthread  -fno-fast-math -DNDEBUG $ENV{GEOS_HOST_FLAGS_CLI} " CACHE STRING "")
# set(CMAKE_C_FLAGS_RELEASE " -O3 -fast -DNDEBUG -L${CUDA_TOOLKIT_LIB_DIR} " CACHE STRING "")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE} " CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g  $ENV{GEOSX_HOST_FLAGS_CLI}  " CACHE STRING "")

set(CMAKE_CXX_COMPILER "$ENV{GEOSX_CXX}" CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_C_FLAGS_RELEASE} " CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g  $ENV{GEOS_FLAGS_CLI} " CACHE STRING "")
# set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

set(CMAKE_Fortran_COMPILER "$ENV{GEOSX_FORT}" CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE " ${CMAKE_C_FLAGS_RELEASE} " CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g ${CMAKE_C_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g  $ENV{GEOS_FLAGS_CLI} " CACHE STRING "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(BLAS_LIBRARIES "/data/saet/software/aarch64/lib/libopenblas.so" CACHE STRING "")
# set(BLAS_LIBRARIES "$ENV{NVHPC_ROOT}/math_libs/lib64/libnvblas.so" CACHE STRING "")
message("## BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

# CHECK : Most recent version does build. Let's wait for an upgrade on our side.
set(ENABLE_HYPRE ON CACHE BOOL "" FORCE)
set(ENABLE_HYPRE_CUDA ON CACHE BOOL "" FORCE)
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)

set(BLAS_LIBRARIES "/data/saet/software/aarch64/lib/libopenblas.so" CACHE STRING "")
# set(BLAS_LIBRARIES "$ENV{NVHPC_ROOT}/math_libs/lib64/libnvblas.so" CACHE STRING "")
message("## BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE )

if(ENABLE_HYPRE_CUDA)
  set(ENABLE_PETSC OFF CACHE BOOL "") 
  set(ENABLE_TRILINOS OFF CACHE BOOL "") 
  set(GEOSX_LA_INTERFACE "Hypre" CACHE STRING "") 
  # set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
endif()


# Hypre ON

if( (ENABLE_HYPRE_CUDA AND ENABLE_TRILINOS_CUDA) OR (NOT ENABLE_TRILINOS_CUDA AND NOT ENABLE_HYPRE_CUDA) )
  MESSAGE(SEND_ERROR "Exactly one of ENABLE_HYPRE and ENABLE_TRILINOS must be defined.")
  MESSAGE(SEND_ERROR "ENABLE_HYPRE = ${ENABLE_HYPRE}.")
  MESSAGE(SEND_ERROR "ENABLE_TRILINOS = ${ENABLE_TRILINOS}.")
endif()

MESSAGE(STATUS "## GEOSX_LA_INTERFACE = ${GEOSX_LA_INTERFACE}")

# enable PAMELA and PVTPackage
set(ENABLE_PAMELA OFF CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

set(ENABLE_VALGRIND OFF CACHE BOOL "")
set(ENABLE_CALIPER ON CACHE BOOL "")

# disable Doxygen
set(ENABLE_DOXYGEN OFF CACHE PATH "")

# enable tests
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

# let GEOSX define some third party libraries information for you
message("## tpls.cmake = ${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake")
include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
message("## CMAKE_CURRENT_LIST_DIR = ${CMAKE_CURRENT_LIST_DIR}")


# Python options
#set(ENABLE_PYLVARRAY ON CACHE BOOL "")
#set(ENABLE_PYGEOSX ON CACHE BOOL "")
#set(PYTHON_DIR "/share/software/user/open/python/3.6.1" CACHE PATH "")
#set(Python3_EXECUTABLE "/share/software/user/open/python/3.6.1/bin/python3" CACHE PATH "")

## Nvidia HPC TDK 

# set(CMAKE_EXE_LINKER_FLAGS " -L${CUDA_TOOLKIT_LIB_DIR} " CACHE STRING "")
# //Flags used by the linker during DEBUG builds.
# set(CMAKE_EXE_LINKER_FLAGS_DEBUG " ${CMAKE_EXE_LINKER_FLAGS} " CACHE STRING "")
# //Flags used by the linker during MINSIZEREL builds.
# CMAKE_EXE_LINKER_FLAGS_MINSIZEREL:STRING=
# //Flags used by the linker during RELEASE builds.
# set(CMAKE_EXE_LINKER_FLAGS_RELEASE " ${CMAKE_EXE_LINKER_FLAGS}  " CACHE STRING "")
# //Flags used by the linker during RELWITHDEBINFO builds.
# set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO " ${CMAKE_EXE_LINKER_FLAGS}  " CACHE STRING "")
# set(CUDA_TOOLKIT_ROOT_DIR "$ENV{NVHPC_ROOT}/cuda/11.7"  CACHE PATH "")
# set(CUDA_TOOLKIT_LIB_DIR "${CUDA_TOOLKIT_ROOT_DIR}/cuda/11.7/lib64:${CUDA_TOOLKIT_ROOT_DIR}/math_libs/lib64"  CACHE PATH "")
# message("## CUDA_TOOLKIT_LIB_DIR = ${CUDA_TOOLKIT_LIB_DIR}")
# set(LAPACK_LIBRARIES "$ENV{NVHPC_ROOT}/math_libs/lib64/liblapack_static.a" CACHE STRING "")
# message("## LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}")
# set(BLAS_LIBRARIES "$ENV{NVHPC_ROOT}/math_libs/lib64/libnvblas.so" CACHE STRING "")
# message("## BLAS_LIBRARIES = ${BLAS_LIBRARIES}")

# CUDA RT libs at 
# /data/saet/mtml/software/aarch64/nvidia/hpc_sdk/Linux_aarch64/22.7/cuda/11.7/lib64/
#   set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations --library-path=${CUDA_TOOLKIT_LIB_DIR} -L${CUDA_TOOLKIT_LIB_DIR} " CACHE STRING "") 
