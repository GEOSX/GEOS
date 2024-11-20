site_name(HOST_NAME)

# Compilers
set(CMAKE_C_COMPILER       "${GCC_ROOT}/bin/gcc"      CACHE PATH "")
set(CMAKE_CXX_COMPILER     "${GCC_ROOT}/bin/g++"      CACHE PATH "")
set(CMAKE_Fortran_COMPILER "${GCC_ROOT}/bin/gfortran" CACHE PATH "")

# OpenMP options
set(ENABLE_OPENMP ON CACHE BOOL "")

# MPI options
set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER       "${MPI_ROOT}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_ROOT}/bin/mpic++"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_ROOT}/bin/mpifort" CACHE PATH "")
set(MPIEXEC              "${MPI_ROOT}/bin/mpirun"  CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# CUDA options
if(ENABLE_CUDA)
  set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "")
  set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
  set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
  set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")
endif()

# Blas/Lapack options
set(BLAS_LIBRARIES "/share/software/user/open/openblas/0.3.10/lib/libblas.so" CACHE STRING "")
set(LAPACK_LIBRARIES "/share/software/user/open/openblas/0.3.10/lib/liblapack.so" CACHE STRING "")

# Python options
#set(ENABLE_PYLVARRAY ON CACHE BOOL "")
#set(ENABLE_PYGEOSX ON CACHE BOOL "")
#set(PYTHON_DIR "/share/software/user/open/python/3.6.1" CACHE PATH "")
#set(Python3_EXECUTABLE "/share/software/user/open/python/3.6.1/bin/python3" CACHE PATH "")

set(ENABLE_VALGRIND OFF CACHE BOOL "")
set(ENABLE_CALIPER ON CACHE BOOL "")

if( ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA" )
  set(ENABLE_PETSC OFF CACHE BOOL "")
  set(ENABLE_TRILINOS OFF CACHE BOOL "")
  set(GEOS_LA_INTERFACE "Hypre" CACHE STRING "")
endif()

set(GEOS_TPL_DIR /home/groups/tchelepi/geosx/thirdPartyLibs/install-${CONFIG_NAME}-release CACHE PATH "")
include(/home/groups/tchelepi/geosx/GEOSX/host-configs/tpls.cmake)
