include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-clang-13-cuda-12.cmake)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran  CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -mtune=power9" CACHE STRING "")
set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")
set(OpenMP_Fortran_FLAGS "-fopenmp" CACHE STRING "")
set(OpenMP_Fortran_LIB_NAMES "" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-13.0.1-gcc-8.3.1/ CACHE PATH "")
set(MPI_Fortran_COMPILER /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpifort CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/lassen-base.cmake)

# Overwrite options set on lassen-base.cmake
set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
set(ENABLE_CUDA_NVTOOLSEXT ON CACHE BOOL "")

# Overwrite ESSL defaults from lassen-base.cmake
# Reason: libesslsmpcuda.so depends on cuda-11
set(ESSL_LIBRARIES ${ESSL_DIR}/lib64/libessl.so
                   ${ESSL_DIR}/lib64/liblapackforessl.so
                   ${ESSL_DIR}/lib64/liblapackforessl_.so
                   CACHE PATH "" FORCE )
