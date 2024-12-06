include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-clang-13-cuda-11.cmake)

# MPI
set(MPI_HOME /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-13.0.1-gcc-8.3.1 CACHE PATH "")
set(MPI_Fortran_COMPILER spectrum-mpi-rolling-release-xl-2023.06.28-cuda-11.8.0/bin/mpifort CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/lassen-base.cmake)

set(ENABLE_CUDA_NVTOOLSEXT ON CACHE BOOL "")

