include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/lassen-gcc-8-cuda-11.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3" CACHE STRING "" FORCE)

# MPI
set(MPI_HOME /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1 CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/lassen-base.cmake)
