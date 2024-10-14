include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/dane-gcc-12.cmake)

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic CACHE PATH "")

# ATS
set(ATS_ARGUMENTS "--machine slurm112"  CACHE STRING "")

# This is here to note the required flags for using valgrind. These will have to be propagated to the TPL's
#set( CMAKE_CXX_FLAGS "-march=x86-64-v2 -mno-avx512f" CACHE STRING "" FORCE)

include(${CMAKE_CURRENT_LIST_DIR}/llnl-cpu-base.cmake)
