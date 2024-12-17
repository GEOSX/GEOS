include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/ruby-gcc-12.cmake)

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic CACHE PATH "")

# ATS
set(ATS_ARGUMENTS "--machine slurm56"  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/llnl-cpu-base.cmake)
