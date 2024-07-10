include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/ruby-clang-14.cmake)

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/llnl-cpu-base.cmake)
