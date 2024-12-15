include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/tioga-cce-18-rocm-6.2.1.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/tioga-base.cmake)

# MPI
set(MPI_HOME /opt/cray/pe/mpich/8.1.31/ofi/crayclang/18.0 CACHE PATH "")
