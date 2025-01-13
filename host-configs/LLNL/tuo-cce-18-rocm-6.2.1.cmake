include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/tuo-cce-18-rocm-6.2.1.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/amdgpu-base.cmake)

# MPI
set(MPI_HOME /opt/cray/pe/mpich/8.1.31/ofi/crayclang/18.0 CACHE PATH "")

# Workaround for linking errors involving missing explicit instantiations of certain functions.
# TODO: we should add the missing explicit instantiations
set( CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -Wl,--allow-shlib-undefined" CACHE STRING "" FORCE )
