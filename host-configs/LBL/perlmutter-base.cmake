###############################################################################
#
# Base configuration for LC Quartz builds
# Calling configuration file must define the following CMAKE variables:
#
# MPI_HOME
#
###############################################################################

# Fortran
set(ENABLE_FORTRAN OFF CACHE BOOL "")

# MPI
set(ENABLE_MPI ON CACHE BOOL "")
find_package(MPI REQUIRED)
#set(MPI_C_COMPILER ${MPI_HOME}/bin/mpicc CACHE PATH "")
#set(MPI_CXX_COMPILER ${MPI_HOME}/bin/mpicxx CACHE PATH "")
#set(MPIEXEC srun CACHE STRING "")
#set(MPIEXEC_NUMPROC_FLAG -N CACHE STRING "")
#set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

# OpenMP
set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)

# CUDA
# LvArray sets this to the CMAKE_CXX_COMPILER.
#set(ENABLE_CUDA ON CACHE BOOL "")
#set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "")
set(ENABLE_CUDA ON CACHE BOOL "")
#set(CUDA_TOOLKIT_ROOT_DIR /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/cuda/11.7 CACHE STRING "")
#set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "")
#set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE STRING "")
set(CUDA_ARCH sm_80 CACHE STRING "")
set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr -Werror cross-execution-space-call,reorder,deprecated-declarations" CACHE STRING "")
set(BLT_CUDA_FLAGS "${BLT_CUDA_FLAGS} ${CUDA_NVCC_FLAGS}" CACHE STRING "")

# ESSL
#set(ENABLE_ESSL ON CACHE BOOL "")
#set(ESSL_INCLUDE_DIRS /usr/tcetmp/packages/essl/essl-6.2.1/include CACHE STRING "")
#set(ESSL_INCLUDE_DIRS /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/targets/x86_64-linux/include CACHE STRING "")
#set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.2.1/lib64/libesslsmpcuda.so
#                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlsmp.so
#                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlfmath.so
#                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxlf90_r.so
#set(ESSL_LIBRARIES /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/targets/x86_64-linux/lib/libcublas.so
#                   /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/targets/x86_64-linux/lib/libcusparse.so
#                   /opt/nvidia/hpc_sdk/Linux_x86_64/22.5/math_libs/11.7/targets/x86_64-linux/lib/libcurand.so
#                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so
#                   /usr/tcetmp/packages/essl/essl-6.2.1/lib64/liblapackforessl.so
#                   /usr/tcetmp/packages/essl/essl-6.2.1/lib64/liblapackforessl_.so
#                   /usr/tce/packages/xl/xl-beta-2019.06.20/alllibs/libxl.a
#                   CACHE PATH "")

# TPL
set(ENABLE_PAPI OFF CACHE BOOL "")
#set(SILO_BUILD_TYPE powerpc64-unknown-linux-gnu CACHE STRING "")

# GEOSX specific options
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_PETSC OFF CACHE BOOL "" FORCE )


set( CUDA_cusparse_LIBRARY "/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/11.7/lib64/libcusparse.so " CACHE PATH "" FORCE )
set( CUDA_cublas_LIBRARY "/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/11.7/lib64/libcublas.so" CACHE PATH "" FORCE )
set( CUDA_curand_LIBRARY "/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/11.7/lib64/libcurand.so" CACHE PATH "" FORCE )

set(ENABLE_HYPRE ON CACHE BOOL "" FORCE )
set( ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "" FORCE )
if( ${ENABLE_HYPRE_DEVICE} STREQUAL "HIP" OR ${ENABLE_HYPRE_DEVICE} STREQUAL "CUDA" )
    set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE )
else()
    set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE )
    set(GEOSX_LA_INTERFACE "Trilinos" CACHE STRING "" FORCE )
endif()


# Documentation
set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "" FORCE)
set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)

# Other
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
