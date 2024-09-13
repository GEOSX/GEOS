include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/ascent-gcc@8.1.1.cmake)

# asmjit doesn't work on PowerPC
set(ENABLE_MATHPRESSO OFF CACHE BOOL "")

# Silo configure script doesn't recognize systype
set(SILO_BUILD_TYPE powerpc64-unknown-linux-gnu CACHE STRING "")

set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_PETSC OFF CACHE BOOL "" FORCE)

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_PAPI OFF CACHE BOOL "")

set(ENABLE_ESSL ON CACHE BOOL "")
set(ESSL_INCLUDE_DIRS /sw/ascent/essl/6.2.0-20190419/essl/6.2/include CACHE STRING "")
set(ESSL_LIBRARIES /sw/ascent/essl/6.2.0-20190419/essl/6.2/lib64/libesslsmpcuda.so
                   /sw/ascent/xl/16.1.1-3/lib/libxlsmp.so.1
                   /sw/ascent/xl/16.1.1-3/lib/libxlfmath.so.1
                   /sw/ascent/xl/16.1.1-3/lib/libxlf90_r.so.1
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas.so
                   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so
		   ${GEOS_TPL_ROOT_DIR}/liblapackforessl.a
                   /sw/ascent/xl/16.1.1-3/xlC/16.1.1/lib/libxl.a
                   CACHE PATH "")

#set(PETSC_OMP_DIR /ccsopen/home/klevtsov/thirdPartyLibs/omp-links-for-petsc CACHE STRING "")

