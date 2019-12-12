site_name(HOST_NAME)

set(CMAKE_C_COMPILER "$ENV{CC}" CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER "$ENV{CXX}" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "$ENV{FC}" CACHE PATH "" FORCE)
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER "$ENV{MPICC}" CACHE PATH "" FORCE)
set(MPI_CXX_COMPILER "$ENV{MPICXX}" CACHE PATH "" FORCE)
set(MPI_Fortran_COMPILER "$ENV{MPIFC}" CACHE PATH "" FORCE)
set(MPIEXEC_EXECUTABLE "$ENV{MPIEXEC}" CACHE PATH "" FORCE)

set( ENABLE_PETSC ON CACHE PATH "" )

set(ENABLE_CUDA ON CACHE PATH "" FORCE)
set(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda CACHE PATH "" FORCE)
set(CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE PATH "" FORCE)
set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc CACHE PATH "" FORCE)
set(CUDA_ARCH sm_70 CACHE PATH "" FORCE)
set(CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations" CACHE STRING "")

set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")
set(gtest_disable_pthreads ON CACHE BOOL "")
