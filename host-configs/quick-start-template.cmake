set( CONFIG_NAME "quick-start" ) 

# Set compilers path
set(CMAKE_C_COMPILER "path-to-gcc/bin/gcc" CACHE PATH "")   # This is typically something like /usr/bin/gcc ... or clang
set(CMAKE_CXX_COMPILER "path-to-gcc/bin/g++" CACHE PATH "") # This is typically something like /usr/bin/g++ ... or clang++
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# Set paths to mpi
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "path-to-mpi/bin/mpicc" CACHE PATH "")    # This is typically something like /usr/bin/mpicc
set(MPI_CXX_COMPILER "path-to-mpi/bin/mpicxx" CACHE PATH "") # This is typically something like /usr/bin/mpicxx
set(MPIEXEC "path-to-mpi/bin/mpirun" CACHE PATH "")          # This is typically something like /usr/bin/mpirun

# Set paths to blas and lapack
set( BLAS_LIBRARIES "path-to-blas" CACHE PATH "" FORCE )   # This is typically something like /usr/lib64/libblas.so 
set( LAPACK_LIBRARIES  CACHE PATH "path-to-lapack" FORCE ) # This is typically something like /usr/lib64/liblapack.so

# Cuda and openMP
set( ENABLE_CUDA OFF CACHE PATH "" FORCE )
set( ENABLE_OPENMP OFF CACHE PATH "" FORCE )

# TPLs
set( ENABLE_TRILINOS OFF CACHE PATH "" FORCE )
set( ENABLE_CALIPER OFF CACHE PATH "" FORCE )
set( ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
set( ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE )

if(NOT ( EXISTS "${GEOS_TPL_DIR}" AND IS_DIRECTORY "${GEOS_TPL_DIR}" ) )
   set(GEOS_TPL_DIR "${CMAKE_SOURCE_DIR}/../../thirdPartyLibs/install-${CONFIG_NAME}-release" CACHE PATH "" FORCE )
endif()

include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
