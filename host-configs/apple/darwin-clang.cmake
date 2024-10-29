site_name(HOST_NAME)
set(CONFIG_NAME "${HOST_NAME}-darwin-x86_64-clang@apple-mp" CACHE PATH "")
message("CONFIG_NAME = ${CONFIG_NAME}")

set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "/usr/local/bin/mpirun" CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)

set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

set(ENABLE_CUDA "OFF" CACHE PATH "" FORCE)
set(ENABLE_OPENMP "OFF" CACHE PATH "" FORCE)

set(ENABLE_CALIPER "OFF" CACHE PATH "" FORCE )

#set(GEOS_BUILD_OBJ_LIBS ON CACHE BOOL "" FORCE)

set( BLAS_LIBRARIES /usr/local/opt/openblas/lib/libblas.dylib CACHE PATH "" FORCE )
set( LAPACK_LIBRARIES /usr/local/opt/openblas/lib/liblapack.dylib CACHE PATH "" FORCE )

set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)

#set( DOXYGEN_EXECUTABLE /usr/local/bin/doxygen CACHE PATH "" FORCE )
#set( SPHINX_EXECUTABLE /usr/local/bin/sphinx-build CACHE PATH "" FORCE )

set(GEOS_TPL_DIR "/usr/local/GEOSX/GEOS_TPL" CACHE PATH "" FORCE )
if(NOT ( EXISTS "${GEOS_TPL_DIR}" AND IS_DIRECTORY "${GEOS_TPL_DIR}" ) )
    set(GEOS_TPL_DIR "${CMAKE_SOURCE_DIR}/../../thirdPartyLibs/install-darwin-clang-release" CACHE PATH "" FORCE )
endif()

include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
