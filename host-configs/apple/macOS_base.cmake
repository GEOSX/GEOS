message( "this hostconfig assumes you are using homebrew")
message( "brew install bison cmake gcc git-lfs open-mpi openblas python3")

message( "CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}" )
message("CONFIG_NAME = ${CONFIG_NAME}")

set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "${HOMEBREW_DIR}/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "${HOMEBREW_DIR}/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "${HOMEBREW_DIR}/bin/mpirun" CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)

set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

set(ENABLE_CUDA "OFF" CACHE PATH "" FORCE)
set(ENABLE_OPENMP "OFF" CACHE PATH "" FORCE)

set(ENABLE_CALIPER "OFF" CACHE PATH "" FORCE )

set( BLAS_LIBRARIES ${HOMEBREW_DIR}/opt/lapack/lib/libblas.dylib CACHE PATH "" FORCE )
set( LAPACK_LIBRARIES ${HOMEBREW_DIR}/opt/lapack/lib/liblapack.dylib CACHE PATH "" FORCE )

set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)
set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE )
#set(GEOSX_BUILD_OBJ_LIBS ON CACHE BOOL "" FORCE)



#set( DOXYGEN_EXECUTABLE /usr/local/bin/doxygen CACHE PATH "" FORCE )
#set( SPHINX_EXECUTABLE /usr/local/bin/sphinx-build CACHE PATH "" FORCE )

if(NOT ( EXISTS "${GEOSX_TPL_DIR}" AND IS_DIRECTORY "${GEOSX_TPL_DIR}" ) )
   set(GEOSX_TPL_DIR "${CMAKE_SOURCE_DIR}/../../thirdPartyLibs/install-${CONFIG_NAME}-release" CACHE PATH "" FORCE )
endif()

# ATS
set(ATS_ARGUMENTS "--machine openmpi --ats openmpi_mpirun=${MPIEXEC}"  CACHE PATH "")

include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
