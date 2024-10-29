message( "\nProcessing GeosxOptions.cmake" )
message( "CMAKE_HOST_SYSTEM_NAME = ${CMAKE_HOST_SYSTEM_NAME}" )
message( "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}" )
message( "CMAKE_HOST_APPLE = ${CMAKE_HOST_APPLE}" )

### OPTIONS ###
option( GEOS_ENABLE_TESTS "" ON )
option( ENABLE_CALIPER "Enables Caliper instrumentation" OFF )

option( ENABLE_MATHPRESSO "" ON )

option( ENABLE_CHAI "Enables CHAI" ON )
option( BUILD_LOCAL_CHAI "Use the local mirrored CHAI" OFF )

option( ENABLE_RAJA "Enables RAJA" ON )
option( BUILD_LOCAL_RAJA "Use the local mirrored RAJA" OFF )
option( RAJA_ENABLE_TBB "" OFF )
option( RAJA_ENABLE_OPENMP "" OFF )
option( RAJA_ENABLE_CUDA "" OFF )
option( RAJA_ENABLE_HIP "" OFF )
option( RAJA_ENABLE_TESTS "" OFF )

option( ENABLE_PVTPackage "" ON )

option( ENABLE_UNCRUSTIFY "" ON )

option( ENABLE_XML_UPDATES "" ON )

option( ENABLE_FORTRAN "Enables Fortran support" OFF )

option( ENABLE_METIS "Enables METIS" ON )
option( ENABLE_PARMETIS "Enables PARMETIS" ON )
option( ENABLE_SCOTCH "Enables SCOTCH" ON )

option( ENABLE_SILO "Enables SILO output" ON )
option( ENABLE_VTK "Enables VTK" ON )

option( ENABLE_TOTALVIEW_OUTPUT "Enables Totalview custom view" OFF )

option( ENABLE_SUPERLU_DIST "Enables SUPERLU_DIST" ON )
option( ENABLE_TRILINOS "Enables TRILINOS" ON )
option( ENABLE_HYPRE "Enables HYPRE" ON )
option( ENABLE_PETSC "Enables PETSC" OFF )
option( ENABLE_SUITESPARSE "Enables SUITESPARSE" ON )

option( ENABLE_HYPRE_MIXINT "Enables mixed int32/int64 local/global" ON )

set( HYPRE_DEVICE_OPTIONS CPU CUDA HIP )
if( NOT ENABLE_HYPRE_DEVICE )
  set( ENABLE_HYPRE_DEVICE CPU )
endif()
if(NOT ${ENABLE_HYPRE_DEVICE} IN_LIST HYPRE_DEVICE_OPTIONS )
    message(FATAL_ERROR "Set ENABLE_HYPRE_DEVICE to CPU, CUDA, or HIP.")
endif()

#if ( "${CMAKE_HOST_APPLE}" )
#  option( ENABLE_PETSC "Enables PETSC" OFF )
#else()
#  option( ENABLE_PETSC "Enables PETSC" ON )
#endif()

### LAI SETUP ###

set( supported_LAI Trilinos Hypre Petsc )
set( GEOS_LA_INTERFACE "Hypre" CACHE STRING "Linear algebra interface to use in solvers" )
message( STATUS "GEOS_LA_INTERFACE = ${GEOS_LA_INTERFACE}" )

if( NOT ( GEOS_LA_INTERFACE IN_LIST supported_LAI ) )
  message( FATAL_ERROR "GEOS_LA_INTERFACE must be one of: ${supported_LAI}" )
endif()

### MPI/OMP/CUDA/HIP SETUP ###

option( ENABLE_MPI "" ON )

option( ENABLE_CUDA "" OFF )

option( ENABLE_HIP "" OFF )

if( CMAKE_HOST_APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang" )
  option( ENABLE_OPENMP "Enables OpenMP compiler support" OFF )
else()
  option( ENABLE_OPENMP "Enables OpenMP compiler support" ON )
endif()

### BUILD & BLT SETUP ###

option( GEOS_INSTALL_SCHEMA "Enables schema generation and installation" ON )

option( GEOS_BUILD_OBJ_LIBS "Builds coreComponent modules as object libraries" OFF )

option( GEOS_BUILD_SHARED_LIBS "Builds geosx_core as a shared library " ON )

set( GEOS_PARALLEL_COMPILE_JOBS "" CACHE STRING "Maximum number of concurrent compilation jobs" )
if( GEOS_PARALLEL_COMPILE_JOBS )
    set_property( GLOBAL APPEND PROPERTY JOB_POOLS compile_job_pool=${GEOS_PARALLEL_COMPILE_JOBS} )
    set( CMAKE_JOB_POOL_COMPILE compile_job_pool )
endif()

set( GEOS_PARALLEL_LINK_JOBS "" CACHE STRING "Maximum number of concurrent link jobs" )
if( GEOS_PARALLEL_LINK_JOBS )
    set_property( GLOBAL APPEND PROPERTY JOB_POOLS link_job_pool=${GEOS_PARALLEL_LINK_JOBS} )
    set( CMAKE_JOB_POOL_LINK link_job_pool )
endif()

#set(CMAKE_POSITION_INDEPENDENT_CODE ON  CACHE BOOL "" FORCE)
#blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -rdynamic)
#set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")

#set(SPHINX_EXECUTABLE "sphinx-build" CACHE PATH "")
#include(cmake/blt/cmake/thirdparty/FindSphinx.cmake)
#message( "SPHINX_FOUND = ${SPHINX_FOUND}" )
#message( "SPHINX_EXECUTABLE = ${SPHINX_EXECUTABLE}" )

if( NOT BLT_CXX_STD STREQUAL c++17 )
    MESSAGE( FATAL_ERROR "c++17 is NOT enabled. GEOSX requires c++17" )
endif( NOT BLT_CXX_STD STREQUAL c++17 )

message( "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}" )

blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT "${OpenMP_CXX_FLAGS}" )
blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS
                                 GNU   "-Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wcast-align -Wcast-qual"
                                 CLANG "-Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wno-cast-align -Wcast-qual"
                               )

blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS_DEBUG
                                 GNU "-Wno-unused-parameter -Wno-unused-variable -Wno-dangling-reference"
                                 CLANG "-Wno-unused-parameter -Wno-unused-variable -fstandalone-debug"
                               )

blt_append_custom_compiler_flag( FLAGS_VAR GEOS_NINJA_FLAGS
                                 DEFAULT " "
                                 GNU     "-fdiagnostics-color=always"
                                 CLANG   "-fcolor-diagnostics"
                               )

# clang-13 and gcc complains about unused-but-set variable.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-Wunused-but-set-variable" CXX_UNUSED_BUT_SET_VAR)
if (ENABLE_GBENCHMARK)
    blt_add_target_compile_flags(TO benchmark
                                FLAGS $<$<AND:$<BOOL:${CXX_UNUSED_BUT_SET_VAR}>,$<COMPILE_LANGUAGE:CXX>>:-Wno-unused-but-set-variable>
                                )
endif()

if( ${CMAKE_MAKE_PROGRAM} STREQUAL "ninja" OR ${CMAKE_MAKE_PROGRAM} MATCHES ".*/ninja$" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GEOS_NINJA_FLAGS}" )
endif()


if( CMAKE_HOST_APPLE )
#    set(GEOS_LINK_PREPEND_FLAG "-Wl,-force_load" CACHE STRING "")
#    set(GEOS_LINK_POSTPEND_FLAG "" CACHE STRING "")
# elseif( ENABLE_CUDA )
#     set( GEOS_LINK_PREPEND_FLAG  "-Xcompiler \\\\\"-Wl,--whole-archive\\\\\""    CACHE STRING "" )
#     set( GEOS_LINK_POSTPEND_FLAG "-Xcompiler \\\\\"-Wl,--no-whole-archive\\\\\"" CACHE STRING "" )
else()
    set( GEOS_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE STRING "" )
    set( GEOS_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE STRING "" )
endif()

set( GEOS_LOCALINDEX_TYPE "int" CACHE STRING "" )
if( ENABLE_HYPRE_MIXINT )
  set( GEOS_GLOBALINDEX_TYPE "long long int" CACHE STRING "" )
else()
  set( GEOS_GLOBALINDEX_TYPE "int" CACHE STRING "" )
endif()

if( GEOS_LOCALINDEX_TYPE STREQUAL "int" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "0" CACHE STRING "" FORCE )
elseif( GEOS_LOCALINDEX_TYPE STREQUAL "long int" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "1" CACHE STRING "" FORCE )
elseif( GEOS_LOCALINDEX_TYPE STREQUAL "long long int" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "2" CACHE STRING "" FORCE )
elseif( GEOS_LOCALINDEX_TYPE STREQUAL "std::ptrdiff_t" )
    set( GEOS_LOCALINDEX_TYPE_FLAG "3" CACHE STRING "" FORCE )
else( TRUE )
    message( FATAL_ERROR "GEOS_LOCALINDEX_TYPE_FLAG not set for ${GEOS_LOCALINDEX_TYPE}" )
endif()



if( GEOS_GLOBALINDEX_TYPE STREQUAL "int" )
    set( GEOS_GLOBALINDEX_TYPE_FLAG "0" CACHE STRING "" FORCE )
elseif( GEOS_GLOBALINDEX_TYPE STREQUAL "long int" )
    set( GEOS_GLOBALINDEX_TYPE_FLAG "1" CACHE STRING "" FORCE )
elseif( GEOS_GLOBALINDEX_TYPE STREQUAL "long long int" )
    set( GEOS_GLOBALINDEX_TYPE_FLAG "2" CACHE STRING "" FORCE )
else( TRUE )
    message( FATAL_ERROR "GEOS_GLOBALINDEX_TYPE_FLAG not set for ${GEOS_GLOBALINDEX_TYPE}" )
endif()

set( GEOS_BLOCK_SIZE 32 )
if( ENABLE_CUDA )
  set( GEOS_BLOCK_SIZE 32 )
endif()
if( ENABLE_HIP )
  set( GEOS_BLOCK_SIZE 64 )
endif()

message( "localIndex is an alias for ${GEOS_LOCALINDEX_TYPE}" )
message( "globalIndex is an alias for ${GEOS_GLOBALINDEX_TYPE}" )
message( "GEOS_LOCALINDEX_TYPE_FLAG = ${GEOS_LOCALINDEX_TYPE_FLAG}" )
message( "GEOS_GLOBALINDEX_TYPE_FLAG = ${GEOS_GLOBALINDEX_TYPE_FLAG}" )


message( "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}" )
message( "GEOS_LINK_PREPEND_FLAG=${GEOS_LINK_PREPEND_FLAG}" )
message( "GEOS_LINK_POSTPEND_FLAG=${GEOS_LINK_POSTPEND_FLAG}" )
message( "Leaving GeosxOptions.cmake\n" )
