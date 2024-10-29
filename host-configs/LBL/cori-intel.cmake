##################################
# cori - KNL build
# Run scripts/cori-setup-environment.sh before building
##################################

set( ENABLE_WARNINGS_AS_ERRORS "OFF" CACHE BOOL "" FORCE)

set(CMAKE_CXX_COMPILER "CC" CACHE PATH "" FORCE)
set(CMAKE_C_COMPILER "cc" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "ftn" CACHE PATH "" FORCE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -axMIC-AVX512,AVX -qoverride-limits" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -axMIC-AVX512,AVX -qoverride-limits" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -axMIC-AVX512,AVX" CACHE STRING "" FORCE)

# NERSC recommends no optimization flags for intel or cray compilers.
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG" CACHE STRING "" FORCE)

set(ENABLE_MPI ON CACHE BOOL "" FORCE)
set(MPI_CXX_COMPILER "${CC}" CACHE PATH "")
set(MPI_C_COMPILER "${cc}" CACHE PATH "")
set(MPI_Fortran_COMPILER "ftn" CACHE PATH "" FORCE)
set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(GEOS_TPL_DIR "/global/project/projectdirs/m1411/GEOSX/tpls/install-cori-intel-release-22-07-20" CACHE PATH "" )

set(GEOS_LINK_PREPEND_FLAG  "-Wl,--whole-archive"    CACHE STRING "" FORCE)
set(GEOS_LINK_POSTPEND_FLAG "-Wl,--no-whole-archive" CACHE STRING "" FORCE)

set(ENABLE_SPHINX_EXECUTABLE OFF CACHE BOOL "")
set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")
set(ENABLE_XML_UPDATES OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_ADIAK ON CACHE BOOL "")
set(ENABLE_PAPI ON CACHE BOOL "")
set(PAPI_PREFIX "/opt/cray/pe/papi/5.7.0.1" CACHE PATH "" FORCE)

set(ENABLE_OPENMP ON CACHE BOOL "")
set(ENABLE_CUDA OFF CACHE BOOL "")

set(ENABLE_TOTALVIEW_OUTPUT OFF CACHE BOOL "Enables Totalview custom view" FORCE)

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

set(ENABLE_PETSC OFF CACHE BOOL "")
