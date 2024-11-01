##################################
# cori - KNL build
# Run scripts/cori-setup-environment.sh before building
##################################
set(CMAKE_CXX_COMPILER "/opt/gcc/8.1.0/bin/g++" CACHE PATH "" FORCE)
set(CMAKE_C_COMPILER "/opt/gcc/8.1.0/bin/gcc" CACHE PATH "" FORCE)
set(CMAKE_Fortran_COMPILER "/opt/gcc/8.1.0/bin/gfortran" CACHE PATH "" FORCE)

# These flags should be set for the TPL build but they break some of our numerically sensative code.
# This needs to be fixed.
set(CMAKE_C_FLAGS "-march=knl -mtune=knl -Wl,-rpath=/opt/gcc/8.1.0/snos/lib64/" CACHE STRING "")
set(CMAKE_CXX_FLAGS "-march=knl -mtune=knl -Wl,-rpath=/opt/gcc/8.1.0/snos/lib64/" CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-march=knl -mtune=knl -Wl,-rpath=/opt/gcc/8.1.0/snos/lib64/" CACHE STRING "")
# set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath=/opt/gcc/8.1.0/snos/lib64/" CACHE STRING "")

set(ENABLE_FORTRAN OFF CACHE BOOL "")

set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_HOME /global/common/software/m3169/cori/openmpi/4.0.2/gnu/ CACHE PATH "")
set(MPI_C_COMPILER       ${MPI_HOME}/bin/mpicc   CACHE PATH "")
set(MPI_CXX_COMPILER     ${MPI_HOME}/bin/mpicxx  CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")

set(MPIEXEC              /usr/bin/srun CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(GEOS_TPL_DIR "/global/project/projectdirs/m1411/GEOSX/tpls/install-cori-gcc\@8.1.0-release-24-07-20" CACHE PATH "" )

set(ENABLE_SPHINX_EXECUTABLE OFF CACHE BOOL "")
set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(ENABLE_PVTPackage ON CACHE BOOL "")

set(ENABLE_CALIPER ON CACHE BOOL "")
set(ENABLE_ADIAK ON CACHE BOOL "")
set(ENABLE_PAPI ON CACHE BOOL "")
set(PAPI_PREFIX "/opt/cray/pe/papi/5.7.0.1" CACHE PATH "" FORCE)

set(ENABLE_OPENMP ON CACHE BOOL "")
set(ENABLE_CUDA OFF CACHE BOOL "")

set(ENABLE_TOTALVIEW_OUTPUT OFF CACHE BOOL "Enables Totalview custom view" FORCE)
set(ENABLE_PETSC OFF CACHE BOOL "")

set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /opt/intel/compilers_and_libraries_2019.3.199/linux/mkl)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

