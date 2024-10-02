set(COMPILER_HOME "/apps/gcc/8.2.0/x86_64")
set(MPI_HOME "/hrtc/apps/mpi/openmpi/4.0.1/RDHPC/gcc/8.2.0")

set(CMAKE_C_COMPILER ${COMPILER_HOME}/bin/gcc CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER ${COMPILER_HOME}/bin/g++ CACHE PATH "" FORCE)
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(ENABLE_OPENMP ON CACHE PATH "" FORCE)

set(MPI_C_COMPILER "${MPI_HOME}/bin/mpicc" CACHE PATH "" FORCE)
set(MPI_CXX_COMPILER "${MPI_HOME}/bin/mpicxx" CACHE PATH "" FORCE)
set(MPIEXEC_EXECUTABLE "${MPI_HOME}/bin/mpirun" CACHE PATH "" FORCE)
#set(MPIEXEC_NUMPROC_FLAG   "-p pecan -n" CACHE STRING "")
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

set(ENABLE_PETSC OFF CACHE BOOL "" FORCE )
set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE )
set(ENABLE_HYPRE ON CACHE BOOL "" FORCE )

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)
set(ENABLE_CALIPER ON CACHE BOOL "")

set(ENABLE_GEOS_PTP ON CACHE BOOL "" FORCE)

set(ENABLE_MKL ON CACHE BOOL "")
set(INTEL_ROOT "/apps/intel/2019/u5/compilers_and_libraries_2019.5.281/linux" )
set(MKL_ROOT "${INTEL_ROOT}/mkl" )
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  ${INTEL_ROOT}/compiler/lib/intel64_lin/libiomp5.so
                  CACHE STRING "")

set(GEOS_TPL_DIR "$ENV{GEOSX_TPL_DIR}" CACHE PATH "" FORCE)
include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
