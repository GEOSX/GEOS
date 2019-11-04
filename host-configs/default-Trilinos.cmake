
site_name(HOST_NAME)

set(CMAKE_C_COMPILER "cc" CACHE PATH "" FORCE)
set(CMAKE_CXX_COMPILER "c++" CACHE PATH "" FORCE)
set(ENABLE_FORTRAN ON CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "" FORCE)
set(MPI_C_COMPILER "mpicc" CACHE PATH "" FORCE)
set(MPI_CXX_COMPILER "mpic++" CACHE PATH "" FORCE)
set(MPIEXEC "mpirun" CACHE PATH "" FORCE)

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )

# PETSc doesn't seem to work correctly with clang.
set(ENABLE_PETSC OFF CACHE BOOL "")

