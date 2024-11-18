###############################################################################
#
# Base configuration for LC Dane builds
# Calling configuration file must define the following CMAKE variables:
#
# MPI_HOME
#
###############################################################################

# Fortran
set(ENABLE_FORTRAN OFF CACHE BOOL "")

# MPI
set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER ${MPI_HOME}/bin/mpicc CACHE PATH "")
set(MPI_CXX_COMPILER  ${MPI_HOME}/bin/mpicxx CACHE PATH "")
set(MPI_Fortran_COMPILER ${MPI_HOME}/bin/mpifort CACHE PATH "")
set(MPIEXEC /usr/bin/srun CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

# PAPI (For TPL caliper builds)
set(ENABLE_PAPI ON CACHE BOOL "")
set(PAPI_PREFIX /usr/tce/packages/papi/papi-5.4.3 CACHE PATH "")

# OpenMP
set(ENABLE_OPENMP ON CACHE BOOL "")

# GEOSX specific options
set(ENABLE_PVTPackage ON CACHE BOOL "")
set(ENABLE_PETSC OFF CACHE BOOL "Enables PETSc." FORCE)

# PYGEOSX
set(ENABLE_PYGEOSX ON CACHE BOOL "")
set(Python3_ROOT_DIR /usr/gapps/GEOSX/thirdPartyLibs/python/dane-gcc-python/python CACHE PATH "")
set(Python3_EXECUTABLE ${Python3_ROOT_DIR}/bin/python3 CACHE PATH "")

# YAPF python formatting
set(YAPF_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/dane-gcc-python/python/bin/yapf CACHE PATH "" FORCE)

# Sphinx
set(SPHINX_EXECUTABLE /usr/gapps/GEOSX/thirdPartyLibs/python/dane-gcc-python/python/bin/sphinx-build CACHE PATH "" FORCE)


set(ENABLE_FESAPI OFF CACHE BOOL "" FORCE)

# caliper 
set(ENABLE_CALIPER ON CACHE BOOL "" FORCE)
set(ENABLE_CALIPER_HYPRE ON CACHE BOOL "" FORCE)

# MKL
set(ENABLE_MKL ON CACHE BOOL "")
set(MKL_ROOT /usr/tce/packages/mkl/mkl-2022.1.0)
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "")
set(MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                  ${MKL_ROOT}/lib/intel64/libmkl_gnu_thread.so
                  ${MKL_ROOT}/lib/intel64/libmkl_core.so
                  CACHE STRING "")

# ATS
set(ATS_ARGUMENTS "--machine slurm36"  CACHE STRING "")
# set(USER $ENV{USER} CACHE STRING "")
# set(ATS_WORKING_DIR "/p/lustre2/${USER}/integratedTests/${CONFIG_NAME}"  CACHE PATH "")
# set(ATS_BASELINE_DIR "/p/lustre2/${USER}/integratedTests/baselines"  CACHE PATH "")

# Temporary argument for python module change testing
# set(GEOS_PYTHON_PACKAGES_BRANCH "feature/sherman/outOfPlaceATS" CACHE STRING "" FORCE)


include(${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake)
