#
# Performance portability
#
message("in tpls.cmake GEOSX_TPL_DIR=${GEOSX_TPL_DIR}")

#
# General TPL Folder verifications
#
if(NOT EXISTS ${GEOSX_TPL_DIR})
  message(WARNING "'GEOSX_TPL_DIR' does not exist.\n")
endif()


if(EXISTS ${GEOSX_TPL_DIR}/raja)
  set(RAJA_DIR ${GEOSX_TPL_DIR}/raja CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/chai)
  set(UMPIRE_DIR ${GEOSX_TPL_DIR}/chai CACHE PATH "" FORCE)
  set(CHAI_DIR ${GEOSX_TPL_DIR}/chai CACHE PATH "" FORCE)
endif()

#
# IO TPLs
#
if(EXISTS ${GEOSX_TPL_DIR}/hdf5)
  set(HDF5_DIR ${GEOSX_TPL_DIR}/hdf5 CACHE PATH "" FORCE)
  message(STATUS "HDF5_DIR = ${HDF5_DIR}")
endif()

if(EXISTS ${GEOSX_TPL_DIR}/conduit)
  set(CONDUIT_DIR ${GEOSX_TPL_DIR}/conduit CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/silo)
  set(SILO_DIR ${GEOSX_TPL_DIR}/silo CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/adiak)
  set(ADIAK_DIR ${GEOSX_TPL_DIR}/adiak CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/caliper)
  set(CALIPER_DIR ${GEOSX_TPL_DIR}/caliper CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/pugixml)
  set(PUGIXML_DIR ${GEOSX_TPL_DIR}/pugixml CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/vtk)
  set(VTK_DIR ${GEOSX_TPL_DIR}/vtk CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/fmt)
  set(FMT_DIR ${GEOSX_TPL_DIR}/fmt CACHE PATH "" FORCE)
endif()

#
# Math TPLs
#
if(EXISTS ${GEOSX_TPL_DIR}/metis)
  set(METIS_DIR ${GEOSX_TPL_DIR}/metis CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/parmetis)
  set(PARMETIS_DIR ${GEOSX_TPL_DIR}/parmetis CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/superlu_dist)
  set(SUPERLU_DIST_DIR ${GEOSX_TPL_DIR}/superlu_dist CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/suitesparse)
  set(SUITESPARSE_DIR ${GEOSX_TPL_DIR}/suitesparse CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/trilinos)
  set(TRILINOS_DIR ${GEOSX_TPL_DIR}/trilinos CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/hypre)
  set(HYPRE_DIR ${GEOSX_TPL_DIR}/hypre CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/scotch)
  set(SCOTCH_DIR ${GEOSX_TPL_DIR}/scotch CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/petsc AND (NOT DEFINED ENABLE_PETSC OR ENABLE_PETSC))
  set(PETSC_DIR ${GEOSX_TPL_DIR}/petsc CACHE PATH "" FORCE)
endif()

#
# Development tools
#
if(EXISTS ${GEOSX_TPL_DIR}/uncrustify/bin/uncrustify)
  set(UNCRUSTIFY_EXECUTABLE ${GEOSX_TPL_DIR}/uncrustify/bin/uncrustify CACHE PATH "" FORCE)
endif()

if(EXISTS ${GEOSX_TPL_DIR}/doxygen/bin/doxygen)
  set(DOXYGEN_EXECUTABLE ${GEOSX_TPL_DIR}/doxygen/bin/doxygen CACHE PATH "" FORCE)
endif()

#
# Other
#
if(EXISTS ${GEOSX_TPL_DIR}/mathpresso)
  set(MATHPRESSO_DIR ${GEOSX_TPL_DIR}/mathpresso CACHE PATH "" FORCE)
endif()
