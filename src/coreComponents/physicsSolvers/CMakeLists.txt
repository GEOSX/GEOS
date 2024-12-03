# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: physicsSolvers

Contains:
 - physics solvers base and manager classes.
 - implementations of different physics models.
 - physics solver wrapper designed for PyGEOSX interface.
#]]

#
# Specify solver headers
#
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     PhysicsSolverBase.hpp
     PhysicsSolverBaseKernels.hpp
     SolverStatistics.hpp
     FieldStatisticsBase.hpp
     LogLevelsInfo.hpp )

#
# Specify solver sources
#
set( physicsSolvers_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverManager.cpp
     PhysicsSolverBase.cpp
     SolverStatistics.cpp )

if( GEOS_ENABLE_CONTACT )
  add_subdirectory( contact )
endif()

if( GEOS_ENABLE_FLUIDFLOW )
  add_subdirectory( fluidFlow )
endif()

if( GEOS_ENABLE_INDUCEDSEISMICITY )
  add_subdirectory( inducedSeismicity )
endif()

if( GEOS_ENABLE_MULTIPHYSICS )
  add_subdirectory( multiphysics )
  include( multiphysics/poromechanicsKernels/PoromechanicsKernels.cmake)
endif()

if( GEOS_ENABLE_SIMPLEPDE )
  add_subdirectory( simplePDE )
endif()

if( GEOS_ENABLE_SOLIDMECHANICS )
  add_subdirectory( solidMechanics )
  include( solidMechanics/kernels/SolidMechanicsKernels.cmake)
endif()

if( GEOS_ENABLE_SURFACEGENERATION )
  add_subdirectory( surfaceGeneration )
endif()

if( GEOS_ENABLE_WAVEPROPAGATION )
  add_subdirectory( wavePropagation )
endif()

set( dependencyList ${parallelDeps} fileIO discretizationMethods events linearAlgebra )

if( ENABLE_PYGEOSX )
  list( APPEND physicsSolvers_headers
  python/PySolver.hpp
  python/PySolverType.hpp )
  list( APPEND physicsSolvers_sources
        python/PySolver.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       physicsSolvers
                 SOURCES    ${physicsSolvers_sources}
                 HEADERS    ${physicsSolvers_headers}
                 DEPENDS_ON ${decoratedDependencies} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS physicsSolvers LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()


