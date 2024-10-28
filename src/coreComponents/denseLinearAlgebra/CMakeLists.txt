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
Package: denseLinearAlgebra

Contains dense linear algebra functions and interfaces to BLAS and LAPACK packages.
#]]

#
# Specify all headers
#
set( denseLinearAlgebra_headers
     common/layouts.hpp
     denseLASolvers.hpp
     interfaces/blaslapack/BlasLapackFunctions.h
     interfaces/blaslapack/BlasLapackLA.hpp )

#
# Specify all sources
#
set( denseLinearAlgebra_sources
     interfaces/blaslapack/BlasLapackLA.cpp )

set( dependencyList codingUtilities blas lapack ${parallelDeps} )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       denseLinearAlgebra
                 SOURCES    ${denseLinearAlgebra_sources}
                 HEADERS    ${denseLinearAlgebra_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( denseLinearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

install( TARGETS denseLinearAlgebra LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

