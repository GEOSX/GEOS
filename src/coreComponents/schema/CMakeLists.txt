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
Package: schema

Contains XML schema definition and functions to generate it.
#]]

#
# Specify all headers
#
set( schema_headers
     schemaUtilities.hpp )

#
# Specify all sources
#
set( schema_sources
     schemaUtilities.cpp )

set( dependencyList ${parallelDeps} dataRepository )

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       schema
                 SOURCES    ${schema_sources}
                 HEADERS    ${schema_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )


target_include_directories( schema PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

install( TARGETS schema LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
