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
Package : functions

Contains classes for storing and computing arbitrary N-dimensional functions.
#]]

#
# Specify all headers
#
set( functions_headers
     FunctionBase.hpp
     FunctionManager.hpp
     TableFunction.hpp
   )

#
# Specify all sources
#
set( functions_sources
     FunctionBase.cpp
     FunctionManager.cpp
     TableFunction.cpp
     MultivariableTableFunction.cpp
   )

if( ENABLE_MATHPRESSO )
  list( APPEND functions_headers
        SymbolicFunction.hpp
        CompositeFunction.hpp )
  list( APPEND functions_sources
        SymbolicFunction.cpp
        CompositeFunction.cpp )
endif()


set( dependencyList ${parallelDeps} dataRepository )

if( ENABLE_MATHPRESSO )
   list( APPEND dependencyList mathpresso )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       functions
                 SOURCES    ${functions_sources}
                 HEADERS    ${functions_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

target_include_directories( functions PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS functions LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif()
