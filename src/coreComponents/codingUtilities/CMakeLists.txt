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
 Package: codingUtilities

 Contains definition of runtime types. Also contains common traits.
#]]

#
# Specify all headers
#
set( codingUtilities_headers
     RTTypes.hpp
     Parsing.hpp
     SFINAE_Macros.hpp
     UnitTestUtilities.hpp
     Utilities.hpp
     traits.hpp
   )

#
# Specify all sources
#
set( codingUtilities_sources
     Parsing.cpp
     RTTypes.cpp
   )

set( dependencyList ${parallelDeps} common fast_float )
geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       codingUtilities
                 SOURCES    ${codingUtilities_sources}
                 HEADERS    ${codingUtilities_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )

# Avoid compiling with nvcc which sometimes crashes on fast_float
set_source_files_properties( Parsing.cpp PROPERTIES LANGUAGE CXX )

target_include_directories( codingUtilities PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS codingUtilities LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )


if( GEOS_ENABLE_TESTS )
  add_subdirectory(tests)
endif( )
